// -*- mode:groovy -*-
// vim: syntax=nextflow

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
*/

// Helper closure to check input files
fileExists = { fn ->
              if (fn.exists())
                  return fn;
              else
                  error("File not found: $fn")
}

getCanonicalFileType = { name ->
    def result
    def filename = name.toString()
    switch(filename) {
         // extension is (traditionally) txt, so simple extraction is not possible
        case ~/.*flag.relatives.txt/:
            result = "relatives"
            break
        case ~/.*annotation.txt/:
            result = "annotation"
            break
        case ~/.*country.eval/:
            result = "country_eval"
            break
        case ~/.*nearly_monomorphic.txt/:
            result = "nearly_monomorphic"
            break
        case ~/.*monomorphic.txt/:
            result = "monomorphic"
            break
        case ~/.*country.pca.evec.*/:
            result = "country_evec"
            break
        case ~/.*.pca.evec.*/:
            result = "evec"
            break
        // catches bim,bed,fam,log,eval,evec and all the other "simple" types
        default:
            result = filename.drop(filename.lastIndexOf('.')+1) // if the last '.' is on index 10, drop the first 11 characters, leaving the rest
            break
    }
//    println "Filing $filename in $result"
    result
}

// Transform a list of file names (strings) into a map with the file extension as key and the filename as value
mapFileList = { fn -> fn.collectEntries {
        [(getCanonicalFileType(it)): it]
    }
}

for_snprelate_prune = Channel.create()
for_snprelate_ann = Channel.create()
for_snprelate_ann_atcg = Channel.create()
for_draw_final_pca_histograms_ds = Channel.create()
for_final_pca_1kg_frauke_ann = Channel.create()
for_twstats_final_pruned_ann = Channel.create()
for_sex_check = Channel.create()
for_det_monomorphics_final = Channel.create()
for_plot_maf = Channel.create()
for_prepare_imputation = Channel.create()

Channel.from([[params.bed, params.bim, params.fam, params.individuals_annotation, params.evec ].collect { fileExists(file(it)) }])
       .separate(for_snprelate_prune, for_draw_final_pca_histograms_ds, for_sex_check, for_det_monomorphics_final,for_plot_maf,for_prepare_imputation) {a -> [a,a,a,a,a,a] }

Channel.from(fileExists(file(params.individuals_annotation)))
       .separate(for_snprelate_ann, for_snprelate_ann_atcg, for_final_pca_1kg_frauke_ann, for_twstats_final_pruned_ann) {a -> [a,a,a,a] }

// part 4
process prune_final {
    validExitStatus 0,128
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true, pattern: '*.prune.{in,out,out.unknown_variants}'
    time '2 h'
    tag "${params.collection_name}"
    input:
    file ds_staged from for_snprelate_prune

    output:
    file "$prefix{.bed,.bim,.fam,.log}" into for_snprelate, for_twstats_final_pruned, for_merge_1kg_pruned_final
    file "${prefix}_with_atcg{.bed,.bim,.fam,.log}" into for_snprelate_atcg
    file "${params.disease_data_set_prefix_release}.prune.in"
    file "${params.disease_data_set_prefix_release}.prune.out"
    file "${params.disease_data_set_prefix_release}.prune.out.unknown_variants"

    shell:
    dataset = mapFileList(ds_staged)
    prefix = params.disease_data_set_prefix_release + "_pruned"
    println "DATASET: ${dataset}"

    if (params.PCA_SNPList != "" && params.PCA_SNPList != "nofileexists") {
'''
    module load "IKMB"
    module load "Plink/1.9"
echo Using PCA SNP List file and sample outliers for variant selection
plink --bfile !{dataset.bed.baseName} --extract "!{params.PCA_SNPList}" --make-bed --out !{prefix}_pruned --allow-no-sex
'''
    } else {
'''
    module load "IKMB"
    module load "Plink/1.9"

<!{dataset.bim} tr -s '\\t ' ' ' | cut -f2 -d' ' | grep ^unk_ >unknowns
plink --bim !{dataset.bim} --bed !{dataset.bed} --fam !{dataset.fam} --exclude unknowns --make-bed --out no-unknowns

plink --bfile no-unknowns --indep-pairwise 50 5 0.2 --out after-indep-pairwise --allow-no-sex
plink --bfile no-unknowns --extract after-indep-pairwise.prune.in --maf 0.05 --make-bed --out after-correlated-remove --allow-no-sex
mv after-indep-pairwise.prune.in "!{params.disease_data_set_prefix_release}.prune.in"
mv after-indep-pairwise.prune.out "!{params.disease_data_set_prefix_release}.prune.out"
mv unknowns "!{params.disease_data_set_prefix_release}.prune.out.unknown_variants"

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("no-unknowns.bim", "include-variants")'
python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noIndels("no-unknowns.bim", "include-variants-with-atcg")'
plink --bfile after-correlated-remove --extract include-variants --make-bed --out "!{prefix}" --allow-no-sex
plink --bfile after-correlated-remove --extract include-variants-with-atcg --make-bed --out "!{prefix}_with_atcg" --allow-no-sex
'''
    }
}

process final_pca_con_projection {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    memory {25.GB * task.attempt}
    errorStrategy 'retry'
    tag "${params.collection_name}"

    input:
    file ds_pruned_staged from for_snprelate
    file ds_final_staged from for_snprelate_ann

    output:
    file "${params.disease_data_set_prefix_release}_final.{pca.evec,eval}" into for_draw_final_pca_histograms, for_twstats_final_pruned_pcaresults
    file "*.png"

    shell:
    dataset = mapFileList(ds_pruned_staged)
    ds_final = mapFileList(ds_final_staged)
    prefix = params.disease_data_set_prefix_release + "_final"
    draw_evec_FLASHPCA2 = SCRIPT_DIR + "/draw_evec_FLASHPCA2.r"
    pcaplot_1KG = SCRIPT_DIR + "/pcaplot_1KG_v2.R"

//    gds_script = SCRIPT_DIR + "/SNPRelate_convert2gds.r"
//    snprelate_script = SCRIPT_DIR + "/SNPRelate_PCA_32PCAs.r"
//    projection_samples = params.disease_data_set_prefix_release + "_pruned_PCAprojection_control_samples.txt"
'''
    module load "IKMB"
    module load "FlashPCA/2.0"
    module load "Eigensoft/4.2"
flashpca2 -d 10 --bfile "!{dataset.bed.baseName}" \
    --outval !{prefix}_eigenvalues_flashpca2 \
    --outvec !{prefix}_eigenvectors_flashpca2 \
    --outpc  !{prefix}_pcs_flashpca2 \
    --numthreads !{task.cpus} \
    --outload !{prefix}_loadings_flashpca2 \
    --outmeansd !{prefix}_meansfd_flashpca2 \
    --memory $((!{task.memory.toMega()}-1000))
#    --memory 64000

echo Adding batch info | ts
python -c 'from SampleQCI_helpers import *; addphenoinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", "!{ds_final.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info | ts
R --slave --args "!{prefix}" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"

echo Adding country info | ts
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.country.pca.evec", "!{prefix}.country.eval", "!{ds_final.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with country info | ts
R --slave --args "!{prefix}.country" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"


'''
}


process final_pca_con_projection_atcg {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    memory {12.GB*task.attempt}
    tag "${params.collection_name}"

    input:
    file ds_pruned_staged from for_snprelate_atcg
    file ds_final_staged from for_snprelate_ann_atcg

    output:
    file "${params.disease_data_set_prefix_release}_final_atcg.{pca.evec,eval}"
    file "*.png"

    shell:
    dataset = mapFileList(ds_pruned_staged)
    ds_final = mapFileList(ds_final_staged)
    prefix = params.disease_data_set_prefix_release + "_final_atcg"
    draw_evec_FLASHPCA2 = SCRIPT_DIR + "/draw_evec_FLASHPCA2.r"
    pcaplot_1KG = SCRIPT_DIR + "/pcaplot_1KG_v2.R"

//    gds_script = SCRIPT_DIR + "/SNPRelate_convert2gds.r"
//    snprelate_script = SCRIPT_DIR + "/SNPRelate_PCA_32PCAs.r"
//    projection_samples = params.disease_data_set_prefix_release + "_pruned_PCAprojection_control_samples.txt"
'''
    module load "IKMB"
    module load "FlashPCA/2.0"
    module load "Eigensoft/4.2"
flashpca2 -d 10 --bfile "!{dataset.bed.baseName}" \
    --outval !{prefix}_eigenvalues_flashpca2 \
    --outvec !{prefix}_eigenvectors_flashpca2 \
    --outpc  !{prefix}_pcs_flashpca2 \
    --numthreads !{task.cpus} \
    --outload !{prefix}_loadings_flashpca2 \
    --outmeansd !{prefix}_meansfd_flashpca2 \
    --memory $((!{task.memory.toMega()}-1000))
#    --memory 64000

echo Adding batch info | ts
python -c 'from SampleQCI_helpers import *; addphenoinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", "!{ds_final.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info | ts
R --slave --args "!{prefix}" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"

echo Adding country info | ts
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.country.pca.evec", "!{prefix}.country.eval", "!{ds_final.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with country info | ts
R --slave --args "!{prefix}.country" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"


'''
}


process draw_final_pca_histograms {
    when:
    params.run_final_snprelate == true

    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    cpus 2
    tag "${params.collection_name}"

    input:
    file ds_staged from for_draw_final_pca_histograms_ds
    file pca_staged from for_draw_final_pca_histograms

    output:
    file "*.png"

    shell:
    ds = mapFileList(ds_staged)
    pca = mapFileList(pca_staged)
    evec_script = SCRIPT_DIR + "/draw_evec_FAMStyle.r"
    con_case_script = SCRIPT_DIR + "/draw_histos_CON_CASE_FAMStyle.r"
    con_all_script = SCRIPT_DIR + "/draw_histos_CON_PS_AS_CD_UC_PSC_FAMStyle.r"
'''
R --slave --args "!{pca.eval.baseName}" < "!{evec_script}"
if [ "!{params.hf_test_CON_only}" == "True" ]; then
    R --slave --args "!{pca.evec}" "!{ds.annotation}" < "!{con_case_script}"
else
    R --slave --args "!{pca.evec}" "!{ds.annotation}" < "!{con_all_script}"
fi
'''
}

process final_merge_pruned_with_1kg {
    tag "${params.collection_name}"

    input:

    file dataset_pruned_staged from for_merge_1kg_pruned_final

    output:

    file "*_1kG.{bed,bim,fam,log}" into for_pca_plot_1KG_frauke_final



    shell:

    dataset = mapFileList(dataset_pruned_staged)
    base_pruned_1kG = "${dataset.bed.baseName}_1kG"
    if(params.PCA_SNPexcludeList == "") {
        snpexclude = ""
    } else {
        snpexclude = "${params.PCA_SNPexcludeList}"
    }
'''
    module load "IKMB"
    module load "Plink/1.9"
echo Merge with 1kG

DONE=0
BASE_PRUNED="!{dataset.bed.baseName}"

while [ $DONE -lt 1 ]
do
    python -c 'from SampleQCI_helpers import *; merge__new_plink_collection_pruned__1kG("'$BASE_PRUNED'", "!{base_pruned_1kG}", "!{snpexclude}", "!{params.preQCIMDS_1kG}")' || true

    if [ -e "!{base_pruned_1kG}-merge.missnp" ]; then
        NEW_PRUNED=${BASE_PRUNED}_clean
        rm -f removelist
        while read f
        do
            CHR=$(echo $f | cut -f1 -d:)
            POS=$(echo $f | cut -f2 -d:)
            grep -E "^$CHR\\\\s.*$POS" ${BASE_PRUNED}.bim | cut -f2 -d$'\\t' >>removelist
            echo "$f" >>removelist
        done <"!{base_pruned_1kG}-merge.missnp"
        plink --bfile "$BASE_PRUNED" --exclude removelist --make-bed --out "$NEW_PRUNED"
        BASE_PRUNED=$NEW_PRUNED
        mv "!{base_pruned_1kG}-merge.missnp" "!{base_pruned_1kG}.missnp.removed"
    else
        DONE=1
    fi
done
'''
}

process pca_plot_1kg_frauke_final {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    memory {25.GB * task.attempt}
    tag "${params.collection_name}"

    input:
    file dataset_staged from for_pca_plot_1KG_frauke_final
    file ds_original from for_final_pca_1kg_frauke_ann

    output:
    file "*.pdf"

    shell:
    dataset = mapFileList(dataset_staged)
    dataset_orig = mapFileList(ds_original)

    pcaplot_1KG = NXF_DIR + "/bin/pcaplot_1KG.R"
    prefix = dataset.bed.baseName
'''
# Rscript $NXF_DIR/bin/pcaplot_1KG.R "!{dataset.bim.baseName}" "!{params.numof_pc}" "!{params.preQCIMDS_1kG_sample}"

# Run PCA

    module load "IKMB"
    module load "FlashPCA/2.0"
    module load "Eigensoft/4.2"
flashpca2 -d 10 --bfile "!{dataset.bed.baseName}" \
    --outval !{prefix}_eigenvalues_flashpca2 \
    --outvec !{prefix}_eigenvectors_flashpca2 \
    --outpc  !{prefix}_pcs_flashpca2 \
    --numthreads !{task.cpus} \
    --outload !{prefix}_loadings_flashpca2 \
    --outmeansd !{prefix}_meansfd_flashpca2 \
    --memory $((!{task.memory.toMega()}-1000))
#    --memory 64000

echo Adding batch info | ts

python -c 'from SampleQCI_helpers import *; addphenoinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", \
 "!{dataset_orig.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info
Rscript "!{pcaplot_1KG}" "!{prefix}" 5 "!{params.preQCIMDS_1kG_sample}"

'''
}


process twstats_final_pruned {
    tag "${params.collection_name}"
    when:
    params.run_final_snprelate == true

    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:

    file dataset_final_staged from for_twstats_final_pruned_ann
    file dataset_pruned_staged from for_twstats_final_pruned
    file pca_staged from for_twstats_final_pruned_pcaresults

    output:

    file "*.tracy_widom"

    shell:

    dataset = mapFileList(dataset_pruned_staged)
    annotation = mapFileList(dataset_final_staged).annotation
    pca = mapFileList(pca_staged)

    //    template "run_tracy_widom_stats.sh"
    '''
    module load 'IKMB'
    module load 'Eigensoft/4.2'
# NUM_CASES=$(grep -P 'After.*\\d+.cases' !{dataset.log} | cut -d' ' -f3)
# NUM_CONTROLS=$(grep -P 'After .*\\d+.controls' !{dataset.log} | cut -d' ' -f5)
NUM_CASES=$(grep controls !{dataset.log} | cut -d' ' -f4)
NUM_CONTROLS=$(grep controls !{dataset.log} | cut -d' ' -f8)
echo Cases: $NUM_CASES, controls: $NUM_CONTROLS
NUM_SAMPLES=$(($NUM_CASES + $NUM_CONTROLS))
NUM_SNPS=$(grep 'pass filters and QC' !{dataset.log} | cut -d' ' -f1)
echo Samples: $NUM_SAMPLES, markers: $NUM_SNPS

head -n 2000 !{pca.eval} >eval.tw

if [ "!{params.projection_on_populations_CON_only}" == "True" ]; then
    if [ "$NUM_CONTROLS" -gt "$NUM_SNPS" ]; then
        let NUM_CONTROLS--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom -n "$NUM_CONTROLS"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom
    fi
else
    if [ "$NUM_SAMPLES" -gt "$NUM_SNPS" ]; then
        let NUM_SAMPLES--
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom -n "$NUM_SAMPLES"
    else
        twstats -t !{params.twtable} -i eval.tw -o !{dataset.bed.baseName}.eval.tracy_widom
    fi
fi
'''
}

process sex_check {
    publishDir params.qc_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file ds_staged from for_sex_check

    output:
    file("*.pdf") optional true
    file("sex-check-not-possible.txt") optional true

    shell:
    ds = mapFileList(ds_staged)
'''
module load IKMB
module load Plink/1.9

# Check if we have data on both X and Y chromosomes
plink --allow-no-sex --bfile "!{ds.bim.baseName}" --missing --out !{ds.bim.baseName}
plink --allow-no-sex --bfile "!{ds.bim.baseName}" --het --out !{ds.bim.baseName}
plink --allow-no-sex --bfile "!{ds.bim.baseName}" --chr X --recode --out !{ds.bim.baseName}_sex.X || true
plink --allow-no-sex --bfile "!{ds.bim.baseName}" --chr Y --recode --out !{ds.bim.baseName}_sex.Y || true

if [ -e "!{ds.bim.baseName}"_sex.X.ped -a -e "!{ds.bim.baseName}"_sex.Y.ped ]; then
    Rscript $NXF_DIR/bin/indivplot.R "!{ds.bed.baseName}" "!{ds.bim.baseName}"_sex.X.ped "!{ds.bim.baseName}"_sex.Y.ped
else
    echo "Sex check requires X and Y chromosomes to be present in the dataset" | tee sex-check-not-possible.txt
fi
rm -f "!{ds.bim.baseName}.log"

'''
}

process det_monomorphics_final {
publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    tag "${params.collection_name}"
    input:

    file dataset_staged from for_det_monomorphics_final // SNPQCII_final

    output:

    //file "${params.collection_name}_SNPQCII.flag.{nearly_,}monomorphic"
    file "${params.collection_name}_QCed_final.flag.{nearly_,}monomorphic"
    // into for_compile_variants_exclusion_monomorphics

    shell:

    dataset = mapFileList(dataset_staged)
    prefix = "${params.collection_name}_QCed_final"
'''
module load IKMB
module load Plink/1.7

echo "PYLIB_DIR: $PYLIB_DIR"
echo "PYTHONPATH: $PYTHONPATH"

# determine monomorphic variants
plink --noweb --bfile "!{dataset.bed.baseName}" --freq --out "!{prefix}_freq" --allow-no-sex
python -c 'from SNPQC_helpers import *; frq =  Frq(frq_file="!{prefix}_freq.frq", write_monomorphic_file="!{prefix}.flag.monomorphic"); frq.write_monomorphic_variants_file(); del frq'

# determine nearly monomorphic variants
if [ "!{params.hf_test_CON_only}" != "True" ]; then
    grep -v -e "Unknown" "!{dataset.annotation}" >clean-annotations
    module switch Plink/1.7 Plink/1.9
    plink --bfile "!{dataset.bed.baseName}" --freq --mwithin 7 --within clean-annotations --out "!{prefix}_freq" --allow-no-sex
    python -c 'from SNPQC_helpers import *; \
        frq = FrqStrat(frq_file="!{prefix}_freq.frq.strat", write_nearly_monomorphic_file="!{prefix}.flag.nearly_monomorphic", maf_thresh=!{params.maf_thresh}); \
        frq.write_nearly_monomorphic_variants_PS_AS_IBD_PSC_file()'
else
    touch "!{prefix}.flag.nearly_monomorphic"
fi

'''
}


process plot_maf {
    publishDir params.qc_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"

    input:
    file ds_staged from for_plot_maf

    output:
    file "*.{frq,png}"

    shell:
    ds = mapFileList(ds_staged)
    logmaf = SCRIPT_DIR + "/logmaf.r"
'''
module load IKMB
module load Plink/1.9




    plink --bfile "!{ds.bim.baseName}" --freq --out "!{ds.bim.baseName}_freq" --allow-no-sex
    R --slave --args "!{ds.bim.baseName}_freq.frq" <"!{logmaf}"
'''
}


process prepare_sanger_imputation {
    publishDir params.qc_dir ?: '.', mode: 'copy'
    tag "${params.collection_name}"
    time {8.h * task.attempt}

    input:
    file ds_staged from for_prepare_imputation

    output:
    file "${ds.bim.baseName}.vcf.refchecked.gz"
    file "${ds.bim.baseName}.vcf.refchecked.gz.tbi"

    shell:
    ds = mapFileList(ds_staged)

    '''
module load IKMB
module load Plink/1.9

ANNOTATION=/ifs/data/nfs_share/sukmb388/human_g1k_v37.fasta.gz

grep D !{ds.bim} | awk '{ print $2 }' >!{ds.bim}.indels.txt
plink --allow-no-sex --bfile !{ds.bim.baseName} --chr 1-22 --exclude !{ds.bim}.indels.txt --make-bed --out final_noindels
plink --allow-no-sex --bfile final_noindels --recode vcf --out !{ds.bim.baseName}_tmp1
bgzip !{ds.bim.baseName}_tmp1.vcf || true

echo "Check the REF allele ...";
bcftools norm --check-ref w -f $ANNOTATION !{ds.bim.baseName}_tmp1.vcf.gz >/dev/null || true
echo "Fix the REF allele ...";
bcftools norm --check-ref s -f $ANNOTATION !{ds.bim.baseName}_tmp1.vcf.gz >!{ds.bim.baseName}.vcf.refchecked || true
bgzip !{ds.bim.baseName}.vcf.refchecked || true

echo "Tabix ...";
tabix -p vcf !{ds.bim.baseName}.vcf.refchecked.gz

    '''
}




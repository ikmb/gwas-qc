// -*- mode:groovy -*-

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


// initialize configuration
params.qc_dir = "."
params.PCA_SNPList = ""

// match auto-generated "no file exists" to actual not-existing files
if (params.PCA_SNPList == "nofileexists") {
    params.PCA_SNPList = ""
}
if (params.PCA_SNPexcludeList == "nofileexists") {
    params.PCA_SNPexcludeList = ""
}

// default value, will be overwritten by config
params.projection_on_populations_CON_only = "False"

// original, "final" SampleQCI dataset
SampleQCI_final = ["${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.bed", 
                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.bim", 
                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.fam", 
                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_flag_relatives.txt", 
                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_annotation.txt"].collect { fileExists(file(it)) }

// based on "original" but without relatives
SampleQCI_final_wr = ["${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.bed", 
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.bim", 
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.fam", 
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.annotation.txt", 
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.pca.evec" ].collect { fileExists(file(it)) }

for_prune = Channel.create()
for_flashpca_1kg_final = Channel.create()
for_final_sample_cleaning_outliers = Channel.create()
for_final_sample_cleaning = Channel.create()

/*
 The number of markers that the HF test will be performed on at once.
 If you use 100 markers out of 10,000 markers, this will queue up
 100 jobs to be run in parallel.
 Default: 100 (as extracted from SNP_QCII_CON_PS_AS_CD_UC_PSC_parallel_part1.py)
*/
hf_test_chunk_size = 1000

/*

// from SNP_QCII_CON_PS_AS_CD_UC_PSC_parallel_part1.py
process hf_test_prepare {

    // It is okay if kill could not successfully kill Rserve, it might have died on its own,
    // so allow an error-indicating return value.
    validExitStatus 0,1


    input:

    // get them staged
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()

    output:
    file 'chunk_*' into hf_test_chunk
    file 'HF_PCA.R' into hf_test_script

    shell:

    dataset = mapFileList(SampleQCI_final_wr_staged)
    println dataset
'''
    module load "IKMB"
    module load "Plink/1.9"
# hf_test(plink, new_plink, pcfile, qced_annotations)

# Generate script file from template
TEMPLATE="!{SCRIPT_DIR + '/template_HF_PCA.CON_PS_AS_CD_UC_PSC.R'}"
if [ "!{params.hf_test_CON_only}" == "True" ]; then
    TEMPLATE="!{SCRIPT_DIR + '/template_HF_PCA.CON.R'}"
fi

# Replace placeholders in template. Note the use of '|' in the regex,
# so we don't need to escape the replacement parameter.
sed 's|NUMOFPCS|!{params.numof_pc}|g' $TEMPLATE >HF_PCA.R
sed 's|INDIVIDUALS_ANNOTATION|!{dataset.annotation}|g' HF_PCA.R | TMPDIR=. sponge HF_PCA.R
sed 's|PCAEVECFILE|!{dataset.evec}|g' HF_PCA.R | TMPDIR=. sponge HF_PCA.R

# Likely not necessary
# chmod u+x HF_PCA.R

# Print batches information
echo Found the following batches:
cut -f7 "!{dataset.annotation}" | tail -n +2 | sort -u

cut -f2 "!{dataset.bim}" | split -l !{hf_test_chunk_size} -a 4 -d - chunk_
'''
}


process hf_test {
    input:
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()
    file r_script from hf_test_script
    file(chunk)from hf_test_chunk.flatten()

    output:
    file "${params.collection_name}_SNPQCII_${chunk}.auto.R" into hf_test_merge_chunks


    tag { chunk } // append current chunk to job name

    shell:
    dataset = mapFileList(SampleQCI_final_wr_staged)
'''
    module load "IKMB"
    module load "Plink/1.9"
# Spawn Rserve. Detaches into background before setting up socket, so wait for socket to appear.
R CMD Rserve --RS-socket /scratch/rserve.sock --no-save --RS-pidfile /scratch/rserve.pid
while [ ! -e /scratch/rserve.sock ]
do
    sleep 0.5
done
RSERVE_PID=$(cat /scratch/rserve.pid)

echo Checking prerequisites...
head -n2 !{dataset.annotation} 2>&1 >/dev/null
echo "!{dataset.annotation} exists: $?  (should be 0)"
echo -n /scratch/rserve.sock exists and is readable:
if [ -S /scratch/rserve.sock ]; then
   echo -n is a socket,
else
   echo -n is not a socket,
fi
if [ -r /scratch/rserve.sock ]; then
   echo is readable
else
   echo is not readable
fi

MYPWD=$(pwd)
echo "setwd('$MYPWD')" >hf-local.r
cat !{r_script} >>hf-local.r

plink --bfile "!{dataset.bim.baseName}" --extract "!{chunk}" --R hf-local.r --R-socket /scratch/rserve.sock --allow-no-sex --out "!{params.collection_name}_SNPQCII_!{chunk}"

# Kill Rserve but give it a chance...
sleep 0.5
if [ -e /scratch/rserve.sock ]; then
    # Might fail as it could have been died by now but that's okay.
    kill $RSERVE_PID
fi
'''
}

// SNP_QCII_CON_PS_AS_CD_UC_PSC_parallel_part2.py starts here

process hf_test_merge {
        publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:

    file chunks from hf_test_merge_chunks.collect()
    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()

    shell:

    dataset = mapFileList(SampleQCI_final_staged)
'''
CHUNKS=$(ls -X !{params.collection_name}_SNPQCII_chunk_*)
OUTFILE="!{params.collection_name}_SNPQCII_hf"
BIMFILE="!{dataset.bim}"

for chunk in $CHUNKS; do
    cat $chunk >>$OUTFILE
done

# Count number of variants in split/merged and original files
LINES_MERGED=$(wc -l $OUTFILE | cut -d " " -f 1)
LINES_ORIG=$(wc -l $BIMFILE} | cut -d " " -f 1)

if [ $LINES_MERGED ne $LINES_ORIG ]; then
    echo Unexpected difference in merged HF file.
    echo Expected $LINES_ORIG variants in $OUTFILE
    echo Observed $LINES_MERGED variants in $OUTFILE
fi

# Count number if NAs in split/merged and original files
NAS_MERGED=$(gawk '{ print $5 }' <$OUTFILE | grep NA | wc -l | cut -d " " -f 1)
NAS_ORIG=$(grep NA <$BIMFILE | wc -l | cut -d " " -f 1)

if [ $NAS_MERGED ne $NAS_ORIG ]; then
   echo Unexpected difference in 'NA' valiues in HD file.
   echo Expected $NAS_ORIG NAs in $OUTFILE
   echo Observed $NAS_MERGED NAs in $OUTFILE
fi
'''
}
*/

// see part2.py: "TODO mit jan kaessens: auskommentieren wenn nur eine kontrollgruppe"

process generate_hf_excludes {
    input:

    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()

    output:

    file 'hf-excludes' into for_exclude_variants

    shell:
    dataset = mapFileList(SampleQCI_final_staged)
'''
touch hf-excludes
'''
}


process exclude_variants {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:

    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()
    file exclude from for_exclude_variants

    output:

    file "*SNPQCII_final{.bed,.bim,.fam,.log,_flag.relatives.txt,_annotation.txt}" into for_final_cleaning, for_det_monomorphics, for_det_unknown_diagnosis
    file "*test_missingness.*" into for_det_diff_missingness

    shell:

    prefix = params.collection_name + "_SNPQCII_final"
    dataset = mapFileList(SampleQCI_final_staged)
//    println dataset
'''
    module load "IKMB"
    module load "Plink/1.7"
plink --noweb --bfile "!{dataset.bed.baseName}" --exclude "!{exclude}" --make-bed --out "!{prefix}" --allow-no-sex
cp "!{dataset.relatives}" "!{prefix}_flag.relatives.txt"
cp "!{dataset.annotation}" "!{prefix}_annotation.txt"

plink --noweb --bfile "!{dataset.bed.baseName}" --test-missing --out !{prefix}_test_missingness

'''
}


process prune {
    input:

    file dataset_staged from for_prune.collect() // output from exclude_variants

    output:

    file "*SNPQCII_final_pruned{.bed,.bim,.fam,.log}" into for_flashpca_pruned, for_draw_histograms_pruned, for_twstats_pruned, for_merge_1kg_pruned



    shell:

    dataset = mapFileList(dataset_staged)
    prefix = params.collection_name + "_SNPQCII_final"
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
    module load "Plink/1.7"

plink --noweb --bfile !{dataset.bed.baseName} --indep-pairwise 50 5 0.2 --out _prune --allow-no-sex
plink --noweb --bfile !{dataset.bed.baseName} --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate --allow-no-sex

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{dataset.bim}", "include_variants")'
plink --noweb --bfile intermediate --extract include_variants --make-bed --out "!{prefix}_pruned" --allow-no-sex
'''
    }
}



process merge_pruned_with_1kg {

    input:

    file dataset_pruned_staged from for_merge_1kg_pruned

    output:

    file "*_1kG.{bed,bim,fam,log}" into for_flashpca_1kg_merged,for_pca_plot_1KG_frauke



    shell:

    dataset = mapFileList(dataset_pruned_staged)
    base_pruned_1kG = "${dataset.bed.baseName}_1kG"
    if(params.PCA_SNPexcludeList == "") {
        snpexclude = ""
    } else {
        snpexclude = BATCH_DIR + "/${params.PCA_SNPexcludeList}"
    }
'''
    module load "IKMB"
    module load "Plink/1.7"
echo Merge with 1kG
python -c 'from SampleQCI_helpers import *; merge__new_plink_collection_pruned__1kG("!{dataset.bed.baseName}", "!{base_pruned_1kG}", "!{snpexclude}", "!{params.preQCIMDS_1kG}")'
'''
}

process flashpca_1kg {

    input:

    file dataset_staged from for_flashpca_1kg_merged // pruned and merged with 1kG
    file dataset_final from for_flashpca_1kg_final // unpruned


    output:

    file "*.{evec,eval}" into for_pca_plot_1kg_pcaresults, for_pca_plot_1KG_frauke_evec



    shell:

    dataset = mapFileList(dataset_staged)
    annotation = mapFileList(dataset_final).annotation
    prefix = params.collection_name + "_SNPQCII_final_pruned_1kG"
    evec = prefix + "_" + params.numof_pc + "PC.pca.evec"
    eval = prefix + "_" + params.numof_pc + "PC.eval"
    draw_evec_FLASHPCA2 = SCRIPT_DIR + "/draw_evec_FLASHPCA2.r"

    //    template "run_flashpca.sh"
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
    --memory !{task.memory.toMega()}
#    --memory 64000

echo Adding batch info | ts
python -c 'from SampleQCI_helpers import *; addphenoinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", "!{annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info | ts
R --slave --args "!{dataset.bed.baseName}" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"

echo Adding country info | ts
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.country.pca.evec", "!{prefix}.country.eval", "!{annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with country info | ts
R --slave --args "!{dataset.bed.baseName}.country" "!{params.preQCIMDS_1kG_sample}" <"!{draw_evec_FLASHPCA2}"

'''
}

process pca_plot_1kg_frauke {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:
    file dataset_staged from for_pca_plot_1KG_frauke
    file pcaresults_staged from for_pca_plot_1KG_frauke_evec

    output:
    file "*.pdf"
    file "outliers" into for_final_sample_cleaning_outliers

    shell:
    dataset = mapFileList(dataset_staged)
    // pcaresults = mapFileList(pcaresults_staged)
    plink_pca_1kG = mapFileList(pcaresults_staged).eval.baseName
    pcaplot_1KG = NXF_DIR + "/bin/pcaplot_1KG.R"
'''
# Rscript $NXF_DIR/bin/pcaplot_1KG.R "!{dataset.bim.baseName}" "!{params.numof_pc}" "!{params.preQCIMDS_1kG_sample}"

echo Drawing FLASHPCA2 eigenvectors for set with batch info
echo Calling R --slave --args "!{plink_pca_1kG}" "!{params.preQCIMDS_1kG_sample}"  lt "!{pcaplot_1KG}"
# R --slave --args "!{plink_pca_1kG}" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"

Rscript "!{pcaplot_1KG}" "!{plink_pca_1kG}" 5 "!{params.preQCIMDS_1kG_sample}"
mv fail-pca-1KG-qc.txt fail-pca-1KG-qc.txt-batch

echo Drawing FLASHPCA2 eigenvectors for set with country info
echo Calling R --slave --args "!{plink_pca_1kG}.country" "!{params.preQCIMDS_1kG_sample}" lt "!{pcaplot_1KG}"
# R --slave --args "!{plink_pca_1kG}.country" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"
Rscript "!{pcaplot_1KG}" "!{plink_pca_1kG}.country" 5 "!{params.preQCIMDS_1kG_sample}"
mv fail-pca-1KG-qc.txt  fail-pca-1KG-qc.txt-country

cat fail-pca-1KG-qc.txt-batch fail-pca-1KG-qc.txt-country | sort | uniq >outliers
'''
}
/*
process pca_plot_1kg {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    file pcaresults_staged from for_pca_plot_1kg_pcaresults

    output:

    file "*.pdf"
    file "*.txt"

    shell:

    plink_pca_1kG = mapFileList(pcaresults_staged).eval.baseName
    pcaplot_1KG = NXF_DIR + "/bin/pcaplot_1KG_v2.R"
'''
echo Drawing FLASHPCA2 eigenvectors for set with batch info
echo Calling R --slave --args "!{plink_pca_1kG}" "!{params.preQCIMDS_1kG_sample}"  lt "!{pcaplot_1KG}"
R --slave --args "!{plink_pca_1kG}" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"

echo Drawing FLASHPCA2 eigenvectors for set with country info
echo Calling R --slave --args "!{plink_pca_1kG}.country" "!{params.preQCIMDS_1kG_sample}" lt "!{pcaplot_1KG}"
R --slave --args "!{plink_pca_1kG}.country" "!{params.preQCIMDS_1kG_sample}" <"!{pcaplot_1KG}"
'''
}
*/

// part 3 begins here


// (1) and (2): determine monomorphic variants and nearly monomorphic variants
process det_monomorphics {
    input:

    file dataset_staged from for_det_monomorphics // SNPQCII_final

    output:

    file "${params.collection_name}_SNPQCII.{nearly_,}monomorphic" into for_compile_variants_exclusion_monomorphics

    shell:

    dataset = mapFileList(dataset_staged)
    prefix = "${params.collection_name}_SNPQCII"
'''
module load IKMB
module load Plink/1.7

echo "PYLIB_DIR: $PYLIB_DIR"
echo "PYTHONPATH: $PYTHONPATH"

# determine monomorphic variants
plink --noweb --bfile "!{dataset.bed.baseName}" --freq --out "!{prefix}_freq" --allow-no-sex
python -c 'from SNPQC_helpers import *; frq =  Frq(frq_file="!{prefix}_freq.frq", write_monomorphic_file="!{prefix}.monomorphic"); frq.write_monomorphic_variants_file(); del frq'

# determine nearly monomorphic variants
if [ "!{params.hf_test_CON_only}" != "True" ]; then
    grep -v -e "Unknown" "!{dataset.annotation}" >clean-annotations
    module switch Plink/1.7 Plink/1.9
    plink --bfile "!{dataset.bed.baseName}" --freq --mwithin 7 --within clean-annotations --out "!{prefix}_freq" --allow-no-sex
    python -c 'from SNPQC_helpers import *; \
        frq = FrqStrat(frq_file="!{prefix}_freq.frq.strat", write_nearly_monomorphic_file="!{prefix}.nearly_monomorphic", maf_thresh=!{params.maf_thresh}); \
        frq.write_nearly_monomorphic_variants_PS_AS_IBD_PSC_file()'
else
    touch "!{prefix}.nearly_monomorphic"
fi
'''
}


process det_diff_missingness {

    prefix = "${params.collection_name}_SNPQCII"

    input:

    file missingness_staged from for_det_diff_missingness // from e

    output:

    file "${prefix}.differential_missingness" into for_compile_variants_exclusion_diff_missingness

    shell:
    prefix = "${params.collection_name}_SNPQCII"

    missingness = mapFileList(missingness_staged)
'''
#!/usr/bin/env python

from SNPQC_helpers import *

test_missing = Test_missing(missing_file="!{missingness.missing}", write_file="!{prefix}.differential_missingness", threshold=!{params.test_missing_thresh})
test_missing.write_variants_file();
del test_missing;

#python -c 'from SNPQC_helpers import *; \
#    test_missing = Test_missing(missing_file="!{missingness.missing}", write_file="!{prefix}.differential_missingness", threshold=!{params.test_missing_thresh}); \
#    test_missing.write_variants_file(); \
#    del test_missing'
'''
}

process det_unknown_diagnosis {
    prefix = "${params.collection_name}_SNPQCII"

    input:

//    file dataset from for_det_unknown_diagnosis

    output:

    file "${prefix}.individuals_remove_final" into for_final_cleaning_individuals

    shell:

    annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"
    prefix = "${params.collection_name}_SNPQCII"
'''
# TODO: Might need to process "!{params.diagnoses}" into an actual list before calling
python -c 'from SNPQC_helpers import *; determine_unknown_diagnosis(annotationfile="!{annotation}", outfile="!{prefix}.unknown_diagnosis", diagnoses="!{params.diagnoses}")'

if [ -e "!{params.individuals_remove_manually}" ]; then
    cp "!{params.individuals_remove_manually}" "!{prefix}.individuals_remove_manually"
else
    touch "!{prefix}.individuals_remove_manually"
fi

gawk '{ print $1, $2 }' "!{prefix}.unknown_diagnosis" "!{prefix}.invididuals_remove_manually" | sort | uniq >"!{prefix}.individuals_remove_final"
'''
}


process compile_variants_exclusion {

    input:

    file differential_missingness from for_compile_variants_exclusion_diff_missingness
    file monomorphics_staged from for_compile_variants_exclusion_monomorphics

    output:

    file "${params.collection_name}.variants_exclude_final" into for_final_cleaning_variants

    shell:
    prefix = "${params.collection_name}_SNPQCII"
    monomorphics = mapFileList(monomorphics_staged)
'''
if [ "!{params.hf_test_CON_only}" == "True" ]; then
    gawk '{ print $1 }' "!{differential_missingness}" | sort | uniq >"!{params.collection_name}.variants_exclude_final"
else
    gawk '{ print $1 }' "!{monomorphics.monomorphic}" "!{monomorphics.nearly_monomorphic}" "!{differential_missingness}" | sort | uniq >"!{params.collection_name}.variants_exclude_final"
fi
'''
}

process final_snp_cleaning {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    // from diff_missingness/det_monomorphics/det_diagnoses to prune/1kg.
    file dataset_staged from for_final_cleaning
    file individuals from for_final_cleaning_individuals
    file variants from for_final_cleaning_variants

    output:

//    file "${params.disease_data_set_prefix_release}{.bed,.bim,.fam,.log,_annotation.txt,_flag.relatives.txt}" into for_snprelate_prune,for_twstats_final_pruned_ann,for_draw_final_pca_histograms_ds,for_plot_maf,for_eigenstrat_convert_ann,for_snprelate_ann,for_twstats_final_pruned_eigenstrat_ann

    file "${params.disease_data_set_prefix_release}{.bed,.bim,.fam,.log,_annotation.txt,_flag.relatives.txt}" into for_prune, for_flashpca_1kg_final, for_final_sample_cleaning
    
    shell:

    annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"
    dataset = mapFileList(dataset_staged)
    prefix = params.disease_data_set_prefix_release
'''
    module load "IKMB"
    module load "Plink/1.7" 

# Remove sample and SNP outliers
plink --noweb --bfile "!{dataset.bed.baseName}" --remove "!{individuals}" --exclude "!{variants}" --make-bed --out "!{prefix}" --allow-no-sex

# Fix annotations
python -c 'from SNPQC_helpers import *; extract_QCsamples_annotationfile_relativesfile( \
    fam="!{prefix}.fam", individuals_annotation_QCed="!{prefix}_annotation.txt", \
    related_samples_file="!{params.collection_name}_SNPQCII_final_flag.relatives.txt", \
    related_samples_file_QCed="!{prefix}_flag.relatives.txt", \
    individuals_annotation="!{annotation}", \
    diagnoses="!{params.diagnoses}")'

touch outliers

'''
}


process final_cleaning {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    // from Frauke's PCA outlier detection
    file sample_outliers from for_final_sample_cleaning_outliers
    file dataset_staged from for_final_sample_cleaning

    output:

    file "${params.disease_data_set_prefix_release}_final{.bed,.bim,.fam,.log,_annotation.txt,_flag.relatives.txt}" into for_snprelate_prune,for_twstats_final_pruned_ann,for_draw_final_pca_histograms_ds,for_plot_maf,for_eigenstrat_convert_ann,for_snprelate_ann,for_twstats_final_pruned_eigenstrat_ann
shell:
    dataset = mapFileList(dataset_staged)
    prefix = params.disease_data_set_prefix_release

'''
module load "IKMB"
module load "Plink/1.9"

# Remove sample outliers
plink --bfile "!{dataset.bed.baseName}" --remove "!{sample_outliers}" --make-bed --out "!{prefix}_final" --allow-no-sex

# Fix annotations
python -c 'from SNPQC_helpers import *; extract_QCsamples_annotationfile_relativesfile( \
    fam="!{prefix}_final.fam", individuals_annotation_QCed="!{prefix}_final_annotation.txt", \
    related_samples_file="!{dataset.relatives}", \
    related_samples_file_QCed="!{prefix}_final_flag.relatives.txt", \
    individuals_annotation="!{dataset.annotation}", \
    diagnoses="!{params.diagnoses}")'
'''
}

// part 4
// TODO: nochmal überlegen, ob diese ganze SNPRelate-Sache überhaupt noch sinnvoll ist

process prune_final {

    input:
    file ds_staged from for_snprelate_prune
    
    output:
    file "$prefix{.bed,.bim,.fam,.log}" into for_snprelate, for_twstats_final_pruned

    shell:
    dataset = mapFileList(ds_staged)
    prefix = params.disease_data_set_prefix_release + "_pruned"

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
    module load "Plink/1.7"

plink --noweb --bfile !{dataset.bed.baseName} --indep-pairwise 50 5 0.2 --out after-indep-pairwise --allow-no-sex
plink --noweb --bfile !{dataset.bed.baseName} --extract after-indep-pairwise.prune.in --maf 0.05 --make-bed --out after-correlated-remove --allow-no-sex

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{dataset.bim}", "include-variants")'
plink --noweb --bfile after-correlated-remove --extract include-variants --make-bed --out "!{prefix}" --allow-no-sex
'''
    }
}

process final_pca_con_projection {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:
    file ds_pruned_staged from for_snprelate
    file ds_final_staged from for_snprelate_ann
    
    output:
    file "${params.disease_data_set_prefix_release}.{pca.evec,eval}" into for_draw_final_pca_histograms, for_twstats_final_pruned_pcaresults
    file "*.png"
    
    shell:
    dataset = mapFileList(ds_pruned_staged)
    ds_final = mapFileList(ds_final_staged)
    prefix = params.disease_data_set_prefix_release
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
    --memory !{task.memory.toMega()}
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


process twstats_final_pruned {
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
NUM_CASES=$(grep -P 'After.*\\d+.cases' !{dataset.log} | cut -d' ' -f3)
NUM_CONTROLS=$(grep -P 'After .*\\d+.controls' !{dataset.log} | cut -d' ' -f5)
echo Cases: $NUM_CASES, controls: $NUM_CONTROLS
NUM_SAMPLES=$(($NUM_CASES + $NUM_CONTROLS))
NUM_SNPS=$(grep -P 'After.*\\d+.SNPs' !{dataset.log} | cut -d' ' -f8)
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



process plot_maf {
    publishDir params.qc_dir ?: '.', mode: 'copy'
    
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

workflow.onComplete {
    println "Generating phase summary..."
    def cmd = ["./generate-phase-summary", "SNPQCII", params.collection_name ?: params.disease_data_set_prefix, workflow.workDir, params.trace_target].join(' ')
    def gensummary = ["bash", "-c", cmd].execute()
    gensummary.waitFor()
}



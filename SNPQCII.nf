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

// from SNP_QCII_CON_PS_AS_CD_UC_PSC_parallel_part1.py
process hf_test_prepare {
   memory 16.GB
    // It is okay if kill could not successfully kill Rserve, it might have died on its own,
    // so allow an error-indicating return value.
    validExitStatus 0,1

    input:

    // get them staged
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()

    output:
    file 'chunk_*' into hf_test_chunk
    file 'HF_PCA.R' into hf_test_script
    file 'cleaned-annotation' into hf_test_anno, hf_test_excludes_anno
    file 'cleaned-evec' into hf_test_evec

    shell:

    dataset = mapFileList(SampleQCI_final_wr_staged)
    println dataset
'''
    module load "IKMB"
    module load "Plink/1.9"
# hf_test(plink, new_plink, pcfile, qced_annotations)
#
# Determine template based on number of batches:
# #cases > 1 and #controls > 1 => HF_CasesControls
# #cases > 1 and #controls =< 1 => HF_Cases
# #cases =< 1 and #controls > 1 => HF_Controls
# #cases =< 1 and #controls =< 1 => skip test

# Fix annotations and evec file to match FAM file
<!{dataset.fam} tr -s '\\t ' ' ' | cut -f2 -d" " >samples_to_be_used
head -n1 !{dataset.annotation} >cleaned-annotation
grep -f samples_to_be_used !{dataset.annotation} >>cleaned-annotation
head -n1 !{dataset.evec} >cleaned-evec
grep -f samples_to_be_used !{dataset.evec} >>cleaned-evec

N_CASES=$(<cleaned-annotation tr -s '\\t ' ' ' | cut -f7,9 -d' ' | uniq | tail -n +2 | cut -f2 -d' ' | grep -v -c Control)
N_CONTROLS=$(<cleaned-annotation tr -s '\\t ' ' ' | cut -f7,9 -d' ' | uniq | tail -n +2 | cut -f2 -d' ' | grep -c Control)

echo "$N_CASES case batches and $N_CONTROLS control batches were detected."

if (( ("$N_CASES" > "4") && ("$N_CONTROLS" > "4") )); then
   TEMPLATE="!{SCRIPT_DIR}/template_HF_PCA.CaseControl.R"
elif [ "$N_CASES" -gt "4" ]; then
   TEMPLATE="!{SCRIPT_DIR}/template_HF_PCA.Case.R"
elif [ "$N_CONTROLS" -gt "4" ]; then
   TEMPLATE="!{SCRIPT_DIR}/template_HF_PCA.CON.R"
else
   echo "At least 5 case or two control batches are required for HF/ANOVA testing."
   touch chunk_0
   touch HF_PCA.R
   exit 0
fi

# Replace placeholders in template.
sed 's|NUMOFPCS|!{params.numof_pc}|g' $TEMPLATE >HF_PCA.R
sed -i 's|INDIVIDUALS_ANNOTATION|cleaned-annotation|g' HF_PCA.R
sed -i 's|PCAEVECFILE|cleaned-evec|g' HF_PCA.R

# Likely not necessary
# chmod u+x HF_PCA.R

# Print batches information
echo Found the following batches:
cut -f7 cleaned-annotation | tail -n +2 | sort -u

cut -f2 "!{dataset.bim}" | split -l !{hf_test_chunk_size} -a 4 -d - chunk_
'''
}


process hf_test {
    // Keine Ahnung, wo das Script "1" produziert. Muss irgendwo beim NA-zaehlen sein
    validExitStatus 0,1
    errorStrategy 'retry'
    memory 4.GB
    input:
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()
    file r_script from hf_test_script
    file(chunk)from hf_test_chunk.flatten()
    file anno from hf_test_anno
    file evec from hf_test_evec

    output:
    file "${params.collection_name}_SNPQCII_${chunk}.auto.R" into hf_test_merge_chunks


    tag { chunk } // append current chunk to job name

    shell:
    dataset = mapFileList(SampleQCI_final_wr_staged)
'''
    module load "IKMB"
    module load "Plink/1.9"

# Empty R script -> no suitable template candidate for HF/ANOVA testing
if [ ! -s "!{r_script}" ]; then
    touch "!{params.collection_name}_SNPQCII_!{chunk}.auto.R"
    exit 0
fi

SOCKET=/scratch/rserve_$$.sock
PIDFILE=/scratch/rserve_$$.pid

RETRY=1
while [ $RETRY -gt 0 ]
do
    echo Start $RETRY
    # Spawn Rserve. Detaches into background before setting up socket, so wait for socket to appear.
    R CMD Rserve --RS-socket $SOCKET --no-save --RS-pidfile $PIDFILE
    while [ ! -e $SOCKET ]
    do
        sleep 0.5
    done
    RSERVE_PID=$(cat $PIDFILE)

    MYPWD=$(pwd)
    echo "setwd('$MYPWD')" >hf-local.r
    cat !{r_script} >>hf-local.r
    echo Script generated $RETRY
    plink --bfile "!{dataset.bim.baseName}" --extract "!{chunk}" --R hf-local.r --R-socket $SOCKET --allow-no-sex --out "!{params.collection_name}_SNPQCII_!{chunk}"
    echo Plink done $RETRY
    # Kill Rserve but give it a chance...
    sleep 3
    if [ -e $SOCKET ]; then
        echo Killing Rserve
        # Might fail as it could have been died by now but that's okay.
        kill -9 $RSERVE_PID || true
        rm -f $SOCKET
    else
        # Socket is closed but maybe the process is still running
        kill -9 $RSERVE_PID || true
        echo Rserve terminated on its own.
    fi

    echo Socket removed or killed $RETRY
    touch !{params.collection_name}_SNPQCII_!{chunk}.auto.R
    echo Result file touched $RETRY
    # Retry if there are NAs (only check score columns, names might (rightfully) contain "NA")
    NAS=$(tr -s '\\t ' ' ' <!{params.collection_name}_SNPQCII_!{chunk}.auto.R | cut -f6,7 -d" " | grep -c NA)
    echo NAs counted $RETRY
    if [ "$NAS" -gt "0" ]; then
       echo Found $NAS NAs
       echo Removing files $RETRY
       rm -f $SOCKET !{params.collection_name}_SNPQCII_!{chunk}.auto.R $PIDFILE || true
    else
       echo Resetting retry counter $RETRY
       RETRY=0
    fi
    echo NAs counted $RETRY
done

'''
}

// SNP_QCII_CON_PS_AS_CD_UC_PSC_parallel_part2.py starts here

process hf_test_merge {
    validExitStatus 0,1

        publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:

    file chunks from hf_test_merge_chunks.collect()
    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()

    output:
    file "${params.collection_name}_SNPQCII_hf" into hf_test_results

    shell:

    dataset = mapFileList(SampleQCI_final_staged)
'''
CHUNKS=$(ls -X !{params.collection_name}_SNPQCII_chunk_*)
OUTFILE="!{params.collection_name}_SNPQCII_hf"
BIMFILE="!{dataset.bim}"

rm -f $OUTFILE

for chunk in $CHUNKS; do
    cat $chunk >>$OUTFILE
done

if [ ! -s "$OUTFILE" ]; then
    touch "!{params.collection_name}_SNPQCII_hf"
    exit 0
fi

# Count number of variants in split/merged and original files
LINES_MERGED=$(wc -l $OUTFILE | cut -d " " -f 1)
LINES_ORIG=$(wc -l $BIMFILE | cut -d " " -f 1)

if [ "$LINES_MERGED" -ne "$LINES_ORIG" ]; then
    echo Unexpected difference in merged HF file.
    echo Expected $LINES_ORIG variants in $OUTFILE
    echo Observed $LINES_MERGED variants in $OUTFILE
    exit 5
fi

# Count number if NAs in split/merged and original files
NAS_MERGED=$(tr -s '\\t ' ' ' <$OUTFILE | cut -d" " -f6,7 | grep -c NA)
NAS_ORIG=$(tr -s '\\t ' ' ' <$BIMFILE | cut -d" " -f1,3,4,5,6 | grep -c NA)

if [ "$NAS_MERGED" -ne "$NAS_ORIG" ]; then
   echo Unexpected difference in 'NA' valiues in HD file.
   echo Expected $NAS_ORIG NAs in $OUTFILE
   echo Observed $NAS_MERGED NAs in $OUTFILE
   exit 5
else
   echo Successful
   exit 0
fi

'''
}

process generate_hf_excludes {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()
    file results from hf_test_results
    file hf_test_excludes_anno

    output:
    file "*.png"
    file 'hf-excludes' into for_exclude_variants

    shell:
        dataset = mapFileList(SampleQCI_final_staged)
    plotscript = SCRIPT_DIR + "/SNP_QCII_draw_FDR_CaseControl.r"
'''
touch hf-excludes
python -c 'from SNPQC_helpers import *; generate_exclude_file_CaseControl("!{results}", "!{hf_test_excludes_anno}", "hf-excludes", "!{params.batches_names}", "!{params.prefix_merged_SNPQCII}", !{params.FDR_index_remove_variants}, "!{plotscript}")'
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
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    file dataset_staged from for_prune.collect() // output from exclude_variants

    output:

//    file "*SNPQCII_final_pruned{.bed,.bim,.fam,.log}" into for_flashpca_pruned, for_draw_histograms_pruned, for_twstats_pruned, for_merge_1kg_pruned
    file "*QCed_pruned{.bed,.bim,.fam,.log}" into for_merge_1kg_pruned
//    file "*SNPQCII_final_pruned_with_atcg{.bed,.bim,.fam,.log}" into for_merge_1kg_pruned_atcg
    file "*.prune.in"
    file "*.prune.out"

    shell:

    dataset = mapFileList(dataset_staged)
    prefix = params.collection_name + "_QCed"
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

<!{dataset.bim} tr -s '\\t ' ' ' | cut -f2 -d' ' | grep ^unk_ >unknowns
plink --bfile !{dataset.bed.baseName} --exclude unknowns --make-bed --out no-unknowns

plink --noweb --bfile no-unknowns --indep-pairwise 50 5 0.2 --out _prune --allow-no-sex
plink --noweb --bfile no-unknowns --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate --allow-no-sex

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("no-unknowns.bim", "include_variants")'
#python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noIndels("no-unknowns.bim", "include_variants_with_atcg")'
plink --noweb --bfile intermediate --extract include_variants --make-bed --out "!{prefix}_pruned" --allow-no-sex
#plink --noweb --bfile intermediate --extract include_variants_with_atcg --make-bed --out "!{prefix}_pruned_with_atcg" --allow-no-sex

mv _prune.prune.in !{dataset.bed.baseName}.prune.in
mv _prune.prune.out !{dataset.bed.baseName}.prune.out

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

process flashpca_1kg {

    input:

    file dataset_staged from for_flashpca_1kg_merged // pruned and merged with 1kG
    file dataset_final from for_flashpca_1kg_final // unpruned


    output:

    file "*.{evec,eval}" into for_pca_plot_1kg_pcaresults, for_pca_plot_1KG_frauke_evec



    shell:

    dataset = mapFileList(dataset_staged)
    annotation = mapFileList(dataset_final).annotation
    prefix = params.collection_name + "_QCed_pruned_1kG"
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
    file "*.sample-excludes" into for_final_sample_cleaning_outliers

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
mv outliers "!{dataset.bim.baseName}.sample-excludes"
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
publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    file dataset_staged from for_det_monomorphics // SNPQCII_final

    output:

    file "${params.collection_name}_SNPQCII.flag.{nearly_,}monomorphic"
    //file "${params.collection_name}_QCed_final.flag.{nearly_,}monomorphic"
    // into for_compile_variants_exclusion_monomorphics

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

#cp !{prefix}.flag.nearly_monomorphic !{params.collection_name}_QCed_final.flag.nearly_monomorphic
#cp !{prefix}.flag.monomorphic !{params.collection_name}_QCed_final.flag.monomorphic
'''
}


process det_diff_missingness {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
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
publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    prefix = "${params.collection_name}_SNPQCII"

    input:

//    file dataset from for_det_unknown_diagnosis

    output:

    file "${prefix}.individuals_remove_final" into for_final_cleaning_individuals
    file "${prefix}.unknown_diagnosis"

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
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true

    input:

    file differential_missingness from for_compile_variants_exclusion_diff_missingness
    //file monomorphics_staged from for_compile_variants_exclusion_monomorphics

    output:

    file "${params.collection_name}.variants_exclude_final" into for_final_cleaning_variants

    shell:
    prefix = "${params.collection_name}_SNPQCII"
    //monomorphics = mapFileList(monomorphics_staged)
'''
    gawk '{ print $1 }' "!{differential_missingness}" | sort | uniq >"!{params.collection_name}.variants_exclude_final"
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

    file "${params.disease_data_set_prefix_release}_final{.bed,.bim,.fam,.log,_annotation.txt,_flag.relatives.txt}" into for_snprelate_prune,for_twstats_final_pruned_ann,for_draw_final_pca_histograms_ds,for_plot_maf,for_eigenstrat_convert_ann,for_snprelate_ann,for_snprelate_ann_atcg,for_twstats_final_pruned_eigenstrat_ann,for_sex_check,for_ibd_verify,for_prepare_imputation,for_final_pca_1kg_frauke_ann,for_det_monomorphics_final

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
    time '2 h'
    input:
    file ds_staged from for_snprelate_prune

    output:
    file "$prefix{.bed,.bim,.fam,.log}" into for_snprelate, for_twstats_final_pruned, for_merge_1kg_pruned_final
    file "${prefix}_with_atcg{.bed,.bim,.fam,.log}" into for_snprelate_atcg

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
    module load "Plink/1.9"

<!{dataset.bim} tr -s '\\t ' ' ' | cut -f2 -d' ' | grep ^unk_ >unknowns
plink --bfile !{dataset.bed.baseName} --exclude unknowns --make-bed --out no-unknowns

plink --bfile no-unknowns --indep-pairwise 50 5 0.2 --out after-indep-pairwise --allow-no-sex
plink --bfile no-unknowns --extract after-indep-pairwise.prune.in --maf 0.05 --make-bed --out after-correlated-remove --allow-no-sex

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("no-unknowns.bim", "include-variants")'
python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noIndels("no-unknowns.bim", "include-variants-with-atcg")'
plink --bfile after-correlated-remove --extract include-variants --make-bed --out "!{prefix}" --allow-no-sex
plink --bfile after-correlated-remove --extract include-variants-with-atcg --make-bed --out "!{prefix}_with_atcg" --allow-no-sex
'''
    }
}

process final_pca_con_projection {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    memory 12.GB

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


process final_pca_con_projection_atcg {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    memory 12.GB

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

process final_merge_pruned_with_1kg {

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
        snpexclude = BATCH_DIR + "/${params.PCA_SNPexcludeList}"
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
    memory "16 G"

    input:
        file dataset_staged from for_pca_plot_1KG_frauke_final
    file ds_original from for_final_pca_1kg_frauke_ann

// need to run 1kg-merge and flashpca manually
//    file pcaresults_staged from for_pca_plot_1KG_frauke_evec_final

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
    --memory !{task.memory.toMega()}
#    --memory 64000

echo Adding batch info | ts

python -c 'from SampleQCI_helpers import *; addphenoinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.pca.evec", "!{prefix}.eval", \
 "!{dataset_orig.annotation}", "!{params.preQCIMDS_1kG_sample}")'


echo Adding country info | ts
python -c 'from SampleQCI_helpers import *; addcountryinfo_10PCs("!{prefix}_pcs_flashpca2", "!{prefix}_eigenvalues_flashpca2", "!{prefix}.country.pca.evec", "!{prefix}.country.eval", "!{dataset_orig.annotation}", "!{params.preQCIMDS_1kG_sample}")'

echo Drawing FLASHPCA2 eigenvectors for set with batch info
Rscript "!{pcaplot_1KG}" "!{prefix}" 5 "!{params.preQCIMDS_1kG_sample}"
mv fail-pca-1KG-qc.txt fail-pca-1KG-qc.txt-batch

echo Drawing FLASHPCA2 eigenvectors for set with country info
Rscript "!{pcaplot_1KG}" "!{prefix}.country" 5 "!{params.preQCIMDS_1kG_sample}"
mv fail-pca-1KG-qc.txt  fail-pca-1KG-qc.txt-country

cat fail-pca-1KG-qc.txt-batch fail-pca-1KG-qc.txt-country | sort | uniq >outliers_final
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

    input:
    file ds_staged from for_sex_check

    output:
    file "*.{pdf,txt}"

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


process ibd_verify {
    publishDir params.qc_dir ?: '.', mode: 'copy'
    cpus 8
    time 8.h
    input:
    file ds_wr from for_ibd_verify

    output:
    file "fail-IBD-percent-qc.txt"
    file "relatedness.png"

    shell:
        ds = mapFileList(ds_wr)
'''
module load IKMB
module load Plink/1.9

Rscript $NXF_DIR/bin/ibdplot.R "!{ds.bed.baseName}" $NXF_DIR/bin/high-LD-regions.txt || true

if [ ! -e relatedness.png ]; then
    echo "No relatives found, no diagrams generated"
    touch relatedness.png
    touch fail-IBD-percent-qc.txt
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


process prepare_sanger_imputation {
    publishDir params.qc_dir ?: '.', mode: 'copy'

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

/* Frauke's version.
process prepare_sanger_imputation {
    publishDir params.qc_dir ?: '.', mode: 'copy'

    input:
    file ds_staged from for_prepare_imputation

    output:
    file "${params.disease_data_set_prefix_release}_final_sanger_{flip.txt,exclude.txt,clean.vcf}"
    shell:
    ds = mapFileList(ds_staged)

'''
module load IKMB
module load Plink/1.9

ANNOTATION=/ifs/data/nfs_share/sukmb388/human_g1k_v37.fasta.gz

awk '($5=="A" && $6=="T") || ($5=="T" && $6=="A")|| ($5=="C" && $6=="G") || ($5=="G" && $6=="C")' "!{ds.bim}" >atcg.bim
sed 's/_/./g' <"!{ds.fam}" >nounderscore.fam

plink --bed "!{ds.bed}" --bim "!{ds.bim}" --fam nounderscore.fam --maf 0.01 --chr 1-22 --exclude atcg.bim --recode-vcf --out "clean_sanger"

bcftools norm -c w  -f $ANNOTATION clean_sanger.vcf &>tmp
<tmp sed 's/REF_MIS/\\nREF_MIS/g' | grep REF_MIS | awk '{print "chr"$2":"$3"\t"$4}' > sanger_force.txt
plink --vcf "clean_sanger.vcf" --reference-allele sanger_force.txt --recode-vcf --out "clean_sanger"

<tmp sed 's/REF_MIS/\\nREF_MIS/g' | grep REF_MIS | awk '{print "chr"$2":"$3}' > sanger_flip.txt

if [ -s sanger_flip.txt ]; then
    comm -12 <(cut -f1 sanger_force.txt | sort) <(comm -3 <(sort sanger_flip.txt) <(cut -f1 sanger_force.txt | sort)) >sanger_force2.txt
    plink --vcf clean_sanger.vcf --reference-allele sanger_force2.txt --flip sanger_flip.txt --recode vcf --out clean_sanger
    <tmp sed 's/REF_MIS/\\nREF_MIS/g' | grep REF_MIS |  awk '{print "chr"$2":"$3}' > sanger_exclude.txt
    if [ -s sanger_exclude.txt ]; then
        plink --vcf clean_sanger.vcf --reference-allele sanger_force2.txt --exclude sanger_exclude.txt --recode vcf --out clean_sanger
    else
        echo "Data prepared. No exclude."
    fi
else
    echo "Data prepared. No flip. No exclude."
fi
echo "Data prepared."
touch sanger_flip.txt
touch sanger_exclude.txt

mv sanger_flip.txt "!{ds.bed.baseName}_sanger_flip.txt"
mv sanger_exclude.txt "!{ds.bed.baseName}_sanger_exclude.txt"
mv clean_sanger.vcf "!{ds.bed.baseName}_sanger_clean.vcf"
'''
}
*/
workflow.onComplete {
    println "Generating phase summary..."
    def cmd = ["./generate-phase-summary", "SNPQCII", params.collection_name ?: params.disease_data_set_prefix, workflow.workDir, params.trace_target].join(' ')
    def gensummary = ["bash", "-c", cmd].execute()
    gensummary.waitFor()
}



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
//SampleQCI_final = ["${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.bed",
//                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.bim",
//                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final.fam",
//                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_flag_relatives.txt",
//                   "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_annotation.txt"].collect { fileExists(file(it)) }

// based on "original" but without relatives
SampleQCI_final_wr = ["${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.bed",
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.bim",
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.fam",
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.annotation.txt",
                      "${params.sampleqci_dir}/${params.collection_name}_SampleQCI_final_withoutRelatives.pca.evec" ].collect { fileExists(file(it)) }
SampleQCI_final = SampleQCI_final_wr

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

    file "*SNPQCII_final{.bed,.bim,.fam,.log,_annotation.txt}" into for_final_cleaning, for_det_monomorphics, for_det_unknown_diagnosis
    file "*test_missingness.*" into for_det_diff_missingness

    shell:

    prefix = params.collection_name + "_SNPQCII_final"
    dataset = mapFileList(SampleQCI_final_staged)
//    println dataset
'''
    module load "IKMB"
    module load "Plink/1.7"
plink --noweb --bfile "!{dataset.bed.baseName}" --exclude "!{exclude}" --make-bed --out "!{prefix}" --allow-no-sex

cp "!{dataset.annotation}" "!{prefix}_annotation.txt"

plink --noweb --bfile "!{dataset.bed.baseName}" --test-missing --out !{prefix}_test_missingness

'''
}



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

    file "${prefix}.variants_exclude_diff_missingness" into for_compile_variants_exclusion_diff_missingness

    shell:
    prefix = "${params.collection_name}_SNPQCII"

    missingness = mapFileList(missingness_staged)
'''
#!/usr/bin/env python

from SNPQC_helpers import *

test_missing = Test_missing(missing_file="!{missingness.missing}", write_file="!{prefix}.variants_exclude_diff_missingness", threshold=!{params.test_missing_thresh})
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

    file ds_staged from for_det_unknown_diagnosis

    output:
    
    file "${prefix}.individuals_remove_final" into for_final_cleaning_individuals
    file "${prefix}.unknown_diagnosis"

    shell:
    dataset = mapFileList(ds_staged)
//    annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"
    annotation = dataset.annotation
    prefix = "${params.collection_name}_SNPQCII"
'''
# TODO: Might need to process "!{params.diagnoses}" into an actual list before calling
python -c 'from SNPQC_helpers import *; determine_unknown_diagnosis(annotationfile="!{annotation}", outfile="!{prefix}.unknown_diagnosis", diagnoses="!{params.diagnoses}")'

if [ -e "!{params.individuals_remove_manually}" ]; then
    cp "!{params.individuals_remove_manually}" "!{prefix}.individuals_remove_manually"
else
    touch "!{prefix}.individuals_remove_manually"
fi

gawk '{ print $1, $2 }' "!{prefix}.unknown_diagnosis" "!{prefix}.individuals_remove_manually" | sort | uniq >"!{prefix}.individuals_remove_final"
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

process final_cleaning {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    input:

    // from diff_missingness/det_monomorphics/det_diagnoses to prune/1kg.
    file dataset_staged from for_final_cleaning
    file individuals from for_final_cleaning_individuals
    file variants from for_final_cleaning_variants

    output:

    file "${params.disease_data_set_prefix_release}{.bed,.bim,.fam,.log,_annotation.txt}" into for_snprelate_prune,for_twstats_final_pruned_ann,for_draw_final_pca_histograms_ds,for_plot_maf,for_eigenstrat_convert_ann,for_snprelate_ann,for_snprelate_ann_atcg,for_twstats_final_pruned_eigenstrat_ann,for_sex_check,for_prepare_imputation,for_final_pca_1kg_frauke_ann,for_det_monomorphics_final
    shell:

    /*annotation = ANNOTATION_DIR + "/${params.individuals_annotation}"*/
    dataset = mapFileList(dataset_staged)
    prefix = params.disease_data_set_prefix_release
'''
    module load "IKMB"
    module load "Plink/1.7"

# Remove sample and SNP outliers
plink --noweb --bfile "!{dataset.bed.baseName}" --remove "!{individuals}" --exclude "!{variants}" --make-bed --out "!{prefix}" --allow-no-sex

touch !{params.collection_name}_SNPQCII_final_flag.relatives.txt

# Fix annotations
python -c 'from SNPQC_helpers import *; extract_QCsamples_annotationfile_relativesfile( \
    fam="!{prefix}.fam", individuals_annotation_QCed="!{prefix}_annotation.txt", \
    related_samples_file="!{params.collection_name}_SNPQCII_final_flag.relatives.txt", \
    related_samples_file_QCed="!{prefix}_flag.relatives.txt", \
    individuals_annotation="!{params.individuals_annotation}", \
    diagnoses="!{params.diagnoses}")'
'''
}

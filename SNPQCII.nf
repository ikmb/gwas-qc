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

process hf_test_prepare {
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
    tag "${params.collection_name}"
    validExitStatus 0,1

    input:
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()

    output:
    file 'chunk_*' into hf_test_chunk
    file 'HF_PCA.R' into hf_test_script
    file 'cleaned-annotation' into hf_test_anno, hf_test_excludes_anno
    file 'cleaned-evec' into hf_test_evec
    file 'diagnoses' into hf_test_diagnoses,hf_generate_excludes_diagnoses

    shell:
    dataset = mapFileList(SampleQCI_final_wr_staged)
'''
    module load "IKMB"
    module load "Plink/1.9"

# Fix annotations and evec file to match FAM file
<!{dataset.fam} tr -s '\\t ' ' ' | cut -f2 -d" " >samples_to_be_used
head -n1 !{dataset.annotation} >cleaned-annotation
head -n1 !{dataset.evec} >cleaned-evec

# NXF sets -ue automatically but we do need +ue for associative arrays in Bash 4.3
set +u
set +e
awk 'NR==FNR{ids[$0];next} {f=($2 in ids)} f' samples_to_be_used !{dataset.annotation} >>cleaned-annotation
awk 'NR==FNR{ids[$0];next} {f=($2 in ids)} f' samples_to_be_used !{dataset.evec} >>cleaned-evec

MIN_BATCH_COUNT=5

declare -A diagnoses

TABLE=$(<cleaned-annotation tail -n+2 | awk '{print $9,$7}' | sort | uniq)

function join_by { local IFS="$1"; shift; echo -n "$*"; }

# Parse table line-by-line
while IFS=  read -r line; do
    # set positional parameters to $line contents
    set -- $line
    DIAG="$1"
    shift

    # make a diagnoses dictionary and count the number of batches for each diagnosis
    for batch; do
        ((diagnoses[$DIAG]++))
    done
done <<<"$TABLE"

LENGTH=${#diagnoses[@]}
echo "Found $LENGTH diagnoses"
ELIGIBLE=()
for diag in "${!diagnoses[@]}";do
    printf "[%s]=%s\\n" "$diag" "${diagnoses[$diag]}"
    if [ "${diagnoses[$diag]}" -ge "$MIN_BATCH_COUNT" ]; then
        ELIGIBLE+=("$diag")
    fi
done

echo "Diagnoses with more than $MIN_BATCH_COUNT batches: "
join_by , ${ELIGIBLE[@]} | tee diagnoses

if [ ! -s diagnoses ]; then
    touch chunk_0
    touch HF_PCA.R
    echo "At least $MIN_BATCH_COUNT batches per diagnosis are required for HF/ANOVA testing."
    exit 0
fi
cp "!{workflow.projectDir}/bin/HF_PCA.R" .
cut -f2 "!{dataset.bim}" | split -l !{hf_test_chunk_size} -a 4 -d - chunk_
'''
}


process hf_test {
    // Keine Ahnung, wo das Script "1" produziert. Muss irgendwo beim NA-zaehlen sein
    validExitStatus 0
    errorStrategy 'retry'
    memory 4.GB
    tag "${params.collection_name}/${chunk}"
    input:
    file SampleQCI_final_wr_staged from Channel.from(SampleQCI_final_wr).collect()
    file r_script from hf_test_script
    file(chunk)from hf_test_chunk.flatten()
    file anno from hf_test_anno
    file evec from hf_test_evec
    file diagnoses from hf_test_diagnoses

    output:
    file "${params.collection_name}_SNPQCII_${chunk}.auto.R" into hf_test_merge_chunks

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

MYPWD=$(pwd)
echo "setwd('$MYPWD')" >hf-local.r
cat !{r_script} >>hf-local.r
plink --bfile "!{dataset.bim.baseName}" --extract "!{chunk}" --allow-no-sex --make-bed --out "!{params.collection_name}_SNPQCII_!{chunk}"
Rscript hf-local.r "!{params.collection_name}_SNPQCII_!{chunk}" 10 "!{anno}" "!{evec}" "$(cat !{diagnoses})" "!{params.collection_name}_SNPQCII_!{chunk}.auto.R"
'''
}

process hf_test_merge {
    validExitStatus 0,1
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    tag "${params.collection_name}"

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
    tag "${params.collection_name}"

    input:
    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()
    file results from hf_test_results
    file annotations from hf_test_excludes_anno
    file diagnoses from hf_generate_excludes_diagnoses

    output:
    file ("*.png") optional
    file 'hf-excludes' into for_exclude_variants

    shell:
    dataset = mapFileList(SampleQCI_final_staged)
    plotscript = SCRIPT_DIR + "/SNP_QCII_draw_FDR_CaseControl.r"
    hf_test_excludes_anno = ""
'''
touch hf-excludes

# make the comma-separated content of 'diagnoses' an array
IFS=, read -r -a diagnoses <<<"$(cat diagnoses)"

MIN_BATCH_COUNT=5
BATCHES_PROCESSED=0

# scan through the array and create an HF-test exclude list for each diagnosis
for i in "${!diagnoses[@]}"
do
    DIAG="${diagnoses[$i]}"

    # prepare annotations and hf_results to contain only one disease
    awk "{if(\\$9==\\"$DIAG\\"||\\$9==\\"diagnosis\\"){print \\$0}}" !{annotations} >"annotations_$DIAG"
    #awk "OFS=\\"\\t\\" { print \\$1,\\$2,\\$3,\\$4,\\$5,\\$6,\\$$((7+$i)),\\$$((7+$i+1)) }" !{results} >"!{results}_$DIAG"
    COL_START=$((7 + $BATCHES_PROCESSED * 2))
    BATCH_COUNT=$(tail -n+2 "annotations_$DIAG" | awk '{print $7}' | sort | uniq | wc -l)
    COL_END=$((7 + $BATCHES_PROCESSED * 2 + $BATCH_COUNT * 2))
    echo "$DIAG: $BATCH_COUNT batches; col start: $COL_START col end: $COL_END; PROCESSED SO FAR: $BATCHES_PROCESSED"
    awk "OFS=\\"\\t\\" { printf(\\"%s\\t%s\\t%s\\t%s\\t%s\\t%s\\",\\$1,\\$2,\\$3,\\$4,\\$5,\\$6); for(x=$COL_START;x<$COL_END;x++) { printf (\\"\\t%s\\", \\$x) }; printf(\\"\\\\n\\");}" !{results} >"!{results}_$DIAG"
    if [ "$BATCH_COUNT" -ge "$MIN_BATCH_COUNT" ]; then
        # only count diseases when there are more than MIN_BATCH_COUNT batches
        continue
    fi
    BATCHES_PROCESSED=$(( $BATCHES_PROCESSED + $BATCH_COUNT))
    # awk "{if(\\$)}"
    python -c "from SNPQC_helpers import *; \\
        generate_exclude_file_for_diagnosis('!{results}_$DIAG',\\
            'annotations_$DIAG', \\
            'hf-excludes-$DIAG', \\
            '_SNPQCII_$DIAG', \\
            !{params.FDR_index_remove_variants}, \\
            '!{plotscript}')"
    cat "hf-excludes-$DIAG" >>hf-excludes-raw
done
touch hf-excludes-raw
sort hf-excludes-raw | uniq >hf-excludes

'''
}



process exclude_variants {
    publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    tag "${params.collection_name}"

    input:
    file SampleQCI_final_staged from Channel.from(SampleQCI_final).collect()
    file exclude from for_exclude_variants

    output:
    file "*SNPQCII_final{.bed,.bim,.fam,.log,_annotation.txt}" into for_final_cleaning, for_det_monomorphics, for_det_unknown_diagnosis
    file "*test_missingness.*" into for_det_diff_missingness

    shell:
    prefix = params.collection_name + "_SNPQCII_final"
    dataset = mapFileList(SampleQCI_final_staged)
'''
    module load "IKMB"
    module load "Plink/1.7"

excludes=$(wc -l <!{exclude})
variants=$(wc -l <!{dataset.bim})

if [ "$variants" -eq "$excludes" ]; then
    echo "All variants have been classified as outliers in HF testing. Check log files and intermediates." | tee "!{prefix}.log" | tee hftest-failed
    ln -s "!{dataset.bed}" "!{prefix}.bed"
    ln -s "!{dataset.bim}" "!{prefix}.bim"
    ln -s "!{dataset.fam}" "!{prefix}.fam"
else
    plink --noweb --bfile "!{dataset.bed.baseName}" --exclude "!{exclude}" --make-bed --out "!{prefix}" --allow-no-sex
fi

cp "!{dataset.annotation}" "!{prefix}_annotation.txt"

plink --noweb --bfile "!{dataset.bed.baseName}" --test-missing --out !{prefix}_test_missingness --allow-no-sex

'''
}



// part 3 begins here


// (1) and (2): determine monomorphic variants and nearly monomorphic variants
process det_monomorphics {
publishDir params.qc_dir ?: '.', mode: 'copy', overwrite: true
    tag "${params.collection_name}"

    input:
    file dataset_staged from for_det_monomorphics // SNPQCII_final

    output:
    file "${params.collection_name}_SNPQCII.flag.{nearly_,}monomorphic"

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
    tag "${params.collection_name}"
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
    tag "${params.collection_name}"
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
    tag "${params.collection_name}"

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
    tag "${params.collection_name}"
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

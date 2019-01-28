// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
*/

def ChipDefinitions = this.class.classLoader.parseClass(new File(params.chip_defs))

// initialize configuration
params.snpqci_dir = "."
hwe_template_script = file(params.hwe_template)

// Lots of indirection layers require lots of backslash escaping
individuals_annotation = file(ANNOTATION_DIR + "/" + params.individuals_annotation)
definetti_r = file(SCRIPT_DIR + "/DeFinetti_hardy.r")
autosomes = file(ANNOTATION_DIR + "/" + params.chip_build + "/" + ChipDefinitions.Producer(params.chip_producer) + "/" + ChipDefinitions.RsAutosomes(params.chip_version))
draw_fdr = file("bin/SNP_QCI_draw_FDR.r")
draw_fdr_allbatches = file("bin/SNP_QCI_draw_FDR_Fail_allbatches.r")

// set up channels
to_calc_hwe_script = Channel.create()
to_calc_hwe = Channel.create()

process merge_batches {

    memory 12.GB

    output:
    file "${params.collection_name}_Rs.bim" into merged_bim, to_split_bim, to_hwe_bim, to_verify_bim, to_exclude_bim, to_miss_bim, to_miss_batch_bim
    file "${params.collection_name}_Rs.bed" into merged_bed, to_hwe_bed, to_verify_bed, to_exclude_bed, to_miss_bed, to_miss_batch_bed
    file "${params.collection_name}_Rs.fam" into merged_fam, to_hwe_fam, to_verify_fam, to_exclude_fam, to_miss_fam, to_miss_batch_fam
    file "${params.collection_name}_Rs.log"
    file "${params.collection_name}.indels" into preqc_hwe_indels,postqc_hwe_indels

    shell:

'''
#!/usr/bin/env bash

module load IKMB
module load Plink/1.9

# Param input: space-separated file names, relative to params.rs_dir
BIMS=(!{params.rs_bims})
BEDS=(!{params.rs_beds})
FAMS=(!{params.rs_fams})
INDELS=(!{params.rs_indels})

IFS=' '

rm -f merge-list
for idx in ${!BIMS[@]}; do
    # All values except the first go into the merge list
    if [ "$idx" -gt 0 ]; then
        echo -n !{params.rs_dir}/${BEDS[$idx]} >>merge-list
        echo -n " " >>merge-list
        echo -n !{params.rs_dir}/${BIMS[$idx]} >>merge-list
        echo -n " " >>merge-list
        echo !{params.rs_dir}/${FAMS[$idx]} >>merge-list
    fi
done

if [ "${#BIMS[*]}" -gt 1 ]; then
    plink --bed !{params.rs_dir}/${BIMS[0]} --bim !{params.rs_dir}/${BIM[0]} --fam !{params.rs_dir}/${FAM[0]} --make-bed --out !{params.collection_name}_Rs --allow-no-sex
else
    echo "No merge necessary, only one batch found." | tee !{params.collection_name}_Rs.log
    ln -s !{params.rs_dir}/${BIMS[0]} !{params.collection_name}_Rs.bim
    ln -s !{params.rs_dir}/${BEDS[0]} !{params.collection_name}_Rs.bed
    ln -s !{params.rs_dir}/${FAMS[0]} !{params.collection_name}_Rs.fam
fi
ln -s !{params.rs_dir}/${INDELS[0]} !{params.collection_name}.indels

'''
}

/*
 Generate HWE tables and draw DeFinetti plots of the whole data set
 */
process hwe_definetti_preqc {
  memory 8192.MB
  cpus 1
  publishDir params.snpqci_dir ?: '.', mode: 'copy', overwrite: true

  input:
    //file plink from to_hwe_diagram
//    file autosomes
    file definetti_r
    file merged_bim
    file merged_bed
    file merged_fam
    file indels from preqc_hwe_indels
  output:
    file "${params.collection_name}_controls_DeFinetti.jpg"
    file "${params.collection_name}_cases_DeFinetti.jpg"
    file "${params.collection_name}_cases_controls_DeFinetti.jpg"
    file "${params.collection_name}_hardy.hwe"
    file "${params.collection_name}_hardy.log"

"""
    module load IKMB
    module load Plink/1.9

plink --bfile ${merged_bim.baseName} --exclude ${indels} --make-bed --out no-indels
plink --bfile no-indels --hardy --out ${params.collection_name}_hardy --hwe 0.0 --chr 1-22 --allow-no-sex
R --slave --args ${params.collection_name}_hardy.hwe ${params.collection_name}_controls_DeFinetti ${params.collection_name}_cases_DeFinetti ${params.collection_name}_cases_controls_DeFinetti <$definetti_r
"""
}

/*
 Split the data set into 1000-SNP chunks so they can be evaluated in parallel.
 */

process split_dataset {
  input:
    file input_bim from to_split_bim

  output:
    file 'chunk_*' into to_calc_hwe

  """
  cut -f 2 $input_bim | split -l 10000 -a 5 -d - chunk_
  """
}


/*
 Run the HWE calculation R script on every 1000-SNP chunk
 */

process calculate_hwe {
  // Should have been zero, but killall -q returns 1 if it didn't find anything
  validExitStatus 0
  errorStrategy 'retry'
  memory 2.GB
  time 1.h

  input:
    each file(chunk) from to_calc_hwe.flatten()
    file input_bim from to_hwe_bim
    file input_bed from to_hwe_bed
    file input_fam from to_hwe_fam
    file individuals_annotation // not directly used in the script below but hwe-script.r expects this to be staged
    file hwe_template_script

  output:
  file "${chunk}-out.auto.R" into from_calc_hwe
//  file "${chunk}-out.nosex"

  tag { chunk }


shell:
'''
#!/usr/bin/env bash

    module load IKMB
    module load Plink/1.9

<!{individuals_annotation} awk ' { if($9 == "Control" || $9 == "diagnosis") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10 }' >individuals_annotation.txt

RETRY=1
while [ $RETRY -eq 1 ]; do
    rm -f /scratch/rserve.sock /scratch/rserve.pid "!{chunk}-out.auto.R"

    # Start Rserve process in background. Keep in mind that --RS-pidfile is an undocumented feature that writes the PID after daemonizing (that is, not the R pid but the Rserve pid) into the specified file. It might be changed in the future
    R CMD Rserve --RS-socket /scratch/rserve.sock --no-save --RS-pidfile /scratch/rserve.pid

    # Wait for R to spawn the Rserve daemon and open its socket. Rule out race conditions with plink not finding the Rserve socket on time
    echo Waiting for Rserve to create the socket
    while [ ! -e /scratch/rserve.sock ]
    do
       sleep 0.5
    done

    # File is there but it does not need to be populated on slow nodes. Wait for that to happen.
    echo Waiting for Rserve to create the pid file
    while [ ! -s /scratch/rserve.pid ]
    do
        sleep 1
    done

    RSERVE_PID=$(cat /scratch/rserve.pid)
    echo Rserve process spawned with pid $RSERVE_PID

    # Cwd to the staging directory where the scripts and chunks are stored
    THEPWD=$(pwd)
    echo "setwd('$THEPWD')" >hwe-script-local.r
    cat !{hwe_template_script} >>hwe-script-local.r
    plink --bfile "!{new File(input_bim.toString()).getBaseName()}" --filter-controls --R-socket /scratch/rserve.sock --R hwe-script-local.r --threads 1 --memory 512 --allow-no-sex --extract !{chunk} --out !{chunk}-out
    sleep 1

    RETRY=0
    kill $RSERVE_PID || {
        echo "Rserve/PLINK failed. Retrying."
        rm -f /scratch/rserve.sock /scratch/rserve.pid "!{chunk}-out.auto.R"
        RETRY=1
    }

    # Did we create a result file? If not, retry.
    if [ ! -e "!{chunk}-out.auto.R" ]; then
        RETRY=1
    else
        # Does our result file contain NAs? If so, retry.
        LINES=$(grep -E '\\sNA\\s?$' "!{chunk}-out.auto.R" | wc -l)
        if [ $LINES -gt 1 ]; then
            RETRY=1
        else
            RETRY=0
        fi
    fi
done
'''
}


/*
 Merge all chunked HWE tables into a single file and screen for obvious errors (i.e. N/A HWE values or wrong SNP counts)
 */

process hwe_fdr_filter {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    input:
    file 'chunk' from from_calc_hwe.collect()
    file input_bim from to_verify_bim
    file input_bed from to_verify_bed
    file input_fam from to_verify_fam

    file individuals_annotation
    file draw_fdr

    output:
    file "${params.collection_name}_exclude-whole-collection-worst-batch-removed" into excludes_whole
    file "${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches" into excludes_perbatch
    file "${params.collection_name}_exclude-whole-collection-all-batches" into excludes_allbatches
    file "${params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt.png"
    file "${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt.png"
    file "${params.collection_name}_exclude-whole-collection-worst-batch-removed.FDRthresholds.SNPQCI.1.txt"
    file "${params.collection_name}_exclude-per-batch-fail-in-2-plus-batches.FDRthresholds.SNPQCI.2.txt"

    shell:
    hwe_result = file("${params.collection_name}_chunks_combined.hwe")
'''
#!/usr/bin/env bash
HWE_RESULTS="!{params.collection_name}_chunks_combined.hwe"
cat chunk[0-9]* >$HWE_RESULTS

combined_hwe_nas=`grep -c NA $HWE_RESULTS`
input_nas=`grep -c NA !{input_bim}`
if [ "$combined_hwe_nas" -ne "$input_nas" ]
then
   echo "There are missing HWE calculations, $combined_hwe_nas NAs in chunked HWE calculations and $input_nas in input."
   exit 1
fi

chunked_snps=`wc -l $HWE_RESULTS | cut -f1 -d" "`
input_snps=`wc -l !{input_bim} | cut -f1 -d" "`
if [ "$chunked_snps" -ne "$input_snps" ]
then
   echo "Some SNPs have gone missing during input splitting. I was expecting $input_snps but got only $chunked_snps."
   exit 1
fi

N_CONTROLS=$(<!{individuals_annotation} tr -s '\\t ' ' ' | cut -f7,9 -d' ' | uniq | tail -n +2 | cut -f2 -d' ' | grep -c Control)

echo Found $N_CONTROLS control batches

if [ "$N_CONTROLS" -gt 2 ]; then
    SNPQCI_fdr_filter.py "$HWE_RESULTS" "!{individuals_annotation}" "!{draw_fdr}" \
       !{params.collection_name}_exclude-whole-collection-worst-batch-removed \
       !{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches \
       !{params.collection_name}_exclude-whole-collection-all-batches \
       !{params.FDR_index_remove_variants}
else
    SNPQCI_fdr_filter.py "$HWE_RESULTS" "!{individuals_annotation}" "!{draw_fdr_allbatches}"\
       !{params.collection_name}_exclude-whole-collection-worst-batch-removed\
       !{params.collection_name}_exclude-per-batch-fail-in-2-plus-batches \
       !{params.collection_name}_exclude-whole-collection-all-batches \
       !{params.FDR_index_remove_variants}
fi

'''
}


/*
 Determine the missingness for the entire collection
 */
process determine_missingness_entire {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    memory { 8.GB * task.attempt }

    input:
    file input_bim from to_miss_bim
    file input_bed from to_miss_bed
    file input_fam from to_miss_fam

    output:
    file 'missingness-excludes-entire' into excludes_miss_entire



"""
module load IKMB
module load Plink/1.7

plink  --bfile "${new File(input_bim.toString()).getBaseName()}" --missing --out missingness_entire --allow-no-sex
SNPQCI_extract_missingness_entire.py missingness_entire.lmiss ${params.geno_entire_collection} missingness-excludes-entire
"""
}


process determine_missingness_per_batch {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    memory { 8.GB * task.attempt }

    input:
    file individuals_annotation
    file input_bim from to_miss_batch_bim
    file input_bed from to_miss_batch_bed
    file input_fam from to_miss_batch_fam

    output:
    file 'missingness-excludes-perbatch' into excludes_miss_perbatch



"""
    module load IKMB
    module load Plink/1.7
awk '{print \$1, \$2, \$7 }' "${individuals_annotation}" | grep -v "familyID" >cluster_file
plink --noweb --bfile "${new File(input_bim.toString()).getBaseName()}" --missing --out missingness_perbatch --missing-phenotype 0 --allow-no-sex --within cluster_file
SNPQCI_extract_missingness_perbatch.py missingness_perbatch.lmiss ${params.geno_batch} "$individuals_annotation" missingness-excludes-perbatch
"""
}

process exclude_bad_variants {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'
    errorStrategy 'retry'
    memory { 8.GB * task.attempt }

    input:
    file individuals_annotation
    file input_bim from to_exclude_bim
    file input_bed from to_exclude_bed
    file input_fam from to_exclude_fam
    file excludes_whole from excludes_whole
    file excludes_perbatch from excludes_perbatch
    file excludes_allbatches from excludes_allbatches
    file missingness_excludes_entire from excludes_miss_entire
    file missingness_excludes_perbatch from excludes_miss_perbatch

    output:
    file "${params.collection_name}_QCI.{bim,bed,fam}" into draw_definetti_after
    file "${params.collection_name}_QCI.log"
    file "variant-excludes"


"""
    module load IKMB
    module load Plink/1.7

NUM_CTRL_BATCHES=\$(tr -s '\\t ' ' ' <${individuals_annotation} | cut -f7,9 -d" " | grep Control | uniq | wc -l)

if [ "\$NUM_CTRL_BATCHES" -gt "4" ]; then
    echo "Found \$NUM_CTRL_BATCHES control batches."
    (tail -n +2 "$excludes_whole" | cut -f1; cat "$excludes_perbatch"; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
else
    echo "Found \$NUM_CTRL_BATCHES control batches. Skipping 'worst batch removed' excludes from HWE testing."
    (tail -n +2 "$excludes_allbatches" | cut -f1; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
fi

plink --noweb --bfile "${new File(input_bim.toString()).getBaseName()}" --exclude variant-excludes --make-bed --out ${params.collection_name}_QCI
"""
}

process hwe_definetti_qci {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'

    def prefix = "${params.collection_name}_QCI_hardy"

    input:
    file autosomes
    file new_plink from draw_definetti_after
    file definetti_r
    file indels from postqc_hwe_indels

    output:
    file prefix+".hwe"
    file prefix+"_{cases,controls,cases_controls}_DeFinetti.jpg"
    file prefix+".log"



"""
    module load IKMB
    module load Plink/1.9
plink --bfile "${new File(new_plink[0].toString()).getBaseName()}" --exclude ${indels} --make-bed --out no-indels
plink --bfile no-indels --hardy --out ${prefix} --hwe 0.0 --chr 1-22 --allow-no-sex
mv no-indels.log ${prefix}.log
R --slave --args ${prefix}.hwe ${prefix}_controls_DeFinetti ${prefix}_cases_DeFinetti ${prefix}_cases_controls_DeFinetti <"$definetti_r"
"""
}

workflow.onComplete {
    println "Generating phase summary..."
    def cmd = ["./generate-phase-summary", "SNPQCI", params.collection_name ?: params.disease_data_set_prefix, workflow.workDir, params.trace_target].join(' ')
    def gensummary = ["bash", "-c", cmd].execute()
    gensummary.waitFor()
}

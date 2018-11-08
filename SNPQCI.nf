// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
*/

def ChipDefinitions = this.class.classLoader.parseClass(new File("config/ChipDefinitions.groovy"))

// initialize configuration
params.rs_dir = "." // Not initialized if we skipped Rs
params.snpqci_dir = "."
//input_basename = params.input
hwe_template_script = file(params.hwe_template)

// Lots of indirection layers require lots of backslash escaping
individuals_annotation = file(ANNOTATION_DIR + "/" + params.individuals_annotation)
definetti_r = file(SCRIPT_DIR + "/DeFinetti_hardy.r")
autosomes = file(ANNOTATION_DIR + "/" + params.chip_build + "/" + ChipDefinitions.Producer(params.chip_producer) + "/" + ChipDefinitions.RsAutosomes(params.chip_version))
draw_fdr = file("bin/SNP_QCI_draw_FDR.r")

// println "Autosomes: " + ANNOTATION_DIR + "/" + params.chip_build + "/" + ChipDefinitions.Producer(params.chip_producer) + "/" + ChipDefinitions.RsAutosomes(params.chip_version)

// set up channels
to_calc_hwe_script = Channel.create()
to_calc_hwe = Channel.create()

process merge_batches {

    memory 12.GB

    output:
    file "${params.collection_name}_Rs.bim" into merged_bim, to_split_bim, to_hwe_bim, to_verify_bim, to_exclude_bim, to_miss_bim, to_miss_batch_bim
    file "${params.collection_name}_Rs.bed" into merged_bed, to_split_bed, to_hwe_bed, to_verify_bed, to_exclude_bed, to_miss_bed, to_miss_batch_bed
    file "${params.collection_name}_Rs.fam" into merged_fam, to_split_fam, to_hwe_fam, to_verify_fam, to_exclude_fam, to_miss_fam, to_miss_batch_fam
    file "${params.collection_name}_Rs.log"
    file "${params.collection_name}.indels" into preqc_hwe_indels,batchqc_hwe_indels,postqc_hwe_indels

    //.getAbsolutePath();
    shell:
//    if(params.skip_rs == 1) {
//    	rsdir = file(".")
//    } else {
//	    rsdir = file(params.rs_dir)
//	}
'''
#!/usr/bin/env bash

IFS=','
NAMES=($(echo "!{params.disease_names}"))
PREFIXES=($(echo "!{params.disease_data_set_prefix_rs}"))

BASEDIR="!{params.rs_dir}"
if [[ ${BASEDIR:0:1} == "/" ]]; then
    echo "$BASEDIR seems to be absolute. No change necessary."
else
    BASEDIR="$NXF_DIR/$BASEDIR"
fi

BASE="$BASEDIR/${PREFIXES[0]}"

module load IKMB
module load Plink/1.9

for idx in ${!PREFIXES[@]}; do
echo idx: $idx
echo name: ${NAMES[$idx]}
echo prefix: ${PREFIXES[$idx]}

cp ${BASEDIR}/!{params.collection_name}.indels !{params.collection_name}.indels

if [ "$idx" -gt 0 ]; then
    echo -n !{params.rs_dir}/${PREFIXES[$idx]}.bed  >>merge-list
    echo -n " " >>merge-list
    echo -n !{params.rs_dir}/${PREFIXES[$idx]}.bim >>merge-list
    echo -n " " >>merge-list
    echo !{params.rs_dir}/${PREFIXES[$idx]}.fam >>merge-list
    cat !{params.rs_dir}/${PREFIXES[$idx]}.indels >>!{params.collection_name}.indels
fi
done

if [ "${#PREFIXES[*]}" -gt 1 ]; then
plink --bfile $BASE --merge-list merge-list --make-bed --out "!{params.collection_name}_Rs" --allow-no-sex --memory 11000
else
cp $BASE.bim !{params.collection_name}_Rs.bim
cp $BASE.bed !{params.collection_name}_Rs.bed
cp $BASE.fam !{params.collection_name}_Rs.fam
echo "No merge perfomed, only one disease name found" >!{params.collection_name}_Rs.log
fi
'''
}

/*
 Generate HWE tables and draw DeFinetti plots of the whole data set
 */
process generate_hwe_diagrams {
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
    file "${params.collection_name}_controls_DeFinetti_preQC.jpg"
    file "${params.collection_name}_cases_DeFinetti_preQC.jpg"
    file "${params.collection_name}_cases_controls_DeFinetti_preQC.jpg"
    file "${params.collection_name}_hardy.hwe"


//    memory '8 GB'
//    cpus 1

//    def basename = new File(merged_bim.toString()).getBaseName()
"""
    module load IKMB
    module load Plink/1.9
    
plink --bfile ${merged_bim.baseName} --exclude ${indels} --make-bed --out no-indels
plink --bfile no-indels --hardy --out ${params.collection_name}_hardy --hwe 0.0 --chr 1-22 --allow-no-sex
R --slave --args ${params.collection_name}_hardy.hwe ${params.collection_name}_controls_DeFinetti_preQC ${params.collection_name}_cases_DeFinetti_preQC ${params.collection_name}_cases_controls_DeFinetti_preQC <$definetti_r
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
 Generate a HWE calculation R script to be used as a Plink slave.
 */

process generate_hwe_script {
  input:
    file individuals_annotation
    file hwe_template_script

  output:
    file "hwe-script.r" into to_calc_hwe_script


  shell:
  '''
  sed "s|INDIVIDUALS_ANNOTATION|!{individuals_annotation}|g" "!{hwe_template_script}" >hwe-script.r
  '''
}

/*
 Run the HWE calculation R script on every 1000-SNP chunk
 */

process calculate_hwe {
  // Should have been zero, but killall -q returns 1 if it didn't find anything
  validExitStatus 0
  errorStrategy 'retry'
    memory { 4.GB * task.attempt }
  input:
    set file(chunk), file('hwe-script.r') from to_calc_hwe.flatten().combine(to_calc_hwe_script)
    file input_bim from to_hwe_bim
    file input_bed from to_hwe_bed
    file input_fam from to_hwe_fam
    file individuals_annotation // not directly used in the script below but hwe-script.r expects this to be staged

  output:
  file "${chunk}-out.auto.R" into from_calc_hwe
//  file "${chunk}-out.nosex"

  tag { chunk }


shell:
'''
#!/usr/bin/env bash

    module load IKMB
    module load Plink/1.9
    
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
	cat hwe-script.r >>hwe-script-local.r
	plink --bfile "!{new File(input_bim.toString()).getBaseName()}" --R-socket /scratch/rserve.sock --R hwe-script-local.r --threads 1 --memory 512 --allow-no-sex --extract !{chunk} --out !{chunk}-out
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

process merge_and_verify_chunked_hwe {
    input:
    file 'chunk' from from_calc_hwe.collect()
    file input_bim from to_verify_bim
    file input_bed from to_verify_bed
    file input_fam from to_verify_fam

    output:
    file "${params.collection_name}_chunks_combined.hwe" into excl_failed_hwe, excl_miss_whole, excl_miss_batch

    shell:
'''
#!/usr/bin/env bash
cat chunk[0-9]* >!{params.collection_name}_chunks_combined.hwe

combined_hwe_nas=`grep -c NA !{params.collection_name}_chunks_combined.hwe`
input_nas=`grep -c NA !{input_bim}`
if [ "$combined_hwe_nas" -ne "$input_nas" ]
then
   echo "There are missing HWE calculations, $combined_hwe_nas NAs in chunked HWE calculations and $input_nas in input."
   exit 1
fi

chunked_snps=`wc -l !{params.collection_name}_chunks_combined.hwe | cut -f1 -d" "`
input_snps=`wc -l !{input_bim} | cut -f1 -d" "`
if [ "$chunked_snps" -ne "$input_snps" ]
then
   echo "Some SNPs have gone missing during input splitting. I was expecting $input_snps but got only $chunked_snps."
   exit 1
fi
'''
}


/*
Make a lists of variants that
 a) fail a HWE test over the entire collection with its worst batch removed
 b) fail at least two tests where the HWE is calculated for each batch

Additionally, FDR values are ploted for both lists
*/

process exclude_lists_for_failed_hwe {

    publishDir params.snpqci_dir ?: '.', mode: 'copy'

    input:
    file hwe_result from excl_failed_hwe
    file individuals_annotation
    file draw_fdr

    output:
    file "${params.collection_name}_exclude-whole-collection" into excludes_whole
    file "${params.collection_name}_exclude-per-batch" into excludes_perbatch
    file "${params.collection_name}_exclude-whole-collection.FDRthresholds.SNPQCI.1.txt.png"
    file "${params.collection_name}_exclude-per-batch.FDRthresholds.SNPQCI.2.txt.png"

//    script:
"""
SNPQCI_fdr_filter.py "$hwe_result" "$individuals_annotation" "$draw_fdr" ${params.collection_name}_exclude-whole-collection ${params.collection_name}_exclude-per-batch ${params.FDR_index_remove_variants}
"""
}

/*
process dummy_hwe_replacement {
output:
    file "${params.collection_name}_exclude-whole-collection" into excludes_whole
    file "${params.collection_name}_exclude-per-batch" into excludes_perbatch
script:
"""
touch "${params.collection_name}_exclude-whole-collection"
touch "${params.collection_name}_exclude-per-batch"
"""
}
*/
/*
 Determine the missingness for the entire collection
 */
process determine_missingness_entire {
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
    file input_bim from to_exclude_bim
    file input_bed from to_exclude_bed
    file input_fam from to_exclude_fam
    file excludes_whole from excludes_whole
    file excludes_perbatch from excludes_perbatch
    file missingness_excludes_entire from excludes_miss_entire
    file missingness_excludes_perbatch from excludes_miss_perbatch

    output:
    file "${params.collection_name}_QCI.{bim,bed,fam}" into draw_definetti_after



"""
    module load IKMB
    module load Plink/1.7
(tail -n +2 "$excludes_whole" | cut -f1; cat "$excludes_perbatch"; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n | uniq >variant-excludes
plink --noweb --bfile "${new File(input_bim.toString()).getBaseName()}" --exclude variant-excludes --make-bed --out ${params.collection_name}_QCI
"""
}

process draw_definetti_after_QCI {
    publishDir params.snpqci_dir ?: '.', mode: 'copy'

    def prefix = "${params.collection_name}"

    input:
    file autosomes
    file new_plink from draw_definetti_after
    file definetti_r
    file indels from postqc_hwe_indels

    output:
    file prefix+".hwe"
    file prefix+"_{cases,controls,cases_controls}_DeFinetti_postQC.jpg"
//    file prefix+"_cases_DeFinetti.jpg"



"""
    module load IKMB
    module load Plink/1.9
plink --bfile "${new File(new_plink[0].toString()).getBaseName()}" --exclude ${indels} --make-bed --out no-indels
plink --bfile no-indels --hardy --out ${prefix} --hwe 0.0 --chr 1-22 --allow-no-sex
R --slave --args ${prefix}.hwe ${prefix}_controls_DeFinetti_postQC ${prefix}_cases_DeFinetti_postQC ${prefix}_cases_controls_DeFinetti_postQC <"$definetti_r"
"""
}

workflow.onComplete {
    println "Generating phase summary..."
    def cmd = ["./generate-phase-summary", "SNPQCI", params.collection_name ?: params.disease_data_set_prefix, workflow.workDir, params.trace_target].join(' ')
    def gensummary = ["bash", "-c", cmd].execute()
    gensummary.waitFor()
}

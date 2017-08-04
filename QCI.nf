// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>

 TODO:



*/


// Should be factored out into some config file
def chip_producer_allowed = ["Illumina" : "Illumina", "Affymetrix" : "Affymetrix"]
def chip_versions_allowed = [
  "Illu300v3" : "unknown",
  "Illu550"   : "unknown",
  "Immunochip"    : "ichip.hg18.hg19.dbsnpID.chr1-26.txt",
  "Exomechipv1"     : "v1_exomearray.hg19.dbsnpID.chr1-26.txt",
  "Exomechipv1-1"   : "v1-1_exomearray.hg19.dbsnpID.chr1-26.txt",
  "HumanCoreExome24v1" : "HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.txt",
  "AgiFranceCustom" : "Agifrance_custom.annot_nofilter.txt",
  "Affy6" : "GenomeWideSNP_6.na24.annot_nofilter.txt",
  "Affy5" : "GenomeWideSNP_5.na24.annot_nofilter.txt",
  "Affy500kSet" : "Mapping250K.na24.annot_nofilter.txt"
]
def chip_rs_autosomes  = [
        "Immunochip" : "ichip.hg19.dbsnpID.chr1-22.rs_id.txt",
        "HumanCoreExome24v1" : "HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-22.rs_id.txt"
        ]


// initialize configuration
params.output = "."
input_basename = params.input
hwe_template_script = file(params.hwe_template)

// Lots of indirection layers require lots of backslash escaping
individuals_annotation = file(BATCH_DIR + "/" + params.individuals_annotation)
definetti_r = file(SCRIPT_DIR + "/DeFinetti_hardy.r")
autosomes = file(ANNOTATION_DIR + "/" + params.chip_build + "/" + chip_producer_allowed[params.chip_producer] + "/" + chip_rs_autosomes[params.chip_version])
draw_fdr = file("bin/SNP_QCI_draw_FDR.r")

// set up channels
to_calc_hwe_script = Channel.create()
to_calc_hwe = Channel.create()

input_bim = file(params.input + ".bim")
input_bed = file(params.input + ".bed")
input_fam = file(params.input + ".fam")

// verify inputs
/* not working well with singularity containers as they are started with the processes, not right here
assert file(hwe_template_script).exists() : "Could not find HWE template script: $hwe_template_script"
assert file(individuals_annotation).exists() : "Could not find individuals annotation file: $individuals_annotation"
assert file(autosomes).exists() : "Could not find autosome annotation file: $autosomes"
assert file(definetti_r).exists() : "Could not find DeFinetti plotting script: $definetti_r"
*/

/*
 Generate HWE tables and draw DeFinetti plots of the whole data set
 */
process generate_hwe_diagrams {
  publishDir params.output ?: '.', mode: 'copy', overwrite: true

  input:
    //file plink from to_hwe_diagram
    file autosomes
    file definetti_r
    file input_bim
    file input_bed
    file input_fam

  output:
    file 'controls_DeFinetti.jpg'
    file 'cases_DeFinetti.jpg'
    file 'cases_controls_DeFinetti.jpg'
    file 'hardy.hwe'

    module 'IKMB'
    module 'Plink/1.9b4.5'
//    memory '8 GB'
//    cpus 1

    def basename = new File(input_bim.toString()).getBaseName()
"""
plink --bfile ${basename} --hardy --out hardy --threads 1 --memory 8192 --hwe 0.0 --extract $autosomes
R --slave --args hardy.hwe controls_DeFinetti cases_DeFinetti cases_controls_DeFinetti <$definetti_r
"""
}

/*
 Split the data set into 1000-SNP chunks so they can be evaluated in parallel.
 */
process split_dataset {
  input:
    file input_bim

  output:
    file 'chunk_*' into to_calc_hwe

  """
  cut -f 2 $input_bim | split -l 1000 -a 3 -d - chunk_
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
  validExitStatus 0,1
  input:
    set file(chunk), file('hwe-script.r') from to_calc_hwe.flatten().combine(to_calc_hwe_script)
    file input_bim
    file input_bed
    file input_fam
    file individuals_annotation // not directly used in the script below but hwe-script.r expects this to be staged

  output:
  file "${chunk}-out.auto.R" into from_calc_hwe
  file "${chunk}-out.nosex"

  tag { chunk }
    module 'IKMB'
    module 'Plink/1.9b4.5'

  def basename = new File(input_bim.toString()).getBaseName()

shell:
'''
#!/usr/bin/env bash

# Start Rserve process in background. Keep in mind that --RS-pidfile is an undocumented feature that writes the PID after daemonizing (that is, not the R pid but the Rserve pid) into the specified file. It might be changed in the future
R CMD Rserve --RS-socket /scratch/rserve.sock --no-save --RS-pidfile /scratch/rserve.pid

# Wait for R to spawn the Rserve daemon and open its socket. Rule out race conditions with plink not finding the Rserve socket on time
while [ ! -e /scratch/rserve.sock ]
do
   sleep 0.5
done
RSERVE_PID=$(cat /scratch/rserve.pid)
echo Rserve process spawned with pid $RSERVE_PID

# Cwd to the staging directory where the scripts and chunks are stored
THEPWD=$(pwd)
echo "setwd('$THEPWD')" >hwe-script-local.r
cat hwe-script.r >>hwe-script-local.r
plink --bfile "!{new File(input_bim.toString()).getBaseName()}" --R-socket /scratch/rserve.sock --R hwe-script-local.r --threads 1 --memory 512 --allow-no-sex --extract !{chunk} --out !{chunk}-out
sleep 1
if [ -e /scratch/rserve.sock ]; then
   echo Killing Rserve with pid $RSERVE_PID
   kill $RSERVE_PID
else
   echo Rserve died on its own
fi
'''
}

/*
 Merge all chunked HWE tables into a single file and screen for obvious errors (i.e. N/A HWE values or wrong SNP counts)
 */
process merge_and_verify_chunked_hwe {
    input:
    file 'chunk' from from_calc_hwe.collect()
    file input_bim
    file input_bed
    file input_fam

    output:
    file 'chunks_combined.hwe' into excl_failed_hwe, excl_miss_whole, excl_miss_batch

    shell:
'''
#!/usr/bin/env bash
cat chunk[0-9]* >chunks_combined.hwe

combined_hwe_nas=`grep -c NA chunks_combined.hwe`
input_nas=`grep -c NA !{input_bim}`
if [ "$combined_hwe_nas" -ne "$input_nas" ]
then
   echo "There are missing HWE calculations, $combined_hwe_nas NAs in chunked HWE calculations and $input_nas in input."
   exit 1
fi

chunked_snps=`wc -l chunks_combined.hwe | cut -f1 -d" "`
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

    publishDir params.output ?: '.', mode: 'move'

    input:
    file hwe_result from excl_failed_hwe
    file individuals_annotation
    file draw_fdr

    output:
    file "exclude-whole-collection" into excludes_whole
    file "exclude-per-batch" into excludes_perbatch
    file "exclude-whole-collection.FDRthresholds.SNPQCI.1.txt.png"
    file "exclude-per-batch.FDRthresholds.SNPQCI.2.txt.png"

//    script:
"""
SNPQCI_fdr_filter.py "$hwe_result" "$individuals_annotation" "$draw_fdr" exclude-whole-collection exclude-per-batch ${params.FDR_index_remove_variants}
"""
}

/*
 Determine the missingness for the entire collection
 */
process determine_missingness_entire {
    input:
    file input_bim
    file input_bed
    file input_fam

    output:
    file 'missingness-excludes-entire' into excludes_miss_entire

    module 'IKMB'
    module 'Plink/1.9b4.5'

"""
plink  --bfile "${new File(input_bim.toString()).getBaseName()}" --missing --out missingness_entire --allow-no-sex
SNPQCI_extract_missingness_entire.py missingness_entire.lmiss ${params.geno_entire_collection} missingness-excludes-entire
"""
}

process determine_missingness_per_batch {
    input:
    file individuals_annotation
    file input_bim
    file input_bed
    file input_fam

    output:
    file 'missingness-excludes-perbatch' into excludes_miss_perbatch

    module 'IKMB'
    module 'Plink/1.9b4.5'

"""
awk '{print \$1, \$2, \$7 }' "${individuals_annotation}" | grep -v "familyID" >cluster_file
plink --bfile "${new File(input_bim.toString()).getBaseName()}" --missing --out missingness_perbatch --allow-no-sex --within cluster_file
SNPQCI_extract_missingness_perbatch.py missingness_perbatch.lmiss ${params.geno_batch} "$individuals_annotation" missingness-excludes-perbatch
"""
}

process exclude_bad_variants {
    publishDir params.output ?: '.', mode: 'copy'

    input:
    file input_bim
    file input_bed
    file input_fam
    file excludes_whole from excludes_whole
    file excludes_perbatch from excludes_perbatch
    file missingness_excludes_entire from excludes_miss_entire
    file missingness_excludes_perbatch from excludes_miss_perbatch

    output:
    file "final.{bim,bed,fam}" into draw_definetti_after

    module 'IKMB'
    module 'Plink/1.9b4.5'

"""
(tail -n +2 "$excludes_whole" | cut -f1; cat "$excludes_perbatch"; cat "$missingness_excludes_entire"; cat "$missingness_excludes_perbatch") | sort -n >variant-excludes
plink --bfile "${new File(input_bim.toString()).getBaseName()}" --exclude variant-excludes --make-bed --out final
"""
}

process draw_definetti_after_QCI {
    publishDir params.output ?: '.', mode: 'move'

    def prefix = 'SNPQCI_entire_collection'

    input:
    file autosomes
    file new_plink from draw_definetti_after
    file definetti_r

    output:
    file prefix+".hwe"
    file prefix+"_{cases,controls,cases_controls}_DeFinetti.jpg"
//    file prefix+"_cases_DeFinetti.jpg"

    module 'IKMB'
    module 'Plink/1.9b4.5'

"""
plink --bfile "${new File(new_plink[0].toString()).getBaseName()}" --hardy --out ${prefix} --hwe 0.0 --extract "$autosomes"
R --slave --args ${prefix}.hwe ${prefix}_controls_DeFinetti ${prefix}_cases_DeFinetti ${prefix}_cases_controls_DeFinetti <"$definetti_r"
"""
}



// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>

 TODO:

 - Turn all absolute paths, i.e. files taken from the PoA hierarchy, into staged
   files. This should reduce friction when running things on the cluster.
   Nextflow (and Slurm, for that matter) will then know what to do.

 - Use module() to select the appropriate Plink version. Reduces configuration
   overhead and maybe eliminates pipeline-specific configuration files
   completely.

 - Find out why an Rserve process starts lingering around, blocking the whole
   pipeline execution on some occasions. In all cases, an OnError handler should
   be installed to kill a leftover Rserve process once the pipeline terminates
   abnormally.

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
  publishDir params.output ?: '.', mode: 'move', overwrite: true

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
    module 'Plink/1.9b4.4'
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
    //file plink from to_split // [bim, bed, fam]
    file input_bim

  output:
    file 'chunk_*' into to_calc_hwe
//    file input_bim into to_calc_hwe

  """
# this is a dummy, no need to split anymore
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
    module 'Plink/1.9b4.4'
    module 'R-3.3.1'

    memory '512 MB'

  def basename = new File(input_bim.toString()).getBaseName()

shell:
'''
#!/usr/bin/env bash

R CMD Rserve --RS-port 12345 --save
THEPWD=$(pwd)
echo "setwd('$THEPWD')" >hwe-script-local.r
cat hwe-script.r >>hwe-script-local.r
plink --bfile "!{new File(input_bim.toString()).getBaseName()}" --R-port 12345 --R hwe-script-local.r --threads 1 --memory 512 --allow-no-sex --extract !{chunk} --out !{chunk}-out
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
    file 'chunks_combined.hwe'

    shell:
'''
#!/usr/bin/env bash
cat chunk* >chunks_combined.hwe

combined_hwe_nas=`grep -c NA chunks_combined.hwe`
input_nas=`grep -c NA !{input_bim}`
if [ $combined_hwe_nas -ne $input_nas ]
then
   echo "There are missing HWE calculations, $combined_hwe_nas NAs in chunked HWE calculations and $input_nas in input."
   exit 1
fi

chunked_snps=`wc -l chunks_combined.hwe`
input_snps=`wc -l !{input_bim}`
if [ $chunked_snps -ne $input_snps ]
then
   echo "Some SNPs have gone missing during input splitting. I was expecting $input_snps but got only $chunked_snps."
   exit 1
fi
'''
}

/*
variant_excludes = Channel.create()

process exclude_lists_for_failed_hwe {
    input:
    file hwe_result

    output:
}

process exclude_lists_for_missingness_whole {
    input:
    file missingness

    output:
    file 'excludes-missingness-whole' into variant_excludes
}

process exclude_lists_for_missingness_per_batch {
    input:
    file missingness // same as in missingness_whole

    output:
    file 'excludes-missingness-per-batch' into variant_excludes
}
*/


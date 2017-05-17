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
evaluate(new File("QCI-Part1.config"))
input_basename = file(params.input)
hwe_template_script = file(params.hwe_template)

// set up channels
input_files = Channel.create()
to_hwe_diagram = Channel.create()
to_split = Channel.create()
to_calc_hwe_script = Channel.create()
to_calc_hwe = Channel.create()

// set up input data set
Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true).separate(to_hwe_diagram, to_split) { a -> [a, a] }

/*
 Generate HWE tables and draw DeFinetti plots of the whole data set
 */
process generate_hwe_diagrams {
  publishDir params.output ?: '.', mode: 'move', overwrite: true

  input:
    file plink from to_hwe_diagram

  output:
    file 'controls_DeFinetti.jpg'
    file 'cases_DeFinetti.jpg'
    file 'cases_controls_DeFinetti.jpg'
    file 'hardy.hwe'

  def autosomes = ANNOTATION_DIR + "/" + params.chip_build + "/" + chip_producer_allowed[params.chip_producer] + "/" + chip_rs_autosomes[params.chip_version]
  def definetti_r = "$SCRIPT_DIR/DeFinetti_hardy.r"

"""
$PLINK --noweb --bfile ${input_basename} --hardy --out hardy --hwe 0.0 --seed 123 --extract $autosomes
R --slave --args hardy.hwe controls_DeFinetti cases_DeFinetti cases_controls_DeFinetti <$definetti_r
"""
}

/*
 Split the data set into 1000-SNP chunks so they can be evaluated in parallel.
 */
process split_dataset {
  input:
  file plink from to_split

  output:
  file 'chunk_*' into to_calc_hwe

  def par = input_basename + ".bim"

  """
  cut -f 2 $par | split -l 1000 -a 3 -d - chunk_
  """
}

/*
 Generate a HWE calculation R script to be used as a Plink slave.
 */
process generate_hwe_script {
  output:
  file "hwe-script.r" into to_calc_hwe_script

  // Lots of indirection layers require lots of backslash escaping
  def annotation_file = (ANNOTATION_DIR + "/" + params.chip_build + "/" + chip_producer_allowed[params.chip_producer] + "/" + params.individuals_annotation).replaceAll('/', '\\\\\\\\/')

  script:
  """
  R CMD Rserve --save
  sed s/INDIVIDUALS_ANNOTATION/${annotation_file}/g ${hwe_template_script} >hwe-script.r
  """
}

/*
 Run the HWE calculation R script on every 1000-SNP chunk
 */
process calculate_hwe {
  input:
  set file(chunk), file('hwe-script.r') from to_calc_hwe.flatten().combine(to_calc_hwe_script)

  output:
  file "${chunk}-out.auto.R" into from_calc_hwe
  file "${chunk}-out.nosex"

"""
$PLINK --noweb --bfile ${input_basename} --R hwe-script.r --missing-phenotype 0 --allow-no-sex --extract ${chunk} --out ${chunk}-out
"""
}

/*
 Merge all chunked HWE tables into a single file and screen for obvious errors (i.e. N/A HWE values or wrong SNP counts)
 */
process merge_and_verify_chunked_hwe {
    input:
    file 'chunk' from from_calc_hwe.collect()

    output:
    file 'chunks_combined.hwe'

    shell:
'''
#!/usr/bin/env bash
cat chunk* >chunks_combined.hwe
combined_hwe_nas=`grep -c NA chunks_combined.hwe`
input_nas=`grep -c NA !{input_basename}.bim`
if [ $combined_hwe_nas -ne $input_nas ]
then
   echo "There are missing HWE calculations, $combined_hwe_nas NAs in chunked HWE calculations and $input_nas in input."
   exit 1
fi

chunked_snps=`wc -l chunks_combined.hwe`
input_snps=`wc -l !{input_basename}.bim`
if [ $chunked_snps -ne $input_snps ]
then
   echo "Some SNPs have gone missing during input splitting. I was expecting $input_snps but got only $chunked_snps."
   exit 1
fi
'''
}


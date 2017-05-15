#!/usr/bin/env nextflow
// -*- mode:groovy -*-

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


params.output = "."
evaluate(new File("QCI-Part1.config"))

input_basename = file(params.input)

hwe_template_script = file(params.hwe_template)

input_files = Channel.create()
to_hwe_diagram = Channel.create()
to_split = Channel.create()
to_calc_hwe_script = Channel.create()

Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true).separate(to_hwe_diagram, to_split) { a -> [a, a] }

process generate_hwe_diagrams {
  input:
  file plink from to_hwe_diagram

  def autosomes = ANNOTATION_DIR + "/" + params.chip_build + "/" + chip_producer_allowed[params.chip_producer] + "/" + chip_rs_autosomes[params.chip_version]
  def definetti_r = "$SCRIPT_DIR/DeFinetti_hardy.r"

"""
$PLINK --noweb --bfile ${params.input} --hardy --out hardy --hwe 0.0 --extract $autosomes
R --slave --args hardy.hwe ${params.input}_controls_DeFinetti ${params.input}_cases_DeFinetti ${params.input}_cases_controls_DeFinetti <$definetti_r
"""
}

to_calc_hwe = Channel.create()
split_flatten_glob = Channel.create()
split_flatten_glob.flatten().combine(to_calc_hwe_script).into(to_calc_hwe)

process split_dataset {
  input:
  file plink from to_split

  output:
  file 'chunk_*' into split_flatten_glob

  def par = params.input + ".bim"

"""
#R CMD Rserve --save
cut -f 2 $par | split -l 1000 -a 3 -d - chunk_
"""
}



process generate_hwe_script {
  output:
  file "hwe-script.r" into to_calc_hwe_script

  def annotation_file = (ANNOTATION_DIR + "/" + params.chip_build + "/" + chip_producer_allowed[params.chip_producer] + "/" + params.individuals_annotation).replaceAll('/', '\\\\\\\\/')

  script:
  """
  R CMD Rserve --save
  sed s/INDIVIDUALS_ANNOTATION/${annotation_file}/g ${hwe_template_script} >hwe-script.r
  """
}

process calculate_hwe {
  input:
  set file('chunk'), file('hwe-script.r') from to_calc_hwe

  output:
  file "chunk-out.auto.R" into from_calc_hwe

  publishDir params.output ?: '.', mode: 'move', overwrite: true

"""
$PLINK --noweb --bfile ${input_basename} --R hwe-script.r --missing-phenotype 0 --allow-no-sex --extract chunk --out chunk-out
"""
}


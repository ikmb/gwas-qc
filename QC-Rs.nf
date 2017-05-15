#!/usr/bin/env nextflow
// -*- mode:groovy -*-

params.output = "."

evaluate(new File("QC-Rs.nf"))


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
def chip_strand_info_allowed = [
  "Immunochip_orig_annotation":"ichip.orig_annotation.hg18.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.hg19.txt",
  "Immunochip_TOP_annotation":"ichip.TOP_annotation.hg18.hg19.dbsnpID.chr1-26.FlipToPlusStrandOnly.rs.hg19.txt",
  "Exomechipv1_orig_annotation":"v1_exomearray.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt",
  "Exomechipv1-1_orig_annotation":"v1-1_exomearray.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt",
  "HumanCoreExome24v1_orig_annotation":"HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.minusStrandOnly.rs.txt"
]
def chip_rs_exclude = [
    "Immunochip":"ichip.hg18.hg19.dbsnpID.chr1-26.exclude.HailiangHuang.chr25.chr26.txt", 
    "HumanCoreExome24v1":"HumanCoreExome24_v1.0.hg19.dbsnpID.chr1-26.duplicates.txt"
]


input_files = Channel.create()
to_flipfile = Channel.create()

Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true).separate(input_files, to_flipfile) { a -> [a, a[1]] }

process generate_annotations {
  def input_files = Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true)

  input:
  file plink from input_files

  output:
  file 'annotations.list' into annotations, to_translate_ann

  def annotation_file = file(ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_versions_allowed.get(params.chip_version)).toAbsolutePath()

"""
perl -ne '@l=split(/\\s+/);print "\$l[3] \$l[4] \$l[7] \$l[8] \$l[5] \$l[1] \$l[2]\\n";' $annotation_file >annotations.list
"""
}

// 
process generate_flipfile {
  input:
  file ann from annotations
  file bim from to_flipfile

  output:
  file 'flipfile' into to_plink_flip

  def source = file(ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_strand_info_allowed.get(params.chip_strand_info)).toAbsolutePath()

"""
if [ -e $source ]; then
  cp $source flipfile
else
  bin/generate_flipfile.pl $bim $ann >flipfile
fi
"""
}

process plink_flip {
  input:
  file plink from input_files
  file flip  from to_plink_flip

  output:
  file 'flipped.{bed,fam}' into to_plink_exclude_plink
  file 'flipped.bim' into to_translate_bim

"""
$PLINK --bed ${plink[1]} --bim ${plink[2]} --fam ${plink[3]} --flip $flip --make-bed --out flipped
"""
}

process translate_ids {
  input:
  file bim from to_translate_bim
  file ann from to_translate_ann

  output:
  file 'translated.bim' into to_find_duplicates, to_find_nn, to_exclude_bim

"""
translate_ichip_to_rs.pl $bim $ann ${params.chip_build} ${params.switch_to_chip_build} >translated.bim
"""
}

process find_duplicates {
  input:
  file bim from to_find_duplicates

  output:
  file 'duplicates' into to_merge_exclude_duplicates

"""
cut -f 2 $bim | sort | uniq -d >duplicates
"""
}

process find_nn {
  input:
  file bim from to_find_nn

  output:
  file 'nn' into to_merge_exclude_nn

"""
grep -P "\\tN\\tN" $bim | cut -f2 >nn
"""
}

process merge_exclude_list {
  input:
  file duplicates from to_merge_exclude_duplicates
  file nn from to_merge_exclude_nn

  output:
  file 'exclude' into to_plink_exclude_list

  def source = ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_rs_exclude.get(params.chip_version)

"""
if [ -e $source ]; then
  cat $duplicates $nn $source >exclude
else
  cat $duplicates $nn >exclude
fi
"""
}

process plink_exclude {

  publishDir params.output ?: '.', mode: 'move', overwrite: true

  input:
  file exclude from to_plink_exclude_list
  file plink from to_plink_exclude_plink
  file bim from to_exclude_bim

  output:
  file 'result.{bim,bed,fam}' into final_results
"""
$PLINK --bed ${plink[0]} --bim $bim --fam ${plink[1]} --exclude $exclude --make-bed --out result
"""
}

// -*- mode:groovy -*-

/*
 Author: Jan Kässens <j.kaessens@ikmb.uni-kiel.de>
*/

// Set default output directory
params.output = "."

// evaluate(new File("config/QC-Rs.config"))

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

// Set up channels between processes

input_files = Channel.create()
to_flipfile = Channel.create()

// Move input files (basename given via --input PLINK) into first-stage channel
Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true).separate(input_files, to_flipfile) { a -> [a, a[1]] }

/*
 Transform a chip-specific annotations file into a format that is easier to process by the following stages.
 */
process generate_annotations {

//    memory "128 MB"
//    cpus 1

    def input_files = Channel.fromFilePairs(params.input + "{.bim,.bed,.fam}", size:3, flat: true)

    input:
    file plink from input_files

    output:
    file 'annotations.list' into annotations, to_translate_ann


    def annotation_file = file(params.annotation_dir+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_versions_allowed.get(params.chip_version)).toAbsolutePath()

"""
perl -ne '@l=split(/\\s+/);print "\$l[3] \$l[4] \$l[7] \$l[8] \$l[5] \$l[1] \$l[2]\\n";' $annotation_file >annotations.list
"""
}

/*
 Generate a Plink flipfile based on annotations and strand info
 */
process generate_flipfile {
    input:
    file ann from annotations
    file bim from to_flipfile

    output:
    file 'flipfile' into to_plink_flip

    module 'perl5.22.0'
 //   memory '128 MB'
//    cpus 1

    def source = file(params.annotation_dir+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_strand_info_allowed.get(params.chip_strand_info)).toAbsolutePath()

"""
if [ -e $source ]; then
  cp $source flipfile
else
  generate_flipfile.pl $bim $ann >flipfile
fi
"""
}

/*
 Call Plink to actually flip alleles based on strand information
 */
process plink_flip {
    input:
    file plink from input_files
    file flip  from to_plink_flip

    output:
    file 'flipped.{bed,fam}' into to_plink_exclude_plink
    file 'flipped.bim' into to_translate_bim

    module 'IKMB'
    module 'Plink/1.9b4.4'
//    cpus 1
//    memory '6 GB'

"""
plink --bed ${plink[1]} --bim ${plink[2]} --fam ${plink[3]} --flip $flip --threads 1 --memory 6144 --make-bed --out flipped
"""
}

/*
 Translate Immunochop IDs to Rs names
 */
process translate_ids {
    input:
    file bim from to_translate_bim
    file ann from to_translate_ann

    output:
    file 'translated.bim' into to_find_duplicates, to_find_nn, to_exclude_bim

    module 'perl5.22.0'
//    cpus 1
//    memory '512 MB'

    println "Switching from ${params.chip_build} to ${params.switch_to_chip_build}"

"""
translate_ichip_to_rs.pl $bim $ann ${params.chip_build} ${params.switch_to_chip_build} >translated.bim
"""
}

/*
 Create a list of duplicate SNPs now that all SNPs have standardized Rs names
 */
process find_duplicates {
    input:
    file bim from to_find_duplicates

    output:
    file 'duplicates' into to_merge_exclude_duplicates

//    cpus 1
//    memory '128 MB'

"""
cut -f 2 $bim | sort | uniq -d >duplicates
"""
}

/*
 Create a list of NN SNPs
 */
process find_nn {
    input:
    file bim from to_find_nn

    output:
    file 'nn' into to_merge_exclude_nn

//    cpus 1
//    memory '128 MB'

"""
grep -P "\\tN\\tN" $bim | cut -f2 >nn
"""
}

/*
 Merge duplicates, NNs and, if available, a chip-specific exclude file into a single file
 */
process merge_exclude_list {
    input:
    file duplicates from to_merge_exclude_duplicates
    file nn from to_merge_exclude_nn

    output:
    file 'exclude' into to_plink_exclude_list

//    cpus 1
//    memory '128 MB'

    def source = params.annotation_dir+'/'+params.switch_to_chip_build+'/'+chip_producer_allowed.get(params.chip_producer)+'/'+chip_rs_exclude.get(params.chip_version)
    println "Using chip exclude list $source"

"""
if [ -e $source ]; then
  cat $duplicates $nn $source >exclude
else
  cat $duplicates $nn >exclude
fi
"""
}

/*
 Apply an exclude list to the translated Plink data set
 */
process plink_exclude {

    publishDir params.output ?: '.', mode: 'move', overwrite: true

    input:
    file exclude from to_plink_exclude_list
    file plink from to_plink_exclude_plink
    file bim from to_exclude_bim

    output:
    file 'result.{bim,bed,fam}'

    module 'IKMB:Plink/1.9b4.4'
//    cpus 1
//    memory '6 GB'
"""
plink --bed ${plink[0]} --bim $bim --fam ${plink[1]} --exclude $exclude --threads 1 --memory 6144 --allow-no-sex --make-bed --out result
"""
}


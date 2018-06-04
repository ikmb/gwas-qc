// -*- mode:groovy -*-

def ChipDefinitions = this.class.classLoader.parseClass(new File("config/ChipDefinitions.groovy"))

// Set default output directory, overwritten by --output=<dir>
params.output = "."

// Set up channels between processe
input_files_flip = Channel.create()
input_files_ann = Channel.create()
to_flipfile = Channel.create()

// Prepare input file pairs from batches
Channel.fromFilePairs(params.input + ".{bim,bed,fam}", size:3, flat: true).separate(input_files_flip, input_files_ann, to_flipfile) { a -> [a, a, a] }

/*
 Transform a chip-specific annotations file into a format that is easier to process by the following stages.
 */
process generate_annotations {
    input:
    // Process results without inputs will not get cached, so force caching of results by providing a predictable dummy parameter
    val dummy from Channel.from(1);

    output:
    file 'annotations.list' into to_flipfile_ann, to_translate_ann


    def annotation_file = ANNOTATION_DIR + "/" + params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.SNPAnnotations(params.chip_version)

    tag { params.disease_data_set_prefix_orig }
"""
echo -n Generating annotations.list from ${annotation_file}
if [ "${params.chip_version}" == "Immunochip" ]; then
    echo " for Immunochip"
    perl -ne '@l=split(/\\s+/);print "\$l[3] \$l[4] \$l[7] \$l[8] \$l[5] \$l[1] \$l[2]\\n";' $annotation_file >annotations.list
else
    echo " for ${params.chip_version}"
    perl -ne '@l=split(/\\s+/);print "\$l[2] \$l[3] \$l[6] \$l[7] \$l[4] \$l[1] \$l[1]\\n";' $annotation_file >annotations.list
fi
"""
}

/*
 Generate a Plink flipfile based on annotations and strand info
 */
process generate_flipfile {
//    echo true

    input:
    file ann from to_flipfile_ann
    file ds from to_flipfile

    output:
    file 'flipfile' into to_plink_flip

    def source = file(ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.StrandInfo(params.chip_strand_info)).toAbsolutePath()
    tag { params.disease_data_set_prefix_orig }
"""
if [ -e $source ]; then
  echo Using ${source} as flipfile
  cp "$source" flipfile
else
  echo Generating flipfile from ${ann}
  generate_flipfile.pl "${ds[0].baseName}.bim" "$ann" >flipfile
fi
"""
}

/*
 Call Plink to actually flip alleles based on strand information
 */
process plink_flip {
//    echo true

    input:
    file pl from input_files_flip
    file flip  from to_plink_flip

    output:
    file "${pl[1].baseName}_flipped.{bed,fam}" into to_plink_exclude_plink
    file "${pl[1].baseName}_flipped.bim" into to_translate_bim

    tag { params.disease_data_set_prefix_orig }

shell:
'''
module load IKMB
module load Plink/1.9
echo Flipping strands for "!{pl}"
plink --bed "!{pl[1]}" --bim "!{pl[2]}" --fam "!{pl[3]}" --flip "!{flip}" --threads 1 --memory 6144 --make-bed --out "!{pl[1].baseName}_flipped" --allow-no-sex
'''
}
//
/*
 Translate Immunochop IDs to Rs names
 */
process translate_ids {
    input:
    file bim from to_translate_bim
    file ann from to_translate_ann

    output:
    file "${bim.baseName}_translated.bim" into to_find_duplicates_nn, to_exclude_bim

    tag { params.disease_data_set_prefix_orig }

"""
echo Translating SNP names for ${params.chip_version}
translate_ichip_to_rs.pl ${params.chip_version} "$bim" "$ann" ${params.chip_build} ${params.switch_to_chip_build} >"${bim.baseName}_translated.bim"
"""
}

/*
 Create a list of duplicate SNPs now that all SNPs have standardized Rs names
 */
process find_duplicates_nn {
    input:
    file bim from to_find_duplicates_nn

    output:
    file 'exclude' into to_plink_exclude_list

    tag { params.disease_data_set_prefix_orig }

    shell:
    source = ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.RsExclude(params.chip_version)
'''
#!/bin/bash
cut -f 2 "!{bim}" | sort | uniq -d >duplicates
echo Found $(wc -l < duplicates) duplicate variants

grep -P "\\tN\\tN" "!{bim}" | cut -f2 >nns
echo Found $(wc -l < nns) N/N variants

if [ -e "!{source}" ]; then
  echo Found $(wc -l < "!{source}") variants in chip exclude list "!{source}"
  cat duplicates nns "!{source}" >exclude
else
  echo Chip exclude list "${source}" not accessible - skipping
  cat duplicates nns >exclude
fi

'''
}


/*
 Apply an exclude list to the translated Plink data set
 */
process plink_exclude {
    publishDir params.output ?: '.', mode: 'copy'

    input:
    file exclude from to_plink_exclude_list
    file plink from to_plink_exclude_plink
    file bim from to_exclude_bim

    output:
    file "${params.target_name}.{bim,bed,fam,log}"

    tag { params.disease_data_set_prefix_orig }
    // Note that Plink 1.07 only excludes the first of duplicates, while 1.9+ removes all duplicates


"""
module load 'IKMB'
module load 'Plink/1.9'
echo Excluding SNP list ${exclude} from ${plink} ${bim}
plink  --bed ${plink[0]} --bim $bim --fam ${plink[1]} --exclude $exclude --make-bed --out ${params.target_name} --allow-no-sex
"""
}

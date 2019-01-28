// -*- mode:groovy -*-

// If you need to add a new chip definition, take a look into the ChipDefinitions.groovy file
def ChipDefinitions = this.class.classLoader.parseClass(new File("config/ChipDefinitions.groovy"))

// Set default output directory, overwritten by --output=<dir>
params.output = "."

// Set up channels between processels
input_files_flip = Channel.create()
input_files_ann = Channel.create()
to_flipfile = Channel.create()

input_files_lift = Channel.create()

// Prepare input file pairs from batches
//Channel.fromFilePairs(params.input + ".{bim,bed,fam}", size:3, flat: true).separate(input_files_flip, input_files_ann, to_flipfile) { a -> [a, a, a] }

//println BATCH_DIR + "/" + params.collection_name + "/" + params.basename + ".{bim,bed,fam}"

// 

params.collection_name = false

/*
 Transform a chip-specific annotations file into a format that is easier to process by the following stages.
 */
// process generate_annotations {
//     input:
//     // Process results without inputs will not get cached, so force caching of results by providing a predictable dummy parameter
//     val dummy from Channel.from(1);

//     output:
//     file 'annotations.list' into to_flipfile_ann, to_translate_ann

//     def annotation_file = ANNOTATION_DIR + "/" + params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.SNPAnnotations(params.chip_version)

//     tag { params.disease_data_set_prefix }
// """
// echo -n Generating annotations.list from ${annotation_file}
// if [ "${params.chip_version}" == "Immunochip" ]; then
//     echo " for Immunochip"
//     perl -ne '@l=split(/\\s+/);print "\$l[3] \$l[4] \$l[7] \$l[8] \$l[5] \$l[1] \$l[2]\\n";' $annotation_file >annotations.list
// else
//     echo " for ${params.chip_version}"
//     perl -ne '@l=split(/\\s+/);print "\$l[2] \$l[3] \$l[6] \$l[7] \$l[4] \$l[1] \$l[1]\\n";' $annotation_file >annotations.list
// fi
// """
// }

// BEWARE: the first item in this list seems to be a "1". I have no clue why.

// input_files_check =  Channel.fromFilePairs(BATCH_DIR + "/" + params.batch_name + "/" + params.basename + ".{bim,bed,fam}", size:3, flat: true)
//input_files_check =  Channel.fromFilePairs(BATCH_DIR + "/" + params.basename + ".{bim,bed,fam}", size:3, flat: true)

Channel.fromFilePairs(BATCH_DIR + "/" + params.batch_name + "/" + params.basename + ".{bim,bed,fam}", size:3, flat: true).into {input_files_check; input_files_lift}

println "Input files: " + BATCH_DIR + "/" + params.batch_name + "/" + params.basename + ".{bim,bed,fam}"

process check_chip_type {
    memory 4.GB
    cpus 2
    publishDir params.rs_dir ?: '.', mode: 'copy'

    input:

    file original from input_files_check

    output:

//    file "*.strand"
//    file "${original[1].baseName}.{bed,bim,fam}" into input_files_lift
    file "${original[1].baseName}.flag_atcg"
    file "${original[1].baseName}.chip_detect.log"

    shell:

'''
chipmatch --verbose --output !{original[1].baseName}.chip_detect.log --threads 2 !{original[1].baseName}.bim /work_beegfs/sukmb388/wayne_strands
<!{original[1].baseName}.bim awk '{if(($5=="A" && $6=="T")||($5=="T" && $6=="A")) { printf("%s %s%s\\n", $2, $5, $6); }}' >!{original[1].baseName}.flag_atcg
<!{original[1].baseName}.bim awk '{if(($5=="C" && $6=="G")||($5=="G" && $6=="C")) { printf("%s %s%s\\n", $2, $5, $6); }}' >>!{original[1].baseName}.flag_atcg
'''
}

process lift_genome_build {
    memory 8.GB
    input:

    file original from input_files_lift
 
    output:
    file "${original[1].baseName}_lift.{bed,bim,fam}" into to_normalize_variants, to_plink_flip_bedfam

    shell:

'''
module load Plink/1.9
TARGETNAME="!{original[1].baseName}_lift"
BASENAME="!{original[1].baseName}"
STRAND_FILE="!{params.lift_to}"

if [ -e "$STRAND_FILE" ]; then
    $NXF_DIR/bin/update_build_PLINK1.9.sh "$BASENAME" "$STRAND_FILE" "$TARGETNAME"
else
    echo "No strand file specified for lifting."
    mv "!{original[1].baseName}.bed" "$TARGETNAME.bed"
    mv "!{original[1].baseName}.bim" "$TARGETNAME.bim"
    mv "!{original[1].baseName}.fam" "$TARGETNAME.fam"
fi
'''
}

process normalize_variant_names {
	time 24.h
    publishDir params.rs_dir ?: '.', mode: 'copy'

    input:
    file source from to_normalize_variants

    output:
    file 'flipfile' into to_plink_flip
    file "${params.collection_name}.indels"
    file "${source[0].baseName}_updated.bim" into to_plink_flip_bim
    file "Rs-${params.collection_name}.RsTranslation.log"

    shell:
'''
#!/usr/bin/env perl

use strict;
use warnings;

use Switch;
use DBI;
use DBD::SQLite::Constants qw/:file_open/;
use File::Copy;
use Data::Dumper;

my $logtarget = "Rs-!{params.collection_name}.RsTranslation.log";

my $scratch_dir = $ENV{'TMPDIR'};

# print "Copying database to $scratch_dir}...\\n";

my $variant_db = "!{params.variant_annotation_db}";
my $db_in_scratch = 0;

if(copy("!{params.variant_annotation_db}", "$scratch_dir/annotation.sqlite") == 0) {
	print "Could not copy annotation DB, using remote DB instead. This will take a long time.\\n";
	# Maybe copy() failed half-way through, so better unlink
	unlink("$scratch_dir/annotation.sqlite");
} else {
	$variant_db = "$scratch_dir/annotation.sqlite";
	print "Copied annotation database to local scratch.\\n";
	$db_in_scratch = 1;
}

my $db = DBI->connect("dbi:SQLite:dbname=$variant_db", "", "", { sqlite_open_flags => SQLITE_OPEN_READONLY });
my $query_plus = $db->prepare('SELECT name FROM annotation WHERE chrom=?1 AND position=?2 AND ((strand="+" AND alleleA=?3 AND alleleB=?4) OR (strand="+" AND alleleB=?3 AND alleleA=?4) OR (strand="-" AND alleleA=?5 AND alleleB=?6) OR (strand="-" AND alleleB=?5 AND alleleA=?6))');

my $query_minus = $db->prepare('SELECT name FROM annotation WHERE chrom=?1 AND position=?2 AND ((strand="+" AND alleleA=?5 AND alleleB=?6) OR (strand="+" AND alleleB=?5 AND alleleA=?6) OR (strand="-" AND alleleA=?3 AND alleleB=?4) OR (strand="-" AND alleleB=?3 AND alleleA=?4))');

# Checks each variant against the variant annotation file, replaces variant names and generates a flip file
my $bim = "!{source[0].baseName}.bim";
my $new_bim = "!{source[0].baseName}_updated.bim";

my $num_match = 0;
my $num_replaced = 0;
my $num_noname = 0;
my $num_flip_match = 0;
my $num_flip_replaced = 0;
my $num_flip_noname = 0;
my $num_mismatch = 0;
my $num_indels = 0;

sub make_complement {
    my $orig = shift;
    switch($orig) {
        case 'A' { return 'T' }
        case 'C' { return 'G' }
        case 'G' { return 'C' }
        case 'T' { return 'A' }
        else { print STDERR "Setting invalid allele $orig to X"; return 'X' }
    }
}

sub find_on_strand {
    my $dbh = shift;
    my $strand = shift;
    my $pos = shift;
    my $alla = shift;
    my $allb = shift;

    my $res = $db->execute();

}

my $indels_file = "!{params.collection_name}.indels";

open my $log, '>', $logtarget or die("Could not open $logtarget: $!");
open my $flip, '>', "flipfile" or die("Could not open flipfile: $!");
open my $newbim, '>', $new_bim or die("Could not open $new_bim: $!");
open my $fh, '<', $bim or die("Could not open $bim: $!");
open my $indels, '>', $indels_file or die("Could not open $indels_file: $!");

print $log "OldID\tNewID\tReason\n";

while(<$fh>) {
   chomp;
   my ($chr, $name, $pos_cm, $pos, $alla, $allb) = split(/\\s+/, $_);

   if($. % 1000 == 0) { print $. . "\\n"; }


   my $alla_c = make_complement($alla);
   my $allb_c = make_complement($allb);

   my $result = $query_plus->execute($chr, $pos, $alla, $allb, $alla_c, $allb_c) or die($!);
   my $rows = $query_plus->fetchrow_arrayref();
   if(!defined $rows) {
       $result = $query_minus->execute($chr, $pos, $alla, $allb, $alla_c, $allb_c) or die($!);
       $rows = $query_minus->fetchrow_arrayref();

       if(!defined $rows) { $num_mismatch++; print $log "$name\t"; $name = "unk_$name"; print $log "$name\tMismatch\n"; }
       elsif(!defined $rows->[0]) { $num_flip_noname++; print $log "$name\t"; $name = "$chr:$pos"; print $log "$name\tNo_name_in_DB\n"; print $flip "$name\\n"; }
       elsif($rows->[0] eq $name) { $num_flip_match++; print $flip "$name\\n"; }
       else { $num_flip_replaced++; print $log "$name\t"; $name = $rows->[0]; print $log "$name\tName_replaced\n"; print $flip "$name\\n"; }
   } else {
       if(!defined $rows->[0]) { $num_noname++; print $log "$name\t"; $name = "$chr:$pos"; print $log "$name\tNo_name_in_DB\n";}
       elsif($rows->[0] eq $name) { $num_match++;}
       else { $num_replaced++; print $log "$name\t"; $name = $rows->[0]; print $log "$name\tName_replaced\n"; }
   }

   if( (length($alla) > 1) || (length($allb) > 1) ) {
      print $indels "$name\\n";
      $num_indels++;
   }

   print $newbim "$chr\\t$name\\t$pos_cm\\t$pos\\t$alla\\t$allb\\n";
}

print "Exact matches: $num_match\\n";
print "Variants with new Rs ID: $num_replaced\\n";
print "Variants without known Rs IDs: $num_noname\\n";
print "Flipped matches: $num_flip_match\\n";
print "Flipped variants with new Rs ID: $num_flip_replaced\\n";
print "Flipped variants without known Rs IDs: $num_flip_noname\\n";
print "Unknown variants: $num_mismatch\\n";
print "...of which are indels: $num_indels\\n";

# copy("$scratch_dir/bim", $new_bim) or die($!);
# copy("$scratch_dir/fliplist", "flipfile") or die($!);

unlink("$scratch_dir/annotation.sqlite");
#unlink("$scratch_dir/bim");
#unlink("$scratch_dir/fliplist");
'''
}

// /*
//  Generate a Plink flipfile based on annotations and strand info
//  */
// process generate_flipfile {
// //    echo true

//     input:
//     file ann from to_flipfile_ann
//     file ds from to_flipfile

//     output:
//     file 'flipfile' into to_plink_flip

//     def source = file(ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.StrandInfo(params.chip_strand_info)).toAbsolutePath()
//     tag { params.disease_data_set_prefix }
// """
// if [ -e $source ]; then
//   echo Using ${source} as flipfile
//   cp "$source" flipfile
// else
//   echo Generating flipfile from ${ann}
//   generate_flipfile.pl "${ds[0].baseName}.bim" "$ann" >flipfile
// fi
// """
// }

/*
 Call Plink to actually flip alleles based on strand information
 */
process plink_flip {
//    echo true

    input:
    file bim from to_plink_flip_bim
    file bedfam from to_plink_flip_bedfam
    file flip  from to_plink_flip

    output:
    file "${bedfam[0].baseName}_flipped.{bed,fam}" into to_plink_exclude_plink
    file "${bedfam[0].baseName}_flipped.bim" into to_exclude_bim, to_find_duplicates_nn
    file "duplicates"

    tag { params.disease_data_set_prefix }

shell:
'''
module load IKMB
module load Plink/1.9

echo Deduplicating "!{bim.baseName}"
cut -f2 "!{bim}" | sort -T . | uniq -d >duplicates
echo Found $(wc -l <duplicates) duplicates.

plink --bed "!{bedfam[0].baseName}.bed" --bim "!{bim}" --fam "!{bedfam[0].baseName}.fam" --exclude duplicates --make-bed --out dedup

echo Flipping strands for "!{bim.baseName}"
plink --bfile dedup --flip "!{flip}" --threads 1 --memory 6144 --make-bed --out "!{bedfam[0].baseName}_flipped" --allow-no-sex
'''
}
//
// /*
//  Translate Immunochop IDs to Rs names
//  */
// process translate_ids {
//     input:
//     file bim from to_translate_bim
//     file ann from to_translate_ann

//     output:
//     file "${bim.baseName}_translated.bim" into to_find_duplicates_nn, to_exclude_bim

//     tag { params.disease_data_set_prefix }

// """
// echo Translating SNP names for ${params.chip_version}
// translate_ichip_to_rs.pl ${params.chip_version} "$bim" "$ann" ${params.chip_build} ${params.switch_to_chip_build} >"${bim.baseName}_translated.bim"
// """
// }

/*
 Create a list of duplicate SNPs now that all SNPs have standardized Rs names
 */
process find_duplicates_nn {
    input:
    file bim from to_find_duplicates_nn

    output:
    file 'exclude' into to_plink_exclude_list

    tag { params.disease_data_set_prefix }

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

individuals_annotation = file(ANNOTATION_DIR + "/" + params.individuals_annotation)

/*
 Apply an exclude list to the translated Plink data set
 */
process plink_exclude {
    publishDir params.rs_dir ?: '.', mode: 'copy'

    input:
//    file individuals_annotation
    file exclude from to_plink_exclude_list
    file plink from to_plink_exclude_plink
    file bim from to_exclude_bim

    output:
    file "${params.disease_data_set_prefix_rs}.{bim,bed,fam,log}"
//    file individuals_annotation

    tag { params.disease_data_set_prefix }
    // Note that Plink 1.07 only excludes the first of duplicates, while 1.9+ removes all duplicates


"""
module load 'IKMB'
module load 'Plink/1.9'
echo Excluding SNP list ${exclude} from ${plink} ${bim}
plink  --bed ${plink[0]} --bim $bim --fam ${plink[1]} --exclude $exclude --make-bed --out ${params.disease_data_set_prefix_rs} --allow-no-sex
"""
}

workflow.onComplete {
    println "Generating phase summary..."
    def cmd = ["./generate-phase-summary", "Rs", params.collection_name ?: params.disease_data_set_prefix, workflow.workDir, params.trace_target].join(' ')
    def gensummary = ["bash", "-c", cmd].execute()
    gensummary.waitFor()
}


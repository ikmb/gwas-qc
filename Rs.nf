// -*- mode:groovy -*-

// If you need to add a new chip definition, take a look into the ChipDefinitions.groovy file
//def ChipDefinitions = this.class.classLoader.parseClass(new File(params.chip_defs))

// Set default output directory, overwritten by --output=<dir>
params.output = "."

// Set up channels between processels
input_files_flip = Channel.create()
input_files_ann = Channel.create()
to_flipfile = Channel.create()

input_files_lift = Channel.create()

Channel.fromFilePairs(params.filebase + ".{bim,bed,fam}", size:3, flat: true).into {input_files_check; input_files_lift}

annotations = file(params.filebase + "_individuals_annotation.txt")

println "Input files: " + params.filebase + ".{bim,bed,fam}"

process batch_statistics {
    tag "${params.ds_name}/${params.batch_name}"
    input:
    file annotations

    output:
    file "${params.ds_name}-${params.batch_name}.stats"
shell:
'''
#!/usr/bin/env perl

use strict;
use warnings;

my $infilename = "!{annotations}";
my $outfilename = "!{params.ds_name}-!{params.batch_name}.stats";
##################################################
open my $fh, '<', $infilename or die($!);
<$fh>; # skip header

my %diagnoses;
my %genders;
my %batches;
my %countries;
my %phenotypes;

while(<$fh>) {
    chomp;
    my @f = split '\\s+';
    $diagnoses{$f[8]}++;
    $genders{$f[4]}++;
    $batches{$f[6]}++;
    $countries{$f[9]}++;
    $phenotypes{$f[5]}++;
}

close $fh;

open my $ofh, '>', $outfilename or die($!);
print $ofh "Genders "; while(my ($key,$val)=each(%genders)) { print $ofh "$key $val "; } print $ofh "\n";
print $ofh "Phenotypes "; while(my ($key,$val)=each(%phenotypes)) { print $ofh "$key $val "; } print $ofh "\n";
print $ofh "Batches "; while(my ($key,$val)=each(%batches)) { print $ofh "$key $val "; } print $ofh "\n";
print $ofh "Diagnoses "; while(my ($key,$val)=each(%diagnoses)) { print $ofh "$key $val "; } print $ofh "\n";
print $ofh "Countries "; while(my ($key,$val)=each(%countries)) { print $ofh "$key $val "; } print $ofh "\n";
close $ofh;
'''
}

process check_chip_type {
    tag "${params.ds_name}/${params.batch_name}"
    cpus 2
    publishDir params.rs_dir ?: '.', mode: 'copy'

    input:

    file original from input_files_check

    output:
    file "${original[1].baseName}.flag_atcg"
    file "${original[1].baseName}.chip_detect.log"
    file "${params.batch_name}.sourcedata.txt"
    shell:
    basename= "${original[1].baseName}"
'''

echo "!{basename}.bim;$(stat -L -c %y !{basename}.bim)" >>sourcedata.txt
echo "!{basename}.bed;$(stat -L -c %y !{basename}.bed)" >>sourcedata.txt
echo "!{basename}.fam;$(stat -L -c %y !{basename}.fam)" >>sourcedata.txt
cp sourcedata.txt !{params.batch_name}.sourcedata.txt

plinkinfo.pl !{basename}.bim !{basename}.fam >info.txt

chipmatch --verbose --output !{original[1].baseName}.chip_detect.log --threads 2 !{original[1].baseName}.bim !{params.wrayner_strands}
<!{original[1].baseName}.bim mawk '{if(($5=="A" && $6=="T")||($5=="T" && $6=="A")) { printf("%s %s%s\\n", $2, $5, $6); }}' >!{original[1].baseName}.flag_atcg
<!{original[1].baseName}.bim mawk '{if(($5=="C" && $6=="G")||($5=="G" && $6=="C")) { printf("%s %s%s\\n", $2, $5, $6); }}' >>!{original[1].baseName}.flag_atcg
'''
}

process lift_genome_build {
    label 'big_mem'
    tag "${params.ds_name}/${params.batch_name}"
    input:

    file original from input_files_lift
 
    output:
    file "${original[1].baseName}_lift.{bed,bim,fam}" into to_fix_par

    shell:

'''
TARGETNAME="!{original[1].baseName}_lift"
BASENAME="!{original[1].baseName}"
STRAND_ZIP="/assets/annotations/wrayner_strands/Source/!{params.liftover}-b37.Source.strand.zip"
STRAND_FILE="!{params.liftover}-b37.Source.strand"

module load Plink/1.9

if [ "!{params.liftover}" == "" ]; then
    echo "No strand file specified for lifting."
    ln -s "!{original[1].baseName}.bed" "$TARGETNAME.bed"
    ln -s "!{original[1].baseName}.bim" "$TARGETNAME.bim"
    ln -s "!{original[1].baseName}.fam" "$TARGETNAME.fam"
    exit 0
fi

if [ ! -f "$STRAND_ZIP" ]; then
    echo "Could not find $STRAND_ZIP. Did you specify the correct strand file?"
    exit 1
fi

unzip "$STRAND_ZIP"

if [ ! -f "$STRAND_FILE" ]; then
    echo "Could not find $STRAND_FILE within $STRAND_ZIP."
    exit 1
fi


# TODO: Add memory requirements, maybe integrate script
!{SCRIPT_DIR}/update_build_PLINK1.9.sh  "!{original[1].baseName}" "$STRAND_FILE" "$TARGETNAME"
'''
}

process fix_par {
    label 'long_running'
    tag "${params.ds_name}/${params.batch_name}"

    input:
    file source from to_fix_par

    output:
    file "*fixed.{bim,bed,fam}" into to_normalize_variants, to_plink_flip_bedfam

    shell:
    '''
module load Plink/1.9

MEM=!{task.memory.toMega()-1000}
plink --memory $MEM --bfile "!{source[1].baseName}" --merge-x no-fail --make-bed --out merged-with-missing
plink --memory $MEM --bfile merged-with-missing --split-x hg19 no-fail --make-bed --out split
plink --memory $MEM --bfile split --set-hh-missing --make-bed --out fixed

# remove intermediates, keep logs for parsing
rm -f merged-with-missing.{bed,bim,fam} merged-no-hh.{bed,bim,fam}

    '''
}

process normalize_variant_names {
	label 'long_running'
    tag "${params.ds_name}/${params.batch_name}"
    publishDir params.rs_dir ?: '.', mode: 'copy'

    input:
    file source from to_normalize_variants

    output:
    file 'flipfile' into to_plink_flip
    file "${params.batch_name}.indels"
    file "${source[0].baseName}_updated.bim" into to_plink_flip_bim
    file "Rs-${params.batch_name}.RsTranslation.log"

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

my $logtarget = "Rs-!{params.batch_name}.RsTranslation.log";

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
my $num_atcg = 0;

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

my $indels_file = "!{params.batch_name}.indels";

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

/*
 Call Plink to actually flip alleles based on strand information
 */
process plink_flip {
    label 'big_mem'

    tag "${params.ds_name}/${params.batch_name}"
    input:
    file bim from to_plink_flip_bim
    file bedfam from to_plink_flip_bedfam
    file flip  from to_plink_flip

    output:
    file "${bedfam[0].baseName}_flipped.{bed,fam}" into to_plink_exclude_plink
    file "${bedfam[0].baseName}_flipped.bim" into to_exclude_bim, to_find_duplicates_nn
    file "duplicates"

shell:
'''
module load IKMB
module load Plink/1.9

MEM=!{task.memory.toMega()-1000}

echo Deduplicating "!{bim.baseName}"
cut -f2 "!{bim}" | sort -T . | uniq -d >duplicates
echo Found $(wc -l <duplicates) duplicates.

plink --bed "!{bedfam[0].baseName}.bed" --bim "!{bim}" --memory $MEM --fam "!{bedfam[0].baseName}.fam" --exclude duplicates --make-bed --out dedup

echo Flipping strands for "!{bim.baseName}"
plink --bfile dedup --flip "!{flip}" --threads 1 --memory $MEM --make-bed --out "!{bedfam[0].baseName}_flipped" --allow-no-sex
'''
}

/*
 Create a list of duplicate SNPs now that all SNPs have standardized Rs names
 */
process find_duplicates_nn {
    label 'short_running'
    label 'small_mem'
    tag "${params.ds_name}/${params.batch_name}"
    input:
    file bim from to_find_duplicates_nn

    output:
    file 'exclude' into to_plink_exclude_list

    shell:
//    source = ANNOTATION_DIR+'/'+params.switch_to_chip_build+'/'+ChipDefinitions.Producer(params.chip_producer)+'/'+ChipDefinitions.RsExclude(params.chip_version)
'''
#!/bin/bash
cut -f 2 "!{bim}" | sort | uniq -d >duplicates
echo Found $(wc -l < duplicates) duplicate variants

grep -P "\\tN\\tN" "!{bim}" | cut -f2 >nns
echo Found $(wc -l < nns) N/N variants

  cat duplicates nns >exclude

'''
}

/*
 Apply an exclude list to the translated Plink data set
 */
process plink_exclude {
    publishDir params.rs_dir ?: '.', mode: 'copy'
    tag "${params.ds_name}/${params.batch_name}"

    input:
    file exclude from to_plink_exclude_list
    file plink from to_plink_exclude_plink
    file bim from to_exclude_bim

    output:
    file "${params.batch_name}_Rs.{bim,bed,fam,log}"
    // Note that Plink 1.07 only excludes the first of duplicates, while 1.9+ removes all duplicates

"""
module load 'IKMB'
module load 'Plink/1.9'
echo Excluding SNP list ${exclude} from ${plink} ${bim}
plink  --bed ${plink[0]} --bim $bim --fam ${plink[1]} --exclude $exclude --make-bed --out ${params.batch_name}_Rs --allow-no-sex


plinkinfo.pl ${params.batch_name}_Rs.bim ${params.batch_name}_Rs.fam >info.txt
"""
}


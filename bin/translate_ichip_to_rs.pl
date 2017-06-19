#!/usr/bin/env perl

use warnings;
use strict;

use 5.020;

use Data::Dumper;

die("Syntax: $0 <file.bim> <file.annotations> <from> <to>") unless @ARGV == 4;

my $bimfile = $ARGV[0];
my $annotfile = $ARGV[1];
my $chip_build = $ARGV[2];
my $switch_to = $ARGV[3];

my %annotations;

# parse annotations
# ichip_id -> (rs_id, alleleA, alleleB, strand, pos_hg18, pos_hg19)

open(my $annot_fh, "<$annotfile") or die("Could not open $annotfile: $!");
while(<$annot_fh>) {
    chomp;
    my @l = split(' ');
    $annotations{$l[0]} = { 'rs' => $l[1], 'a1' => $l[2], 'a2' => $l[3], 'strand' => $l[4], 'hg18' => $l[5], 'hg19' => $l[6] };
}
close $annot_fh;


open(my $bim_fh, "<$bimfile") or die("Could not open $bimfile: $!");
while(<$bim_fh>) {
    chomp;
    my @l = split(/\s+/);
    my $ichip = $l[1];

    if(defined $annotations{$ichip}) {
        my %ann = %{$annotations{$ichip}};

        if ($ann{hg19} eq 'NA') {
            $ann{hg19} = '-9999999';
        }

        if ($chip_build eq 'hg18' and $switch_to eq 'hg19') {
            say "$l[0]\t$ann{rs}\t$l[2]\t$ann{hg19}\t$l[4]\t$l[5]";
        } else {
            say "$l[0]\t$ann{rs}\t$l[2]\t$l[3]\t$l[4]\t$l[5]";
        }
    } else {
        say STDERR "Could not find SNP $ichip in annotation file!";
    }

}
close $bim_fh;


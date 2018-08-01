package Report::SNPQCII;

use strict;
use warnings;

use Report::SampleQC;
use Data::Dumper;
use File::Slurp::Tiny qw/ read_file /;

our $VERSION = '1.00';


sub exclude_variants {
    my $workdir = shift;
    my $tag = shift;

    my $s = '\section{SNP QC}';
    $s = '\subsection{Huynh-Feldt Correction for Repeated Measure ANOVA}';

    my $excludes = `wc -l $workdir/hf-excludes`;
    chomp($excludes);
    $s .= "$excludes SNPs were found to violate the sphericity assumption and are removed from the dataset. ";


    return $s;
}

sub prune {
    my $workdir = shift;
    my $tag = shift;

    return Report::SampleQC::prune($workdir, $tag, 'BZZZZT')
}

sub flashpca_pruned {
    my $workdir = shift;
    my $tag = shift;

    return Report::SampleQC::onekg_flashpca($workdir, $tag)
}

1

package NXFQC::PlinkLog;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(parse_plink);


our $VERSION = "1.00";

sub parse_19 {
    my $fh = shift;
    my $l = {
        'variants-same-position' => 0
    };

    while(<$fh>) {
        if(/--indep-pairwise (\S+) (\S+) (\S+)/) {
            $l->{'indep-pairwise-windowsize'} = $1;
            $l->{'indep-pairwise-stepsize'} = $2;
            $l->{'indep-pairwise-r2'} = $3;
            next;
        }

        if(/--maf (\S+)/) {
            $l->{'maf'} = $1;
            next;
        }

        if(/^(\d+) variants loaded from \.bim file/) {
            $l->{'loaded-variants'} = $1;
            next;
        }

        if(/^(\d+) phenotype values loaded from .fam/) {
            $l->{'loaded-phenotypes'} = $1;
            next;
        }

        if(/^\d+ people \((\d+) males, (\d+) females, (\d+) ambiguous\) loaded from \.fam/) {
            $l->{'loaded-males'} = $1;
            $l->{'loaded-females'} = $2;
            $l->{'loaded-ambig'} = $3;
            next;
        }

        if(/^Total genotyping rate is (\S+)\.$/) {
            $l->{'loaded-gt-rate'} = $1;
            next;
        }

        if(/^Warning: Variants.*and.* have the same position/) {
            $l->{'variants-same-position'} += 1;
            next;
        }

        if(/^Total genotyping rate in remaining samples is (\S+)\./) {
            $l->{'final-gt-rate'} = $1;
            next;
        }

        if(/^--exclude: (\d+) variants remaining/) {
            $l->{'variants-after-exclude'} = $1;
            next;
        }

        if(/^--extract: (\d+) variants remaining/) {
            $l->{'variants-after-extract'} = $1;
            next;
        }

        if(/^--remove: (\d+) people remaining/) {
            $l->{'samples-after-remove'} = $1;
            next;
        }

        if(/^--flip: (\d+) SNPs flipped/) {
            $l->{'snps-flipped'} = $1;
            next;
        }

        if(/^Pruning complete.\s+(\d+) of \d+ variants removed/) {
            $l->{'pruned-variants'} = $1;
            next;
        }

        if(/^(\d+) variants removed due to minor allele threshold/) {
            $l->{'maf-variants-removed'} = $1;
            next;
        }

        if(/^Among remaining phenotypes, (\d+) are cases and (\d+) are controls/) {
            $l->{'final-cases'} = $1;
            $l->{'final-controls'} = $2;
            next;
        }
    }
    return $l;
}

sub parse_102_107 {
    my $fh = shift;
    my $l = {};

    while(<$fh>) {
        $l->{'in-variants'} = $1 if /^(\d+) markers to be included from/;
        $l->{'in-samples'} = $1 if /^(\d+) individuals read from/;
        $l->{'final-variants'} = $1 if /^After frequency and genotyping pruning, there are (\d+) SNPs/;
        if(/^(\d+) cases, (\d+) controls and (\d+) missing/) {
            $l->{'in-cases'} = $1;
            $l->{'in-controls'} = $2;
            $l->{'in-missing'} = $3;
        }
        if(/^After filtering, (\d+) cases, (\d+) controls and (\d+) missing/) {
            $l->{'out-cases'} = $1;
            $l->{'out-controls'} = $2;
            $l->{'out-missing'} = $3;
        }
        $l->{'removed'} = $1 if /^(\d+) individuals removed with --remove option/;
    }
    $l->{'out-samples'} = $l->{'out-cases'} + $l->{'out-controls'} + $l->{'out-missing'};

    return $l;
}

sub parse_plink {
    my $log = shift;

    open my $fh, '<', $log or die("$!: $log");

    my $line = <$fh>;
    my $data = {};

    if ($line =~ /PLINK v(.*)$/) {
        $data = parse_19($fh);
        $data->{'version'} = $1;
    } else {
        # first line is empty
        # second line is only graphics
        # third line contains actual information
        <$fh>;
        $line = <$fh>;
        if ($line =~ /PLINK!.*v(\S+)\s+\S\s+(\S+)\s+/) {
            $data = parse_102_107($fh);
            $data->{'version'} = "$1 $2";
        } else {
            print STDERR "Unknown PLINK logfile format: $log\n";
        }
    }

    return $data;
}

1

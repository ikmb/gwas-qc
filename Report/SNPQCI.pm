package Report::SNPQCI;

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use File::Slurp::Tiny qw/ read_file /;

our $VERSION = '1.00';

sub merge_batches {
    my $workdir = shift;
    my $tag = shift;

    my @batches = ();

    my $cases = 0;
    my $controls = 0;
    my $variants = 0;
    my $samples = 0;

    open my $log, '<', $workdir . '/.command.log' or die($!);
    while(<$log>) {
        chomp;
        if (/^name: (.*)/) {
            push @batches, $1;
            next;
        }

        if (/^(\d+) variants and \d+ people pass filters and QC/) {
            $variants = $1;
            next;
        }

        if (/^Among remaining phenotypes, (\d+) are cases and (\d+) are controls/) {
            $cases = $1;
            $controls = $2;
            $samples = $cases + $controls;
        }
    }

    my $s = "\\section{Merged Batches}";

    map { s/\_/\\\_/ } @batches;
    $s .= "After merging the batches " . join(", ", map { "\\texttt{$_}" } @batches) . ", we have $variants variants and $samples samples, of which $cases are cases and $controls are controls.";
    return $s;
}

sub generate_hwe_diagrams {
    my $workdir = shift;
    my $tag = shift;

    my $minipage;

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        chomp;
        if (/hardy.hwe (\S+) (\S+) (\S+) /) {

            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $1 . '.jpg}\hfill';
            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $2 . '.jpg}\\\\';
            $minipage .= '\hfill\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $3 . '.jpg}\hfill\\\\';
            last;
        }
    }

    my $s = "\n" . '\subsection{DeFinetti Diagrams}' . "\n";
    return $s . $minipage;

}

sub exclude_lists_for_failed_hwe {
    my $workdir = shift;
    my $tag = shift;
    my $name;
    open my $fh, '<', $workdir . '/.command.sh' or die($!);
    while(<$fh>) {
        if (/(\S+)_exclude-whole-collection/) {
            $name = $1;
        }
    }
    close $fh;


    # Read whole-collection table
    my $filename = $workdir . '/' . $name . '_exclude-whole-collection.FDRthresholds.SNPQCI.1.txt';
    open my $fh1, '<', $filename or die("Couldn't open $filename for reading: $!");
    <$fh1>; # Skip first line
    my %snpcounts;
    while(<$fh1>) {
        chomp;
        print STDERR "Read line: $_\n";
        if (/^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)$/) {
            print STDERR "Match!\n";
            $snpcounts{$1} = { 'all_fail' => $2, 'all_hwe' => $3, 'oneless_fail' => $4, 'oneless_hwe' => $5 };
        }
    }
    close $fh1;

    # Read per-batch table
    open my $fh2, '<', $workdir . '/' . $name . '_exclude-per-batch.FDRthresholds.SNPQCI.2.txt' or die($!);
    <$fh2>;
    while(<$fh2>) {
        chomp;
        if (/^(\S+)\s(\S+)\s(\S+)$/) {
            $snpcounts{$1}->{'1plus_fail'} = $2;
            $snpcounts{$1}->{'2plus_fail'} = $3;
        }
    }
    close $fh2;

    my $s = <<'HWE_TABLE'
\subsection{HWE Tests}'
    \begin{tabular}{r|rr|rr|r|r}
%    \toprule
    & \multicolumn{2}{r|}{\emph{All Batches}} & \multicolumn{2}{r|}{\emph{Worst Batch Removed}} & & \\
    FDR & \# Failed & HWE p-Val & \# Failed & HWE p-Val & Failed in 1+ batch & Failed in 2+ batches\\
    \hline

HWE_TABLE
        ;

    foreach my $key (sort { $b <=> $a} keys %snpcounts) {
        my $v = $snpcounts{$key};
        $s .= $key . "&" . $v->{'all_fail'} . "&" . $v->{'all_hwe'} . "&" . $v->{"oneless_fail"} . "&" . $v->{"oneless_hwe"} . "&" . $v->{"1plus_fail"} . "&" . $v->{"2plus_fail"} . '\\\\' . "\n";
    }

    $s .= '\end{tabular}\\\\';

    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.FDRthresholds.SNPQCI.1.txt.png,read=.FDRthresholds.SNPQCI.1.txt.png]{' . "$workdir/$name" . '_exclude-whole-collection}';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.FDRthresholds.SNPQCI.2.txt.png,read=.FDRthresholds.SNPQCI.2.txt.png]{' . "$workdir/$name" . '_exclude-per-batch}';


    return $s;
}

sub exclude_bad_variants {
    my $workdir = shift;
    my $tag = shift;

    my $name;
    open my $fh, '<', $workdir . '/.command.sh' or die($!);
    while(<$fh>) {
        if (/(\S+)_QCI/) {
            $name = $1;
        }
    }
    close $fh;

    my $miss_entire = `wc -l < $workdir/missingness-excludes-entire`;
    chomp($miss_entire);
    my $miss_batch = `wc -l < $workdir/missingness-excludes-perbatch`;
    chomp($miss_batch);
    my $hwe_entire = `wc -l < $workdir/${name}_exclude-whole-collection`;
    chomp($hwe_entire);
    my $hwe_batch = `wc -l < $workdir/${name}_exclude-per-batch`;
    chomp($hwe_batch);
    my $all = `wc -l < $workdir/variant-excludes`;
    chomp($all);

    my $miss_entire_thres;
    my $miss_perbatch_thres;

    (undef, my $miss_entire_dir, undef) = fileparse(readlink("$workdir/missingness-excludes-entire"));
    open my $miss_entire_fh, '<', "$miss_entire_dir/.command.sh" or die("$!");
    while(<$miss_entire_fh>) {
        if (/missingness_entire.lmiss (\d+\.\d+) missingness/) {
            $miss_entire_thres = $1;
        }
    }
    close $miss_entire_fh;

    (undef, my $miss_batch_dir, undef) = fileparse(readlink("$workdir/missingness-excludes-perbatch"));
    open my $miss_batch_fh, '<', "$miss_batch_dir/.command.sh" or die("$!");
    while(<$miss_batch_fh>) {
        if (/perbatch.lmiss (\d+\.\d+) /) {
            $miss_perbatch_thres = $1;
        }
    }
    close $miss_batch_fh;

    my $s = '\subsection{Removed Bad Variants}';
    $s .= "$hwe_entire variant(s) have been found that fail the HWE test over the whole collection and $hwe_batch variant(s) that fail the HWE test with the worst batch removed.\\\\";
    $s .= "Additionally, $miss_entire variants with too high missingness " . '($\\geq ' . $miss_entire_thres . '$' . ") over the entire collection";
    $s .= "and $miss_batch variants with batch-wise high missingess " . '($\\geq ' . $miss_perbatch_thres . '$' . ") where identified. ";
    $s .= "$all unique variants are removed from the dataset.";

    return $s;
}

sub draw_definetti_after_QCI {
    my $workdir = shift;
    my $tag = shift;

    my $minipage;

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        chomp;
        if (/collection.hwe (\S+) (\S+) (\S+) /) {

            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $1 . '.jpg}\hfill';
            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $2 . '.jpg}\\\\';
            $minipage .= '\hfill\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $3 . '.jpg}\hfill\\\\';
            last;
        }
    }

    my $s = "\n" . '\subsection{DeFinetti Diagrams after first QC stage}' . "\n";
    return $s . $minipage;
}


1

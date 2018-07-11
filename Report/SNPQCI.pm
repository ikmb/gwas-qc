package Report::SNPQCI;

use strict;
use warnings;

use Data::Dumper;
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
            $minipage = '\begin{figure}\centering';
            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $1 . '.jpg}\hfill';
            $minipage .= '\includegraphics[width=0.5\columnwidth]{' . $workdir . '/' . $2 . '.jpg}\end{figure}';
            $minipage .= '\begin{figure}\centering\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $3 . '.jpg}\end{figure}';
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
        if (/^(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)$/) {
            $snpcounts{$1} = { 'all_fail' => $2, 'all_hwe' => $3, 'oneless_fail' => $4, 'oneless_hwe' => $5 };
        }
    }
    close $fh1;

    # Read per-batch table
    open my $fh2, '<', $workdir . '/' . $name . '_exclude-per-batch.FDRthresholds.SNPQCI.2.txt' or die($!);
    <$fh2>;
    while(<$fh2>) {
        chomp;
        if (/^(\d)\t(\d+)\t(\d+)$/) {
            $snpcounts{$1}->{'1plus_fail'} = $2;
            $snpcounts{$1}->{'2plus_fail'} = $2;
        }
    }
    close $fh2;

    my $s = <<'HWE_TABLE'
\subsection{HWE Tests}'
\begin{table}
    \centering
    \begin{tabular}{rrrrrrr}
    \toprule
    & \multicolumn{2}{r}{All Batches} & \multicolumn{2}{r}{Worst Batch Removed} & & \\
    FDR & \# Failed & HWE & \# Failed & HWE & Failed in 1+ batch & Failed in 2+ batches\\
    \midrule

HWE_TABLE
        ;

    print Dumper(%snpcounts);

    keys %snpcounts; # reset the internal iterator so a prior each() doesn't affect the loop
    while((my $k, my $v) = each %snpcounts) {
        $s .= $k . "&" . $v->{'all_fail'} . "&" . $v->{'all_hwe'} . "&" . $v->{"oneless_fail"} . "&" . $v->{"oneless_hwe"} . "&" . $v->{"1plus_fail"} . "&" . $v->{"2plus_fail"} . '\\';
    }

    $s .= '\bottomrule\end{tabular}\end{tabular}';
    return $s;
}

sub exclude_bad_variants {
    my $workdir = shift;
    my $tag = shift;
    return "";
}

sub draw_definetti_after_QCI {
    my $workdir = shift;
    my $tag = shift;
    return "";
}


1

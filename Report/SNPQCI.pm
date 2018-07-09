package Report::SNPQCI;

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

    my $s = "\\section{Batch Merge}";

    map { s/\_/\\\_/ } @batches;
    $s .= "After merging the batches " . join(", ", map { "\\texttt{$_}" } @batches) . ", we have $variants variants and $samples samples, of which $cases are cases and $controls are controls.";
    return $s;
}

sub generate_hwe_diagrams {
    my $workdir = shift;
    my $tag = shift;

    my @images = ();

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        chomp;
        print STDERR "$_\n";
        if (/hardy.hwe (\S+) (\S+) (\S+) /) {
            push @images, '\begin{figure}\centering\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $1 . '.jpg}\end{figure}';
            push @images, '\begin{figure}\centering\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $2 . '.jpg}\end{figure}';
            push @images, '\begin{figure}\centering\includegraphics[height=0.45\textheight]{' . $workdir . '/' . $3 . '.jpg}\end{figure}';
            last;
        }
    }

    print STDERR @images;

    my $s = '\subsection{DeFinetti Diagrams}';
    return $s . join('', @images);

}

sub exclude_lists_for_failed_hwe {
    my $workdir = shift;
    my $tag = shift;

}

sub exclude_bad_variants {
    my $workdir = shift;
    my $tag = shift;

}

sub draw_definetti_after_QCI {
    my $workdir = shift;
    my $tag = shift;

}


1

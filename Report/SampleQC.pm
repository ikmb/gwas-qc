package Report::SampleQC;

use strict;
use warnings;

use Data::Dumper;
use File::Slurp::Tiny qw/ read_file /;

our $VERSION = '1.00';

sub apply_precalc_remove_list {
    my $workdir = shift;
    my $tag = shift;

    my $s = '\newpage\section{Sample QC}';
    open my $sh, '<', $workdir . '/manually-removed.log' or die($!);
    my $pre = 0;
    my $post = 0;
    my $cases = 0;
    my $controls = 0;

    while(<$sh>) {
        if (/^(\d+) people/) {
            $pre = $1;
            next;
        }

        if (/^\d+ variants and (\d+) people pass filters and QC/) {
            $post = $1;
            next;
        }

        if (/^Among remaining phenotypes, (\d+) are cases and (\d+) are controls/) {
            $cases = $1;
            $controls = $2;
        }
    }
    close $sh;

    $s .= '\subsection{Prefiltering}';
    if ($post == $pre) {
        $s .= 'No (valid) samples where specified in the exclusion list.';
    } else {
        $s .= ($pre-$post) . ' samples where removed because they were specified in the exclusion list. ';
        $s .= $post . ' samples remain in the dataset, where ' . $cases . ' are cases and ' . $controls . ' are controls.';
    }
}

sub determine_miss_het {
    my $workdir = shift;
    my $tag = shift;

    my $s;
    my $name;

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        if (/out (\w+)_miss/) {
            $name = $1;
        }
    }
    close $sh;

    $name = $workdir . "/$name";

    $s .= '\subsection{Missingness}';
#    $s .= '\begin{figure}[H!]\centering';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.lmiss.1.png,read=.lmiss.1.png]{' . $name . '_miss}';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.lmiss.2.png,read=.lmiss.2.png]{' . $name . '_miss}\\\\';
#    $s .= '\caption{Missingness with threshold}';
#    $s .= '\end{figure}';
#    $s .= '\begin{figure}[H!]\centering';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.lmiss.logscale.1.png,read=.lmiss.logscale.1.png]{' . $name . '_miss}';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.lmiss.logscale.2.png,read=.lmiss.logscale.2.png]{' . $name . '_miss}\\\\';
#    $s .= '\caption{Missingness with threshold, log scale}';
#    $s .= '\end{figure}';

    $s .= '\subsection{Heterozygosity}';
#    $s .= '\begin{figure}[H!]\centering';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.het.1.png,read=.het.1.png]{' . $name . '_het}';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.het.2.png,read=.het.2.png]{' . $name . '_het}\\\\';
#    $s .= '\caption{Heterozygosity with threshold}';
#    $s .= '\end{figure}';
#    $s .= '\begin{figure}[H!]\centering';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.het.logscale.1.png,read=.het.logscale.1.png]{' . $name . '_het}';
    $s .= '\includegraphics[width=0.5\columnwidth,type=png,ext=.het.logscale.2.png,read=.het.logscale.2.png]{' . $name . '_het}';
#    $s .= '\caption{Heterozygosity with threshold, log scale}';
#    $s .= '\end{figure}';
    return $s;
}

sub prune {
    my $workdir = shift;
    my $tag = shift;

    my $name;

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        if (/out "(\w+_pruned)"/) {
            $name = $1;
        }
    }
    close $sh;

    my $s = '\subsection{Pruning}';

    # Parameters
    my $window_size = 0;
    my $step_size = 0;
    my $rsq = 0;
    my $maf = 0;

    # (s)ample and (v)ariant counts
    my $s_initial = 0;
    my $s_missings_removed = 0; # number of removed samples
    my $s_final = 0;

    my $s_final_cases = 0;
    my $s_final_controls = 0;
    my $v_initial = 0;
    my $v_pruned = 0; # number of variants removed through indep-pairwise
    my $v_maf_removed = 0; # number of variants removed through MAF filter
    my $v_after_intermediate = 0; # number of variants after intermediate filter
    my $v_region_filtered = 0; # number of variants removed in autosomes/LDregion filter
    my $v_final = 0; # final number of variants

    my $rate_before = 0;
    my $rate_after = 0;

    open my $prune_fh, '<', $workdir . '/_prune.log' or die($!);
    while(<$prune_fh>) {
        if (/indep-pairwise (\d+) (\d+) ([0-9\.]+)/) {
            $window_size = $1;
            $step_size = $2;
            $rsq = $3;
        }

        if (/Total genotyping rate is ([0-9\.]+)./) {
            $rate_before = $1;
        }
    }
    close $prune_fh;

    open my $extract_fh, '<', $workdir . '/intermediate.log' or die($!);
    while(<$extract_fh>) {
        if (/maf ([0-9\.]+)/) {
            $maf = $1; next;
        }

        if (/^(\d+) variants loaded from/) {
            $v_initial = $1; next;
        }

        if (/^(\d+) people/) {
            $s_initial = $1; next;
        }

        if (/^--extract.(\d+) variants remaining/) {
            $v_pruned = $v_intial - $1; next;
        }

        if (/^--remove.(\d+) people remaining/) {
            $s_missings_removed = $s_initial - $1; next;
        }

        if (/^(\d+) variants removed due to minor allele threshold/) {
            $v_maf_removed = $1; next;
        }

        if (/^(\d+) variants and (\d+) people pass filters and QC/) {
            $s_final = $2;
            $v_after_intermediate = $1;
            next;
        }

        if (/^Among remaining phenotypes, (\d+) are cases and (\d+) are controls./) {
            $s_final_cases = $1;
            $s_final_controls = $2;
            next;
        }
    }

    open my $final_fh, '<', $workdir . "/" . $name . ".log" or die($!);
    while(<$final_fh>) {
        if (/Total genotpying rate is ([0-9\.]+)./) {
            $rate_after = $1;
        }

        if (/^(\d+) variants and/) {
            $v_region_filtered = $v_after_intermediate - $1;
            $v_final = $1;
        }
    }

    # Jetzt noch alles in eine Tabelle gießen

}

# DKTNF_POPGEN_GS_SampleQCI_het.het                 DKTNF_POPGEN_GS_SampleQCI_miss.lmiss.1.png
# DKTNF_POPGEN_GS_SampleQCI_het.het.1.png           DKTNF_POPGEN_GS_SampleQCI_miss.lmiss.2.png
# DKTNF_POPGEN_GS_SampleQCI_het.het.2.png           DKTNF_POPGEN_GS_SampleQCI_miss.lmiss.logscale.1.png
# DKTNF_POPGEN_GS_SampleQCI_het.het.logscale.1.png  DKTNF_POPGEN_GS_SampleQCI_miss.lmiss.logscale.2.png
# DKTNF_POPGEN_GS_SampleQCI_het.het.logscale.2.png  DKTNF_POPGEN_GS_SampleQCI_miss.log
# DKTNF_POPGEN_GS_SampleQCI_het.het.outlier.txt     DKTNF_POPGEN_GS_SampleQCI_miss.nosex
# DKTNF_POPGEN_GS_SampleQCI_het.log                 DKTNF_POPGEN_GS_SampleQCI_miss.outlier.txt
# DKTNF_POPGEN_GS_SampleQCI_het.nosex               manually-removed.bed
# DKTNF_POPGEN_GS_SampleQCI_miss.imiss              manually-removed.bim
# DKTNF_POPGEN_GS_SampleQCI_miss.lmiss              manually-removed.fam


1

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
    $s .= '\subsection{Huynh-Feldt Correction for Repeated Measure ANOVA}';

    my $excludes = `wc -l < $workdir/hf-excludes`;
    chomp($excludes);
    $s .= "$excludes SNPs were found to violate the sphericity assumption and are removed from the dataset. ";


    return $s;
}


# Parser Regexes are different. Sample QC uses Plink 1.9, this uses 1.07
sub prune {
    my $workdir = shift;
    my $tag = shift;
    my $intro = shift // '\subsection{Sample Outlier Detection}';

    my $name;

    open my $sh, '<', $workdir . '/.command.sh' or die($!);
    while(<$sh>) {
        if (/out "(\w+_pruned)"/) {
            $name = $1;
        }
    }
    close $sh;

    my $s = $intro . '\subsubsection{Pruning}';

    # Parameters
    my $window_size = 0;
    my $step_size = 0;
    my $rsq = 0;
    my $maf = 0;

    # (s)ample and (v)ariant counts
    my $s_initial = 0;
#    my $s_missings_removed = 0; # number of removed samples
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

        if (/Total genotyping rate in remaining individuals is ([0-9\.]+)./) {
            $rate_before = $1;
        }
    }
    close $prune_fh;

    open my $extract_fh, '<', $workdir . '/intermediate.log' or die($!);
    while(<$extract_fh>) {
        if (/maf ([0-9\.]+)/) {
            $maf = $1; next;
        }

        if (/^(\d+) markers to be included from/) {
            $v_initial = $1; next;
        }

        if (/^(\d+) individuals read from/) {
            $s_initial = $1; next;
        }

        if (/^Before frequency and genotyping pruning, there are (\d+) SNPs/) {
            $v_pruned = $v_initial - $1; next;
        }

#        if (/^--remove..(\d+) people remaining/) {
 #           $s_missings_removed = $s_initial - $1; next;
  #      }

        if (/^(\d+) SNPs failed frequency test/) {
            $v_maf_removed = $1; next;
        }

        if (/^After frequency and genotyping pruning, there are (\d+) SNPs/) {#}(\d+) variants and (\d+) people pass filters and QC/) {
#            $s_final = $2;
            $v_after_intermediate = $1;
            next;
        }

        if (/^After filtering, (\d+) cases, (\d+) controls and 0 missing/) {
            $s_final = $1 + $2;
            $s_final_cases = $1;
            $s_final_controls = $2;
            next;
        }
    }

    open my $final_fh, '<', $workdir . "/" . $name . ".log" or die($!);
    while(<$final_fh>) {
        if (/Total genotyping rate in remaining individuals is ([0-9\.]+)./) {
            $rate_after = $1;
        }

        if (/^After frequency and genotyping pruning, there are (\d+) SNPs/) {
            $v_region_filtered = $v_after_intermediate - $1;
            $v_final = $1;
        }
    }

    # Jetzt noch alles in eine Tabelle gießen
    $s .= 'For sample outlier detection, a separate pruned dataset is created to carry out principal component analyses and IBD/IBS detection. The initial dataset is prepared as follows:\\\\';
    $s .= '\begin{tabularx}{\textwidth}{rrX}\toprule{}';
    $s .= 'Samples & Variants & Process and parameters\\\\\midrule{}';
    $s .= "$s_initial & $v_initial & \\emph{Initial dataset, genotyping rate $rate_before}\\\\";
#    $s .= ($s_initial - $s_missings_removed) . " & $v_initial & Samples removed with high missingness (see \\ref{sec:sampleqc-missingness})\\\\";
    # $s .= "$s_final & $v_initial & Missing-sex samples removed\\\\";
    $s .= "$s_final & " . ($v_initial - $v_pruned) . " & LD pruning (" . 'r$^2$' . "=$rsq, window size $window_size, step size $step_size)\\\\";
    $s .= "$s_final & " . ($v_initial - $v_pruned - $v_maf_removed) . " & MAF filtering (threshold $maf)\\\\";
    $s .= "$s_final & " . ($v_after_intermediate - $v_region_filtered) . " & Remove SNPs from non-autosomes, SNPs from xMHC regions, A/T and C/G SNPs, and D/I SNPs\\\\";
    $s .= "$s_final & $v_final & \\emph{Pruned dataset, genotyping rate $rate_after}\\\\\\bottomrule\\end{tabularx}\\\\[1ex]";
    $s .= "After pruning, the dataset consists of $s_final_cases cases and $s_final_controls controls.";
    return $s;
}


sub flashpca_pruned {
    my $workdir = shift;
    my $tag = shift;

    my $basename;
    my $alt_basename;
    my $num_pcs;

    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {
#        print $_;
        chomp;

        if (/merge__new_plink_collection_pruned__1kG\(\S+, "(\S+)", .*\)/) {
            $basename = $1;
            next;
        }

        if (/flashpca2 -d (\d+) --bfile "(\S+?)"/) {
            $alt_basename = $2;
            $num_pcs = $1;
            last;
        }

        if (/flashpca2 -d (\d+)/) {
            $num_pcs = $1;
            next;
        }
    }

    $basename = $basename // $alt_basename;

    my $s = '\subsubsection{PCA of cases and controls only}';


    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.2PC.png,read=.2PC.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.country.2PC.png,read=.country.2PC.png]{$workdir/$basename}";

    return $s;
}


sub draw_histograms {
    my $workdir = shift;
    my $tag = shift;

    my $basename;

    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {

        if (/^cp (\S+).pca.evec \S+.evec.pca.evec$/) {
            $basename = $1;
            last;
        }
    }
    close $sh;

    my $s;
    $s = '\subsubsection{PCA Histograms}Phenotype-wise histograms of the first 4 PCs:\\\\';
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC1.png,read=.pca.evec.histPC1.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC2.png,read=.pca.evec.histPC2.png]{$workdir/$basename}\\\\";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC3.png,read=.pca.evec.histPC3.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC4.png,read=.pca.evec.histPC4.png]{$workdir/$basename}\\\\";
    $s .= 'Country-wise histograms of the first 4 PCs:\\\\';
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC1.png,read=.pca.evec.histPC1.png]{$workdir/$basename.country}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC2.png,read=.pca.evec.histPC2.png]{$workdir/$basename.country}\\\\";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC3.png,read=.pca.evec.histPC3.png]{$workdir/$basename.country}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.pca.evec.histPC4.png,read=.pca.evec.histPC4.png]{$workdir/$basename.country}\\\\";

    return $s;
}

sub compile_variants_exclusion {
    my $workdir = shift;
    my $tag = shift;

    my $s = '\subsection{Final Cleaning}';
    my $extra = '';
    my $basename;
    my $final_list;


    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {
        if (/gawk.* "(\S+).monomorphic" ".*.nearly_monomorphic" ".*.differential_missingness".*>"(\S+)"/) {
            $basename = $1;
            $final_list = $2;
            next;
        }

        if (/"True" == "True"/) {
            $extra = '(ignored in controls-only setting)';
        }
    }
    close $sh;


    my $monomorphics = `wc -l < $workdir/$basename.monomorphic`; chomp($monomorphics);
    my $nearly_monomorphics = `wc -l < $workdir/$basename.nearly_monomorphic`; chomp($nearly_monomorphics);
    my $diff_missingness = `wc -l < $workdir/$basename.differential_missingness`; chomp($diff_missingness);
    my $total = `wc -l < $workdir/$final_list`; chomp($total);

    $s .= "The following SNPs have been removed:\\\\[1ex]";
    $s .= '\begin{tabular}{r|l}';
    $s .= '\toprule{}';
    $s .= $monomorphics . ' & Monomorphic SNPs ' . $extra . '\\\\';
    $s .= $nearly_monomorphics . ' & Nearly monomorphic SNPs ' . $extra . '\\\\';
    $s .= $diff_missingness . ' & SNPs with high differential missingness\\\\\\midrule';
    $s .= $total . ' & Total unique SNPs to be removed\\\\\\bottomrule\end{tabular}\\\\[1ex]';


    return $s;
}

sub final_cleaning {
    my $workdir = shift;
    my $tag = shift;


    my $cases;
    my $controls;
    my $males;
    my $females;
    my $gt_rate;
    my $ind_removed;
    my $snps;

    open my $log, '<', "$workdir/.command.log" or die($!);
    while(<$log>) {
        if (/^(\d+) individuals removed with/) {
            $ind_removed = $1;
            next;
        }

        if (/^Total genotyping rate in remaining individuals is ([0-9\.]+)/) {
            $gt_rate = $1;
            next;
        }

        if (/After filtering, (\d+) cases, (\d+) controls and 0 missing/) {
            $cases = $1;
            $controls = $2;
            next;
        }

        if (/After filtering, (\d+) males, (\d+) females, and \d+ of unspecified sex/) {
            $males = $1;
            $females = $2;
            next;
        }

        if (/After frequency and genotyping pruning, there are (\d+) SNPs/) {
            $snps = $1;
            next;
        }
    }

    my $s = "After removing $ind_removed sample(s) because of unknown diagnoses or manual exclusion, the final dataset consists of $cases cases and $controls controls ($males males, $females females, ";
    $s .= ($males + $females) . " in total) and $snps SNPs with a genotyping rate of $gt_rate.";
}

1

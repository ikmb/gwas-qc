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

    $s .= '\subsection{Missingness}\label{sec:sampleqc-missingness}';
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

        if (/^--extract..(\d+) variants remaining/) {
            $v_pruned = $v_initial - $1; next;
        }

        if (/^--remove..(\d+) people remaining/) {
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
        if (/Total genotyping rate is ([0-9\.]+)./) {
            $rate_after = $1;
        }

        if (/^(\d+) variants and/) {
            $v_region_filtered = $v_after_intermediate - $1;
            $v_final = $1;
        }
    }

    # Jetzt noch alles in eine Tabelle gießen
    $s .= 'For sample outlier detection, a separate pruned dataset is created to carry out principal component analyses and IBD/IBS detection. The initial dataset is prepared as follows:\\\\';
    $s .= '\begin{tabularx}{\textwidth}{rrX}\toprule{}';
    $s .= 'Samples & Variants & Process and parameters\\\\\midrule{}';
    $s .= "$s_initial & $v_initial & \\emph{Initial dataset, genotyping rate $rate_before}\\\\";
    $s .= ($s_initial - $s_missings_removed) . " & $v_initial & Samples removed with high missingness (see \\ref{sec:sampleqc-missingness})\\\\";
    # $s .= "$s_final & $v_initial & Missing-sex samples removed\\\\";
    $s .= "$s_final & " . ($v_initial - $v_pruned) . " & LD pruning (" . 'r$^2$' . "=$rsq, window size $window_size, step size $step_size)\\\\";
    $s .= "$s_final & " . ($v_initial - $v_pruned - $v_maf_removed) . " & MAF filtering (threshold $maf)\\\\";
    $s .= "$s_final & " . ($v_after_intermediate - $v_region_filtered) . " & Remove SNPs from non-autosomes, SNPs from xMHC regions, A/T and C/G SNPs, and D/I SNPs\\\\";
    $s .= "$s_final & $v_final & \\emph{Pruned dataset, genotyping rate $rate_after}\\\\\\bottomrule\\end{tabularx}\\\\[1ex]";
    $s .= "After pruning, the dataset consists of $s_final_cases cases and $s_final_controls controls.";
    return $s;
}

sub hapmap_eigenstrat {
    my $workdir = shift;
    my $tag = shift;

    my $basename;
    my $sigma_threshold;
    my $num_pcs;
    
    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {
        chomp;
        if (/pca_run\("(\S+?)", ([0-9\.]+), \S+?, (\d+), \d+, .*\)/) {
            $basename = $1;
            $sigma_threshold = $2;
            $num_pcs = $3;
            last;
        }
    }
    $basename .= "_$num_pcs" . 'PC';

    my $s = '\subsubsection{PCA with merged HapMap2 samples}';
    $s .= 'The images on the left side are with, the others are without projection on HapMap samples.\\\\';

    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.2PC.png,read=.2PC.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.2PC.png,read=.2PC.png]{$workdir/$basename.withoutProjection}\\\\";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.4PCpairs.png,read=.4PCpairs.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.4PCpairs.png,read=.4PCpairs.png]{$workdir/$basename.withoutProjection}\\\\";
    return $s;
}

# Note: This function is also used from SNPQCII, with num_pcs = 1000
sub onekg_flashpca {
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

    my $s;
    my $type;
    if( $num_pcs < 1000 ){
        $type = "pdf";
        $basename .= '_' . $num_pcs . 'PC';
        $s = '\subsubsection{PCA with merged 1000 Genomes samples}';
    } else {
        $type = "png";
        $s = '\subsubsection{PCA of cases and controls only}';
    }

    $s .= "\\includegraphics[width=0.5\\textwidth,type=pdf,ext=.2PC.$type,read=.2PC.$type]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=pdf,ext=.country.2PC.$type,read=.country.2PC.$type]{$workdir/$basename}";

    if (-e "$workdir/$basename.fail-pca-1KG-qc.txt") {
        my $count_fail = `wc -l <$workdir/$basename.fail-pca-1KG-qc.txt` or die($!);
        chomp($count_fail);
        my $count_fail_country = `wc -l <$workdir/$basename.country.fail-pca-1KG-qc.txt` or die($!);
        chomp($count_fail_country);
        $s .= "\\\\The batch- and country-wise PCA yielded $count_fail and $count_fail_country outliers, respectively.";
    }


    return $s;
}

sub remove_bad_samples {
    my $workdir = shift;
    my $tag = shift;

    my $s = '\subsubsection{Outlier Removal}';
    $s .= 'The following outliers are scheduled for removal:\\\\[1ex]';
    $s .= '\begin{tabular}{rl}\toprule{}';

    my @files;
    open my $fh, '<', "$workdir/.command.sh" or die($!);
    while(<$fh>) {
        if (/^cat (\S+)/) {
            push @files, $1;
            next;
        }

        if (/cut -f 1,2 (\S+)/) {
            push @files, $1;
            next;
        }

        if (/cut .* -f 1,2 (\S+)/) {
            push @files, $1;
            next;
        }
    }
    close $fh;

    my $precalc = "0";
    if (-e "$workdir/$files[0]") {
        $precalc = `wc -l <$workdir/$files[0]`;
        chomp($precalc);
    }

    my $het = `wc -l <$workdir/$files[1]` or die("$files[1]: $!");
    chomp($het);
    my $miss = `wc -l <$workdir/$files[2]` or die("files[2]: $!");
    chomp($miss);
    my $duplicates = `wc -l <$workdir/$files[3]` or die("files[3]: $!");
    chomp($duplicates);
    my $es_outliers = `wc -l <$workdir/$files[4]` or die("files[4]: $!");
    chomp($es_outliers);
    my $flash_outliers = `wc -l <$workdir/$files[5]` or die("files[5]: $!");
    chomp($flash_outliers);
    my $all = `wc -l <$workdir/remove-samples` or die("remove-samples: $!");

    $s .= 'Count & Process\\\\\midrule{}';
    $s .= "$precalc & Pre-calculated sample exclusion list\\\\";
    $s .= "$het & Heterozygosity outliers\\\\";
    $s .= "$miss & Missingness outliers\\\\";
    $s .= "$duplicates & Duplicates (through IBD analysis)\\\\";
    $s .= "$es_outliers & PCA Outliers (EIGENSTRAT)\\\\";
    $s .= "$flash_outliers & PCA Outliers (FlashPCA)\\\\" . '\midrule{}';
    $s .= ($precalc + $het + $miss + $duplicates + $es_outliers + $flash_outliers) . " & Total\\\\";
    $s .= "$all & Total (unique)\\\\" . '\bottomrule{}';
    $s .= "\\end{tabular}\\\\";
    return $s;


}

sub draw_histograms {
    my $workdir = shift;
    my $tag = shift;

    my $basename;

    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {

        if (/^cp (\S+).evec \S+.evec.pca.evec$/) {
            $basename = $1;
            last;
        }
    }
    close $sh;

    my $s;
    $s .= '\subsubsection{PCA Histograms}Phenotype-wise histograms of the first 4 PCs:\\\\';
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.evec.histPC1.png,read=.evec.histPC1.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.evec.histPC2.png,read=.evec.histPC2.png]{$workdir/$basename}\\\\";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.evec.histPC3.png,read=.evec.histPC3.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.evec.histPC4.png,read=.evec.histPC4.png]{$workdir/$basename}\\\\";
    $s = 'Country-wise histograms of the first 4 PCs:\\\\';
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.country.evec.histPC1.png,read=.country.evec.histPC1.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.country.evec.histPC2.png,read=.country.evec.histPC2.png]{$workdir/$basename}\\\\";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.country.evec.histPC3.png,read=.country.evec.histPC3.png]{$workdir/$basename}";
    $s .= "\\includegraphics[width=0.5\\textwidth,type=png,ext=.country.evec.histPC4.png,read=.country.evec.histPC4.png]{$workdir/$basename}\\\\";

    return $s;
}

sub prune_related {
    my $workdir = shift;
    my $tag = shift;

    my $cases_before;
    my $controls_before;
    my $cases_after;
    my $controls_after;
    my $samples_before;
    my $samples_removed;

    my $log;

    open $log, '<', "$workdir/.command.log" or die("$log: $!");

    while(<$log>) {
        chomp;

        if(/^(\d+) individuals read from/) {
            $samples_before = $1;
            next;
        }

        if(/^(\d+) cases, (\d+) controls and \d+ missing$/) {
            $cases_before = $1;
            $controls_before = $2;
            next;
        }

        if(/^(\d+) individuals removed with/) {
            $samples_removed = $1;
            next;
        }

        if(/^After filtering, (\d+) cases, (\d+) controls and \d+ missing/) {
            $cases_after = $1;
            $controls_after = $2;
            next;
        }
    }
    close $log;

    my $s = '\subsection{Dataset Without Relatives}';
    $s .= 'For subsequent analyses, another secondary dataset is created, where the relatives detected during IBD/IBS analysis are removed. ';
    $s .= "From the original dataset ($cases_before cases, $controls_before controls), $samples_removed individuals are considered to be related and were removed, leaving $cases_after cases and $controls_after controls. ";

    return $s;
}

sub prune_outliers_without_related {
    my $workdir = shift;
    my $tag = shift;
    return prune($workdir, $tag, '');
}

sub ibs_merge_and_verify_wr {
    my $workdir = shift;
    my $tag = shift;

    my $s = '\subsubsection{IBD Z0/Z1 Plot}';
    my $basename;
    open my $sh, '<', "$workdir/.command.sh" or die($!);
    while(<$sh>) {
        if (/^OUTFILE="(\S+)"$/) {
            $basename = $1;
        }
    }

    $s .= "\\includegraphics[width=\\textwidth,type=png,ext=.Z0.Z1.IBD-plot.png,read=.Z0.Z1.IBD-plot.png]{$workdir/$basename}";
    return $s;
}

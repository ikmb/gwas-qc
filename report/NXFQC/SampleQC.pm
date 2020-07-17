package NXFQC::SampleQC;

our $VERSION = "1.00";

use NXFQC::Process;
use NXFQC::PlinkLog;
use NXFQC::PlinkInfo;

# From the Perl core distribution
use Data::Dumper;
#use File::Slurp::Tiny qw(read_file);
use File::Basename;
#use List::MoreUtils qw(uniq);

use strict;
use warnings;

sub sanitize {
  my $s = shift;
  $s =~ s/_/\\_/g;
  $s
}

sub sanitize_img {
  my $f = shift;
  $f =~ /(.*)(\.[a-zA-Z]+)$/;
  "{$1}$2";
}

# Taken from List::MoreUtils
sub uniq (@) {
    my %seen = ();
    my $k;
    my $seen_undef;
    grep { defined $_ ? not $seen{$k = $_}++ : not $seen_undef++  } @_;
}

sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or die("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}
sub new {
    my $class = shift;
    my $trace = shift;

    my $self = {'trace' => $trace};

    bless $self, $class;
    return $self;
}

sub countlines {
  my $filename = shift;
  open my $fh, '<', $filename or die("$!: $filename");
  my $count = 0;
  $count += tr/\n/\n/ while sysread($fh, $_, 2**16);
  $count;
}

sub count_unique {
    my @files = @_;
    my @lines;
    foreach(@files) {
        open my $handle, '<', $_ or die($!);
        push @lines, <$handle>;
        close $handle;
    }
    my @uniquelines = uniq @lines;
    return scalar @uniquelines;
}

sub miss_het {
    my $self = shift;
    my $dir = $self->{'trace'}->{'determine_miss_het'};
    my $dat = {};
    my $miss_basename;
    my $het_basename;

    open my $sh, '<', "$dir/.command.sh";
    while(<$sh>) {
        if(/plink.*--out (\S+) --missing/) {
            $miss_basename = "$dir/$1";
            next;
        }
        if(/plink.*--out (\S+) --het/) {
            $het_basename = "$dir/$1";
            next;
        }
        if(/^R .*imiss (\S+) < .*logimiss\.r/) {
            $dat->{'het_threshold'} = $1;
            next;
        }

        if(/^R .*lmiss (\S+) < .*lmiss\.r/) {
            $dat->{'miss_threshold'} = $1;
        }
    }

    $dat->{'het_outliers'} = countlines("${het_basename}.het.outlier.txt");
    $dat->{'miss_outliers'} = countlines("${miss_basename}.outlier.txt");
    $dat->{'img_het1'} = $het_basename . ".het.1.png";
    $dat->{'img_het2'} = $het_basename . ".het.2.png";
    $dat->{'img_het1log'} = $het_basename . ".het.logscale.1.png";
    $dat->{'img_het2log'} = $het_basename . ".het.logscale.2.png";
    $dat->{'img_miss1'} = $miss_basename . ".lmiss.1.png";
    $dat->{'img_miss2'} = $miss_basename . ".lmiss.2.png";
    $dat->{'img_miss1log'} = $miss_basename . ".lmiss.logscale.1.png";
    $dat->{'img_miss2log'} = $miss_basename . ".lmiss.logscale.2.png";
    $dat;

}

sub proc_ibs_dir {
    my $dir = shift;
    my $basename;
    my $dat = {};

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if (/^OUTFILE="(.*)"$/) {
            $basename = "$dir/$1";
            last;
        }
    }
    $dat->{'ibdplot'} = $basename . ".Z0.Z1.IBD-plot.png";
    $dat;
}

sub pca_run {
    my $self = shift;
    my $dir = $self->{'trace'}->{'pca_with_hapmap'};
    my $dat = {};

    my $num_pcs;
    my $basename;

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        chomp;
        if (/^python.*pca_run\("(\S+)",.*?,.".*?",.(\d+)/) {
            $num_pcs = $2;
            $basename = $1;
            last;
        }
    }

    $dat->{'pcs'} = $num_pcs;
    $dat->{'2pc'} = "$dir/${basename}_${num_pcs}PC.2PC.png";
    $dat->{'4pc'} = "$dir/${basename}_${num_pcs}PC.4PCpairs.png";
    $dat->{'2pc_wp'} = "$dir/${basename}_${num_pcs}PC.withoutProjection.2PC.png";
    $dat->{'4pc_wp'} = "$dir/${basename}_${num_pcs}PC.withoutProjection.4PCpairs.png";
    $dat->{'plot'} = "$dir/${basename}_${num_pcs}.plot.pdf";
    $dat;
}

sub flashpca1kg {
    my $self = shift;
    my $dir = $self->{'trace'}->{'flashpca2_pruned_1kG'};
    my $dat = {};
    my $basename;

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if(/^flashpca2 -d (\d+) --bfile "(\S+)" --outval/) {
            $dat->{'pcs'} = $1;
            $basename = $2;
        }
    }

    $dat->{'batch'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.2PC.pdf";
    $dat->{'country'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.country.2PC.pdf";
    $dat->{'batch_fail'} = countlines("$dir/${basename}_" . $dat->{'pcs'} . "PC.fail-pca-1KG-qc.txt");
    $dat->{'country_fail'} = countlines("$dir/${basename}_" . $dat->{'pcs'} . "PC.country.fail-pca-1KG-qc.txt");
    $dat->{'unique'} = count_unique(("$dir/${basename}_" . $dat->{'pcs'} . "PC.fail-pca-1KG-qc.txt", "$dir/${basename}_" . $dat->{'pcs'} . "PC.country.fail-pca-1KG-qc.txt"));
    $dat;
}

sub flashpca {
    my $self = shift;
    my $dir = $self->{'trace'}->{'flashpca2_pruned'};
    my $dat = {};
    my $basename;

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if(/^flashpca2 -d (\d+) --bfile "(\S+)" --outval/) {
            $dat->{'pcs'} = $1;
            $basename = $2;
        }
    }

    $dat->{'batch'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.2PC.png";
    $dat->{'country'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.country.2PC.png";
    $dat->{'batch_4pc'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.4PCpairs.png";
    $dat->{'country_4pc'} = "$dir/${basename}_" . $dat->{'pcs'} . "PC.country.4PCpairs.png";

    $dat;

}

sub remove_bad_samples {
    my $self = shift;
    my $dir = $self->{'trace'}->{'remove_bad_samples'};
    my $dat = {};

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if(/-e (.*) \]/) {
            if(-e $1) {
                $dat->{'manual'} = countlines("$1");
            } else {
                $dat->{'manual'} = 0;
            }
            next;
        }

        if(/^cat (.*.het.outlier.txt)/) {
            $dat->{'het'} = countlines("$dir/$1");
            next;
        }

        if(/^\s+cat (.*_miss.outlier.txt)/) {
            $dat->{'miss'} = countlines("$dir/$1");
            next;
        }

        if(/^cut -d.*1,2 (\S+_duplicates.txt)\s/) {
            $dat->{'dup'} = countlines("$dir/$1");
            next;
        }

        if(/^cat (\S+PC.outlier.txt)/) {
            $dat->{'eigenstrat'} = countlines("$dir/$1");
            next;
        }

        if(/^cat (\S+1KG-qc.txt)/) {
            $dat->{'flashpca'} = countlines("$dir/$1");
            next;
        }

        if(/^plink.*--out (\S+) --allow/) {
            $dat->{'plink'} = parse_plink("$dir/$1.log");
            next;
        }
    }
    $dat->{'unique'} = countlines("$dir/remove-samples");

    $dat->{'info'} = plink_table("$dir/info.txt");

    $dat;
}

sub draw_histograms {
    my $self = shift;
    my $dir = $self->{'trace'}->{'draw_histograms'};
    my $dat = {};
    my $max_pc = 20;
    my $got_batch_flag = 0;

    my $batch_basename;
    my $country_basename;

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if(/^cp (\S+)\s/ and not $got_batch_flag) {
            $batch_basename ="$dir/$1";
            $got_batch_flag = 1;
            next;
        }

        if(/^cp (\S+)\s/ and $got_batch_flag) {
            $country_basename = "$dir/$1";
            last;
        }
    }

    $dat->{'batch'} = ();
    $dat->{'country'} = ();
    for(my $i = 1; $i < ($max_pc+1); $i++) {
        if(-e "${batch_basename}.histPC$i.png") {
            push(@{$dat->{'batch'}}, "${batch_basename}.histPC$i.png");
            push(@{$dat->{'country'}}, "${country_basename}.histPC$i.png");
        }
    }
    $dat;

}

sub detect_duplicates_related {
    my $self = shift;
    my $dir = $self->{'trace'}->{'detect_duplicates_related'};
    my $dat = {};

    my ($a,$b,$relatives,$duplicates);

    open my $fh, '<', "$dir/.command.sh" or die($!);
    while(<$fh>) {
        if(/python.*imiss", "(\d+\.\d+)", "(\d+\.\d+)", .*"(\S+_relatives.txt)".*"(\S+_duplicates.txt)"/) {
            $a = $1;
            $b = $2;
            $relatives = "$dir/$3";
            $duplicates = "$dir/$4";
        }
    }

    $dat->{'thres_dup'} = $a;
    $dat->{'thres_rel'} = $b;
    $dat->{'dup'} = countlines($duplicates);
    $dat->{'rel'} = countlines($relatives);
    $dat;
}

sub remove_relatives {
    my $self = shift;
    my $dir_ibs = $self->{'trace'}->{'ibs_merge_and_verify'};
    my $dir_ibs_wr = $self->{'trace'}->{'ibs_merge_and_verify_withoutRelatives'};
    my $dir_rem = $self->{'trace'}->{'remove_relatives'};
    my $dat = {};

    print Dumper($self);
    open my $fh, '<', "$dir_ibs/.command.sh" or die($!);
    $dat->{'ibs_img'} = (map { /^OUTFILE="(\S+)"/ ? "$dir_ibs/$1.png" : () } <$fh>)[0];
    close $fh;

    open $fh, '<', "$dir_ibs_wr/.command.sh" or die($!);
    $dat->{'ibs_img_wr'} = (map { /^OUTFILE="(\S+)"/ ? "$dir_ibs_wr/$1.png" : () } <$fh>)[0];
    close $fh;

    open $fh, '<', "$dir_rem/.command.sh" or die($!);
    $dat->{'plink'} = parse_plink((map { /^plink .*--out "(\S+)"/ ? "$dir_rem/$1.log" : ()  } <$fh>)[0]);

    $dat->{'info'} = plink_table("$dir_rem/info.txt");
    $dat;
}

sub add_images {
    my @images = @_;
    my $mboxwidth = 1.0 / @images;
    my $s = '';
    foreach(@images) {
        $s .= '\begin{minipage}{' . $mboxwidth . '\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($_) . '}\end{minipage}';
    }
    $s .= '\\\\';
    $s;
}

sub build_report_chunk {
    my $self = shift;

    my $s = '\section{Sample QC}';
    my $misshet = $self->miss_het();

    $s .= '\subsection{Missingness}';
    $s .= add_images(($misshet->{'img_miss1'}, $misshet->{'img_miss2'}));
    $s .= add_images(($misshet->{'img_miss1log'}, $misshet->{'img_miss2log'}));
    $s .= "In the sample missingness analysis, $misshet->{miss_outliers} individuals have been classified as outliers with respect to a threshold of $misshet->{'miss_threshold'}.";
    $s .= '\subsection{Heterozygosity}';
    $s .= add_images(($misshet->{'img_het1'}, $misshet->{'img_het2'}));
    $s .= add_images(($misshet->{'img_het1log'}, $misshet->{'img_het2log'}));
    $s .= "In the heterozygosity analysis, $misshet->{het_outliers} individuals have been classified as outliers, parting more than 5*SD from the mean.";

    $s .= '\subsection{PCA with HapMap2 Projection}';
    my $pca = $self->pca_run();
    $s .= $pca->{'pcs'} . ' principal components have been analyzed.\\\\';
    $s .= add_images(($pca->{'2pc'}, $pca->{'2pc_wp'}));
    $s .= add_images(($pca->{'4pc'}, $pca->{'4pc_wp'}));

    $s .= '\subsection{PCA with 1000 Genomes Projection}';
    $pca = $self->flashpca1kg();
    $s .= add_images(($pca->{'batch'}, $pca->{'country'}));
    my $out_b = $pca->{'batch_fail'};
    my $out_c = $pca->{'country_fail'};
    my $out_total = $out_b + $out_c;
    my $out_uniq = $pca->{'unique'};
    $s .= "$out_b samples for batch-annotated data and $out_c samples for country-annotated data were classified as outliers (total: $out_total, unique total: $out_uniq).";
    $s .= '\subsection{PCA without Projection}';
    $pca = $self->flashpca();
    $s .= '\subsubsection{Batch-annotated PCA}';
    $s .= add_images(($pca->{'batch'})) . '\\\\' .  add_images(($pca->{'batch_4pc'}));
    $s .= '\subsubsection{Country-annotated PCA}';
    $s .= add_images(($pca->{'country'})) . '\\\\' . add_images(($pca->{'country_4pc'}));

    $s .= '\subsection{Principal Component Histograms}';
    my $hist = $self->draw_histograms();
    $s .= '\subsubsection{By Diagnosis}';
    $s .= add_images(@{$hist->{'batch'}}[0..1]);
    $s .= add_images(@{$hist->{'batch'}}[2..3]);
    $s .= add_images(@{$hist->{'batch'}}[4..5]);
    $s .= add_images(@{$hist->{'batch'}}[6..7]);
    $s .= add_images(@{$hist->{'batch'}}[8..9]);
    $s .= '\subsubsection{By Country}';
    $s .= add_images(@{$hist->{'country'}}[0..1]);
    $s .= add_images(@{$hist->{'country'}}[2..3]);
    $s .= add_images(@{$hist->{'country'}}[4..5]);
    $s .= add_images(@{$hist->{'country'}}[6..7]);
    $s .= add_images(@{$hist->{'country'}}[8..9]);
    my $bad = $self->remove_bad_samples();

    my $ibs = $self->detect_duplicates_related();
    my $relatives = $self->remove_relatives();
    $s .= '\subsection{IBS Relatives Detection}';
    $s .= add_images(($relatives->{'ibs_img'}));
    $s .= ($ibs->{'dup'}) . " samples were identified as duplicates or identical twins (\$\\hat\\pi\\geq\$" . $ibs->{'thres_dup'} . "). ";
    $s .= "Additionally, " . $ibs->{'rel'} . " samples are closely related or inbred but not identical (\$\\hat\\pi\\geq\$ " . $ibs->{'thres_rel'} . "). ";

    $s .= '\subsection{Summary: Sample Removal}';
    $s .= '\begin{tabular}{lr}\toprule{}';
    $s .= "Samples before QC:&".$bad->{'plink'}->{'loaded-phenotypes'}." total\\\\";
    $s .= '\midrule{}';
    $s .= "Missingness outliers:&$bad->{miss}\\\\";
    $s .= "Heterozygosity outliers:&$bad->{het}\\\\";
    $s .= "PCA outliers:&" . ($bad->{'flashpca'}) . "\\\\";
    $s .= "Duplicates:&$ibs->{dup}\\\\";
    $s .= "Manual removals:&$bad->{manual}\\\\\\midrule{}";
    $s .= "Total:&" . ($bad->{'miss'} + $bad->{het} + $bad->{'flashpca'} + $bad->{'dup'} + $bad->{'manual'}) . "\\\\";
    $s .= "Unique:&" . ($bad->{'unique'}) . "\\\\\\midrule{}";
    $s .= "Samples after QC:&" . $bad->{'plink'}->{'final-cases'} . " cases, " . $bad->{'plink'}->{'final-controls'} . " controls, " . ($bad->{'plink'}->{'samples-after-remove'}) . " total\\\\";
    $s .= "\\midrule{}Relatives:&" . $ibs->{'rel'} . "\\\\";
    $s .= "Samples after QC, no relatives:& " . $relatives->{'plink'}->{'final-cases'} . " cases, " . $relatives->{'plink'}->{'final-controls'} . " controls, " . ($relatives->{'plink'}->{'samples-after-remove'}) . " total\\\\\\bottomrule{}";
    $s .= "\\end{tabular}\\\\";
    my $removed = $bad->{'unique'} + $ibs->{'rel'};
    $s .= "Of initially " . $bad->{'plink'}->{'loaded-phenotypes'} . " samples, $removed were removed (" . sprintf("%.2f\\,\\%%). ", 100.0*($removed/$bad->{'plink'}->{'loaded-phenotypes'})) . "\\\\";
    $s .= '\subsection{Phase Summary}';
    $s .= '\begin{minipage}{0.5\textwidth}%' . "\n";
    $s .= '\centering\textbf{With Related Individuals}\\\\';
    $s .= $bad->{'info'};
    $s .= '\end{minipage}\begin{minipage}{0.5\textwidth}%' . "\n";
    $s .= '\centering\textbf{Without Related Individuals}\\\\';
    $s .= $relatives->{'info'};
    $s .= '\end{minipage}';
    $s

}

1


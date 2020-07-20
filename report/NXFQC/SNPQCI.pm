# -*- mode:cperl -*-
package NXFQC::SNPQCI;

our $VERSION = "1.00";

use NXFQC::Process;
use NXFQC::PlinkLog;
use NXFQC::PlinkInfo;

# From the Perl core distribution
use Data::Dumper;
#use File::Slurp::Tiny qw(read_file);
use File::Basename;

use strict;
use warnings;

sub sanitize {
  my $s = shift;
  $s =~ s/_/\\_/g;
  $s
}

sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or die("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}
sub sanitize_img {
  my $f = shift;
  $f =~ /(.*)(\.[a-zA-Z]+)$/;
  "{$1}$2";
}

sub countlines {
  my $filename = shift;
  open my $fh, '<', $filename or die("$!: $filename");
  my $count = 0;
  $count += tr/\n/\n/ while sysread($fh, $_, 2**16);
  $count;
}

sub merge_datasets {
    my $self = shift;
    my $dir = $self->{'trace'}->{'merge_batches'};
    my $dat = {};
    $dat->{'multiallelic'} = 0;

    # Collect multi-allelic (>2) count
    open my $logfh, '<', "$dir/.command.log" or die($!);
    while(<$logfh>) {
        if(/Removing (\d+) variants/) {
            $dat->{'multiallelic'} = $1;
            last;
        }
    }
    close $logfh;

    open my $shfh, '<', "$dir/.command.sh" or die($!);
    while(<$shfh>) {
        chomp;
        if(/^BIMS=\((.*)\)/) {
            my @bims = split /\s+/, $1;
            print "### bims: @bims\n";
            my @ds = map { s/_Rs\.bim//g; $_ } @bims;
            print "### datasets: @ds\n";
            $dat->{'datasets'} = \@ds;
            last;
        }
    }
    close $shfh;

    # Collect post-merge info
    $dat->{'info'} = plink_table("$dir/info.txt");
    return $dat;

}

sub definetti_preqc {
  my $self = shift;
  my $dir = $self->{'trace'}->{'hwe_definetti_preqc'};
  my $dat = {};

  open my $fh, '<', "$dir/.command.sh";
  while(<$fh>) {
    if(/R --slave.*hwe (\S+) (\S+) (\S+) </) {
      $dat->{'controls'} = "$dir/$1.jpg";
      $dat->{'cases'} = "$dir/$2.jpg";
      $dat->{'cases/controls'} = "$dir/$3.jpg";
    }
  }
  $dat;
}

sub definetti_postqc {
  my $self = shift;
  my $dir = $self->{'trace'}->{'hwe_definetti_qci'};
  my $dat = {};

  open my $fh, '<', "$dir/.command.sh";
  while(<$fh>) {
    if(/R --slave.*hwe (\S+) (\S+) (\S+) </) {
      $dat->{'controls'} = "$dir/$1.jpg";
      $dat->{'cases'} = "$dir/$2.jpg";
      $dat->{'cases/controls'} = "$dir/$3.jpg";
    }
  }
  $dat;
}

sub hwe_fdr_filter {
    my $self = shift;
    my $dir = $self->{'trace'}->{'hwe_fdr_filter'};

    my $worstbatch;
    my $twoplus;
    my $all;

    my $dat = {};
    open my $cmdfh, '<', "$dir/.command.sh" or die("$!: $dir/.command.sh");
    while(<$cmdfh>) {
      chomp;
      if(/^fdrfilter\.pl .*RESULTS" (\d+).*\s+(\S+all-batches)\s+(\S+batch-removed)\s+\S+\s+(\S+2-plus-batches)\s*$/) {
          #      if(/^fdrfilter\.pl .*RESULTS" (\d+).*\s(\S+exclude-whole-collection-worst-batch-removed)\s+(\S+exclude-per-batch-fail-in-2-plus-batches)/) {
      #      if(/^fdrfilter\s+(\S+FDRthresholds.SNPQCI.1.txt)\s+(\S+FDRthresholds.SNPQCI.2.txt).+(\d+)/) {
        $all = "$dir/$2";
        $dat->{'twoplus-img'} = "$dir/$4.FDRthresholds.SNPQCI.2.txt.png";
        $worstbatch = "$dir/$3";
        $dat->{'worstbatch-img'} = "$dir/$3.FDRthresholds.SNPQCI.1.txt.png";
        $twoplus = "$dir/$4";
        $dat->{'threshold'} = $1;
        last;
      }
    }



    $dat->{'twoplus-count'} = countlines($twoplus);
    $dat->{'worstbatch-count'} = countlines($worstbatch);
	$dat->{'all-count'} = countlines($all);
    $dat;
  }

sub exclude_missingness {
  my $self = shift;
  my $perbatch = $self->{'trace'}->{'determine_missingness_per_batch'};
  my $whole = $self->{'trace'}->{'determine_missingness_entire'};

  my $dat = {};

  my $outfile;
  open my $cmd, '<', "$whole/.command.log" or die($!);
  while(<$cmd>) {
    chomp;
    $dat->{'whole-thres'} = $1 if /Threshold.*(\d\.\d+)$/;
    $outfile = $1 if /Outfile:  (\S+)/;
  }

  $dat->{'whole-count'} = countlines("$whole/$outfile");

  open $cmd, '<', "$perbatch/.command.log" or die($!);
  while(<$cmd>) {
    chomp;
    $dat->{'perbatch-thres'} = $1 if /^Threshold:\s+(\S+)$/;
    $outfile = $1 if /^Outfile:\s+(\S+)$/;
  }

  $dat->{'perbatch-count'} = countlines("$perbatch/$outfile");

  $dat;
}


sub exclude_bad_variants {
  my $self = shift;
  my $dir = $self->{'trace'}->{'exclude_bad_variants'};
  my $dat = {};

  my $logfile;
  open my $fh, '<', "$dir/.command.log" or die($!);
  while(<$fh>) {
    $logfile = $1 if /^Logging to (\S+)\.$/;
  }

  $dat = parse_plink("$dir/$logfile");
  $dat->{'info'} = plink_table("$dir/info.txt");

  $dat;
}

sub new {

    my $class = shift;
    my $trace = shift;

    my $self = {'trace' => $trace };

    bless $self, $class;
    return $self;
}

sub build_report_chunk {
    my $self = shift;

    my $s = "";

    $s = '\section{SNP QC Phase I -- Missingness and Hardy-Weinberg}';
    my $miss = $self->exclude_missingness();
    my $hwe = $self->hwe_fdr_filter();
    my $bad = $self->exclude_bad_variants();

    my $merge = $self->merge_datasets();

    print Dumper($merge);
    $s .= '\subsection{Dataset Merge}' . "%\n";
    my @ds = $merge->{'datasets'};
    if(@ds > 1) {
        $s .= 'The following datasets were merged: ' . sanitize(join(', ', @{$merge->{'datasets'}})) . '. ';

        if($merge->{'multiallelic'} > 0) {
            $s .= $merge->{'multiallelics'} . " variants were removed that would present 3 or more alleles after the dataset merge. ";
        }
    } else {
        $s .= 'Only one input dataset is specified, dataset merge has been skipped. ';
    }
    $s .= '\\\\';
    $s .= $merge->{'info'};


    #    $s .= '\subsection{Missingness}';
    #    $s .= 'Variants have been checked for too high missingness with respect to the call rate. The first missingness test has been performed with a threshold of ' . $miss->{'whole-thres'} . ' on the whole batch collection. ';
    #    $s .= $miss->{'whole-count'} . ' variants with missingness rates exceeding the threshold were found. The second test has been performed on the same dataset with the worst-performing batch excluded. The threshold was ' . $miss->{'perbatch-thres'} . ' and ' . $miss->{'perbatch-count'} . ' variants were found to exceed the threshold.\\\\';
    $s .= '\subsection{Hardy-Weinberg}';
    my $fdr_thres = '10$^{\text{-' . ($hwe->{'threshold'} + 1) . '}}$';

    if($bad->{'final-controls'} == 0) {
        $s .= 'Hardy-Weinberg tests were skipped because no control samples are available.';
    } else {
        #        my $fdr_thres =        #$s .= 'Variants were tested for Hardy-Weinberg-Equilibrium, corrected for false discovery rates with a threshold of ' . $fdr_thres . ' (Benjamini and Hochberg, 1995). ';
        # $s .= $hwe->{'all-count'} . ' variants where rejected from the whole collection, and ' . $hwe->{'worstbatch-count'} . ' variants when the worst-performing batch is was removed.\\\\[1em]';
        $s .= '\begin{minipage}{1\textwidth}\includegraphics[width=0.9\textwidth]{' . sanitize_img($hwe->{'worstbatch-img'}) . '}\end{minipage}\\\\';
        $s .= '\begin{minipage}{1\textwidth}\includegraphics[width=0.9\textwidth]{' . sanitize_img($hwe->{'twoplus-img'}) . '}\end{minipage}';
    }
    $s .= "\\\\\n";

    $s .= '\subsection{Exclusion Summary}';
    my $percentage=0;
    if($bad->{'loaded-variants'} != 0) {
        $percentage = sprintf("%.1f", ($bad->{'loaded-variants'} - $bad->{'variants-after-exclude'}) * 100.0 / $bad->{'loaded-variants'});
    }
    $s .= 'Due to failing missingness and/or HWE tests, ' . ($bad->{'loaded-variants'} - $bad->{'variants-after-exclude'}) . " variants have been removed ($percentage\\,\\%):\\\\[1em]";
    $s .= '\begin{tabular}{lr}\toprule{}';
    $s .= 'Variants in whole collection failing missingness test (threshold ' . $miss->{'whole-thres'} . '): & ' . $miss->{'whole-count'} . '\\\\';
    $s .= 'Variants failing missingness test with worst batch removed (threshold ' . $miss->{'perbatch-thres'} . '): & ' . $miss->{'perbatch-count'} . '\\\\';
    $s .= 'Variants in whole collection failing HWE test (FDR <' . $fdr_thres . '): & ' . $hwe->{'all-count'} . '\\\\';
    $s .= 'Variants failing HWE test with worst batch removed (FDR <' . $fdr_thres . '): & ' . $hwe->{'worstbatch-count'} . '\\\\';
    $s .= 'Variants failing HWE test in 2+ batches (FDR <' . $fdr_thres . '): & ' . $hwe->{'twoplus-count'} . '\\\\';
    $s .= '\midrule{}';
    $s .= 'Total variants to be excluded: & ' . ($miss->{'whole-count'}+$miss->{'perbatch-count'}+$hwe->{'all-count'}+$hwe->{'worstbatch-count'}+$hwe->{'twoplus-count'}) . '\\\\';
    $s .= 'Unique variants excluded: & ' . ($bad->{'loaded-variants'} - $bad->{'variants-after-exclude'}) . '\\\\\bottomrule{}\end{tabular}';

    my $def_pre = $self->definetti_preqc();
    my $def_post = $self->definetti_postqc();
    $s .= '\subsection{DeFinetti Diagrams}';
    $s .= '\begin{minipage}{0.5\textwidth}{\centering \textbf{Pre-QC}\\\\}\includegraphics[width=\textwidth]{' . sanitize_img($def_pre->{'controls'}) . '}\end{minipage}';
    $s .= '\begin{minipage}{0.5\textwidth}{\centering \textbf{Post-QC}\\\\}\includegraphics[width=\textwidth]{' . sanitize_img($def_post->{'controls'}) . '}\end{minipage}\\\\';
    $s .= '\begin{minipage}{0.5\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($def_pre->{'cases'}) . '}\end{minipage}';
    $s .= '\begin{minipage}{0.5\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($def_post->{'cases'}) . '}\end{minipage}\\\\';
    $s .= '\begin{minipage}{0.5\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($def_pre->{'cases/controls'}) . '}\end{minipage}';
    $s .= '\begin{minipage}{0.5\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($def_post->{'cases/controls'}) . '}\end{minipage}\\\\';

    $s .= '\subsection{Phase Summary}';
    $s .= $bad->{'info'};
    $s
}

1


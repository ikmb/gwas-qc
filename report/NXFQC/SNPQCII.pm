package NXFQC::SNPQCII;

use NXFQC::Process;
use NXFQC::PlinkLog;

# From the Perl core distribution
use Data::Dumper;
use File::Basename;

use strict;
our $VERSION = "1.00";

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

# Taken from List::MoreUtils
sub uniq (@) {
    my %seen = ();
    my $k;
    my $seen_undef;
    grep { defined $_ ? not $seen{$k = $_}++ : not $seen_undef++  } @_;
}

sub sanitize_img {
  my $f = shift;
  $f =~ /(.*)(\.[a-zA-Z]+)$/;
  "{$1}$2";
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

sub hftest {
    my $self = shift;
    my $dir_hfprep = $self->{'trace'}->{'hf_test_prepare'};
    my $dir_exclude = $self->{'trace'}->{'generate_hf_excludes'};
    my $basename;
    my $dat = {};

    open my $fh, '<', "$dir_hfprep/.command.sh" or die("$dir_hfprep/.command.sh: $!");
    while(<$fh>) {
        $dat->{'min'} = $1 if /^MIN_BATCH_COUNT=(\d+)/;
    }
    close $fh;

    open my $fh, '<', "$dir_exclude/.command.sh" or die($!);
    while(<$fh>) {

    }
    close $fh;

    $dat->{'excludes-total'} = countlines("$dir_exclude/hf-excludes");
    $dat->{'diagnoses'} = read_file("$dir_exclude/diagnoses");
    $dat;
}


sub det_monomorphics {
    my $self = shift;
    my $dir = $self->{'trace'}->{'det_monomorphics'};
    my $dat = {};

    open my $fh, '<', "$dir/.command.sh" or die($!);
    my $file = (map { /write_monomorphic_file="(\S+)"/ ? "$dir/$1" : () } <$fh>)[0];
    $dat->{'count'} = countlines($file);
    $dat;
}

sub det_diff_missingness {
    my $self = shift;
    my $dir = $self->{'trace'}->{'det_diff_missingness'};
    my $dat = {};

    open my $fh, '<', "$dir/.command.sh" or die($!);
    while(<$fh>) {
        if (/write_file="(\S+)", threshold=(\S+)\)/) {
            $dat->{'count'} = countlines("$dir/$1");
            $dat->{'threshold'} = $2;
            last;
        }
    }
    $dat;
}


sub det_unknown_diagnosis {
    my $self = shift;
    my $dir = $self->{'trace'}->{'det_unknown_diagnosis'};
    my $dat = {};

    open my $fh, '<', "$dir/.command.sh" or die($!);
    my $file = (map { /outfile="(\S+)"/ ? "$dir/$1" : () } <$fh>)[0];
    $dat->{'count'} = countlines($file);
    $dat;
}

sub build_report_chunk {
    my $self = shift;

    my $s = '\section{SNP QC: HF/ANOVA-Testing and Monomorphics}';
    my $hf = $self->hftest();
    $s .= '\subsection{HF/ANOVA-Testing}';
    if($hf->{'diagnoses'} eq '') {
        $s .= '\begin{note}A HF test has not been performed as no contained diagnosis reached the minimum batch count of ' . $hf->{'min'} . '\end{note}';
    } else {
        $s .='\begin{warning}TODO: include graphics, diagnoses and exclusion counts\end{warning}';
    }
    $s;

    my $diff = $self->det_diff_missingness();
    my $diag = $self->det_unknown_diagnosis();
    my $mono = $self->det_monomorphics();

    my $total = $diff->{'count'} + $hf->{'excludes-total'};
    my $unique = 0;

    $s .= '\subsection{Exclusion Summary}';
    $s .= 'Several smaller tests identify more variants to be excluded:\\\\';
    $s .= '\begin{tabular}{lr}\toprule ';
    $s .= 'Differential Missingness: & ' . $diff->{'count'} . '\\\\';
    $s .= 'HF/ANOVA-Test: & ' . $hf->{'excludes-total'} . '\\\\\\midrule{}';
    $s .= 'Total: & ' . $total . '\\\\';
    $s .= 'Unique: & ' . $unique . '\\\\\\bottomrule{}\end{tabular}\\\\';
    $s .= 'Additionally, ' . $mono->{'count'} . ' variants have been flagged as monomorphic. ' . $diag->{'count'} . ' samples have been removed due to unspecified or undeclared diagnoses. ';
    $s;
}
1


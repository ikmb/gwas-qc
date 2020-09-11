package NXFQC::FinalAnalysis;

use NXFQC::Process;
use NXFQC::PlinkLog;
use NXFQC::PlinkInfo;

# From the Perl core distribution
use Cwd;
use Data::Dumper;
#use File::Slurp::Tiny qw(read_file);
use File::Basename;
#use List::MoreUtils qw(uniq);

use strict;
our $VERSION = "1.00";
sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or die("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}

sub sanitize {
  my $s = shift;
  $s =~ s/_/\\_/g;
  $s
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
        $s .= '\begin{minipage}{' . $mboxwidth . '\textwidth}\includegraphics[width=\textwidth]{' . sanitize_img($_) . '}\end{minipage}' . "\n";
    }
    $s .= '\\\\';
    $s;
}



sub pca {
    my $self = shift;
    my $trace = shift;
    my $dir = $self->{'trace'}->{$trace};
    my $dat = {};

    my $batch;
    my $country;

    open my $fh, '<', "$dir/.command.sh" or die($!);
    while(<$fh>) {
        if(/^R --slave --args "(\S+?)"/) {
            $batch = "$dir/$1";
            $country = "$dir/$1.country";
            last;
        }
    }

    $dat->{'batch_2pc'} = "$batch.2PC.png";
    $dat->{'batch_4pc'} = "$batch.4PCpairs.png";
    $dat->{'country_2pc'} = "$country.2PC.png";
    $dat->{'country_4pc'} = "$country.4PCpairs.png";
    $dat;
}

sub pca_1kg {
    my $self = shift;
    my $trace = shift;
    my $dir = $self->{'trace'}->{'pca_plot_1kg_frauke_final'};
    my $cwd = cwd();
    my $basename;
    my $dat = {};

    open my $fh, '<', "$dir/.command.sh" or die($!);
    while(<$fh>) {
        if(/^flashpca2 -d \d+ --bfile "(\S+?)"/) {
            $basename = $1;
            last;
        }
    }

    chdir($dir) or die($!);
    system("convert -density 300 -trim ${basename}_pca.pdf ${basename}_pca.png") unless -e "${basename}_pca.png";
    system("convert -density 300 -trim ${basename}_pca_PCs.pdf ${basename}_pca_PCs.png") unless -e "${basename}_pca_PCs-0.png";
    $dat->{'overview'} = "$dir/${basename}_pca.pdf";
    $dat->{'pairs'} = ();
    @{$dat->{'pairs'}} = map {"$dir/$_"} sort(glob "${basename}_pca_PCs-*");
    chdir($cwd);
    $dat;
}

sub draw_histograms {
    my $self = shift;
    my $dir = $self->{'trace'}->{'draw_final_pca_histograms'};
    my $dat = {};
    my $max_pc = 20;
    my $basename;

    open my $sh, '<', "$dir/.command.sh" or die($!);
    while(<$sh>) {
        if(/R --slave --args "(\S+.pca.evec)"/) {
            $basename ="$dir/$1";
            next;
        }
    }
    $dat->{'img'} = ();
    for(my $i = 1; $i < ($max_pc+1); $i++) {
        if(-e "${basename}.histPC$i.png") {
            push(@{$dat->{'img'}}, "${basename}.histPC$i.png");
        }
    }
    $dat;

}

sub plot_maf {
    my $self = shift;
    my $dir = $self->{'trace'}->{'plot_maf'};
    my $dat = {};

    open my $fh, '<', "$dir/.command.sh" or die($!);
    while(<$fh>) {
        $dat->{'img'} = "$dir/$1.logscale.2.png" if /R --slave --args "(\S+.frq)"/;
    }

    $dat;
}

sub get_plink_info {

    my $self = shift;
    my $dir = $self->{'trace'}->{'prune_final'};
    my $dat = {};


    $dat->{'info'} = plink_table("$dir/info.txt");
    $dat;
}

sub build_report_chunk {
    my $self = shift;

    my $s = '\section{Final Analysis Results}';
    my $info = $self->get_plink_info();
    my $pca = $self->pca('final_pca_con_projection');
    my $pca_a = $self->pca('final_pca_con_projection_atcg');
    my $pca_1kg = $self->pca_1kg();

    $s .= '\subsection{Phase Summary}';
    $s .= $info->{'info'};

    $s .= '\subsection{PCA without Projection}';
    $s .= '\subsubsection{Batch-annotated}';
    $s .= 'Principal Component Analyses were conducted excluding (left) and including (right) AT/CG variants:\\\\';
    $s .= add_images(($pca->{'batch_2pc'}, $pca_a->{'batch_2pc'}));
    $s .= add_images(($pca->{'batch_4pc'}, $pca_a->{'batch_4pc'}));
    $s .= '\subsubsection{Country-annotated}';
    $s .= add_images(($pca->{'country_2pc'}, $pca_a->{'country_2pc'}));
    $s .= add_images(($pca->{'country_4pc'}, $pca_a->{'country_4pc'}));
    $s .= '\subsection{PCA with 1000 Genomes Projection}' . "\n";
    $s .= add_images(($pca_1kg->{'overview'}));
    my $num_pairs = scalar(@{$pca_1kg->{'pairs'}});
    for(my $i = 1; $i < $num_pairs; $i += 2) {
        $s .= add_images(($pca_1kg->{'pairs'}-> [$i-1], $pca_1kg->{'pairs'}->[$i]));
    }
    #    $s .= add_images(($pca_1kg->{'pairs'}->[0], $pca_1kg->{'pairs'}->[1]));
    #$s .= add_images(($pca_1kg->{'pairs'}->[2], $pca_1kg->{'pairs'}->[3]));
    #$s .= add_images(($pca_1kg->{'pairs'}->[4]));
    $s .= '\subsection{PCA Histograms}';

    my $hist = $self->draw_histograms();
    $s .= add_images(@{$hist->{'img'}}[0..1]);
    $s .= add_images(@{$hist->{'img'}}[2..3]);
    $s .= add_images(@{$hist->{'img'}}[4..5]);
    $s .= '\subsection{MAF Plot}';
    my $maf = $self->plot_maf();
    $s .= add_images(($maf->{'img'}));
    print Dumper($s);
    $s;
}
1


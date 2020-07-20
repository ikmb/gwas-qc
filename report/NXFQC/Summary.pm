package NXFQC::Summary;

use strict;
use warnings;


use NXFQC::Process;
use NXFQC::PlinkLog;
use NXFQC::PlinkInfo;

# From the Perl core distribution
use Data::Dumper;
#use File::Slurp::Tiny qw(read_file);
use File::Basename;

use Carp;

$Carp::Verbose=1;

our $VERSION = "1.00";

use strict;
use warnings;
sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or croak("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}

sub sanitize {
    my $s = shift;
    $s =~ s/_/\\_/g;
    $s
}

sub new {
    my $class = shift;
    my $summary = shift;

    my $self = {};

    bless $self, $class;

    $self->{'summary_file'} = $summary;

    return $self;
}

sub tool_versions {
    my $p = shift;
    my $f = $p->{'summary_file'};
    my @dat;
    open my $fh, '<', $f or die($!);
    while(<$fh>) {
        chomp;
        my @parts = split /;/;
        push(@dat, \@parts);
    }

    return \@dat;
}

sub build_report_chunk {
    my $self = shift;

    my $s = '\appendix\section{Software Details}';
    my $tools = $self->tool_versions();

    $s .= '\begin{tabularx}{\textwidth}{lX}';
    $s .= '\textbf{Component} & \textbf{Self-reported Version}\\\\\midrule' . "\n";
    foreach (@$tools) {
        my @v = @$_;
        if($v[0] eq '') {
            $s .= '\midrule{}';
        } else {
            $s .= sanitize($v[0]).'&' . sanitize($v[1]) . "\\\\\n";

        }
    }
    $s .= '\end{tabularx}';

    $s
}

1


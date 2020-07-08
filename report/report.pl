#!/usr/bin/env perl

use NXFQC::Process;
use NXFQC::Rs;
use NXFQC::SNPQCI;
use NXFQC::SampleQC;
use NXFQC::SNPQCII;
use NXFQC::FinalAnalysis;

use Data::Dumper;
# use File::Slurp::Tiny qw(read_file);

use strict;
use warnings;

# Command-line Args

sub Usage {
    print STDERR "Usage: $0 <NXF-BASEDIR> <preamble.tex> <Rs-trace0> [<Rs-trace1> [<Rs-trace2> ...]] <SNPQCI-trace> <SampleQC-trace> <SNPQCII-trace> <Final-trace>\n";
    print STDERR "Example: $0 \$NXF_WORK preamble.tex Rs-BATCH1.txt Rs-BATCH2.txt SNPQCI.txt SampleQC.txt SNPQCII.txt Final.txt\n";
    print STDERR "\tPrints a LaTeX Report to STDOUT.\n";
    exit 1;
}

Usage if @ARGV < 6;

sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or die("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}

# Expand a given NXF Unit path (i.e. '8e/1234567' and a basename to a valid path within the NXF_WORK directory
sub expand_path {
    my $base = shift;
    my $hash = shift;

    if ($hash =~ m|(..)/(.*)|) {
        my $hash_base = $1;
        my $hash_remainder = $2;
        opendir(my $dh, "$base/$hash_base/") or die("$!: $base/$hash_base/");
        my @results = grep { /^$hash_remainder.*/ } readdir($dh);
        closedir $dh;

        if (@results > 0) {
            return $base . "/" . $hash_base . "/" . $results[0];
        } else {
            return "/dev/null";
        }
    }
}

# Creates a hash with $process_name => $workdir k/v pairs
sub parse_trace {
    my $basedir = shift;
    my $filename = shift;

    my $result = {};
    open my $fh, '<', $filename or die("Could not open trace $filename");
    <$fh>; # Skip header line
    while(<$fh>) {
        chomp;
        my ($process, $procname, $tag, $workdir) = split /;/;
        $result->{$process} = expand_path($basedir, $workdir);
    }
    $result
}

my $nxf_work_basedir = shift @ARGV;
my $preamble = shift @ARGV;
my $final_trace_file = pop @ARGV;
my $snpqcii_trace_file = pop @ARGV;
my $sampleqc_trace_file = pop @ARGV;
my $snpqci_trace_file = pop @ARGV;
my @rs_trace_files = @ARGV;

my @rs_traces = ();
foreach(@rs_trace_files) {
    print "Processing 'Rs' trace $_\n";
    my $trace = parse_trace($nxf_work_basedir, $_);
    print Dumper($trace);
    push(@rs_traces, new NXFQC::Rs($trace));
}
print "Processing 'SNP QC I' trace $snpqci_trace_file\n";
my $snpqci = new NXFQC::SNPQCI(parse_trace($nxf_work_basedir, $snpqci_trace_file));
#print Dumper($snpqci);
print "Processing 'SampleQC' trace $sampleqc_trace_file\n";
my $sampleqci = new NXFQC::SampleQC(parse_trace($nxf_work_basedir, $sampleqc_trace_file));
#print Dumper($sampleqci);
print "Processing 'SNP QC II' trace $snpqcii_trace_file\n";
my $snpqcii = new NXFQC::SNPQCII(parse_trace($nxf_work_basedir, $snpqcii_trace_file));
#print Dumper($snpqcii);
print "Processing 'FinalAnalysis' trace $final_trace_file\n";
my $final = new NXFQC::FinalAnalysis(parse_trace($nxf_work_basedir, $final_trace_file));

open my $texfh, '>report.tex' or die($!);
print $texfh read_file($preamble);
if($ENV{'WARN_SNPQC'} eq '1' || $ENV{'WARN_RELATED'} eq 'true' || $ENV{'WARN_SAMPLEQC'} eq '1') {
    print $texfh '\begin{warning}\begin{itemize}';
    print $texfh '\item \texttt{--skip\_snpqc} set, not removing SNPs during QC' if $ENV{'WARN_SNPQC'} eq '1';
    print $texfh '\item \texttt{--skip\_sampleqc} set, not removing samples during QC' if $ENV{'WARN_SAMPLEQC'} eq '1';
    print $texfh '\item \texttt{--keep\_related} set, not removing relatives during QC' if $ENV{'WARN_RELATED'} eq 'true';
    print $texfh '\end{itemize}\end{warning}';
}
foreach(@rs_traces) {
    print $texfh $_->build_report_chunk();
}
print $texfh $snpqci->build_report_chunk();
print $texfh $sampleqci->build_report_chunk();
print $texfh $snpqcii->build_report_chunk();
# print Dumper($snpqcii->build_report_chunk());
print $texfh $final->build_report_chunk();

print $texfh "\\end{document}";
close $texfh;


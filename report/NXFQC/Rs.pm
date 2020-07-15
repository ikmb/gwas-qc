package NXFQC::Rs;

use strict;
use warnings;


use NXFQC::Process;
use NXFQC::PlinkLog;

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

sub batch_statistics {
    my $self = shift;
    my $dir = $self->{'trace'}{'batch_statistics'};
    my $dat = {};

    my $statsfile;
    open my $basefh, '<', "$dir/.command.sh" or die($!);
    while(<$basefh>) {
        chomp;
        if(/outfilename.=."(.*)";$/) {
            $statsfile = "$dir/$1";
            last;
        }
    }
    close $basefh;
    print "Statsfile: $statsfile\n";
    open my $statsfh, '<', $statsfile or die($!);
    my @cols;
    while(<$statsfh>) {
        chomp;
        @cols = split /\s+/;
        $dat->{$cols[0]} = {};
        for(my $i=1; $i < @cols; $i+= 2) {
            $dat->{$cols[0]}->{$cols[$i]} = $cols[$i+1];
        }
    }
    $dat
}

sub check_chip_type {
    my $self = shift;
    my $dir = $self->{'trace'}{'check_chip_type'};
    my $basename;
    open my $basefh, '<', "$dir/.command.sh" or croak($!);
    while(<$basefh>) {

        $basename = $1 if /output\s(\S+)\.chip_detect\.log\s/;

    }

    my $dat = {};

    # Parse AT/CG file
    open my $atcg, '<', "$dir/${basename}.flag_atcg" or croak("$!: $dir/${basename}.flag_atcg");
    $dat->{'at'} = 0;
    $dat->{'cg'} = 0;
    while(<$atcg>) {
        chomp;
        $dat->{'at'} += 1 if /AT$/ or /TA$/;
        $dat->{'cg'} += 1 if /CG$/ or /GC$/;
    }



    # Parse chip detection log
    open my $chipdetect, '<', "$dir/${basename}.chip_detect.log" or die("$!: $dir/${basename}.chip_detect.log");
    <$chipdetect>; # skip header

    $dat->{'chips'} = [];
    my $build = 0;

    while(<$chipdetect>) {
        my @fields = split /\s+/, $_;
        if (@{$dat->{'chips'}} == 0 or ($dat->{'name_match_rate'} eq $fields[1])) {
            $dat->{'name_match_rate'} = $fields[1];
            $dat->{'pos_match_rate'} = $fields[2];
            $dat->{'original_match_rate'} = $fields[3];
            $dat->{'plus_match_rate'} = $fields[4];
            $dat->{'atcg_match_rate'} = $fields[5];
            $build = $1 if ($build == 0 and $fields[0] =~ m/-b(\d\d).(Ilmn|Source)/);
            $fields[0] =~ s/-b(\d\d)\.(Ilmn|Source)\.strand//g;
            if($dat->{'name_match_rate'} ne '0') {
                push(@{$dat->{'chips'}}, $fields[0]);

            }
        } else {
            last;
        }
    }

    $dat->{'tag'} = $basename;
    $dat->{'build'} = $build;

    # Parse sourcedata
    $dat->{'source'} = {};
    open my $sourcedata, '<', "$dir/sourcedata.txt" or die("$!: $dir/sourcedata.txt");
    while(<$sourcedata>) {
        my @fields = split /;/, $_;
        $dat->{'source'}->{$fields[0]} = $fields[1];
    }

    $dat
}

sub lift_genome_build {
    my $self = shift;
    my $dir = $self->{'trace'}{'lift_genome_build'};

    my $target_name;
    my $strand_name;

    my $dat = {};

    # Parse file name targets
    open my $fh, "<$dir/.command.sh" or die("$1: ${dir}/.command.sh");
    while(<$fh>) {
        if(/TARGETNAME="(.*)"$/) { $target_name = $1; }
        if(/STRAND_FILE="(.*)"$/) { $strand_name = $1; }
    }
    close $fh;

    # Extract name portion from strand
    if ( -e $strand_name) {
        my($name, $path, $suffix) = fileparse($strand_name, qw(.strand));
        $dat->{'strand-file-name'} = $name;
        $dat->{'lift-result'} = parse_plink("$dir/${target_name}.log");
    } else {
        $dat->{'strand-file-name'} = '';
    }
    return $dat;
}

sub fix_par {
    my $self = shift;
    my $dir = $self->{'trace'}{'fix_par'};

    my $merged = 0;
    my $fixed = 0;
    my $hh = 0;

    open my $fh, "<$dir/merged-with-missing.log" or carp("$1: ${dir}/merged-with-missing.log");
    while(<$fh>) {
        $merged = $1 if /merge-x: (\d+) chromosome codes changed/;
    }

    open my $fh2, "<$dir/split.log" or carp("$1: ${dir}/fixed.log");
    while(<$fh2>) {
        $fixed = $1 if /split-x: (\d+) chromosome codes changed/;
    }

    open my $fh3, "<$dir/fixed.log" or carp("$1: ${dir}/fixed.log");
    while(<$fh3>) {
        $hh = $1 if /Warning: (\d+) het. haploid/;
    }

    my $dat = {};
    $dat->{'hh'} = $hh;
    $dat->{'merged'} = $merged;
    $dat->{'fixed'} = $fixed;
    return $dat;
}

sub normalize_variant_names {
    my $self = shift;
    my $dir = $self->{'trace'}{'normalize_variant_names'};

    my $l = {};

    open my $fh, "<$dir/.command.log" or die("$1: $dir/.command.log");
    while(<$fh>) {
        $l->{'matches'} = $1 if /^Exact matches: (\d+)$/;
        $l->{'new-ids'} = $1 if /^Variants with new Rs ID: (\d+)$/;
        $l->{'no-id'} = $1 if /^Variants without known Rs IDs: (\d+)$/;
        $l->{'flip-matches'} = $1 if /^Flipped matches: (\d+)$/;
        $l->{'flip-new-ids'} = $1 if /^Flipped variants with new Rs ID: (\d+)$/;
        $l->{'flip-no-id'} = $1 if /^Flipped variants without known Rs IDs: (\d+)$/;
        $l->{'unknown'} = $1 if /^Unknown variants: (\d+)$/;
        $l->{'indels'} = $1 if /^\.\.\.of which are indels: (\d+)$/;
    }

    return $l;
}

sub plink_flip {
    my $self = shift;
    my $dir = $self->{'trace'}{'plink_flip'};

    my $l = {};

    # Find target files
    open my $fh, "<$dir/.command.sh" or die("$!: $dir/.command.sh");
    while(<$fh>) {
        $l->{'flip'} = parse_plink("$dir/$1.log") if /^plink --bfile dedup --flip.*--out "(.*)" --allow-no-sex$/;
    }
    $l->{'dedup'} = parse_plink("$dir/dedup.log");
    return $l;
}


sub new {
    my $class = shift;
    my $trace = shift;

    my $self = {};

    bless $self, $class;

    $self->{'trace'} = $trace;

    return $self;
}

sub compile_stats {
    my $p = shift;
    my $s = '';
    my $count = 0;
    foreach my $key (sort keys %$p) {
        $count += $p->{$key};
    }

    my @findings;
    foreach my $key (sort keys %$p) {
        push(@findings, sanitize($key) . " (" . sanitize($p->{$key}) . ", " . sprintf("%.1f\\%%)", $p->{$key}*100.0/$count));
    }
    return join(', ', @findings);
}

sub build_report_chunk {
    my $self = shift;

    my $chip = $self->check_chip_type();

    my $s = "\\section{Annotation Check for " . sanitize($chip->{'tag'}) . "}";
    $s .= '\subsection{Input Files}';
    $s .= '\begin{tabular}{l@{\hskip 1cm}l}';
    $s .= '\textbf{Filename} & \textbf{Last modified}\toprule' . "\\\\\n";
    foreach my $f (keys $chip->{'source'}) {
        $s .= sanitize($f) . ' & ' . sanitize($chip->{'source'}->{$f}) . "\\\\\n";
    }
    $s .= '\end{tabular}';
    ############################## Batch Overview Report ###################################
    $s .= '\subsection{Annotation Overview}';
    my $stats = $self->batch_statistics();
    $s .= "\n" . '\\begin{tabularx}{\textwidth}{lX}\\toprule{}' . "%\n";
    $s .= "Genders: & " . compile_stats($stats->{'Genders'}) . "\\\\\n";
    $s .= "Phenotypes: & " . compile_stats($stats->{'Phenotypes'}) . "\\\\\n";
    $s .= "Batches: & " . compile_stats($stats->{'Batches'}) . "\\\\\n";
    $s .= "Diagnoses: & " . compile_stats($stats->{'Diagnoses'}) . "\\\\\n";
    $s .= "Countries: & " . compile_stats($stats->{'Countries'}) . "\\\\\\bottomrule{}\n";
    $s .= "\\end{tabularx}\\\\\n";

    ############################## Chip type detection report ##############################

    if(@{$chip->{chips}} > 1) {
        $s .= "The detected chip type is one of " . sanitize(join(', ', map { "\\texttt{$_}"} @{$chip->{'chips'}} ));
    } elsif(@{$chip->{chips}} == 1) {
        $s .= "The detected chip type is " . sanitize($chip->{'chips'}->[0]);
    }

    my $strand_rate;

    if(@{$chip->{chips}} > 0) {
        $s .= ". The genome build ist most likely GRCh$chip->{build}. ";
        $s .= "Test characteristics for chip type and genome build detection:\\\\[1em]\n";
        $s .= '\begin{tabular}{lr}\\toprule{}';
        $s .= 'Variants matched by name: & ' . sprintf("%.2f\\,\\%%\\\\", $chip->{'name_match_rate'} * 100.0);
        $s .= 'Variants matched by name and position: & ' . sprintf("%.2f\\,\\%%\\\\", $chip->{'pos_match_rate'} * 100.0);
        $s .= '\midrule{}Variants matched on the original annotation (without AT/CG): & ' . sprintf("%.2f\\,\\%%\\\\", $chip->{'original_match_rate'} * 100.0);
        $s .= 'Variants matched on the "+" annotation (without AT/CG): & ' . sprintf("%.2f\\,\\%%\\\\", $chip->{'plus_match_rate'} * 100.0);
        $s .= 'Variants ratio of A/T or C/G alleles: & ' . sprintf("%.2f\\,\\%%\\\\", $chip->{'atcg_match_rate'} * 100.0);
        $s .= '\midrule{}Variants with A/T alleles: & ' . $chip->{'at'} . '\\\\';
        $s .= 'Variants with C/G alleles: & ' . $chip->{'cg'} . '\\\\\\bottomrule{}' . "\n";
        $s .= '\end{tabular}' . "\n";

        my $rate = $chip->{'original_match_rate'};
        if($chip->{'plus_match_rate'} > $chip->{'original_match_rate'}) {
            $s .= '\begin{note}';
            $s .= 'The variant definitions fit better to the + strand than to the original annotation file.\end{note}';
            $rate = $chip->{'plus_match_rate'};
        }

        if($chip->{'name_match_rate'} < 0.85) {
            $s .= '\begin{warning}';
            $s .= 'Only ' . sprintf("%.2f", $chip->{'name_match_rate'} * 100.0) .'\,\% of defined variant names could be matched against the known chip database (i.e. original chip annotation files from companies such as Illumina or Affymetrix). It is likely that either a custom chip build has been used or the dataset is incomplete or the original SNP identifiers have been manipulated.';
            $s .= '\end{warning}';
        }

        if($rate < 0.95) {
            $s .= '\begin{warning}';
            $s .= 'Only ' . sprintf("%.2f", $rate * 100.0) . '\,\% of the variants could be matched against the chip strand information. It is likely that some variants have been flipped to a different strand and are inconsistent with the company-supplied chip manifest.\end{warning}';
        }
    } else {
        $s .= '\begin{note}The chip definition database does not contain matches for this dataset. This mostly happens for incomplete databases, malformed or modified variant names and (currently unsupported) Affymetrix chips.\end{note}';
    }

    ############################## Genome Build Conversion (lift-over) ##############################

    $s .= '\subsection{Genome Build Conversion}';
    my $lift = $self->lift_genome_build();

    if($lift->{'strand-file-name'} eq '') {
        $s .= "No genome build lift-over has been requested. ";
    } else {
        $lift->{'strand-file-name'} =~ /b(\d+)$/;
        my $new_build = $1;
        $s .= "A lift-over from build $chip->{build} to $new_build has been requested using the strand information from \\texttt{" . sanitize($lift->{'strand-file-name'}) . "}. ";
        $s .= "\\begin{note}Although the genome build stays the same, the lift-over is still being performed as requested.\\end{note} " if $new_build == $chip->{'build'};
        my $target_chip = $lift->{'strand-file-name'} =~ s/\-b(\d+)$//r;
        my @source_chips = map { s/\-b(\d+)$//r } @{$chip->{'chips'}};
        if( grep($target_chip, @source_chips) == 0) {
            $s .= "\\begin{warning}The target chip does not match the detected chip type. The lift-over has been performed, but you should double-check the results or consider setting a different target chip for lift-over.\\end{warning}";
        }

        my $orig = $lift->{'lift-result'}->{'loaded-variants'};
        my $afterlift = $lift->{'lift-result'}->{'variants-after-extract'};
        my $percent_removed = ($orig-$afterlift) * 100.0 / ($orig);
        $s .= "Of the original $orig variants, $afterlift could be updated. " . ($orig-$afterlift) . " variants were removed (" . sprintf("%.2f", $percent_removed ) . "\\,\\%). ";
        if($percent_removed > 20.0) {
            $s .= "\\begin{warning}" . sprintf("%.3f", $percent_removed) . "\\,\\% of the original variants were removed during the lift-over process. This means that either a wrong target chip type has been chosen or the original variant names have not been retained. The lift-over process requires the input dataset to match the original annotation database.\\end{warning}";
        }
    }

    $s .= '\subsection{Chromosome X/Y PAR Check}';
    my $par = $self->fix_par();
    if($par->{'merged'} == 0 && $par->{'fixed'} == 0) {
        $s .= 'No PAR regions are available in the input dataset.';
    } elsif($par->{'merged'} == $par->{'fixed'}) {
        $s .= "Found $par->{merged} variants in PAR regions.";
    } elsif($par->{'merged'} == 0 && $par->{'fixed'} != 0) {
        $s .= "No PAR regions were defined in the input dataset but $par->{fixed} variants were found that would have belonged in PAR1 or PAR2. They have been moved to their respective regions.";
    } else {
        $s .= '\begin{warning}Inconsistencies were found in the PAR region contents. '.$par->{'merged'}.' variants were previously defined in PAR regions but '.$par->{'fixed'}.' variants actually belonged there. This is fixed now.\end{warning}';
    }
    if($par->{'hh'} != 0) {
        $s .= '\begin{warning}'.$par->{'hh'}.' heterozygous haploid genotype calls were found. These will be erased.\end{warning}';
    } else {
        $s .= 'No heterozygous haploid calls were found.';
    }
    ############################## Flipping and Dedup ##############################

    $s .= '\subsection{Variant ID and Strand Correction}';
    my $flip = $self->plink_flip();
    my $norm = $self->normalize_variant_names();

    $s .= 'The variant names have been normalized from chip-specific IDs to the NCBI dbSNP IDs using a database compiled from the HRC1.1 reference panel, the UK10k panel, the UK10k + 1000G Phase 3 merge and the dbSNP150 database, in that order and priority.\\\\[1em]';
    $s .= '\begin{tabular}{lr}\toprule{}';
    $s .= 'Variants with correct names: & ' . $norm->{'matches'} . '\\\\';
    $s .= 'Variants having their name replaced: & ' . $norm->{'new-ids'} . '\\\\';
    $s .= 'Known variants with ambiguous naming: & ' . $norm->{'no-id'} . '\\\\\midrule{}';
    $s .= 'Variants with correct names, rev strand: & ' . $norm->{'flip-matches'} . '\\\\';
    $s .= 'Variants having their name replaced, rev strand: & ' . $norm->{'flip-new-ids'} . '\\\\';
    $s .= 'Known variants with ambiguous naming, rev strand: & ' . $norm->{'flip-no-id'} . '\\\\\midrule{}';
    $s .= 'Variants not found in database: & ' . $norm->{'unknown'} . '\\\\';
    $s .= '\ldots{} of which are indels and non-ATCG variants: & ' . $norm->{'indels'} . '\\\\\bottomrule{}';
    $s .= '\end{tabular}\\\\[1em]';

    $s .= 'After name normalization, ' . ($flip->{'dedup'}->{'loaded-variants'} - $flip->{'dedup'}->{'variants-after-exclude'}) . ' duplicate variants were removed. ';
    $s .= $flip->{'flip'}->{'snps-flipped'} . ' variants were flipped to the plus strand. ';
    $flip = $flip->{'flip'};

    $s .= 'The batch dataset now consists of '
       . $flip->{'loaded-variants'} . ' variants and '
       . ($flip->{'loaded-phenotypes'}) .' samples ('
       . $flip->{'final-cases'} . ' cases, ' . $flip->{'final-controls'} . ' controls, ' . ($flip->{'loaded-phenotypes'} - ($flip->{'final-cases'} + $flip->{'final-controls'})) . ' unknown). ';


    $s
}

1


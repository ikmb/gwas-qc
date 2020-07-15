package NXFQC::PlinkInfo;

use strict;
use warnings;

require Exporter;

our @ISA = qw(Exporter);
our @EXPORT = qw(plink_table);


our $VERSION = "1.00";

sub plink_table {
    my $log = shift;

    my %d;

    open my $fh, '<', $log or die("$!: $log");
    while(<$fh>) {
        chomp;
        my @p = split /\s+/, $_;
        $d{$p[0]} = $p[1];
    }
    close $fh;

    my $s = '';

    $s .= '\begin{tabular}{lrrr|r}\\\\\toprule{}' . "%\n";
    $s .= '& \textbf{Male} & \textbf{Female} & \textbf{Unknown} & \textbf{Total}\\\\\midrule{}' . "%\n";
    $s .= '\textbf{Case} & ' . $d{'male-cases'} . '&' . $d{'female-cases'} . '&' . $d{'unk-cases'} . '&' . $d{'cases'} . '\\\\' . "%\n";
    $s .= '\textbf{Control} & ' . $d{'male-controls'} . '&' . $d{'female-controls'} . '&' . $d{'unk-controls'} . '&' . $d{'controls'} . '\\\\' . "%\n";
    $s .= '\textbf{Unknown} & ' . $d{'male-unk'} . '&' . $d{'female-unk'} . '&' . $d{'unk-unk'} . '&' . $d{'unknown-pheno'} . '\\\\\midrule{}' . "%\n";
    $s .= '\textbf{Total} & ' . $d{'males'} . '&' . $d{'females'} . '&' . $d{'unknown-sex'} . '&' . $d{'samples'} . '\\\\\bottomrule{}' . "%\n";
    $s .= '\end{tabular}\\\\The dataset contains ' . $d{'variants'} . ' variants. ';
    return $s;
}

1

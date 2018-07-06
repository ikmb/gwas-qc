package Report::Rs;

use File::Slurp;

our $VERSION = '1.00';

sub plink_flip {
    my $workdir = shift;
    my $tag = shift;

    my $s = "\\section{$tag}";
    my $prefix = read_file("$workdir/input.1");
    my $logfile = $workdir . '/' . $prefix . '_flipped.log';

    open my $log, '<', $logfile or die("Could not open $logfile: " . $!);

    my $variants = 0;
    my $samples = 0;
    my $flipped = 0;

    while(<$log>) {
        chomp;
        if (/^(\d+) variants and (\d+) people pass filters and QC/) {
            $variants = $1;
            $samples = $2;
            next;
        }

        if (/^--flip: (\d+) SNPs flipped/) {
            $flipped = $1;
        }
    }

    close $log;

    $s . "The $tag dataset initially consisted of $samples samples and $variants variants, of which $flipped had their strand flipped to the \\texttt{+}-strand.\\ "
}

sub plink_exclude {
    my $workdir = shift;
    my $tag = shift;

    my $logfile = $workdir . '/' .  $tag . '_rs.log';
    open my $log, '<', $logfile or die("Could not open $logfile: " . $!);

    my $variants = 0;
    my $remaining = 0;
    while(<$log>) {
        chomp;
        if (/^(\d+) variants loaded from/) {
            $variants = $1;
            next;
        }

        if (/^--exclude: (\d+) variants remaining/) {
            $remaining = $1;
        }
    }
    close $log;

    if($variants == $remaining) {
        "No duplicates, N/Ns or chip-specific excludes where removed."
    } else {
        ($variants - $remaining) . " variants where removed because they were duplicates, had N/N alleles or were specified in the chip-specific exclude list. This leaves the $tag dataset with $remaining variants after the first QC stage.\\ "
    }

}

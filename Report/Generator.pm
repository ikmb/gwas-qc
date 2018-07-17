package Report::Generator;

use strict;
use warnings;

use Data::Dumper;

our $VERSION = '1.00';

sub new {
    my ($class, %args) = @_;

    my $self = bless ({}, $class);

    my $logname = `whoami`;
    my $ent = `getent passwd $logname`;
    my @fields = split /:/, $ent;
    $fields[4] =~ s/,+$//; # GECOS field contains tailing commas on weirdly configured systems (i.e. Ubuntu)
    my @gecos = split /,/, $fields[4];
    
    $self->{author} = $args{author} // $gecos[0]; # `getent passwd $logname | cut -d: -f5`;
    $self->{title}  = $args{title} // "No Title";
    $self->{parts} = [];
    $self
}

sub append {
    my $self = shift;
    my $arg = shift;
    push @{$self->{parts}}, $arg;
}


sub add_header {
    my $self = shift;

    my $preamble = <<'HEADER_PREAMBLE'
\documentclass{scrartcl}
\usepackage{lmodern}
\usepackage{graphicx}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{float}
\usepackage[scale=0.85]{geometry}
HEADER_PREAMBLE
        ;
    my $content = '\author{' . $self->{author} . '}'
    . '\title{' . $self->{title} . '}';


   $preamble . $content . '\begin{document}\maketitle\tableofcontents' . "\n"
}

sub add_footer {
    "\\end{document}\n"
}

sub to_string {
    my $self = shift;
    my $s = $self->add_header();
    $s .= join '', @{$self->{parts}};
    $s .= $self->add_footer();
    $s
}


1;

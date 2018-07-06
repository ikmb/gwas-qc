package Report::Generator;

use strict;
use warnings;

use Data::Dumper;

our $VERSION = '1.00';

sub new {
    my ($class, %args) = @_;

    my $self = bless ({}, $class);

    my $logname = `whoami`;
    $self->{author} = $args{author} // $logname; # `getent passwd $logname | cut -d: -f5`;
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
HEADER_PREAMBLE
        ;
    my $content = '\author{' . $self->{author} . '}'
    . '\title{' . $self->{title} . '}';


   $preamble . $content . '\begin{document}\maketitle' . "\n"
}

sub add_footer {
    "\\end{document}\n"
}

sub to_string {
    my $self = shift;
    my $s = $self->add_header();
    for (@{$self->{parts}}) {
        $s .= $_;
    }

    $s .= $self->add_footer();
    $s
}


1;

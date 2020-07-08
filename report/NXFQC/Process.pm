package NXFQC::Process;

# use File::Slurp::Tiny qw/ read_file /;
sub read_file {
    my $filename = shift;
    open my $fh, '<', $filename or die("Cannot open $filename: $!");
    read $fh, my $file_content, -s $fh;
    return $file_content;
}


our $VERSION = '1.00';

sub new {
    my $class = shift;

    my ($working_dir, $process_name) = @_;

    die("Could not access $working_dir for process $process_name.") unless -d $working_dir;

    # Add a trailing slash, if necessary
    $working_dir =~ s!/*$!/!;

    my $self = {
        'working_dir' => $working_dir,
        'process_name' => $process_name,
            'staged_files' => (),
            'generated_files' => ()
    };

    bless $self, $class;
    return $self;
}

sub get_staged {
    my $self = shift;

    if(@{$self->staged_files} == 0) {
        open my $fh, "$self->{working_dir}.command.run", '<' or die("Could not open $self->{working_dir}.command.run");
        while(<$fh>) {
            if( /ln -s \S+ (\S+)$/) {
                push(@{$self->staged_files}, $1);
            }
        }
    }

    $self->staged_files
}

sub get_generated {
    my $self = shift;

    opendir(DIR, $self->{'working_dir'}) or die("Could not open working directory $self->{working_dir}");
    while(my $file = readdir(DIR)) {
        # add only files that are plain files and not symbolic links and not hidden
        push(@{$self->generated_files}, $file) if -f $file && !-l $file && $file[0] != '.';
    }
}

1;

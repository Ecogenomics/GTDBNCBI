#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use File::Spec;
use Carp;

my $options = checkParams();

print STDERR Dumper($options);

my %info_dict;
opendir(my $info_dh, $options->{'info_dir'}) or croak("Can't open info dir");
while (my $entry = readdir($info_dh)) {
    my $full_path = File::Spec->join($options->{'info_dir'}, $entry);
    if (-f $full_path && $full_path =~ /(TIGR([0-9]*))\./) {
        my ($name, $accession, $description) = getInfoFromInfoFile($full_path);
        if ($accession ne $1) {
            croak "Accession in info file doesn't match name of info file: ($accession) $full_path";
        }
        $info_dict{$accession} = {'name' => $name, 'description' => $description};
    }
}

opendir(my $hmm_dh, $options->{'hmm_dir'}) or die;
while (my $entry = readdir($hmm_dh)) {
    my $full_path = File::Spec->join($options->{'hmm_dir'}, $entry);
    if (-f $full_path && $full_path =~ /(TIGR([0-9]*))\./) {
        my $accession = $1;
        
        if (! defined($info_dict{$accession})) {
            croak "No info for HMM file under accession: ($accession) $full_path";
        }
        
        print join("\t", (
            $full_path,
            $info_dict{$accession}->{"name"},
            $info_dict{$accession}->{"description"},
            "TIGRFAM",
            $accession
        )), "\n";
    }
}


sub getInfoFromInfoFile {
    my $filepath = shift;
    
    my ($name, $accession, $description);
    open(my $fh, $filepath) or croak("Can't open file: $filepath");
    
    while (my $line = <$fh>) {
        chomp $line;
        
        # Ignore blank lines
        if (! $line) {
            next;
        }
        if ($line !~ /^([A-Z]{2})(\ {2})(.*)$/) {
            #carp "Unexpected line format in the info file. $line";
        }
        
        my $prefix = $1;
        if ($prefix eq "ID") {
            $name = $3;
        } elsif ($prefix eq "AC") {
            $accession = $3;
        } elsif ($prefix eq "DE") {
            $description = $3;
        }
    }
    
    return ($name, $accession, $description);
}

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("help+", "hmm_dir|d:s", "info_dir|i:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'hmm_dir'} ) { printParamError ("We need a tigrfam hmm directory!"); }
    if(!exists $options{'info_dir'} ) { printParamError ("We need a column number!"); }
    
    return \%options;
}


sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print STDERR "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}



__DATA__

=head1 NAME

tigrfam_to_batchfile.pl

=head1 COPYRIGHT


=head1 DESCRIPTION

=head1 SYNOPSIS

Usage: tigrfam_to_batchfile.pl -d <HMM_DIR> -i <INFO_DIR>

    e.g. tigrfam_to_batchfile.pl -d TIGRFAMs_15.0_HMM/individual_hmms/ -i TIGRFAMs_15.0_INFO/

=cut

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

opendir(my $hmm_dh, $options->{'hmm_dir'}) or die;
while (my $entry = readdir($hmm_dh)) {
    my $full_path = File::Spec->join($options->{'hmm_dir'}, $entry);
    if (-f $full_path && $full_path =~ /(PF(.*))(\.hmm)/) {
        
        my ($name, $accession, $description) = getInfoFromHMMFile($full_path);
        
        print join("\t", (
            $full_path,
            $name,
            $description,
            "PFAM",
            $accession
        )), "\n";
    }
}


sub getInfoFromHMMFile {
    my $filepath = shift;
    
    my ($name, $accession, $description);
    open(my $fh, $filepath) or croak("Can't open file: $filepath");
    
    while (my $line = <$fh>) {
        chomp $line;
        
        # Ignore blank lines
        if (! $line) {
            next;
        }
        if ($line !~ /^([A-Z]+)(\ +)([^\ ].*)$/) {
            next;
            #carp "Unexpected line format in the info file. $line";
        }
        
        my $prefix = $1;
        if ($prefix eq "NAME") {
            $name = $3;
        } elsif ($prefix eq "ACC") {
            $accession = $3;
        } elsif ($prefix eq "DESC") {
            $description = $3;
        }
    }
    
    return ($name, $accession, $description);
}

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("help+", "hmm_dir|d:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'hmm_dir'} ) { printParamError ("We need a tigrfam hmm directory!"); }
    
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

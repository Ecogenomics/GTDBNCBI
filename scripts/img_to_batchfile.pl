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


my %img_id_to_metadata;

open(my $metadata_fh, $options->{'metadata_tsv'}) or croak("Can't open file: " . $options->{'metadata_tsv'});
my %header_to_col;
my $header_line = <$metadata_fh>;
chomp $header_line;
my %map_header_to_col = %{create_headers_to_col_pos_map($header_line)};

if (! defined($map_header_to_col{'taxon_oid'})) {
   croak("IMG metadata file doesn't have a 'taxon_oid' column.");
}

if (! defined($map_header_to_col{'Genome Name / Sample Name'})) {
   croak("IMG metadata file doesn't have a 'Genome Name / Sample Name' column.");
}

while (my $line = <$metadata_fh>) {
    chomp $line;
    my @splitline = split /\t/, $line;
    my $metadata = {
        'name' => $splitline[$map_header_to_col{'Genome Name / Sample Name'}]
    };
    my $img_id = $splitline[$map_header_to_col{'taxon_oid'}];
    $img_id_to_metadata{$img_id} = $metadata;
};

close($metadata_fh);

open(my $img_path_fh, $options->{'genome_path_listfile'}) or croak("Can't open file: " . $options->{'genome_path_listfile'});
while (my $line = <$img_path_fh>) {
    chomp $line;
    if (! -e $line) {
        print STDERR "File: $line doesn't exist, skipping....\n";
        next;
    }

    if ($line !~ /\/([0-9]*).fna$/) {
        print STDERR "Filepath: $line doesn't conform to expected IMG filename ([0-9]*.fna). Skipping....\n";
        next;
    }
    
    my $img_id = $1;
    if (! defined ($img_id_to_metadata{$img_id})) {
        print STDERR "IMG id $img_id not in metadata file. Filepath: $line.  Skipping....\n";
        next;
    }   

    my $metadata = $img_id_to_metadata{$img_id};
    print join("\t", (
        $line, 
        $metadata->{'name'},
        '',
        'IMG',
        $img_id
    )) . "\n";
}
close($img_path_fh);


sub create_headers_to_col_pos_map {
    my $header_line = shift;
    my %map_header_to_col;
    my @splitline = split /\t/, $header_line;
    for (my $i = 0; $i < scalar(@splitline); $i++) {
        my $col_header = $splitline[$i];
        if ($col_header) {
            $map_header_to_col{$col_header} = $i;
        }
    }
    return \%map_header_to_col;
}

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("help+", "genome_path_listfile|f:s", "metadata_tsv|m:s" );
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'genome_path_listfile'} ) { printParamError ("We need a file that lists the paths of the img genomes!"); }
    if(!exists $options{'metadata_tsv'} ) { printParamError ("We need a img metadata tsv file!"); }    

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

img_to_batchfile.pl

=head1 COPYRIGHT


=head1 DESCRIPTION

=head1 SYNOPSIS

Usage: img_to_batchfile.pl -f <genome_paths_file> -m <IMG_metadata_TSV_file>

    e.g. img_to_batchfile.pl -f img_filepaths.list -m metadata.tsv

=cut

#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util qw(min max);
use Carp qw(croak carp);
use File::Spec;

my $options = checkParams();

print STDERR Dumper($options);

opendir(my $dh, $options->{'dir'}) or croak("Cannot open directory: " . $options->{'dir'});
while (my $localpath = readdir($dh)) {
    my $path = File::Spec->join($options->{'dir'}, $localpath);
    if (-f $path) {
        my $genus_info = parse_genus_html($path);
        #print Dumper($genus_info);
        if ($genus_info) {
            my $is_type_genus = ($genus_info->{genus}->{info} =~ /Type genus of the/ ? 1 : 0);
            print join("\t", ("genus", $is_type_genus, $genus_info->{genus}->{name}, $genus_info->{genus}->{info})) . "\n";
            my %species_hash = %{$genus_info->{species}};
            foreach my $species (keys %species_hash) {
                if (ref($species_hash{$species})) {
                    if (! defined $species_hash{$species}->{info}) {
                        print STDERR Dumper($species_hash{$species});
                        print STDERR $genus_info->{genus}->{name}, " ", $species;
                    }
                    print join("\t", (
                        "species",
                        ($species_hash{$species}->{info} =~ /Type species of the genus\./ ? 1 : 0),
                        $genus_info->{genus}->{name} . " " . $species,
                        $species_hash{$species}->{info},
                        @{$species_hash{$species}->{type_strain_synonyms}}
                    )) . "\n";
                }
            }
        }
    }
}

closedir($dh);


sub parse_genus_html {
    my $path = shift;
    my $fh;
    if (! open($fh, $path)) {
        carp "Ignoring " . $path . ". Cannot open file.";
        return;
    }
    my $genus;
    my $genus_info;

    my $current_species_name;
    my $current_species_info;

    my %species_hash;

    my $temp_match;

    while (my $line = <$fh>) {

        # Carriage returns everywhere!
        $line =~ s/\r//g;
        chomp $line;

        # Change all single quotes to double quotes
        $line =~ s/\'/\"/g;

        #$line = '<span class="genusspecies">Saccharomonospora</span> Nonomura and Ohara 1971, <i>genus</i>.';

        # Genus line looks like:
        # <span class="genusspecies">Desulfomusa</span> Finster <i>et al.</i> 2001, gen. nov.

        # Species line looks like
        #<span class="genusspecies"><a name="lacus" id="lacus"></a>Methylohalomonas</span> <span class="specificepithet">lacus</span> Sorokin <i>et al.</i> 2007, sp. nov. (Type species of the genus.)


        my $genusspecies_span_match = '<span class="genusspecies">([^<]*)</span>(.*)';
        my $species_span_match = '<span class="specificepithet">([^<]*)</span>(.*)';
        my $subspecies_span_match = '<span class="subspecificepithet">([^<]*)</span>';

        #my $genusspecies_span_start = '(<span class="genusspecies(-subspecies)?">)([^<]+)';
        my $genusspecies_span_end = '</span>';

        my $type_strain_span = quotemeta('<span class="taxon-subhead">Type strain:</span>');

        $line = clean_LSPN_html_line($line);

        if ($line =~ /$genusspecies_span_match/) {
            my @matches = ($1);
            my $genus_residual_match = $2;

            # Skip empty matches
            if ($matches[0] =~ /^(\s)*$/) {
                next;
            }

            if ($line =~ /$species_span_match/) {
                push(@matches, $1);
                $current_species_info = $2;

                if ($line =~ /$subspecies_span_match/) {
                    push(@matches, $1);
                }
            }

            # Strip leading/trailing whitespace
            for (my $i = 0; $i < scalar @matches; ++$i) {
                $matches[$i] =~ s/^(\s*)([^\s]*)(\s*)$/$2/;
            }

            # Break any cases where the genus and species werent split;
            if (scalar @matches == 1) {
                @matches = split(/\ /, $matches[0]);
                $current_species_info = $genus_residual_match;
            }

            if (scalar @matches == 1) {
                if (defined $genus) {
                    carp "Multiple definitions of genus in $path. Taking first encountered definition";
                    print STDERR "\n" . $temp_match . "\n";
                    print STDERR $line . "\n\n\n";
                } else {
                    $genus = $matches[0];
                    $genus_info = $genus_residual_match;
                    $temp_match = $line;
                }
            }

            # Catch any lines with spaces
            for (my $i = 0; $i < scalar @matches; ++$i) {
                if ($matches[$i] =~ /\ /) {
                    carp "Space found in name, unable to handle. Skipping line: \n" . $line;
                };
            }

            # Don't bother with subspecies
            if (scalar @matches > 3) {
                next;
            }

            if (scalar @matches == 2) {
                $current_species_name = $matches[1];
            }

        } elsif ($line =~ /$type_strain_span([^<]*)/) {

            if (! $current_species_name) {
                next;
            }

            my $strain_info = $1;

            $strain_info =~ s/\(see also StrainInfo\.net\) (strain)*//;

            # Remove leading/trailing whitespace/dots.
            $strain_info =~ s/^[\s]*//;
            $strain_info =~ s/\.[\s]*$//;

            my @strain_synonyms = split /=/, $strain_info;

            # Remove leading/trailing whitespace
            map {$_ =~ s/^[\s]*(.*?)[\s]*$/$1/} @strain_synonyms;

            if (defined $species_hash{$current_species_name}) {
                carp "$current_species_name already previously defined in $path. Taking first defintion";
            } else {
                $species_hash{$current_species_name} = {
                    'info' => $current_species_info,
                    'type_strain_synonyms' => \@strain_synonyms
                };
            }

            $current_species_name = "";
            $current_species_info = "";
        }
    }

    if (! defined $genus) {
        return;
    }

    return {
        'genus' => {
            'name' => $genus,
            'info' => $genus_info
        },
        'species' => \%species_hash
    };

}

sub clean_LSPN_html_line {
    my $line = shift;

    #$line = 'Bob &ampid; <i>Bobette</i>  <i>&eacute;</i> <a href="img.src" class="zzzzzz">Something is <b>here</b></a>';

    # remove HTML markup
    $line =~ s/\&.{1,10}\;/?/g;

    # remove italics
    $line =~ s#<i>([^<]*)<\/i>#$1#g;

    # remove hyperlinks
    $line =~ s#<a[^>]*>(.*?)<\/a>#$1#g;

    # print $line . "\n";
    return $line;
}


sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("help+", "dir|d:s");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options );

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    if(!exists $options{'dir'} ) { printParamError ("We need a directory of genus HTMLs to parse!"); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}


__DATA__

=head1 NAME

LPSN_parser.pl

=head1 COPYRIGHT

Adam Skarshewski

=head1 DESCRIPTION

Parse the LSPN information

=head1 SYNOPSIS


Usage: LPSN_parser.pl -d HTML_directory



    -help   Display this help
    -man    Display extensive help
    -d      Directory containing HTML files


=cut

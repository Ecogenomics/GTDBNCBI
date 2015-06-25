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
        if ($genus_info) {
            print join("\t", ("genus", 0, $genus_info->{genus}->{name}, $genus_info->{genus}->{info})) . "\n";
            my %species_hash = %{$genus_info->{species}};
            foreach my $species (keys %species_hash) {
                if (ref($species_hash{$species})) {
                    print join("\t", (
                        "species",
                        ($species_hash{$species}->{info} =~ /Type species of the genus\./ ? 1 : 0),
                        $species,
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
    
    while (my $line = <$fh>) {
        
        # Carriage returns everywhere!
        $line =~ s/\r//g;        
        chomp $line;
        
        # Change all single quotes to double quotes
        $line =~ s/\'/\"/g;        
        
        #$line = '<span class="genusspecies">Saccharomonospora</span> Nonomura and Ohara 1971, <i>genus</i>.';
        
        my $genusspecies_span_start = '(<span class="genusspecies(-subspecies)?">)([^<]+)';
        my $genusspecies_span_end = '</span>';
        
        my $type_strain_span = quotemeta('<span class="taxon-subhead">Type strain:</span>');
        
        $line = clean_LSPN_html_line($line);
        
        if (my @c = $line =~ /$genusspecies_span_start/g) {
            
            # Divide by 3 as we make two matches
            my $count = (scalar @c)/3;
            
            # Check that the genus isnt both genus and species
            if ($count == 1) {
                
                my $full_putative_genus = $3;
                my $putative_genus = $full_putative_genus;
                #Remove leading/trailing whitespace
                $putative_genus =~ s/^[\s]*(.*?)[\s]*$/$1/;
                
                # Check that the genus isnt genus and species
                my @split_genus = split (/\s+/,  $putative_genus);
    
                $count = scalar @split_genus;
                
                if ($count > 1) {
                    my $regex_start = quotemeta($c[0]);
                    my $regex_end = quotemeta($full_putative_genus);
                    my $corrected_html = join('</span><span class="genusspecies">', @split_genus);
                    $line =~ s/($regex_start)($regex_end)/$1$corrected_html/;
                    
                }
            }
            
            # if you only see the <span class="genusspecies"> tag once it's a species.
            if ($count == 1) {
                if (defined $genus) {
                    carp "Multiple definitions of genus in $path. Ignoring file.";
                    return;
                }
                if ($line =~ /($genusspecies_span_start)$genusspecies_span_end(.*)/) {
                    $genus = $4;
                    $genus_info = $5;
                }
                
            } else {
                $current_species_name = "";
                $current_species_info = "";
                my @name_parts;

                while ($line =~ /($genusspecies_span_start)$genusspecies_span_end/g) {
                    push @name_parts, $4;
                }
                $current_species_name = join(" ", @name_parts);
                
                # Remove the species name and preceding chars
                $line =~ s/(.*)($genusspecies_span_start)$genusspecies_span_end//g;
                

                
                $current_species_info = $line;
                
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
                carp "$current_species_name already previously defined in $path";
                $species_hash{$current_species_name} = 'multiples';
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
    
    #print $line . "\n";
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




#!/usr/bin/env perl

###############################################################################
#
# __Script__Name__
#
# <one line to give the program's name and a brief idea of what it does.>

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
use Data::Dumper;

#CPAN modules
use Bio::SearchIO;
use Bio::SeqIO;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $options = checkParams();

######################################################################
# CODE HERE
######################################################################

opendir(my $dh, $options->{'d'});


while (my $path = readdir($dh)) {
    if ($path =~ /\.blast$/) {
        my $blast_obj = new Bio::SearchIO(-format => 'blast', 
                                            -file   => $path);
        foreach my $result (@{parse_blast($blast_obj)}) {
            print join ("\t", ($path, $result->{'subject'}, $result->{'query'}, $result->{'match'}, $result->{'length'}, $result->{'score'})) . "\n";
        }
    }
}


sub parse_blast {
    my $blast_obj = shift;
    my @results;
    while (my $result = $blast_obj->next_result()) {
        my $best_score = 0;
        my $best_result;
        while (my $hit = $result->next_hit()) {
            while (my $hsp = $hit->next_hsp()) {
                if ($hsp->length('total') > 300) {
                    if ($hsp->score > $best_score) {
                        $best_score = $hsp->score;
                        $best_result = {
                            'query' => $result->query_name,
                            'subject' => $hit->name,
                            'score' => $hsp->score,
                            'percent' => $hsp->frac_identical('total'),
                            'length' => $hsp->length('total'),
                            'match' => $hsp->frac_identical('total') * $hsp->length('total')
                        };
                    }
                }
            }
        }
        if (defined($best_result)) {
            push @results, $best_result;
        }
    }
    return \@results;
}


######################################################################
# CUSTOM SUBS
######################################################################

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options, "d:s");

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'d'} ) { printParamError ("Need to supply a BLAST directory"); }
    
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

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
# -a => ""
# },
# WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC


__DATA__

=head1 NAME

__Script__Name__

=head1 DESCRIPTION

Insert detailed description here

=head1 SYNOPSIS

__Script__Name__ [-help|h]

[-help -h] Displays basic usage information
=cut


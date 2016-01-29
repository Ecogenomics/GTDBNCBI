#!/usr/bin/perl -w

use strict;

use Bio::SeqIO;
use Text::Wrap;

if (@ARGV == 0) { die "
\tExtracts all the CDS nucleotide sequences from a GenBank file (*.gbk). 
\tNB must be a 'full' gbk file, with the source sequence at the end.\n
\tUSAGE: genbank2nt_fasta.pl <myGenbankFile> > myFastaFile.fna\n
"}

## wraps sequence string every 80 characters:
$Text::Wrap::columns = 80;
$Text::Wrap::separator = "\n";

my $in = Bio::SeqIO -> new( -file => $ARGV[0], -format => 'genbank');

my $i;

while (my $seq = $in->next_seq) {
    for my $feat ($seq -> get_SeqFeatures) {
        if ($feat->primary_tag eq "CDS") {
            print "\>";
            ## construct fasta header based on information present for each CDS:
            if ( $feat->has_tag('locus_tag') ) {
                print $feat->get_tag_values('locus_tag')," ";    
            } 
            if ( $feat->has_tag('db_xref') ) {
                print $feat->get_tag_values('db_xref') ,"\|";    
            } 
            if ( $feat->has_tag('protein_id') ) {
                print $feat->get_tag_values('protein_id'),"\|";    
            } 
            if ( $feat->has_tag('gene') ) {
                print $feat->get_tag_values('gene'),"\|";            
            } 
            if ( $feat->has_tag('product') ) {
                print $feat->get_tag_values('product');
            }
            ## print sequence:
            print "\n",wrap("","",$feat->spliced_seq->seq),"\n";
            $i++;
        }
    }
}

print STDERR "\n\tExtracted $i sequences from file $ARGV[0]\n\n";

__END__

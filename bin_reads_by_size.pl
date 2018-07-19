#!/usr/lib/perl5
use warnings;
use strict;

=head DESCRIPTION

bin reads by size.
command line arguments should be:
(1) the fasta file containing the reads
(2) the min and max of the division in bp
(3) the name of the division

=cut

use Bio::SeqIO;

my $fasta_file = $ARGV[0];
my $min_len = $ARGV[1];
my $max_len = $ARGV[2];
my $bin_name = $ARGV[3];
$bin_name = "_".$bin_name;

my $fasta_out = $fasta_file;
$fasta_out =~ s/\.fasta/$bin_name\.fasta/;

my $seq_in = Bio::SeqIO->new(-file => "$fasta_file",
							-format => "fasta");
my $seq_out = Bio::SeqIO->new(-file => ">$fasta_out",
							-format => "fasta");
my $i=0;
while (my $seq = $seq_in->next_seq) {
	if($seq->length > $min_len && $seq->length < $max_len) {
		$seq_out->write_seq($seq);
	}
	$i++;
	if($i % 10000==0) {print $i/10000, "e4 reads\n"; }
}

exit;













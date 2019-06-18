#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long::Descriptive;
use GenOO::Data::File::SAM;
use GenOO::Data::File::BED;

my ($opt, $usage) = describe_options(
	"Usage: %c %o\n" .
	"Convert a SAM file into a BED file.",
	['ifile|i=s', 'SAM input file. If empty it reads from STDIN'],
	['help|h', 'print usage message and exit'],
	['verbose|v', 'show progress'],
);
print($usage->text), exit if $opt->help;

warn "opening SAM file\n" if $opt->verbose;
my $fp = GenOO::Data::File::SAM->new(file => $opt->ifile);

warn "converting to BED\n" if $opt->verbose;
while (my $r = $fp->next_record) {
	next if $r->is_unmapped;

	my $start = $r->start;
	my $stop = $r->stop;
	
	print join("\t", $r->rname, $start, $stop+1, $r->qname, '1', $r->strand_symbol)."\n";
}

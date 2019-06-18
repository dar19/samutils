#!/usr/bin/env perl

# Load modules
use Modern::Perl;
use Getopt::Long::Descriptive;

# Load GenOO library
use GenOO::Data::File::SAM;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Measure value distribution of the provided tag in a SAM file."],
	["Measures how many reads have each particular value for the provided tag."],
	["The output is a tab delimited file"],
	[],
	['input|i=s', 'input SAM file. If not set use STDIN.'],
	['tag|t=s', 'tag for which to measure distribution.', {required => 1}],
	['verbose|v', "print progress"],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

warn "opening input file\n" if $opt->verbose;
my $fp = GenOO::Data::File::SAM->new(file => $opt->input);

warn "measuring distribution for tag ".$opt->tag."\n" if $opt->verbose;
my $notag = 0;
my %counts;
while (my $rec = $fp->next_record) {
	my $val = $rec->tag($opt->tag);

	if (!defined $val) {
		$notag++;
		next;
	}
	$counts{$val}++;
}

warn "writting output table\n" if $opt->verbose;
say join("\t", "tag", "count", "freq");
foreach my $k (keys %counts) {
	say join("\t", $k, $counts{$k}, $counts{$k} / $fp->records_read_count);
}

warn "info: tag not found in $notag reads\n" if $opt->verbose;

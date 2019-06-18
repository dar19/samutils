#!/usr/bin/env perl

=head1 NAME

sam-to-collapsed - Collapse the records of a sorted SAM file

=head1 SYNOPSIS

sam-to-collapsed [options/parameters] <sam>

Collapse the records of a sorted SAM file. Reads that have exactly the same
coordinates are collapsed into a single entry. A new tag that reports the
number of records that have been collapsed is added to the new collapsed
record. The quality field is set to * and the name of the first record is used
for the new collapsed record. The input SAM file must be sorted by coordinates
(rname and pos).

=cut

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use IO::Interactive;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Collapse the records of a sorted SAM file. Reads that have exactly the
	same coordinates are collapsed into a single entry. A new tag that
	reports the number of records that have been collapsed is added to the
	new collapsed record. The quality field is set to * and the name of
	the first record is used for the new collapsed record. The input SAM
	file must be sorted by coordinates (rname and pos)."],
	[],
	['tag|t=s',
		'new tag with the number of collapsed records (Default: XC:i).',
		{default => "XC:i"}],
	['verbose|v', "print progress"],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

# If sam option is not given and no data are coming from a pipe.
if (IO::Interactive::is_interactive()) {
	say "Error: No input data.";
	print($usage->text);
	exit 1;
}

my $new_tag = $opt->tag;

warn "collapsing records\n" if $opt->verbose;
my %coords_counts;
my ($current_pos, $current_rname);
while (my $line = <>) {
	chomp $line;

	# Skip the header
	if ($line =~ /^\@/) {
		say $line;
		next;
	};

	# Split the SAM line and get required fields
	my @fields = split(/\t/, $line);
	my ($flag, $rname, $pos, $seq) = @fields[1, 2, 3, 9];
	$fields[10] = '*'; # delete quality

	# Skip unmapped reads
	next if $flag & 4; #unmapped

	# Get the strand
	my $strand = $flag & 16 ? -1 : 1;

	# Check if a similar record has been read
	if (defined $current_pos and defined $current_rname and
		($rname eq $current_rname) and ($pos == $current_pos)) {

		if (!exists $coords_counts{$strand}{$seq}) {
			$coords_counts{$strand}{$seq} = [@fields, 0];
		}

		$coords_counts{$strand}{$seq}->[-1]++;
	}
	else {
		if (defined $current_pos and defined $current_rname) {
			print_records_on_previous_pos(
				$current_pos, $current_rname, \%coords_counts);
		}
		%coords_counts = ();
		$current_rname = $rname;
		$current_pos = $pos;

		$coords_counts{$strand}{$seq} = [@fields, 1];
	}
}
print_records_on_previous_pos($current_pos, $current_rname, \%coords_counts);

exit;

###########################################
sub print_records_on_previous_pos {
	my ($pos, $rname, $stored_records_hashref) = @_;

	foreach my $strand (keys %{$stored_records_hashref}) {
		foreach my $seq (keys %{$stored_records_hashref->{$strand}}) {
			my $fields = $stored_records_hashref->{$strand}->{$seq};
			$fields->[-1] = join(':', $new_tag, $fields->[-1]);
			say join("\t", @$fields);
		}
	}
}

#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use IO::Uncompress::Gunzip qw($GunzipError) ;

# Define and read command line options
my ($opt, $usage) = describe_options(
	"Usage: %c %o",
	["Add a new tag to each SAM record populated by a corresponding sequence in a FASTQ file. SAM and FASTQ records are combined based on the qname and header respectively."],
	[],
	['input=s', "input SAM file; reads from STDIN if \"-\"", {required => 1}],
	['fastq=s', "input FASTQ file; reads from STDIN if \"-\"", {required => 1}],
	['tag|t=s', "new tag that is populated by FASTQ sequences", {required => 1}],
	['tag-val=s', "if set, the new tag is populated with this value instead of FASTQ sequences"],
	['verbose|v', "print progress"],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

my $new_tag = $opt->tag;
my $new_tag_val = $opt->tag_val;

if ($opt->input eq '-' and $opt->fastq eq '-') {
	die "Cannot read from STDIN for both \'input\' and \'fastq\'\n";
}

# read the sequences from the FASTQ file.
my %seqs;
my $FASTQ = filehandle_for($opt->fastq);
while (my $name = $FASTQ->getline) {
	if ($name !~ /^\@/) {
		next;
	}
	chomp $name;
	$name = substr($name, 1);
	my $seq = $FASTQ->getline;
	chomp $seq;
	$FASTQ->getline;
	$FASTQ->getline;
	$seqs{$name} = $seq;
}
close($FASTQ);

# append the new tags to SAM records.
my $SAM = filehandle_for($opt->input);
while (my $line = $SAM->getline) {
	chomp $line;

	# Skip the header
	if ($line =~ /^\@/) {
		say $line;
		next;
	};

	# Split the SAM line and get required fields
	my @fields = split(/\t/, $line);
	my ($qname) = $fields[0];

	if (exists $seqs{$qname}) {
		if (defined $new_tag_val) {
			push @fields, join(':', $new_tag, $new_tag_val);
		} else {
			push @fields, join(':', $new_tag, $seqs{$qname});
		}
	}

	say join("\t", @fields);
}
close($SAM);

exit;

sub filehandle_for {
	my ($file) = @_;

	if (!defined $file) {
		die "Undefined file\n";
	}
	if ($file eq '-'){
		return IO::File->new("<-");
	}
	if ($file =~ /\.gz$/) {
		my $z = new IO::Uncompress::Gunzip $file or die "gunzip failed: $GunzipError\n";
		return $z;
	}
	return IO::File->new($file, "<");
}


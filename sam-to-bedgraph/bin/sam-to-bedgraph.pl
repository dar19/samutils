#!/usr/bin/env perl

use Modern::Perl;
use autodie;
use Getopt::Long::Descriptive;
use PDL::Lite; $PDL::BIGPDL = 0; $PDL::BIGPDL++; # enable huge pdls
use GenOO::Data::File::SAM;

# Define and read command line options
my ($opt, $usage) = describe_options(
	'Usage: %c %o', 
	['Convert a SAM file into a bedgraph file.'],
	[],
	['sam=s', 'SAM file. If empty read from STDIN'],
	['chr_sizes=s', 'File with chromosome sizes. (<chrom><TAB><size>)',
		{required =>1 }],
	['split|s', 'skip empty fasta lines'],
	['verbose|v', 'show progress'],
	['help|h', 'print usage and exit', {shortcircuit => 1}],
);
print($usage->text), exit if $opt->help;

# If sam option is not given and no data are coming from a pipe.
if (!defined $opt->sam and -t STDIN) {
	say "Error: No input data.";
	print($usage->text), exit 1 if $opt->help;
}

warn "reading chromosome sizes\n" if $opt->verbose;
my %chr_sizes = read_rname_sizes($opt->chr_sizes);

warn "opening input file\n" if $opt->verbose;
my $fp = GenOO::Data::File::SAM->new(file => $opt->sam);

warn "converting to bedgraph\n" if $opt->verbose;
my $pdl;
my $prev_rname = '';
my %already_read;
while (my $r = $fp->next_record) {
	if ($r->rname ne $prev_rname) {
		if (exists $already_read{$r->rname}) {
			die "Error: Input file is not sorted by chromosome\n";
		}
		print_bedgraph($pdl, $prev_rname) if defined $pdl;
		$already_read{$r->rname} = 1;
		$prev_rname = $r->rname;
		$pdl = PDL->zeros(PDL::long(), $chr_sizes{$r->rname});
	}

	if ($opt->split) {
		my ($starts, $stops) = $r->aligned_areas_starts_and_stops();
		for (my $i=0; $i<@$starts; $i++) {
			$pdl->slice([$starts->[$i], $stops->[$i]]) += 1;
		}
	}
	else {
		$pdl->slice([$r->start, $r->stop]) += 1;
	}
}
print_bedgraph($pdl, $prev_rname) if defined $pdl;

######################################################################
sub print_bedgraph {
	my ($pdl, $chr) = @_;

	my ($factor, $minS, $initS) = (2, 16, 512);
	my ($i, $pos, $prev_score) = (0, 0, 0);
	my $len = $pdl->getdim(0);
	while ($i < $len) {
		#jump if no values found within $S
		my $S = $initS;
		my $found = 0;
		while ($i < $len and $S >= $minS) {
			if ($i + $S - 1 >= $len) {
				$S = $S / $factor;
			}
			elsif ($pdl->slice([$i, $i+$S-1])->sum() != 0) {
				$S = $S / $factor;
				$found = 1;
			}
			else {
				$i += $S;
				if (not $found) {
					$S = $S * $factor;
				}
			}
		}

		while ($i < $len) {
			my $score = $pdl->at($i);
			if ($score != $prev_score) {
				if ($prev_score != 0) {
					say join("\t", $chr, $pos, $i, $prev_score);
				}
				$pos = $i;
				$prev_score = $score;
				if ($score == 0) {
					$i++;
					last;
				}
			}
			$i++;
		}
	}
	if ($prev_score != 0) {
		say join("\t", $chr, $pos, $pdl->getdim(0), $prev_score);
	}
}

sub read_rname_sizes {
	my ($file) = @_;

	my %rname_size;
	open (my $CHRSIZE, '<', $file);
	while (my $line = <$CHRSIZE>) {
		chomp $line;
		my ($chr, $size) = split(/\t/, $line);
		$rname_size{$chr} = $size;
	}
	close $CHRSIZE;
	return %rname_size;
}

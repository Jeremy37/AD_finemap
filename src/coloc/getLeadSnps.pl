#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Common qw(open_for_read2 open_for_write2 verbose checkFileExists checkDirExists);

sub usage($);
my $usage = "$0 - Get the lead SNP and the p value difference to the second-best SNP.
    usage:
       $0 [OPTIONS] <inFile|stdin>
    options:
       -verbose=XXX - set verbosity level
       -XXX\n";

## Read options from CL
my $FvalField = 11;
my ($inFile, $statcol, $pcol, $genecol, $insep);
GetOptions(	'statcol=i' => \$statcol,
			'genecol=i' => \$genecol,
			'pcol=i' => \$pcol,
			'file|f=s' => \$inFile,
			'sep=s' => \$insep);

(!defined $inFile) and usage("Error: required argument --file is missing.\n");
(!defined $pcol and !defined $statcol) and usage("Error: either --pcol or --statcol argument is required.\n");
(!defined $genecol) and usage("Error: --genecol argument is required.\n");

(defined $genecol) and $genecol -= 1;
(defined $pcol) and $pcol -= 1;
(defined $statcol) and $statcol -= 1;
my $sep = "\t";
(defined $insep) and $sep = $insep;

## Check CL args
#&checkFileExists($inFile) unless $inFile eq 'stdin';

## Read input files
my @bestLine;
my @secondBestLine;
my $curGene = 'NA';
my $numTests = 1;
$" = "\t";
*IN = &open_for_read2($inFile);
while(<IN>) {
	chomp();
	my @l = split($sep, $_, -1);
	if ($l[$pcol] eq "") {
		next; # Ignore lines with empty p value
	}
	if ($curGene ne 'NA') {
		if ($curGene eq $l[$genecol]) {
			$numTests += 1;
			if (defined $pcol) {
				if ($l[$pcol] < $bestLine[$pcol]) {
					@secondBestLine = ();
					@secondBestLine = @bestLine;
					@bestLine = @l;
				} elsif (scalar(@secondBestLine) < 2 or $l[$pcol] < $secondBestLine[$pcol]) {
					@secondBestLine = @l;
				}
			} else {
				if (abs($l[$statcol]) > abs($bestLine[$statcol])) {
					@secondBestLine = ();
					@secondBestLine = @bestLine;
					@bestLine = @l;
				} elsif (scalar(@secondBestLine) < 2 or abs($l[$statcol]) > abs($secondBestLine[$statcol])) {
					@secondBestLine = @l;
				}
			}
		} else {
			print STDERR ($l[$genecol]."\n");
			my $secondBest = 1;
			if (scalar(@secondBestLine) > 1) {
				if (defined $pcol) {
					$secondBest = $secondBestLine[$pcol];
				} else {
					$secondBest = $secondBestLine[$statcol];
				}
			}
			print "@bestLine\t$secondBest\t$numTests\n";
			@bestLine = ();
			@bestLine = @l;
			@secondBestLine = ();
			$curGene = $l[$genecol];
			$numTests = 1;
		}
	} else {
		print STDERR ($l[$genecol]."\n");
		@bestLine = @l;
		@secondBestLine = ();
		$curGene = $l[$genecol];
		$numTests = 1;
	}
}
close(IN) unless $inFile eq 'stdin'; ## Avoids perl warning if we try and close STDIN;

my $secondBest = 1;
if (scalar(@secondBestLine) > 1) {
	if (defined $pcol) {
		$secondBest = $secondBestLine[$pcol];
	} else {
		$secondBest = $secondBestLine[$statcol];
	}
}
print "@bestLine\t$secondBest\t$numTests\n";

###############################################################################
## Functions
sub usage($) {
	print $_[0]."\n";
    print STDERR $usage;
    exit;
}


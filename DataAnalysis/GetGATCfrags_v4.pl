# Author: Nicolas Altemose
# Date: 2/25/19

use warnings;
use strict;
my $tic = time;
print "\n\n";

my $usage = 'perl GetGATCfrags_v4.pl <input fasta file> <output bed file>';
my @chrs = qw(
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
);
my %chrkeep;
foreach my $chr(@chrs){
	$chrkeep{$chr}=0;
}

#initialise input variables
my $infile = '0_reference/GRCh38_DamLMNB1_mTracer.fasta';
my $outfile = '5_fragcov/GRCh38_ALL_DpnI_Fragments_v4.bed';

#obtain input variable values from command line
if(defined $ARGV[0]){
	$infile = $ARGV[0];
	chomp($infile); #remove any newlines from the end of the input value
}
else{
	die "$usage\n\n";
}
if(defined $ARGV[1]){
	$outfile = $ARGV[1];
	chomp($outfile);
}
else{
	die "$usage\n\n";
}

#open input and output files and process

open(OUT,'>'.$outfile); #open output file for writing
local $/ = ">"; #change the record separator (locally)-i.e. split files by '>' instead of newlines
open(FASTA, $infile) or die "ERROR: unable to open $infile\n"; #open input file for reading

my $ct=0; #initialise a counter to count the number of substrings printed
my $junk = <FASTA>; #remove first > from beginning of file
while(my $frecord = <FASTA>){ #read in input file one fasta record at a time
	chomp($frecord); #remove any newlines from the end of the record
	if($frecord=~/^(\S+).+?\n(.*)$/s){ #check whether each record matches the desired chromosome
		if(exists $chrkeep{$1}){
			$chrkeep{$1}=1;
			my $chr = $1;
			my $seq = uc($2); #store the entire chromosome sequence in a string
			$seq=~s/\s//g; #eliminate any newline characters from the chromosome sequence
			my $len = length($seq);
			print "$chr\t$len\n";
			my $prevpos=0;
			while($seq =~ /GATC/gso){
				my $idx = $-[0];
				print OUT "$chr\t$prevpos\t$idx\n";
				$prevpos = $idx + 4;
				$ct++;
			}#closes while
			print OUT "$chr\t$prevpos\t$len\n";
			$ct++;
		}#closes if
	}#closes if	
}#closes while
close FASTA;
close OUT;


foreach my $chr(@chrs){
	if($chrkeep{$chr}==0){
		print "ERROR: chromosome $chr not found in input file $infile\n";
	}
}

print "printed $ct fragments to $outfile\n";


#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));

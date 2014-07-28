#!/usr/bin/perl
use strict;
#require "/mnt/crick/asab/data/local/subroutines.pl";

# will make a summary of the star mapping to genome
# input the star_assembly folder  


my $indir = $ARGV[0];
unless ($indir =~ /\/$/) { $indir .= "/"; }

my @samples = glob($indir."*");

my $nread = "Number of input reads";
my $len = "Average input read length";
my $uniq = "Uniquely mapped reads number";
my $mlen = "Average mapped length";
my $nsplice = "Number of splices: Total";
my $noncanonical = "Number of splices: Non-canonical";
my $multimap = "Number of reads mapped to multiple loci";
my $multimap2 = "Number of reads mapped to too many loci";
my $u1 = "% of reads unmapped: too many mismatches";
my $u2 = "% of reads unmapped: too short";
my $u3 = "% of reads unmapped: other";

print (join "\t", ("Sample","Reads","Readlength","Uniqely mapping (%)","Multimapping reads (%)","Unmapped reads (%)","Splices","Non-canonical splices (%)","Coverage first 5%","Mid 5%", "Last 5%","Junctions","Known","Partial","Novel\n"));

sub sum_array { #input is an array, sum of all the elements is returned
    my $sum = 0;
    foreach my $a (@_){
	$sum += $a;
#	print ":$a:";
    }
    return ($sum);
}

foreach my $sample (@samples) { 
    my ($sname)=($sample =~ /$indir(.*)$/);
    #read summary from star alignment
    my $infile = $sample."/Log.final.out"; 
    unless (-e $infile) {  print "$sname no file Log.final.out\n";  next; }
    open(IN, $infile) || die "no file $infile\n";
    my %hash;
    while (<IN>) { 
	chomp;
	unless (/\|/) { next; }
	my ($name,$val)=split(/\s\|\s+/,$_);
	$name =~ s/^\s+//;
	$val =~ s/\%//;
	$hash{$name}=$val;
#	print ":$name: :$val: $_\n";
    }
    close IN;
    if (!defined $hash{$nread} || $hash{$nread}==0) { print "\tno data for $sname\n"; next;}
    my $nsplice = $hash{$nsplice};
    my $splice_perc;
    if($nsplice == 0){
	$splice_perc = 0;
    }else{
	$splice_perc = $hash{$noncanonical} / $nsplice * 100;
    }
    printf ("%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\t%.4f\t",$sname,$hash{$nread},$hash{$len},$hash{$uniq}/$hash{$nread}*100,($hash{$multimap}+$hash{$multimap2})/$hash{$nread}*100,$hash{$u1}+$hash{$u2}+$hash{$u3},$hash{$nsplice}, $splice_perc);
    # then read gene coverage stats (100 bins, report diff to mean coverage in first 5, mid 5 and last 5 bins)
    my $infile = $sample."/qc/geneBody_coverage.geneBodyCoverage.txt";
    unless (-e $infile) {  print "-\t-\t-\t"; }
    else {
   open(IN, $infile) || die "no file $infile\n";
    my @coverage; my $sum=0;
    while (<IN>) {
	unless (/^\d/) { next; }
	chomp;
	my ($quant,$count)=split;
	$coverage[$quant]=$count;
	$sum += $count;
    }
    close IN;
    my $mean = $sum/100;
   if ($mean == 0) { print "-\t-\t-\t"; }
   else {
    my $first = &sum_array(@coverage[0..4])/5/$mean;
    my $last = &sum_array(@coverage[95..99])/5/$mean;
    my $mid =  &sum_array(@coverage[48..52])/5/$mean;
    printf ("%.4f\t%.4f\t%.4f\t",$first,$mid,$last);
    }
    }
    my $infile = $sample."/qc/jxn_annotation_ensgene.txt";
    unless (-e $infile) {  print "-\t-\t-\t-\n"; }
    else {
    open(IN, $infile) || die "no file $infile\n";
    my $total=0; my $known=0; my $partial=0; my $novel=0;
    while (<IN>) {
        unless (/Junctions/) { next; }
        chomp;
	my ($num) = ($_ =~ /(\d+)/);
	if (/Total/) { $total = $num; }
	elsif (/Known/) { $known = $num; }
	elsif (/Partial/) { $partial = $num; }
	elsif (/Novel/) { $novel = $num; }
    }
    close IN;
    if ($total == 0) {  print "-\t-\t-\t-\n"; }
    else {    printf ("%d\t%.4f\t%.4f\t%.4f\n",$total,$known/$total*100, $partial/$total*100, $novel/$total*100); }
    }

}

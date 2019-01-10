use Getopt::Long;
use Data::Dumper;
use vcfCompare;


mkdir("Comparison_results");
open(CSUM,">>Comparison_results\\countSum.txt");
print CSUM "Standard\tScriptOutput\tTP\tFN\tFP\n";
close(CSUM);

my ($file1, $file2)="";
my $result = GetOptions(
"file1=s" => \$file1,
"file2=s" => \$file2,
);

if($file1 eq "" or $file2 eq ""){
	print "Please Specify Correct FILE NAME !!!";
	exit();
}

if($file1=~/\.txt?$/){
	print"Processing excel reports !!";
	vcfCompare::compareExcels($file1,$file2);
}
elsif($file1=~/\.vcf?$/){
	print"Processing vcf !!";
	vcfCompare::compareVcfs($file1,$file2);
}
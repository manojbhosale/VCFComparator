use Getopt::Long;
use Data::Dumper;
use vcfCompare;



my ($file1, $file2)="";
my $result = GetOptions(
"F1=s" => \$file1,
"F2=s" => \$file2,
"D=s" => \$destinationFolderPath,
) or die("Error in commandline arguments !!");

my $resFolderPath = "";

if(defined $destinationFolderPath){
	$resFolderPath = File::Spec->catfile( $destinationFolderPath, "Comparison_results" );
}else{
	$resFolderPath = File::Spec->catfile( File::Spec->curdir(), "Comparison_results" );
}	

mkdir($resFolderPath) or die($!);
my $countSumPath = File::Spec->catfile( $resFolderPath, "countSummary.txt" );
open(CSUM,">>","$countSumPath");
print CSUM "Standard\tScriptOutput\tTP\tFN\tFP\n";
close(CSUM);

if($file1 eq "" or $file2 eq ""){
	print "Please Specify Correct FILE NAME !!!\n";
	exit();
}

vcfCompare::compareJobFoldersOrFile($file1,$file2,$resFolderPath);



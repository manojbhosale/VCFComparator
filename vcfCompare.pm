package vcfCompare;
use Data::Dumper;
use strict;
use File::Spec;
use File::Basename;


sub compareVcfs{
	my $file1 = shift;
	my $file2 = shift;
	my $flag = 0;
	
	my $file1_hash = {};
	my $combined_hash = {};
	open(vcf1,$file1) or die($!);
	open(vcf2,$file2) or die($!);

	while(my $line = <vcf1>){
	
		if($line=~/^\#\#/){
			next;
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t/;
			
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//g;
			$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FN";
		}
		
	}
	close(vcf1);
	my $flag = 0;
	while(my $line = <vcf2>){
		if($line=~/^\#\#/){
			next
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+/;
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//;
			if(defined($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall})){
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="TP";
			}
			else{
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FP";
			}
		}
		
	}
	my $fn=0;
	my $tp=0;
	my $fp=0;
	my (@fname1,@fname2) = ();
	 $fname1[0] = $file1;
	 $fname2[0] = $file2;
	if($file1=~/\\/g){
		 @fname1 = ($file1=~/\\([\w\s\.]+)$/);
	}
	if($file2=~/\\/g){
		 @fname2 = ($file2=~/\\([\w\s\.]+)$/);
	}
	#print "$fname1[0]\_$fname2[0]\_common.txt";
	#print "@fname1\n";
	#print "Manoj\n";
	#print "$file1 $file2";
	open(CSUM,">>Comparison_results\\countSum.txt");
	open(TP,">Comparison_results\\$fname1[0]\_$fname2[0]\_common.txt") or die("$!");
	open(FP,">Comparison_results\\$fname2[0]\_exclusive.txt") or die("$!");
	open(FN,">Comparison_results\\$fname1[0]\_exclusive.txt") or die("$!");
	print TP "CHROM\tPOS\tREF\tALT\n";
	print FP "CHROM\tPOS\tREF\tALT\n";
	print FN "CHROM\tPOS\tREF\tALT\n";
	foreach my $chrom(keys %{$file1_hash}){
		foreach my $pos(keys %{$file1_hash->{$chrom}}){
			foreach my $refall(keys %{$file1_hash->{$chrom}->{$pos}}){
				foreach my $altall(keys %{$file1_hash->{$chrom}->{$pos}->{$refall}}){
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "TP"){
						print TP "$chrom\t$pos\t$refall\t$altall\n" ;
						$tp++;
					}
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FP"){
						print FP "$chrom\t$pos\t$refall\t$altall\n";
						$fp++;
					}
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FN"){
						print FN "$chrom\t$pos\t$refall\t$altall\n" ;
						$fn++;
					}
				}
			}
		}
	}
	print CSUM "@fname1\t@fname2\t$tp\t$fn\t$fp\n";
	close(TP);
	close(FP);
	close(CSUM);
	close(FN);
}

sub compareExcels{
	my $file1 = shift;
	my $file2 = shift;
	
	#Variable declare and initialize
	my $line_cnt = 0;
	my @req_headers = ("Chrom","Pos","Ref Allele","Alt Allele","Quality","Allele Frequency","Read Depth(Raw)","Forward Ref Alleles","Reverse Ref Alleles","Forward Alt Alleles","Reverse Alt Alleles","Mapping Quality");
	my (@record,%headers_hash,@headers,%req_hed_hash) = ();
	my ($std_hash,$std_hash_extra_first,$std_hash_extra_second) ={};
	$"= "\t";

	#GENERATE REQUIRED HEADERS HASH
	foreach(@req_headers){
		$req_hed_hash{$_} ="";
	}


	#Input file opening
	open(IN,$file1) or die("$!");
	open(IN_REPORT,$file2) or die("$!");;

	#Intermediate OUTPUT file opening
	while(<IN>){
		my $line = $_;
		chomp($line);
		if($line_cnt == 0 && $line=~/^.+$/){
		@headers= split(/\t/,$line);
		$line_cnt++;
			foreach(@headers){
			$headers_hash{$_}=0;
			}
		next;
		}
		
		#ASSIGN VALUES OF THE SELECTED HEADERS FROM THE RECORDS ACCORDINGLY
		@record = split(/\t/,$line);
		my $head_count = 0;
		foreach my $head(@headers){
		
			if(exists($req_hed_hash{$head})){
				$req_hed_hash{$head}=$record[$head_count];
				}
		$head_count++;
		}
		
		#TO GET RID OF MULTIPLE ALLELES 
		my @major_allele_arr = split(/\,/,$req_hed_hash{"Alt Allele"});
		$std_hash->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]}="FN";
		
		#GENERATING EXTRA FIELDS
		my $extra_fields=();
		
		foreach("Quality","Allele Frequency","Read Depth(Raw)","Forward Ref Alleles","Reverse Ref Alleles","Forward Alt Alleles","Reverse Alt Alleles","Mapping Quality"){
			$extra_fields.=$req_hed_hash{$_}."\t";
		}
		
		$std_hash_extra_first->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]}=$extra_fields;
				
		#EMPTYING HASH SO AS TO USE IN NEXT ITERATION WITH NEW VALUES.
		empty_hash(\%req_hed_hash);
		
	}

	close(IN);
	# Again empty variables to clear up the assigned values in previous file reading
	$line_cnt = 0;
	@headers = ();
	empty_hash(\%req_hed_hash);

	while(<IN_REPORT>){
			
			my $line = $_;
			chomp($line);
			if($line_cnt == 0 && $line=~/^.+$/){
			@headers= split(/\t/,$line);
			$line_cnt++;
				foreach(@headers){
				$headers_hash{$_}=0;
				}
			next;
			}
			
			#ASSIGN VALUES OF THE SELECTED HEADERS FROM THE RECORDS ACCORDINGLY
			@record = split(/\t/,$line);
			
			my $head_count = 0;
			foreach my $head(@headers){
					if(exists($req_hed_hash{$head})){
					$req_hed_hash{$head}=$record[$head_count];
					}
			$head_count++;
			}
			
			
			# TO GET RID OF MULTIPLE ALLELES 
			my @major_allele_arr = split(/\,/,$req_hed_hash{"Alt Allele"});
			
			# GENERATING EXTRA FIELDS
			my $extra_fields=();
			
			foreach("Quality","Allele Frequency","Read Depth(Raw)","Forward Ref Alleles","Reverse Ref Alleles","Forward Alt Alleles","Reverse Alt Alleles","Mapping Quality"){
			$extra_fields.=$req_hed_hash{$_}."\t";
			}
			
			$std_hash_extra_second->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]}=$extra_fields;
			if(defined($std_hash->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]})){
					$std_hash->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]} = "TP";
			}
			else{
				$std_hash->{$req_hed_hash{"Chrom"}}->{$req_hed_hash{"Pos"}}->{$req_hed_hash{"Ref Allele"}}->{$major_allele_arr[0]}="FP";
			}
			

			#EMPTYING HASH SO AS TO USE IN NEXT ITERATION WITH NEW VALUES.
			empty_hash(\%req_hed_hash);
			
	}
		
	close(IN_REPORT);
	my (@fname1,@fname2) = ();
	 $fname1[0] = $file1;
	 $fname2[0] = $file2;

	if($file1=~/\\/){
		 @fname1 = ($file1=~/\\([\w\s\.]+)$/);
	}
	if($file2=~/\\/){
		 @fname2 = ($file2=~/\\([\w\s\.]+)$/);
	}
	
	open(CSUM,">>Comparison_results\\countSum.txt");
	open(TP,">Comparison_results\\$fname1[0]\_$fname2[0]\_common.txt") or die("$!");
	open(FP,">Comparison_results\\$fname2[0]\_exclusive.txt") or die("$!");
	open(FN,">Comparison_results\\$fname1[0]\_exclusive.txt") or die("$!");
	print FN "@req_headers","\n";
	print TP "@req_headers","\n";
	print FP "@req_headers","\n";

	my $fn=0;
	my $tp=0;
	my $fp=0;
	foreach my $chr(keys %{$std_hash}){
		foreach my $pos(keys %{$std_hash->{$chr}}){
			foreach my $ref_allele(keys %{$std_hash->{$chr}->{$pos}}){
				foreach my $alt_allele(keys %{$std_hash->{$chr}->{$pos}->{$ref_allele}}){
							if($std_hash->{$chr}->{$pos}->{$ref_allele}->{$alt_allele} eq "FN"){
								print FN "$chr\t$pos\t$ref_allele\t$alt_allele\t",$std_hash_extra_first->{$chr}->{$pos}->{$ref_allele}->{$alt_allele},"\n" ;
								$fn++;
							}
							if($std_hash->{$chr}->{$pos}->{$ref_allele}->{$alt_allele} eq "TP"){
								print TP "$chr\t$pos\t$ref_allele\t$alt_allele\t",$std_hash_extra_first->{$chr}->{$pos}->{$ref_allele}->{$alt_allele},"\n" ;
								$tp++;
							}
							if($std_hash->{$chr}->{$pos}->{$ref_allele}->{$alt_allele} eq "FP"){
								print FP "$chr\t$pos\t$ref_allele\t$alt_allele\t",$std_hash_extra_second->{$chr}->{$pos}->{$ref_allele}->{$alt_allele},"\n";
								$fp++;
							}
				}
			}	
		}
	}
	print CSUM "@fname1\t@fname2\t$tp\t$fn\t$fp\n";
	close(CSUM);
	close(FN);
	close(TP);
	close(FP);
}


sub empty_hash{
	my $hash = shift;
	foreach my $key(keys %{$hash}){
	${$hash}{$key} = "";
	}
}


sub getVcfPathFromFolder{
	my $folderPath = shift;
	opendir(DIR, $folderPath) or die "Couldn't open directory, $!";
	my $abs_path = "";
	my @paths = File::Spec->no_upwards( readdir DIR);
	foreach my $file (@paths) {
		if($file=~/_4_2/){
		#print("*** ",$file,"\n");
			$abs_path = File::Spec->rel2abs( $file, $folderPath ) ;
			#print(basename($abs_path),"\n");
		}
	}
	closedir DIR;
	
	if($abs_path eq ""){
	
		die("Couldn't find the VCF with \"_4_2\" in name !!")
	
	}
	
	return $abs_path;
}

sub getLoggingTime {

    my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
    my $nice_timestamp = sprintf ( "%04d%02d%02d_%02d%02d%02d",
                                   $year+1900,$mon+1,$mday,$hour,$min,$sec);
    return $nice_timestamp;
}

sub compareJobFolders{
	my $folder1 = shift;
	my $folder2 = shift;
	my $resFolderPath = shift;
	#my $resultFolderName = "comparisonResults";
	#my $resFolderPath = File::Spec->catfile( $destinationFolderPath, $resultFolderName );
	my $file1 = getVcfPathFromFolder($folder1);
	my $file2 = getVcfPathFromFolder($folder2);
	
	#print($file1," ",$file2);
	my $fileName1 = basename($file1);
	my $fileName2 = basename($file2);
	
	my $modFile1 = $fileName1=~ s/\./_/rg;
	my $modFile2 = $fileName2=~ s/\./_/rg;
	#print($modFile1," ",$modFile2);

	#Create the output file paths
	my $csumPath = File::Spec->catfile( $resFolderPath, "countSummary.txt" );
	#my $tpPath = File::Spec->catfile( $resFolderPath, $modFile1."_".$modFile2."_common.txt" );
	my $tpPath = File::Spec->catfile( $resFolderPath, "common.txt" );
	my $fpPath = File::Spec->catfile( $resFolderPath, $modFile2."_exclusive.txt" );
	my $fnPath = File::Spec->catfile( $resFolderPath, $modFile1."_exclusive.txt");
	print($csumPath,"\n",$tpPath,"\n",$fpPath,"\n",$fnPath);
	
	#mkdir($resFolderPath) or die ($!);
	
	my $flag = 0;
	
	
	my $file1_hash = {};
	my $combined_hash = {};
	open(vcf1,$file1) or die($!);
	open(vcf2,$file2) or die($!);

	while(my $line = <vcf1>){
	
		if($line=~/^\#\#/){
			next;
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t/;
			
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//g;
			$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FN";
		}
		
	}
	close(vcf1);
	my $flag = 0;
	while(my $line = <vcf2>){
		if($line=~/^\#\#/){
			next
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+/;
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//;
			if(defined($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall})){
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="TP";
			}
			else{
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FP";
			}
		}
		
	}
	my $fn=0;
	my $tp=0;
	my $fp=0;
	
	#print "$fname1[0]\_$fname2[0]\_common.txt";
	#print "@fname1\n";
	#print "Manoj\n";
	#print "$file1 $file2";
	open(CSUM,">>","$csumPath") or die("$!");
	open(TP,">","$tpPath") or die("$!");
	open(FP,">","$fpPath") or die("$!");
	open(FN,">","$fnPath") or die("$!");
	print TP "CHROM\tPOS\tREF\tALT\n";
	print FP "CHROM\tPOS\tREF\tALT\n";
	print FN "CHROM\tPOS\tREF\tALT\n";
	foreach my $chrom(keys %{$file1_hash}){
		foreach my $pos(keys %{$file1_hash->{$chrom}}){
			foreach my $refall(keys %{$file1_hash->{$chrom}->{$pos}}){
				foreach my $altall(keys %{$file1_hash->{$chrom}->{$pos}->{$refall}}){
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "TP"){
						print TP "$chrom\t$pos\t$refall\t$altall\n" ;
						$tp++;
					}
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FP"){
						print FP "$chrom\t$pos\t$refall\t$altall\n";
						$fp++;
					}
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FN"){
						print FN "$chrom\t$pos\t$refall\t$altall\n" ;
						$fn++;
					}
				}
			}
		}
	}
	print CSUM "$fileName1\t$fileName2\t$tp\t$fn\t$fp\n";
	close(TP);
	close(FP);
	close(CSUM);
	close(FN);
}

sub compareJobFoldersOrFile{
	my $folder1 = shift;
	my $folder2 = shift;
	my $resFolderPath = shift;
	my $file1 = "";
	my $file2 = "";
	
	if(-d $folder1 and -e $folder1){
		 $file1 = getVcfPathFromFolder($folder1);
	}else{
		$file1  = $folder1;
	}
	if(-d $folder2 and -e $folder2){
		 $file2 = getVcfPathFromFolder($folder2);
	}else{
		$file2  = $folder2;
	}

	my $fileName1 = basename($file1);
	my $fileName2 = basename($file2);
	
	my $modFile1 = $fileName1=~ s/\./_/rg;
	my $modFile2 = $fileName2=~ s/\./_/rg;

	#Create the output file paths
	my $csumPath = File::Spec->catfile( $resFolderPath, "countSummary.txt" );
	#my $tpPath = File::Spec->catfile( $resFolderPath, $modFile1."_".$modFile2."_common.txt" );
	my $tpPath = File::Spec->catfile( $resFolderPath, "common.txt" );
	my $fpPath = File::Spec->catfile( $resFolderPath, $modFile2."_exclusive.txt" );
	my $fnPath = File::Spec->catfile( $resFolderPath, $modFile1."_exclusive.txt");
	
	

	my $flag = 0;
	
	
	my $file1_hash = {};
	my $combined_hash = {};
	open(vcf1,$file1) or die($!);
	open(vcf2,$file2) or die($!);

	while(my $line = <vcf1>){
	
		if($line=~/^\#\#/){
			next;
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t/;
			
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//g;
			$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FN";
		}
		
	}
	close(vcf1);
	my $flag = 0;
	while(my $line = <vcf2>){
		if($line=~/^\#\#/){
			next
		}
		if($line=~/^\#CHROM/){
			$flag = 1;
			next;
		}
		if($flag == 1){
			$line =~/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+/;
			my $chrom = $1;
			my $pos = $2;
			my $refall = $3;
			my $altall = $4;
			$chrom =~s/chr//;
			if(defined($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall})){
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="TP";
			}
			else{
				$file1_hash->{$chrom}->{$pos}->{$refall}->{$altall}="FP";
			}
		}
		
	}
	my $fn=0;
	my $tp=0;
	my $fp=0;
	
	#print "$fname1[0]\_$fname2[0]\_common.txt";
	#print "@fname1\n";
	#print "Manoj\n";
	#print "$file1 $file2";
	open(CSUM,">>","$csumPath") or die("$!");
	open(TP,">","$tpPath") or die("$!");
	open(FP,">","$fpPath") or die("$!");
	open(FN,">","$fnPath") or die("$!");
	print TP "CHROM\tPOS\tREF\tALT\n";
	print FP "CHROM\tPOS\tREF\tALT\n";
	print FN "CHROM\tPOS\tREF\tALT\n";
	foreach my $chrom(keys %{$file1_hash}){
		foreach my $pos(keys %{$file1_hash->{$chrom}}){
			foreach my $refall(keys %{$file1_hash->{$chrom}->{$pos}}){
				foreach my $altall(keys %{$file1_hash->{$chrom}->{$pos}->{$refall}}){
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "TP"){
						print TP "$chrom\t$pos\t$refall\t$altall\n" ;
						$tp++;
					}
					
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FP"){
						print FP "$chrom\t$pos\t$refall\t$altall\n";
						$fp++;
					}
					if($file1_hash->{$chrom}->{$pos}->{$refall}->{$altall} eq "FN"){
						print FN "$chrom\t$pos\t$refall\t$altall\n" ;
						$fn++;
					}
				}
			}
		}
	}
	print CSUM "$fileName1\t$fileName2\t$tp\t$fn\t$fp\n";
	close(TP);
	close(FP);
	close(CSUM);
	close(FN);
}




1;
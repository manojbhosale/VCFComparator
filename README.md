# VCFComparator
NGS VCF file comparison(diff) script

## Usage

> perl compare_cmdline.pl -file1 vcfFile1.vcf -file2 vcfFile2.vcf

### Results
The utility will generate below listed file in **/_Comparison_results_** directory
+ countSum.txt
+ <file1>_<file2>_common.txt
+ <file1>_exclusive.txt
+ <file2>_exclusive.txt


Note: The **_.pl_** file and the **_.pm_** file should be in same directory.

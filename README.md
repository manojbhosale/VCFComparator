# VCFComparator
NGS VCF file comparison(diff) script

## Usage

>perl compare_cmdline.pl -F1 <folder or file path> -F2 <folder or file path> -D <optional destination directory>

If -D is not specified then current working directory will be used as destination. -D will create folder with name “Comparison_results” at the destination location in which the result files will be generated.
Also, make sure that the “Comparison_results” is not already present in the destination directory otherwise it will fail.

### Results
The utility will generate below listed file in **/_Comparison_results_** directory
+ countSum.txt
+ <file1>_<file2>_common.txt
+ <file1>_exclusive.txt
+ <file2>_exclusive.txt


Note: The **_.pl_** file and the **_.pm_** file should be in same directory.

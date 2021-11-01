# Hairpin Analyzer
The Hairpin Analyzer takes the output files of the [BiQ Analyzer HT](http://biq-analyzer.bioinf.mpi-inf.mpg.de) and scans them. The following output files are produces:

 - One heatmap as a Portable network graphic (png) for each input file
 - One text file with mapped CpG positions and further data for each input file
 - One text file with result summary for each set of files having the same 
   amplicon type.

**ATTENTION:**
When mapping the CpG files with the non-cpg files, the ids in the non-cpg files do not contain the string "CG" but instead the string "NN". The Hairpin Analyzer contains a dirty work-around in which "NN" in the ids from the non-cpg files is replaced by "CG" to match the ids in the CpG files. The correct solution would be to replace the "CG"s in the CpG files by "NN"s and compare these with the non-cpg files. But since we need the CpG file ids for comparison with the linker and SNP-files, this simple work-around was implemented. Feel free to fix this if you need to.

The file "[Reference]_[Sample]_results.txt" contains the processed information  of the "results.tsv" file in the subfolder "[Reference]/[Sample]HP/". The following columns are contained:

### ID
The column "Id" of the analyzed sequence copied from the CpG data file

### Sequence_Identity
The column "Sequence Identity" copied from the CpG data file

### Mapped CpG Methylation Pattern
The mapped CpG methylation pattern. The Hairpin Analyzer takes the column CG_Methylation_Pattern from the CpG data file, cuts specified positions, folds the methylation pattern in the middle, such that always two positions map on each other. The resulting information is coded as follows:
```
0 - both sides are not methylated (00)
1 - the left  side of the original data is methylated, the other one not (10)
3 - the right side of the original data is methylated, the other one not (01)
4 - both sides are methylated (11)
5 - one side is unmethylated, the other one is mutated (x0 or 0x)
6 - the left  side is methylated, the other one is mutated (1x)
8 - the right side is methylated, the other one is mutated (x1)
9 - both sides are mutated (xx)
```

### Conversion Read	
The column "Conversion" copied from the CpG data file

### Conversion Linker
This value is calculatd by 1 - the column "Mean_CG_Methylation" from the  linker data file

### Reference
The column "Reference" copied from the CpG data file (e.g. "9.5dpc")

### Sample
The column "Sample" copied from the CpG data file (e.g. "L1HP")

### l1, l2, ...
The column "CG_Methylation_Pattern" from the linker file, one column for each character. These columns are not present if there was no linker file found.

### c1, c2, ...
The column "CG_Methylation_Pattern" from the "non CpGs" file, one column for each character. These columns are not present if there was no "non CpGs" file found.

### snp1, snp2, ...
The column "CG_Methylation_Pattern" from the SNPs file, one column for each character. These columns are not present it there was no SNP file found.

The file "achieved_results_[Sample].txt" contains a summary of the processed 
information of all "results.tsv" having [Sample] as sample. The 
following columns are contained:

### Reference
The column "Reference" copied from the CpG data file.

### Amplicon
The column "Sample" copied from the CpG data file.

### reads
How many lines have been in the CpG data file and remained after all filtering (all these lines can be found in the heatmap).

### CpG
How many mapped CpG positions are there in the result file having as character 4, 3, 1 or 0.

### 4
How many mapped CpG positions with 4 as result character (methylated at both sides) are there in the result file.

### 3
How many mapped CpG positions with 3 as result character (no methylation at left side, methylated at right side) are there in the result file.

### 1
How many mapped CpG positions with 1 as result character (methylated at left side, no methylation at right side) are there in the result file.

### 0
How many mapped CpG positions with 0 as result character (no methylation at both sides) are there in the result file.

### na
How many mapped CpG positions are there in the result file having as character 5, 6, 8 or 9.

### 5
How many mapped CpG positions with 5 as result character (one side not methylated, the other side mutated) are there in the result file.

### 6
How many mapped CpG positions with 6 as result character (left side methylated, right side mutated) are there in the result file.

### 8
How many mapped CpG positions with 8 as result character (left side mutated, right side methylated) are there in the result file.

### 9
How many mapped CpG positions with 8 as result character (left side mutated, right side methylated) are there in the result file.

### m04
How many fully methylated lines (only 4 as result character) are there in the result file. Mutations (x) are ignored.


### m043_041
How many fully and halfly methylated lines (4 and 3 or 4 and 1 as result character) are there in the result file. Mutations (x) are ignored.

### m01
How many halfly left methylated lines (1 or 0 as result character) are there in the result file. Mutations (x) are ignored.

### m03
How many halfly right methylated lines (3 or 0 as result character) are there in the result file. Mutations (x) are ignored.

### m013
How many halfly left and halfly right methylated lines (1, 3 and 0 as result character) are there in the result file. Mutations (x) are ignored.

### m0134
How many fully, halfly left and halfly right methylated lines (4, 1, 3 and 0 as result character) are there in the result file. Mutations (x) are ignored.

### m0
How many completely unmethylated lines (only 0 as result character) are there in the result file. Mutations (x) are ignored.

### only_one_methylation
How many lines are there which have exactly one methylation on one side and zero or one on the other (e.g. 040, 030, 010, 031). Mutations are ignored.

### mosaic_pattern
How many lines are there which have a mosaic methylation pattern on at least one side (e.g. 404, 101, 303, 104, 401, 304 and 403, but not 301 or 103 (number of zeros can be more than one)). Mutations are ignored.

### continuous_methylation
How many lines are there which have a continuous methylation pattern on at least one side (e.g. 0440, 0110, 0330, 0140, 0410, 0340 and 0430 but not 0130 or 0310 (number of follow-up methylations has to be at least 2)). Mutations are ignored.

### Data<->Linker: Conversion rate linker >= 80%)
How many lines which are present in the CpG data as well as in the linker data have a conversion rate bigger than or equal to 80%?

### Conversion rate linker
The column "Conversion" copied from the linker file.

### Linker under threshold
From all lines in the CpG data file that have a matching line in the linker file, the percentage that was under the threshold (80%).

### Linker not available
From all lines in the CpG data file, the percentage that doesn't have a matching line in the linker file.

### SNPs perfect match
The percentage of the lines where the matching SNPs file shows a perfect match when mapping the SNP positions in the middle.

### l1, l2, ...
For each linker position, the percentage of ones out of the set {0,1} (ignoring the mutations).

### c1, c2, ...
For each "non CpG" position, the percentage of ones out of the set {0,1} (ignoring the mutations).

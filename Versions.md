## Hairpin Analyzer version 0.3   (23. August 2011)
 - changed order of results in results_summary.txt (m0 column moved 
   3 positions to the right)
 - made configurable: how many digits after the decimal point should be output 
   for percentages
 - added TestSuite with Unit Tests (classes/HairpinTest.py)
 - seperate summary file is created for each amplicon type
 - scaning linker files, only takes lines for results which have a mach
   in data files. for each line added the following information:
   - conversion rate of linker
   - linker methylation pattern (l0, l1, l2, ...)
   In the summary file, added the folling information
   - conversion rate for all mapped lines where Mean methy. level
     is above threshold (currently 0.8)
   - conversion rate linker (from summary.dat file)
   - Percentage of all paired lines (data <--> linker) that were 
     under the threshold
   - Percentage of all data lines for which no linker lines were available
   - Percentage of lines with perfect SNP match
 - scan files with "non CpG" positions
 - scan files with SNPs, map SNPs in the middle, calculate percentage 
   of correct mappings

## Hairpin Analyzer version 0.2   (8. August 2011)
 - showing process of analysis with progress bar
 - recognize wether running on windows or unix system
 - calculate the following additional statistics:
   - number of rows with exactly one methylation on one strand and
     zero or one on the other
   - number of rows with mosaic methylation pattern on one of the strands
   - number of rows with continuous methylation pattern and without 
     mosaic pattern

## Hairpin Analyzer version 0.1   (8. August 2011)
 - map CpG positions and calculate result values
 - sort result values according to specific order
 - produce text output for each input file
 - produce heatmap as png file for each input file
   - each heatmap with same size (result columns are scaled)
   - each heatmap with additional summary column on left side (1/5 image width)
 - produce results_summary with statistics on all input files
 - configurable:
   - path to data
   - heatmap size
   - filename for results_summary.txt
   - amplicons where middle position should be removed before mapping
   - amplicons where special positions should be removed before mapping
   - amplicoins where special positions should be removed after mapping

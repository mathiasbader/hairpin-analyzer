# How to
The Hairpin Analyzer is a python based script which takes the output files of the [BiQ Analyzer HT](http://biq-analyzer.bioinf.mpi-inf.mpg.de) and restores the double strand information. To run HairpinAnalyzer, the user has to install python and provide the following output folders of BiQHT:

**Mandatory folder:**
* `CpGs`

**Optional folders:**
* `nonCpGs`    provides nonCpG methylation data
* `Conversion` provides conversion rate of cytosine
* `SNPs`       outputs the base at a given SNP

A detailed description, on how the individual folders are generated, can be found in *Giehr at al. 2018*. Make sure that all folders are stored in one main folder:

    BiQHToutputProjectXY/CpGs
    BiQHToutputProjectXY/nonCpGs
    BiQHToutputProjectXY/Conversion
    BiQHToutputProjectXY/SNPs

We recommend using the given folder names to minimize the adjustments on the configuration file `configuration.py`.

1. Copy the folder HairpinAnalyzer to the desired location
2. Open `configuration.py`. Go to the section `Additional data paths` and give the path to the main folder under the point `data_path_main_folder = `
   
    Example:
```
data_path_main_folder = /data/BiQHToutputProjectXY/
```
3. Save and close the configuration file.
4. Open command line and access the HairpinAnalyzer folder:

Example:
```
cd tools/HairpinAnalyzer
```

5. Run HairpinAnalyzer:
```
python HarpinAnalyzer.py
```

6. The results will be stored as individual files (txt and png) in the folder `Hairpinanalyzer/results/`. In addition, all results from one reference sequence will be stored in form of a summary file ‘achieved_results.txt’ For a detailed description of the output files see `documentation.txt`.

The HairpinAnalyzer always expects the same number of CpGs at upper and lower DNA strand. However, depending on the used restriction enzyme or primer placement the number of CpGs from upper and lower DNA strand can differ. In such a case, the user can declare that specific CpG position should be ignored or rather omitted when restoring the double strand information.

For this, go to the section ’**Deleting CpG positions before mapping**’ in the file `configuration.py`. Three option are possible:

1. *remove_middle_position*: Use in case the restriction site contains CpG positions i.e. TaqI T|CGA. Here the CpG position will later only be present at one DNA strand.
```
NN-CG-NN-NN-NN-TG-Linker-CG-NN-NN-CG-NN-NN

remove_middle_position = [‘name amplicon 1’, ‘name amplicon 2’] 
```

2. *remove_special_positions*: Individual CpG positions can be removed i.e. in case of SNPs or asymmetrical CPG distribution due to primer placement.

```
NN-CG-NN-CG-NN-CG-NN-Linker-NN-NN-CG-NN-CG-NN    --> remove first CpG for even numbers of CpGs (example 1)
NN-CG-NN-CG-NN-NN-Linker-NN-NN-CG-NN-CG-NN-CG-CG --> remove position four and five for even numbers of CpGs (example 2)

remove_special_positions = [(‘name amplicon example 1’, [1]), (‘name amplicon example 2’, [4,5])]
```

3. *delete_mapped_columns*: Removes !CpG dyads! from the pattern map and also from the data files; meaning the corresponding CpGs from upper and lower strand. This function might be used in case of SNPs resulting in TpGs or TpAs instead of CpGs.

```
NN-CG-NN-CG-NN-NN-Linker-NN-NN-CG-NN-CG-NN 

delete_mapped_columns = [(‘name amplicon example 3’, [1]), (‘name amplicon example 4’, [2])]
```
Example 3: Removes the first `!CpG dyad!` from the Analysis, meaning `CpG 1` and `CpG 4`

Example 4: Removes the first `!CpG dyad!` from the Analysis, meaning `CpG 2` and `CpG 3`

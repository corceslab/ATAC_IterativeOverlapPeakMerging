[![DOI](https://zenodo.org/badge/416004561.svg)](https://zenodo.org/badge/latestdoi/416004561)

# ATAC_IterativeOverlapPeakMerging

This repository is home to a single R script used for performing the iterative overlap peak merging procedure. This was originally described in [Corces & Granja et. al. Science 2018](https://pubmed.ncbi.nlm.nih.gov/30361341/) and details can be found in the Supplementary Materials in the "Materials and Methods" section under the heading "Peak Calling".

The provided script was written by Jeff Granja and Ryan Corces

## Description
This R script was written to generate merged peak sets for large-scale ATAC-seq experiments. Given metadata
describing how samples are grouped, this script creates a merged peakset for each group and then a single
merged peak set for all samples. Even though the script is designed for use on large-scale projects, it can and should
also be used for small-scale projects. Inputs are described below. Script is meant to be run from the
command line on a UNIX-based system like so:

```
Rscript createIterativeOverlapPeakSet.R
  --metadata /path/to/metadata/file.txt
  --macs2dir /path/to/macs2/dir/with/summit_files/
  --outdir /path/to/the/desired/output/dir/
  --suffix _summits.bed
  --blacklist /path/to/the/blacklist/file.bed
  --genome hg38
  --spm 5
  --rule "(n+1)/2"
  --extend 250
```

In creating this iterative overlap merged peak set, this script also (i) filters for peak significance, (ii) removes "cliffed" peaks that extend past the ends of chromosomes, (iii) removes peaks that map to "chrY", and (iv) removes peaks that span a region in the reference genome where an "N" (unknown) base is present.

Because of some assumptions made by the code, reference genomes other than mouse and human may encounter issues. As we dont work with data from these organisms, we cannot support code tailored to their analysis.

## Requirements
This script requires the following R packages to be installed:
1. optparse
2. dplyr
3. yaml
4. rtracklayer
5. Biostrings
6. SummarizedExperiment
7. GenomeInfoDb
8. GenomicRanges >1.44.0
9. readr
10. edgeR
11. BSgenome

Additionally, you will need to install the BSGenome package for your genome of interest. For example for the UCSC hg38 reference genome:
```
library(BiocInstaller)
biocLite("BSgenome.Hsapiens.UCSC.hg38")
```

## Input Files and Parameters
##### MACS2 Summit Files
This script used the "summits" files that are provided by MACS2 when the `--call-summits` parameter is used. For ATAC-seq, BAM files used as input to MACS to should be adjusted for the Tn5 offset as described in our [Nature Protocols]() publication. We recommend running MACS2 with the following parameters:
`--shift -75 --extsize 150 --nomodel --call-summits --nolambda --keep-dup all -p 0.01`
All input summit files are expected to be located in the same directory.

##### Metadata file
A 2-column tab-delimited text file containing the metadata for the project as "_sample_name_ __tab__ _group_" 
where `sample_name` is expected to be the prefix of a summit file in the format "_prefixsuffix_" and
group is an identifier that can be shared among multiple sample_names and will be used to create group-specific peak sets. See below for a description of `suffix`.
Duplicate sample_name entries are not allowed.
The first line of the metadata file is expected to have the column names. The second column should be "Group". This is case sensitive. The first column name does not matter.

##### Blacklist file
A blacklist is a comprehensive set of regions in the genome that have anomalous, unstructured, or high signal in next-generation sequencing experiments independent of cell line or experiment. You can read more about blacklists and download one for your reference genome [here](https://github.com/Boyle-Lab/Blacklist). Please note that these are gzipped bed files and this script expects an uncompressed BED file as input.

#### Parameters
1. `metadata` - This is the full path to the metadata file. For example `/path/to/metadata/file.txt`. See above for formating guidelines. No default value.
2. `macs2dir` - This is the full path to the directory containing the summit files and is expected to have "/" on the end. For example `/path/to/macs2/dir/with/summit_files/`. No default value.
3. `outdir` - This is the full path to the desired output directory and is expected to have "/" on the end. For example `/path/to/the/desired/output/dir/`. No default value.
4. `suffix` - This is a character string of the suffix of MACS2 summit files (typically `_summits.bed`). Peak
     callers other than MACS2 are not intended to work with this script. No default value.
5. `blacklist` - This is the full path to a blacklist bed file. For example `/path/to/the/blacklist/file.bed`. No default value.
6. `genome` - This is the shorthand genome name from BSGenome (i.e. "hg38"). Can be determined using the "provider_version"
    column from the output of "available.genomes(splitNameParts=TRUE)". Due to changes in BSGenome, these exact commands may not work for you depending on your version. No default value.
7. `spm` - The score-per-million cutoff to be used in determining which peak calls represent true peaks. The higher the `spm` cutoff, the more significant a peak call has to be to be retained. We recommended setting `spm` to between 2 and 5 as a starting point. Default value is 5.
8. `rule` - This is a character string denoting the selection rule for peak filtering. For example, "2" means at least 2 samples
    must have a peak overlapping the given summit location. Any overlap (i.e. 1 bp) is sufficient. Similarly, "(n+1)/2" means that a majority of samples must have a peak overlapping the given summit location. This rule must be in quotes. Default value is 2. 
9. `extend` - An integer of how many base pairs on either side of the summit to extend the peak. For example, if set to 250, then you will create 501-bp peaks containing the summit position and 250 bp on either side of the summit. Default value is 250.

## Example Usage
To assist with usage, we provide an example metadata file that accompanies the below description.

Lets say we have 3 cell types, named `CellTypeA`, `CellTypeB`, and `CellTypeC`. Each cell type has three replicates named `Rep1`, `Rep2`, and `Rep3`. This would give us 9 summit files from MACS2. We will pretend these files are located in `/macs2/peakCalls/` like so:

```
CellTypeA_Rep1_summits.bed
CellTypeA_Rep2_summits.bed
CellTypeA_Rep3_summits.bed

CellTypeB_Rep1_summits.bed
CellTypeB_Rep2_summits.bed
CellTypeB_Rep3_summits.bed

CellTypeC_Rep1_summits.bed
CellTypeC_Rep2_summits.bed
CellTypeC_Rep3_summits.bed
```

In this example, the `suffix` of our files is `_summits.bed` because each of these files ends with `_summits.bed`. Thus, our prefixes are everything before `_summits.bed`. For example "CellTypeA_Rep1".

From this, we would make our metadata file. Because we have three cell types, we will create three groups. This metadata file is tab-delimited and has two columns. The first line of the second column must be "Group" (case sensitive). The first column contains the sample prefixes (i.e. the file names minus the suffix).
```
Sample	Group
CellTypeA_Rep1	CellTypeA
CellTypeA_Rep2	CellTypeA
CellTypeA_Rep3	CellTypeA
CellTypeB_Rep1	CellTypeB
CellTypeB_Rep2	CellTypeB
CellTypeB_Rep3	CellTypeB
CellTypeC_Rep1	CellTypeC
CellTypeC_Rep2	CellTypeC
CellTypeC_Rep3	CellTypeC
```

If we set `rule = 2` in this example, then two of the three replicates will be required to have a peak call overlapping a given summit. If we set `rule = "(n+1)/2"`, then the outcome will be exactly the same because the majority here is 2. However, you can imagine how these different values for `rule` would have very different outcomes in experiments with more replicates. 

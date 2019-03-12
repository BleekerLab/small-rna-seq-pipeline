# small_rna_landscape
A pipeline to annotate small RNA loci (miRNAs, phased siRNAs) using a reference genome and experimental sRNA-Seq datasets.

This [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) pipeline is meant to characterize in details small RNA sequencing files.
*  It  will generate statistics on the number of small RNAs that aligns to different databases (phytoviruses, RFAM, repeats, transposons, etc.)
*  It will characterize the type of small RNAs that can be found in the samples (microRNAs, siRNAs, phased siRNAs, etc.) using specialized programs like [ShortStack](https://github.com/MikeAxtell/ShortStack)

# How to use it

## Install all softwares and packages needed with the [Conda package manager](https://conda.io/docs/using/envs.html)

### Create a virtual environment named "srna" from the `environment.yaml` file
conda env create --name srna --file environment.yaml

### Activate this virtual environment
source activate srna

## Usage
If Snakemake is properly set and if the pipeline is called `Snakefile` then just type `Snakemake` and it will run.
You can try a dry run (nothing happens, it checks the different steps and print the shell command) by typing `Snakemake -pn`


# Main outputs:
*  Number of reads mapped at each step (before and after trimming, spikes, virus, rRNA, etc.)  
*  Length distribution of small RNAs  
*  Log files of Bowtie mapping steps  
*  A complete ShortStack analysis that includes:  
    * a fasta file of non-redundant miRNAs
    * a *de novo* analysis of small RNA loci (miRNAs, phased siRNAs)
*  BigWig coverage files (indexed binary format of bedGraph files)

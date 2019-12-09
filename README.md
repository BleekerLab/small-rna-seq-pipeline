# The small RNA-Seq pipeline 

- [Summary](#summary)
- |[Installation](#installation)
  * [Create a Conda environment](#create-a-conda-environment)
  * [Dependencies](#dependencies)
- [Usage](#usage)
  * [Example](#example)
  * [Samples](#samples)
  * [Configuration](#configuration)
  * [Genomic references](#genomic-references)
- [Authors](#authors)
  * [Contributors](#contributors)
  * [Maintainers](#maintainers)
- [Citation](#citation)
- [License](#license)
- [Versioning](#versioning)
- [Acknowledgments](#acknowledgments)
- [References](#references)

## Summary
The small RNA-Seq description pipeline is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to annotate small RNA loci (miRNAs, phased siRNAs) using one or more reference genomes and based on experimental small RNA-Seq datasets.  
This pipeline heavily relies on the [ShortStack](https://github.com/MikeAxtell/ShortStack) software that annotates and quantifies small RNAs using a reference genome.  

Upon completion, several outputs will be generated for each sample:
- One Shortstack result file called `Results.txt`. See the description of this file in the [Shortstack manual](https://github.com/MikeAxtell/ShortStack).
- Two fasta files for each sample: one fasta file containing the predicted hairpins and one containing the predicted mature microRNAs.
- Two blast result files (in tabular format) based on the blast of predicted hairpins and mature miRNAs against mirbase (the version of miRBase is specified in the config file). See the [miRBase website](http://www.mirbase.org/) for releases.

## Installation

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Create a Conda environment
This Snakemake pipeline make use of the [conda package manager](https://docs.conda.io/en/latest/) to install softwares and dependencies.
1. First, make sure you have conda installed on your system. Use [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) and follow the [installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).  
2. Using `conda`, create a virtual environment called `snakemake` to install Snakemake (version 5.4.3 or higher) by executing the following code in a Shell window: `conda env create -f envs/snakemake.yaml`. This will install `snakemake version 5.4.3` and `pandas version 0.25.0` in a new environment called __snakemake__.
3. Activate this environment using: `source activate snakemake`
4. You can now run the pipeline (see below) since snakemake will use conda to install softwares and packages for each rule.  

If you have set up `conda` and created the `snakemake` environment, that's all you need to do! Snakemake will take care of the rest of the software and package installation specified in the _yaml_ files in the `envs/` folder.

### Dependencies

* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - The Snakemake workflow management system is a tool to create reproducible and scalable data analyses.
* [NCBI blast+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) - A program to perform sequence similarity search. See [NCBI Blast webpage for more info](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
* [ShortStack](https://github.com/MikeAxtell/ShortStack) - Small RNA loci annotation and quantification.
* [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) - Read trimming for NGS data. 
* [bioawk](https://github.com/lh3/bioawk) -  Bioawk is an extension to Brian Kernighan's awk, adding the support of several common biological data formats, including optionally gzip'ed BED, GFF, SAM, VCF, FASTA/Q and TAB-delimited formats with column names. 

A series of custom Python functions are also used and can be found in the `helpers.py` file.  
Versions of softwares and packages can be seen in their respective environment `.YAML` file in the `envs/` folder.


## Usage

### Example 
A small dataset is available in `test/` to run some tests rapidly. It will use the genome and miRBase reference fasta files stored in `refs/`.  
To run the test, open a new Shell window and:
1. Activate your working environment: `source activate snakemake`
2. Type `snakemake --use-conda -np` for a dry run. No analysis is run but it checks that the Directed Acyclic Graph of jobs is OK (input and output from each rule chained to each other).
3. For the real run, type `snakemake --use-conda`

### Samples
A `samples.tsv` file can be used to specify sample names, their corresponding genomic reference to use and the location of their sequencing file.

### Configuration
Configuration settings can be changed in the `config.yaml` file. For instance, one could modify the minimal coverage required by Shorstack to discover sRNA loci.  

### Genomic references
Different genomic references can be used for each sample. Simply provide a genomic reference corresponding to your sample.

## Authors

### Contributors
* **Marc Galland** - *Initial work* - [Github profile](https://github.com/mgalland)
* **Michelle van der Gragt** - *Initial work* - [Github profile](https://github.com/MvanderGragt)

### Maintainers
* **Marc Galland** - *Initial work* - [Github profile](https://github.com/mgalland)

## Citation
...as soon as we have published this software!

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Versioning

[SemVer](http://semver.org/) is used for versioning. For the versions available, see the [releases on this repository](https://github.com/BleekerLab/small-rna-seq-pipeline/releases).

## Acknowledgments

* [RNA Biology & Applied Bioinformatics group](http://sils.uva.nl/content/research-groups/rna-biology--applied-bioinformatics/rna-biology--applied-bioinformatics.html)
* [Mike Axtell at the Penn State University](https://bio.psu.edu/directory/mja18)

## References
* Bioawk tutorial: https://isugenomics.github.io/bioinformatics-workbook/Appendix/bioawk-basics

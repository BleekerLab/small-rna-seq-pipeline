# The small RNA landscape pipeline
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline to annotate small RNA loci (miRNAs, phased siRNAs) using a reference genome and experimental sRNA-Seq datasets.  
This pipeline heavily relies on the [ShortStack](https://github.com/MikeAxtell/ShortStack) software that annotates and quantifies small RNAs using a reference genome.  
Upon completion, several outputs will be generated for each sample:
- Number of small RNA clusters per predominant sRNA length ("Dicer Call")
- Abundance of small RNA clusters per predominant sRNA length ("Dicer Call")
- Proportion of sRNA classes per sample (pie chart)
- Number of MIR genes per family (barplot)


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites
This Snakemake pipeline make use of the [conda package manager](https://docs.conda.io/en/latest/) to install softwares and dependencies.
1. First, make sure you have conda installed on your system. Use [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) and follow the [installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).  
2. Using `conda`, create a virtual environment called `snakemake` to install Snakemake (version 5.4.3 or higher) by executing the following code in a Shell window: `conda create --name snakemake -c bioconda snakemake=5.4.3 pandas=0.24.2`. This will install `snakemake version 5.4.3` and `pandas version 0.24.2` in a new environment called __snakemake__.
3. Activate this environment using: `source activate snakemake`
4. You can now run the pipeline (see below) since snakemake will use conda to install softwares and packages for each rule.  


### Installing
If you have set up `conda` and installed `snakemake` in your environment, that's all you need to do! Snakemake will take care of the rest of the software and package installation specified in the _yaml_ files in the `envs/` folder.


## Running the tests
A small dataset is available in `test/` to run some tests rapidly.   

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment
You can deploy Snakemake on your local machine or on more powerful cluster systems.
On the SURFsara LISA cluster, ...


Add additional notes about how to deploy this on a live system

## Built With

* [Snakemake](https://snakemake.readthedocs.io/en/stable/) - The Snakemake workflow management system is a tool to create reproducible and scalable data analyses.
* [ShortStack](https://github.com/MikeAxtell/ShortStack) - Small RNA loci annotation and quantification.
* [R](https://www.r-project.org/) - Used to generate plots.

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/BleekerLab/small_rna_landscape).

## Authors

* **Marc Galland** - *Initial work* - [Github profile](https://github.com/mgalland)
* **Michelle van der Gragt** - *Initial work* - [Github profile](https://github.com/MvanderGragt)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [RNA Biology & Applied Bioinformatics group](http://sils.uva.nl/content/research-groups/rna-biology--applied-bioinformatics/rna-biology--applied-bioinformatics.html)

## References and link
* Bioawk tutorial: https://isugenomics.github.io/bioinformatics-workbook/Appendix/bioawk-basics

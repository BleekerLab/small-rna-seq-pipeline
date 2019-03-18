import os
import subprocess
from snakemake.utils import validate, min_version
#import pandas as pd


##### set minimum snakemake version #####
min_version("5.4.3")

#########################
## Pipeline configuration
##########################
configfile: "config.yaml"

wildcard_constraints:
  dataset="[Aa-Zz0-9]+"

# directories
FQ_DIR = config["fastqdir"]
WORKING_DIR = config["temp_dir"]
RES_DIR = config["result_dir"]

# reference genome
GENOME = config["refs"]["genome"]

# ShortStack parameters
SHORTSTACK_PARAMS = " ".join(config["shortstack"].values())

####################
## Desired outputs
####################

# ShortStack
SHORTSTACK = expand(RES_DIR + "shortstack/{sample}/Results.txt",sample=config["samples"])
PLOTS = [RES_DIR + "plots/cluster_abundance_per_dicercall.png",RES_DIR + "plots/abundance_of_clusters_per_dicer_call.png"]

rule all:
    input:
        SHORTSTACK,
        PLOTS
    output:
    message:"All done! Removing intermediate files"
    shell:
        "rm -rf {WORKING_DIR}"

rule make_plots:
    input:
        expand(RES_DIR + "shortstack/{sample}/Results.txt",sample=config["samples"])
    output:
        RES_DIR + "plots/cluster_abundance_per_dicercall.png",
        RES_DIR + "plots/abundance_of_clusters_per_dicer_call.png"
    message: "making plots based on ShortStack results"
    conda:
        "envs/plots.yaml"
    params:
        shortstack_resdir = RES_DIR + "shortstack/"
    shell:
        "Rscript --vanilla scripts/number_srna_clusters_per_dicer_call.R {params} {RES_DIR} " # plot 1
        "Rscript --vanilla scripts/abundance_of_clusters_per_dicer_call.R {params} {RES_DIR} " # plot 2

rule shortstack:
    input:
        reads =  WORKING_DIR + "trim/{sample}.trimmed.size.fastq",
        genome = GENOME
    output:
        RES_DIR + "shortstack/{sample}/Results.txt",
        RES_DIR + "shortstack/{sample}/{sample}.trimmed.size.bam"
    message:"Shortstack analysis of {wildcards.sample} using {input.genome} reference"
    params:
        RES_DIR + "shortstack/{sample}/"
    threads: 10
    conda:
        "envs/shortstack.yaml"
    shell:
        "perl bin/ShortStack.pl "
        "--outdir {wildcards.sample} "
        "--bowtie_cores {threads} "
        "--sort_mem 4G "
        "{SHORTSTACK_PARAMS} "
        "--readfile {input.reads} "
        "--genome {input.genome};"
        "cp -r {wildcards.sample}/* {params};"
        "rm -r {wildcards.sample};"

#############################
## Trim reads for all samples
#############################
rule keep_reads_shorter_than:
    input:
        WORKING_DIR + "trimmed/{sample}.trimmed.fastq"
    output:
        WORKING_DIR + "trim/{sample}.trimmed.size.fastq"
    message: "Discarding reads longer than {params.max_length} nucleotides"
    params:
        max_length = config["length"]["max_length"]
    conda:
        "envs/bioawk.yaml"
    shell:
        """
        bioawk -c fastx '{{ length($seq) <= {params.max_length} }} {{print "@"$name; print $seq ;print "+";print $qual}}' {input} > {output}
        """


rule trimmomatic:
    input:
        FQ_DIR + "{sample}.fastq"
    output:
        WORKING_DIR + "trimmed/{sample}.trimmed.fastq",
    message: "trimming {wildcards.sample} on quality and length"
    log:
        RES_DIR + "logs/trimmomatic/{sample}.log"
    params :
        LeadMinTrimQual =           str(config['trimmomatic']['LeadMinTrimQual']),
        TrailMinTrimQual =          str(config['trimmomatic']['TrailMinTrimQual']),
        windowSize =                str(config['trimmomatic']['windowSize']),
        avgMinQual =                str(config['trimmomatic']['avgMinQual']),
        minReadLen =                str(config['length']['min_length']),
        phred = 		            str(config["trimmomatic"]["phred"])
    threads: 10
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "trimmomatic SE {params.phred} -threads {threads} "
        "{input} "
        "{output} "
        "LEADING:{params.LeadMinTrimQual} "
        "TRAILING:{params.TrailMinTrimQual} "
        "SLIDINGWINDOW:{params.windowSize}:{params.avgMinQual} "
        "MINLEN:{params.minReadLen} &>{log}"

#####
## QC
#####

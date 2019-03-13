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

GENOME = config["refs"]["genome"]

###########
## Programs
###########

# ShortStack parameters
SHORTSTACK_PARAMS = " ".join(config["shortstack"].values())

###################
## Other parameters
###################

####################
## Desired outputs
####################

# ShortStack
SHORTSTACK = expand(RES_DIR + "shortstack/{sample}/Results.txt",sample=config["samples"])


rule all:
     input:
         SHORTSTACK
    #message:"All done! Removing intermediate files"
    #shell: "rm -rf {WORKING_DIR}"

rule shortstack:
    input:
        reads =  WORKING_DIR + "trim/{sample}.trimmed.size.fastq",
        genome = GENOME
    output:
        RES_DIR + "shortstack/{sample}/Results.txt",
        RES_DIR + "shortstack/{sample}/{sample}_aligned.bam"
    message:"Shortstack analysis of {wildcards.sample} using {input.genome} reference"
    params:
        RES_DIR + "shortstack/{sample}/"
    threads: 10
    shell:
        "ShortStack "
        "--outdir {wildcards.sample} "
        "--bowtie_cores {threads} "
        "--sort_mem 4G "
        "{SHORTSTACK_PARAMS} "
        "--readfile {input.reads} "
        "--genome {input.genome};"
        "cp -r {wildcards.sample}/* {params};"
        "mv {params}{wildcards.sample}_unaligned.bam {params}{wildcards.sample}_aligned.bam ;"
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

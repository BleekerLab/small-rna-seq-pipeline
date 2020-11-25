import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
from datetime import datetime 

from helpers import extract_hairpin_name_and_sequence
from helpers import collect_clusterfiles_path
from helpers import converts_list_of_sequence_dictionary_to_fasta
from helpers import add_blast_header_to_file
from helpers import add_sample_name_and_hairpin_seq_to_shortstack
from helpers import concatenate_shortstacks_and_assign_unique_cluster_ids
from helpers import extract_mature_micrornas_from_concatenated_shortstack_file
from helpers import extract_hairpins_from_concatenated_shortstack_file
from helpers import extract_mature_mirna_fasta_file_from_shortstack_file
from helpers import create_df_of_seq_length_distributions
from helpers import check_samples_tsv_file

###############################
# OS and related configurations
###############################

##### set minimum snakemake version #####
min_version("5.4.3")

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

#########################
## Pipeline configuration
##########################
configfile: "config/config.yaml"

wildcard_constraints:
  dataset="[Aa-Zz0-9]+"

# directories
WORKING_DIR = config["temp_dir"]
RES_DIR = config["result_dir"]
CURRENT_TIME = datetime.now().strftime("%Y-%m-%d_%H-%M")  # time when pipeline started (will be used to rename result directory)


# Samples: verify 
# get list of samples
# The first functions verify that the provided sample file 
samples_df = check_samples_tsv_file(sample_tsv_file = "config/samples.tsv")
SAMPLES = samples_df.index.values.tolist()

# get fastq file
def get_fastq_file(wildcards):
    fastq_file = samples_df.loc[wildcards.sample,"fastq"]
    return fastq_file

# ShortStack parameters
SHORTSTACK_PARAMS = " ".join(config["shortstack"].values())

####################
## Desired outputs
####################
QC = RES_DIR + "qc/multiqc_report.html"

SEQ_DISTRI = RES_DIR + "seq_length_distribution.tsv"

SHORTSTACK = expand(RES_DIR + "shortstack/{sample}/Results.with_sample_name_and_hairpins.tsv",sample = SAMPLES)
SHORTSTACK_CONCAT = RES_DIR + "concatenated_shortstacks.tsv"

MIRNAS = [expand(RES_DIR + "fasta/{sample}.mature_mirnas.fasta",sample = SAMPLES), RES_DIR + "mature_mirnas.fasta"]
HAIRPINS = [expand(RES_DIR + "fasta/{sample}.hairpin.fasta",sample = SAMPLES), RES_DIR + "hairpins.fasta"]

MFEs = expand(RES_DIR + "hairpins.mfe", sample=SAMPLES) # minimal free energy secondary structures of hairpins

BLAST = expand(RES_DIR + "blast/{sample}.{type}_mirbase.header.txt",sample = SAMPLES, type = ["mature","hairpin"])

rule all:
    input:
        QC,
        SEQ_DISTRI,
        SHORTSTACK,
        SHORTSTACK_CONCAT,
        MIRNAS,
        HAIRPINS,
        BLAST, 
        MFEs
    message:
        "All done! Removing intermediate files in {WORKING_DIR}. Adding date and current time to {RES_DIR} folder name"
    params:
        new_result_dir_name =  CURRENT_TIME + "_" + RES_DIR
    shell:
        "rm -rf {WORKING_DIR};" # removes unwanted intermediate files
        "mv {RES_DIR} {params.new_result_dir_name};"


#######
# Rules
#######

#########################
# Folding of RNA hairpins
#########################

rule rna_fold: 
    input:
        hairpins = RES_DIR + "hairpins.fasta"
    output: 
        mfe = RES_DIR + "hairpins.mfe"
    message: "Calculate minimum free energy secondary structures of hairpins"
    params:
        temp_name = "hairpins.mfe"
    threads: 10
    shell:
        "RNAfold --jobs={threads} --infile={input.hairpins} --outfile={params.temp_name};"
        "mv {params.temp_name} {RES_DIR}{params.temp_name};"
        "rm cluster*"

#############################################
# Produce a concatenated Shortstack dataframe
#############################################

rule extract_fasta_files_for_hairpins_and_mature_miRNAs_from_concatenated_shortstack: 
    input:
        RES_DIR + "concatenated_shortstacks.tsv" 
    output:
        hairpins = RES_DIR + "hairpins.fasta",
        mature = RES_DIR + "mature_mirnas.fasta"
    message: "extract hairpins and mature miRNAs from {input}"
    run:
        extract_hairpins_from_concatenated_shortstack_file(input[0], output[0])
        extract_mature_micrornas_from_concatenated_shortstack_file(input[0], output[1])


rule concatenate_shorstacks_and_assign_unique_cluster_ids:
    input:
        expand(RES_DIR + "shortstack/{sample}/Results.with_sample_name_and_hairpins.tsv", sample=SAMPLES)
    output:
        RES_DIR + "concatenated_shortstacks.tsv"
    message: "Row-bind all Shortstacks and assign a unique id to each sRNA cluster"
    run: 
        dfs = [pd.read_csv(f,sep="\t") for f in input]
        df = pd.concat(dfs)
        df["cluster_unique_id"] = ["cluster_" + str(i+1).zfill(10) for i in range(0,df.shape[0],1)]  
        df.to_csv(output[0], sep="\t", index=False, header=True, na_rep = "NaN")


rule add_sample_name_and_hairpin_seq_to_shortstack:
    input:
        RES_DIR + "shortstack/{sample}/Results.txt", 
        RES_DIR + "fasta/{sample}.hairpin.fasta"
    output:
        RES_DIR + "shortstack/{sample}/Results.with_sample_name_and_hairpins.tsv"
    message: "Add sample name and discovered hairpin sequences to {wildcards.sample} Shortstack dataframe"
    run:
        add_sample_name_and_hairpin_seq_to_shortstack(
            path_to_shortstack_results = input[0],
            sample_name = wildcards.sample,
            hairpin_fasta_file = input[1],
            outfile = output[0]
            )

##################
# mirbase analysis
##################
rule add_blast_header:
    input:
        hairpin = WORKING_DIR + "blast/{sample}.hairpin_mirbase.txt",
        mature = WORKING_DIR + "blast/{sample}.mature_mirbase.txt"
    output:
        hairpin = RES_DIR + "blast/{sample}.hairpin_mirbase.header.txt",
        mature = RES_DIR + "blast/{sample}.mature_mirbase.header.txt"
    message: "adding blast header for {wildcards.sample}"
    params:
        blast_header = "qseqid \t subject_id \t pct_identity \t aln_length \t n_of_mismatches gap_openings \t q_start \t q_end \t s_start \t s_end \t e_value \t bit_score"
    run:
        add_blast_header_to_file(input[0],output[0]) # hairpin
        add_blast_header_to_file(input[1],output[1]) # mature

rule blast_hairpin_against_mirbase:
    input:
        db = config["refs"]["mirbase"]["hairpin"] + ".nhr",
        fasta = RES_DIR + "fasta/{sample}.hairpin.fasta"
    output:
        WORKING_DIR + "blast/{sample}.hairpin_mirbase.txt"
    message:"blasting {wildcards.sample} hairpins against mirbase"
    params:
        dbname = config["refs"]["mirbase"]["hairpin"],
        max_target_seqs = config["blastn"]["hairpin"]["max_target_seqs"]
    shell:
        "blastn -db {params.dbname} "
        "-max_target_seqs {params.max_target_seqs} "
        "-outfmt 6 "                                  # tabular output format
        "-query {input.fasta} "
        "-out {output}"

rule extract_hairpin_fasta_file:
    input:
        RES_DIR + "shortstack/{sample}/Results.txt" # not used in the actual rule but necessary to chain this rule to the shortstack rule
    output:
        RES_DIR + "fasta/{sample}.hairpin.fasta"
    message: "extracting hairpin sequences for {wildcards.sample} clusters annotated as true MIRNAs by shortstack"
    params:
        mirna_clusterpath = RES_DIR + "shortstack/{sample}/MIRNAs/",
        samples = list(config["samples"])
    run:
        for sample in params[1]:
            # make a list of sequence dictionaries (clusterName:hairpinSequence)
            l = [extract_hairpin_name_and_sequence(cluster,sample) for cluster in collect_clusterfiles_path(params[0])]
            # writes this dictionary to a fasta file
            converts_list_of_sequence_dictionary_to_fasta(l,output[0])

rule blast_mature_mirna_against_mirbase:
    input:
        db = config["refs"]["mirbase"]["mature"] + ".nhr",
        fasta = RES_DIR + "fasta/{sample}.mature_mirnas.fasta"
    output:
        WORKING_DIR + "blast/{sample}.mature_mirbase.txt"
    message:"blasting {wildcards.sample} mature miRNAs against mirbase"
    params:
        dbname = config["refs"]["mirbase"]["mature"],
        qcov_hsp_perc = config["blastn"]["mature"]["qcov_hsp_perc"],
        max_target_seqs = config["blastn"]["mature"]["max_target_seqs"]
    shell:
        "blastn -db {params.dbname} "
        "-task blastn-short "                         # BLASTN program optimized for sequences shorter than 50 bases
        "-qcov_hsp_perc {params.qcov_hsp_perc} "
        "-max_target_seqs {params.max_target_seqs} "
        "-outfmt 6 "                                  # tabular output format
        "-query {input.fasta} "
        "-out {output}"

rule make_mirbase_blastdb:
    input:
        mature = config["refs"]["mirbase"]["mature"],
        hairpin = config["refs"]["mirbase"]["hairpin"]
    output:
        mature = config["refs"]["mirbase"]["mature"] + ".nhr",
        hairpin = config["refs"]["mirbase"]["hairpin"] + ".nhr",
    message: "creating blastdb databases for mature miRNA and hairpins"
    shell:
        "makeblastdb -in {input.mature} -dbtype nucl;"
        "makeblastdb -in {input.hairpin} -dbtype nucl"

rule extract_mature_mirna_fasta_file_from_shortstack_file:
    input:
        RES_DIR + "shortstack/{sample}/Results.txt"
    output:
        RES_DIR + "fasta/{sample}.mature_mirnas.fasta"
    message: "extracting mature miRNA fasta file of {wildcards.sample}"
    run:
        extract_mature_mirna_fasta_file_from_shortstack_file(
            path_to_shortstack_file = input[0],
            out_fasta_file = output[0])



######################
## Shortstack analysis
######################

rule shortstack:
    input:
        reads =  WORKING_DIR + "trimmed/{sample}.trimmed.fastq"
    output:
        RES_DIR + "shortstack/{sample}/Results.txt"
    message:"Shortstack analysis of {wildcards.sample} using {params.genome} reference"
    params:
        resdir = RES_DIR + "shortstack/{sample}/",
        genome = lambda wildcards: samples_df.loc[wildcards.sample,"genome"]
    threads: 10
    shell:
        "ShortStack "
        "--outdir {wildcards.sample} "
        "--bowtie_cores {threads} "
        "--sort_mem 4G "
        "{SHORTSTACK_PARAMS} "
        "--readfile {input.reads} "
        "--genome {params.genome};"
        "cp -r {wildcards.sample}/* {params.resdir};"
        "rm -r {wildcards.sample};"


###########################################################
## Get read length distribution (before and after trimming)
##########################################################

rule read_length_distribution:
    input: 
        expand(WORKING_DIR + "trimmed/{sample}.trimmed.fastq", sample = SAMPLES)
    output:
        RES_DIR + "seq_length_distribution.tsv"
    message: 
        "Computing sequence length distribution for all samples"
    params:
        path_to_fastq_files = WORKING_DIR + "trimmed/"
    run:
        create_df_of_seq_length_distributions(path_to_fastq_files =  params.path_to_fastq_files, outfile = output[0])


###########################################
## Trim reads for all samples and QC report
###########################################
rule multiqc_report:
    input:
        expand(WORKING_DIR + "trimmed/{sample}_fastp.json", sample = SAMPLES)
    output:
        RES_DIR + "qc/multiqc_report.html"
    message:
        "Compiling QC reports from fastp"
    params:
        input_directory = WORKING_DIR + "trimmed/",
        output_directory = RES_DIR + "qc/"
    shell:
        "multiqc "
        "--force " # force directory to be created
        "--outdir {params.output_directory} "
        "{params.input_directory}"

rule fastp:
    input:
        get_fastq_file
    output:
        fastq = WORKING_DIR + "trimmed/{sample}.trimmed.fastq",
        json = WORKING_DIR + "trimmed/{sample}_fastp.json",
        html = WORKING_DIR + "trimmed/{sample}_fastp.html"
    message:
        "trimming {wildcards.sample} reads on quality and adapter presence"
    params:
        adapters_fasta =                config["fasta_adapters"],
        qualified_quality_phred =       config["fastp"]["qualified_quality_phred"],
        average_quality =               config["fastp"]["average_quality"],
        min_length =                    config["fastp"]["min_length"],
        max_length =                    config["fastp"]["max_length"]
    shell:
        "fastp -i {input} "
        "--stdout "
        "--json {output.json} "
        "--html {output.html} "
        "--qualified_quality_phred {params.qualified_quality_phred} "
        "--average_qual {params.average_quality} "
        "--length_required {params.min_length} "
        "--length_limit {params.max_length} "
        "--adapter_fasta {params.adapters_fasta} > {output.fastq}"


      
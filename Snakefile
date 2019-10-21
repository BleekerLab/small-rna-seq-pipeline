import os
import subprocess
from snakemake.utils import min_version
import pandas as pd
from helpers import extract_hairpin_name_and_sequence
from helpers import collect_clusterfiles_path
from helpers import converts_list_of_sequence_dictionary_to_fasta
from helpers import add_blast_header_to_file

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
SHORTSTACK = expand(RES_DIR + "shortstack/{sample}/Results.txt",sample=config["samples"])



MIRNAS = expand(RES_DIR + "fasta/{sample}.mature_mirnas.fasta",sample=config["samples"])
HAIRPINS = expand(RES_DIR + "fasta/{sample}.hairpin.fasta",sample=config["samples"])
BLAST = expand(RES_DIR + "blast/{sample}.{type}_mirbase.header.txt",sample=config["samples"],type=["mature","hairpin"])
PLOTS = [RES_DIR + "plots/n_clusters_per_dicercall.png",
         expand(RES_DIR + "plots/piecharts/{sample}.piechart.png",sample=config["samples"]),
         expand(RES_DIR + "plots/MIR/{sample}.mirgenes.png",sample=config["samples"])
         ]

rule all:
    input:
        SHORTSTACK,
        MIRNAS,
        HAIRPINS,
        BLAST,
        PLOTS
    message:"All done! Removing intermediate files"
    shell:
        "rm -rf {WORKING_DIR}" # removes unwanted intermediate files

#######
# Rules
#######

########
## Plots
########
rule mir_gene_families_barplot:
    input:
        RES_DIR + "blast/{sample}.hairpin_mirbase.header.txt"
    output:
        png = RES_DIR + "plots/MIR/{sample}.mirgenes.png",
        svg = RES_DIR + "plots/MIR/{sample}.mirgenes.svg"
    message:"Counting and plotting the number of MIR genes per family"
    conda:
        "envs/plots.yaml"
    shell:
        "Rscript --vanilla scripts/mir_gene_families.R {input} {output.png} {output.svg}"

rule pie_chart_srna_classes:
    input:
        RES_DIR + "shortstack/{sample}/Results.txt"
    output:
        RES_DIR + "plots/piecharts/{sample}.piechart.png"
    message: "making a pie chart of {wildcards.sample} small RNA classes based on ShortStack results"
    conda:
        "envs/plots.yaml"
    shell:
        "Rscript --vanilla scripts/piechart.R {input} {output} "

rule plot_number_srna_clusters_per_dicer_call:
    input:
        expand(RES_DIR + "shortstack/{sample}/Results.txt",sample=config["samples"])
    output:
        png = RES_DIR + "plots/n_clusters_per_dicercall.png",
        svg = RES_DIR + "plots/n_clusters_per_dicercall.svg"
    message: "Making the 'number of sRNA cluster = f(DicerCall)' plot based on ShortStack results"
    conda:
        "envs/plots.yaml"
    params:
        shortstack = RES_DIR + "shortstack/"
    shell:
        "Rscript --vanilla scripts/number_srna_clusters_per_dicer_call.R {params.shortstack} {output.png} {output.svg} "

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
    conda:
        "envs/blast.yaml"
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
    conda:
        "envs/blast.yaml"
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
    conda:
        "envs/blast.yaml"
    shell:
        "makeblastdb -in {input.mature} -dbtype nucl;"
        "makeblastdb -in {input.hairpin} -dbtype nucl"

rule extract_mature_mirna_fasta_file:
    input:
        RES_DIR + "shortstack/{sample}/Results.txt"
    output:
        RES_DIR + "fasta/{sample}.mature_mirnas.fasta"
    message: "extracting mature miRNA fasta file of {wildcards.sample}"
    run:
        df = pd.read_csv(input[0],sep="\t",index_col=1)
        df_mirnas = df[df["MIRNA"] == "Y"]
        mirnas_dict = df_mirnas["MajorRNA"].to_dict()
        # convert to fasta format
        with open(output[0],"w") as fileout:
            for name,sequence in mirnas_dict.items():
                fileout.write(">" + name + "\n" + sequence + "\n")

######################
## Shortstack analysis
######################

rule shortstack:
    input:
        reads =  WORKING_DIR + "trim/{sample}.trimmed.size.fastq",
        genome = GENOME
    output:
        RES_DIR + "shortstack/{sample}/Results.txt"
    message:"Shortstack analysis (v3.6.0) of {wildcards.sample} using {input.genome} reference"
    params:
        RES_DIR + "shortstack/{sample}/"
    threads: 10
    conda:
        "envs/shortstack.yaml"
    shell:
        "./scr/ShortStack-3.6/ShortStack "
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

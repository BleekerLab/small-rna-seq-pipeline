FROM continuumio/miniconda3:latest


LABEL author="l.fokkens@uva.nl" \
      description="A Docker image used to build a container that is used to run the small-rna-seq-pipeline" \
      url="https://github.com/BleekerLab/small-rna-seq-pipeline/envs/Dockerfile" \
      usage="docker run mgalland/small-rna-seq-pipeline" 

# Create the environment:
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment and run Snakemake
ENTRYPOINT ["conda", "run", "-n", "small", "/bin/bash", "-c", "snakemake"]




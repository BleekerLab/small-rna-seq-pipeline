FROM continuumio/miniconda:4.7.12

WORKDIR /home/snakemake/

COPY ["environment.yml", "./"]

# mamba is a faster C++ re-implemenbtation of conda
# name of the environment is rnaseq
RUN conda install -c conda-forge mamba --yes \
  && mamba env create -f environment.yaml \
  && conda clean --all

RUN echo "source activate small" > ~/.bashrc
ENV PATH /opt/conda/envs/small/bin:$PATH

ENTRYPOINT ["snakemake"]
CMD ["--version"]
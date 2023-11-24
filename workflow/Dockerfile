FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

## Dockerfile to build Docker image for this Snakemake workflow.
FROM ubuntu:20.04

# Set environment variables for non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe


# Step 1: Install essential dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        git \
        python3-pip \
        curl

# Install Miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Set conda on the PATH
ENV PATH /opt/conda/bin:$PATH
RUN conda init --all

# Step 2: Install R without prompts
RUN apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/salmon.yaml
#   name: salmon
#   channels:
#       - bioconda
#   dependencies:
#       - salmon

RUN conda install -c bioconda salmon
RUN conda install -c conda-forge r-essentials

# Install required R packages
RUN Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com")' &&\
               'install.packages("RColorBrewer2", repos="https://cran.rstudio.com")' &&\
	       'install.packages("magrittr", repos="https://cran.rstudio.com")' &&\
	       'install.packages("stringr", repos="https://cran.rstudio.com")' &&\
	       'install.packages("ggplot2", repos="https://cran.rstudio.com")' && \
	       'install.packages("ggfortify", repos="https://cran.rstudio.com")' && \
	       'install.packages("ggthemes", repos="https://cran.rstudio.com")' && \
	       'install.packages("ggrepel", repos="https://cran.rstudio.com")' && \
	       'install.packages("stringi", repos="https://cran.rstudio.com")' && \
	       'install.packages("tidyverse", repos="https://cran.rstudio.com")' && \
	       'install.packages("ggplot2", repos="https://cran.rstudio.com")' &&\
	       'install.packages("limma", repos="https://cran.rstudio.com")' &&\
	       'install.packages("voom", repos="https://cran.rstudio.com")' &&\
	       'install.packages("edgeR", repos="https://cran.rstudio.com")'

# Step 3: Install and test Snakemake

RUN conda install wheel &&\
    conda install snakemake

RUN git clone https://github.com/JasonCharamis/QuteS.git

RUN curl -O -L https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/run_DE_analysis.pl &&\
    curl -O -L https://raw.githubusercontent.com/trinityrnaseq/trinityrnaseq/master/Analysis/DifferentialExpression/analyze_diff_expr.pl

RUN mv *.pl QuteS/workflow/scripts

ENV PATH="~/.local:$PATH"
ENV PATH="QuteS:$PATH"
ENV PATH="QuteS/workflow/scripts:$PATH"

CMD ["Everything installed!"]

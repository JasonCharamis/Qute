FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

## Dockerfile to build Docker image for this Snakemake workflow.

FROM r-base:latest

# Step 1: Install the essentials
RUN apt-get update && \
    apt-get install -y apt-utils &&\
    apt-get install -y build-essential &&\
    apt-get install -y libcurl4-openssl-dev &&\
    apt-get install -y libssl-dev &&\
    apt-get install -y libxml2-dev &&\
    apt-get install -y git &&\
    apt-get install -y python3-pip &&\
    apt-get install -y curl

RUN apt-get install -y -f

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/salmon.yaml
#   name: salmon
#   channels:
#       - bioconda
#   dependencies:
#       - salmon

RUN conda install -c bioconda salmon

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

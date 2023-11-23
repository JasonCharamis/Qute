FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

## Dockerfile to build Docker image for this Snakemake workflow.

# Step 1: Install the essentials

RUN apt-get update &&\
    apt-get install -f -y &&\
    apt-get install build-essential -y &&\
    apt-get install git -y &&\
    apt-get install python3-pip -y &&\
    apt-get install r-base-core 

# Step 2: Retrieve conda environments

# Conda environment:
#   source: envs/salmon.yaml
#   name: salmon
#   channels:
#       - bioconda
#   dependencies:
#       - salmon

RUN conda install -c bioconda salmon

# Install R dependencies

RUN Rscript -e 'install.packages("dplyr", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("RColorBrewer2", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("magrittr", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("stringr", repos="https://cran.rstudio.com")' &&\
    Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("ggfortify", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("ggthemes", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("ggrepel", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("stringi", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("tidyverse", repos="https://cran.rstudio.com")' && \
    Rscript -e 'install.packages("ggplot2", repos="https://cran.rstudio.com")'

# Step 3: Install and test Snakemake

RUN conda install wheel &&\
    conda install snakemake

RUN git clone https://github.com/JasonCharamis/QuteS.git

ENV PATH="~/.local:$PATH"
ENV PATH="QuteS:$PATH"
ENV PATH="QuteS/workflow/scripts:$PATH"

CMD ["Everything installed!"]

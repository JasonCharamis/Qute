# QuteS (Quantification of transcript expression using Salmon)
Fully automated Snakemake pipeline for streamlining end-to-end salmon-based RNAseq analysis, from fastq reads to Differential Expression (DE) analysis.

Salmon and all other dependencies is automatically installed through conda.

To use as a Docker container, run the following commands:

1. git clone https://github.com/JasonCharamis/QuteS.git
2. cd QuteS/workflow/ && sudo docker build -t automated_transcript_quantification_using_salmon:latest .
3. Travel to the directory where your data live
4. sudo docker run -it -v $(pwd):/mnt/workdir -w /mnt/workdir automated_rnaseq_analysis:latest snakemake --snakefile ARuS/workflow/Snakefile --cores 1 --use-conda --conda-frontend mamba

Of course, to customize the run edit the config/config.yaml file. 
That's it! The pipeline will run automatically.

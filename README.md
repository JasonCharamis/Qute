# Quantification of transcript expression using Salmon (QuteS)
Automated Snakemake workflow for quantification of transcript expression using Salmon (QuteS) and running DE analysis using edgeR.

Salmon and all other dependencies is automatically installed and downloaded through conda.
Moreover, the pipeline is available as a Docker container.

To run as a Docker container, run the following commands:

1. git clone https://github.com/JasonCharamis/QuteS.git

2. cd QuteS/workflow/ && sudo docker build -t automated_transcript_quantification_using_salmon:latest .

3. Travel to the directory where your data live

4. sudo docker run -it -v $(pwd):/workflow -w /workflow automated_transcript_quantification_using_salmon:latest snakemake --cores 20 --use-conda --snakefile Snakefile

Of course, to customize the run edit the config/config.yaml file. 
That's it! The pipeline will run automatically.

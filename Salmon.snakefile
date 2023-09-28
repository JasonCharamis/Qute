import os
import re

samples = { samples[:-19] for samples in os.listdir (".") if samples.endswith(".fastq.gz") }
species = { sample[:-1] for sample in samples }

rule all:
     input:
	  expand ( '{spec}.salmon.index', spec=species ),
	  expand ( '{sample}.salmon.quants', sample = samples )

rule index:
      input: "{species}.fasta"
      output: directory("{species}.salmon.index")
      shell: """ salmon index -t {input} -i {output} --keepDuplicates """
      

rule selective_alignment:
      input:  "{samples}_1.trimmed.fastq.gz", "{samples}_2.trimmed.fastq.gz"
      output: directory ( "{samples}.salmon.quants" )
      params: species = lambda wildcards: wildcards.samples[:-1]
      shell: """ salmon quant -i {params.species}.trinity_assembly.salmon.index -l A -1 {input[0]} -2 {input[1]} --threads 20 --validateMappings -o {output} """


rule quantmerge:
     input: expand ( "{sample}.salmon.quants/quant.sf", sample = samples )
     output:  "{params.prefix}.all.genes.counts"
     params: prefix = c[sample]
     shell: """ salmon quantmerge --quants {input} --column numreads --output {params.prefix}.all.genes.counts && {output} """

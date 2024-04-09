                                        # Principal components analysis using log2TPM values

.load_packages <- function(tools) {
  tmp <- as.data.frame(installed.packages()) 
  max_version <- max(as.numeric(substr(tmp$Built, 1, 1)))
  tmp <- tmp[as.numeric(substr(tmp$Built, 1, 1)) == max_version,]

  for (pkg in tools) {
    if (pkg %in% tmp$Package) {
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    } else {
      print(sprintf("%s %s", pkg, "is not installed. Installing it!"))
      
      if (pkg %in% BiocManager::available(pkg)) {
        BiocManager::install(pkg, dependencies = TRUE, update = TRUE)
      } else {
        install.packages(pkg, dependencies = TRUE, ask = FALSE)
      }
    }
  }
}

# Load required packages or install them if necessary
dependencies <- c("edgeR", 
                  "magrittr",
                  "stringr",
                  "ggfortify",
                  "ggrepel",
                  "ggthemes",
                  "stringi",
                  "dplyr",
                  "tidyverse")

.load_packages(dependencies)

## Get input files from standard input ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))
f_counts <- args

# Open input files: a gene counts table and a file declaring sample grouping of replicates
tpm <- read.table(args[1], header=TRUE)
samples_list <- read.table(args[2], header=TRUE, sep="\t")


## Manipulate data and convert TPM to log2TPM values ------------------------------------------------------------------

# Set first column as rownames
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]

# Transform colnames to same syntax
colnames(tpm) <- str_replace(colnames(tpm), "results.|.s.bam|_", "")

# Convert dataframe to matrix, to allow for log2 transformation
tpm <- as.matrix(tpm)

# Function to convert tpm to log2tpm values
logTPM <- function(tpm, dividebyten=TRUE) {
  if(dividebyten) {
    logtpm <- log(tpm/10+1, 2)}
  else if(!dividebyten) {
    logtpm <- log(tpm+1, 2)}
  return(logtpm)
}

logtpms<-logTPM(tpm, dividebyten = FALSE)


## Run PCA analysis ------------------------------------------------------------------

# Prepare data structures
xt = t(logtpms)
xt <- as.data.frame(xt)
groups <- samples_list$V1 # get group column

# Add column with sample name to group replicates
xtl <- xt %>% add_column(Sample = groups)

# Run PCA analysis 
pca_res = prcomp(xt, center=T, scale.=F)

svg("PCA.svg")

# Draw auto-PCA plot with color mappings for groups
print ( autoplot(pca_res,
                 data = xtl,
                 colour = 'Sample',
                 legend.size = 5) +
        labs(title = "Principal Component Analysis using log2TPM values" ) +
        theme_bw() +
        geom_text_repel( aes( label = rownames(xt), color = Sample ), show.legend = FALSE  ) +
        theme( plot.title = element_text(face = "bold", hjust = 0.5),
              legend.title = element_text( size = 12),
              legend.text = element_text( size = 12 ),
              axis.title = element_text( size = 12 )
              )
       )

dev.off()

# Before running this code:

## Install R Studio

### Download: https://rstudio.com/products/rstudio/download/

### Tutorial: https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu

## Install tidyverse package

### Instructions: https://tidyverse.tidyverse.org/

## Please also install the following packages: csv, dada2, phyloseq

#### dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

#### phyloseq
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")

#### csv
install.packages("csv")

#### check your package installation:

library(csv)

library(dada2)

library(phyloseq)

library(tidyverse)

#### The code for each of these should run without any errors

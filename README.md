# microbiomeQTL
A pipeline for QTL analysis of host microbiomes

The input required is a OTU table and a mapping table

## Installation
```{r]
Various thinsg were required to install on the NIOO servers
install.packages("digest")
install.packages("devtools")
install.packages("curl")
devtools::install_github("klutometis/roxygen")
getwd()
setwd("/mnt/nfs/bioinfdata/home/NIOO/beno/microbiomeQTL")
create("microbiomeQTL")
library(roxygen2)
library("devtools")
create("microbiomeQTL")


# Prerequisites 
install_github("cytoscape/RCy3")

```

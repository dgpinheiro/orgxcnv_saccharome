## Diretório
setwd("/work/rcsilva/projects/saccharome-wgs/")

## Matrizes de entrada

#### Taxonomias
table_kraken <- read.delim(file = "./analises/2019-04-29_kraken_all_official/tabletax.noplants.txt")
colnames(table_kraken)[1] <- "OTU_ID"

table_16S <- read.delim(file = "../saccharome-16s-rizo/analises/outputs/P01_A01_zotus_work_in_R/output/L1_OTUs.txt")
table_ITS <- read.delim(file = "../saccharome-its-rizo/analises/E14_work_in_R/P01_A01_work_in_R/output/L1_OTUs.txt")

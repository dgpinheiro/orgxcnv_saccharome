## Download NR via biomaRt

## setwd
setwd("/work/rcsilva/projects/saccharome-wgs")

## Bibliotecas
library(biomartr)

## DownloadDB
download.database.all(db = "nr", path = "/data/db/metagenomics/db-NR")

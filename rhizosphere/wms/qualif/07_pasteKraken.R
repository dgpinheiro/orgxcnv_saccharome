## Juntar taxonomia do Kraken com Anotações
setwd('/work/rcsilva/projects/saccharome-wgs/')

# Lista de anotações
View(masterlist)

# Objetivo: colar taxonomia do Kraken
annot.kraken <- read.delim(file = './analises/2019-04-30_kraken_contigs/contig_name.tsv',
                           stringsAsFactors = F)

colnames(annot.kraken)

## Juntar as listas
masterlist <- merge(annot.bkp, annot.kraken, by = "contig", all.x = T)
masterlist.p <- subset(masterlist, ! is.na(pvalue)) 

## Adicionar contagens
masterlist.p.cpm <- merge(masterlist.p, cpm_table, by = "contig")

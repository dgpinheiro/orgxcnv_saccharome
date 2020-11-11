## Cruzando core OTUs

## preproc
setwd('/work/rcsilva/projects/saccharome-wgs')

## Core OTUs (gênero)
rs.tab <- read.delim(file = './dados/core_microbiome/genus.tsv',
                     stringsAsFactors = F)

## Tabela de taxons
tax.de.tab <- read.delim(file = './dados/subsets/deseq2-kaiju-taxonomy/de_table.tsv',
                         stringsAsFactors = F)

## Table de localiza;'ao
rs.locations <- read.delim(file = './dados/core_microbiome/otus_per_location.tsv',
                           stringsAsFactors = F)

## Cruzam os gêneros
tax.intersec <- subset(tax.de.tab, tax.de.tab$genus %in% rs.tab$genus)

## Adiciona os OTUs
tax.inter.2 <- merge(tax.intersec, rs.tab, by = "genus")

## Cola a localização
tax.inter.2$is.rhizo <- tax.inter.2$OTU_id %in% rs.locations$RZO
tax.inter.2$is.bulk.soil <- tax.inter.2$OTU_id %in% rs.locations$BKS
tax.inter.2$is.endo.root <- tax.inter.2$OTU_id %in% rs.locations$ERA

## Tabela final
core.tax.final <- tax.inter.2

write.table(core.tax.final,
            file="./dados/subsets/deseq2-core-microbiome/tax.tsv",
            quote = F,
            sep = '\t',
            row.names = F)




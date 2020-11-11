## Tables for paper
## Objetivo: criar uma tabela completa, massiva de anotação para todas as proteínas!
## 1. Usar a tabela de contigs
## 2. Colar os EC no
## 3. Colar os KEGG
## 4. Colar os GO
## 5. Colar os COGs
## 6. Colar as anotações do Prokka
## 7. Colar as anotações do EGGNOG
## 8. Colar l2fc
## 9. Colar a taxonomia do Kaiju (prots)

## Bibliotecas
library(data.table)
library(tidyr)
setwd('/work/rcsilva/projects/saccharome-wgs/')

## 1. Tabela de contigs
head(contig_X_prokka_table)
dim(contig_X_prokka_table) # 4.962.721 proteínas

## 2. Adicionar EC number
annot_table <- data.table(contig_X_prokka_table)
annot_table <- merge(annot_table, ec_table,
                     by = "prk_id",
                     all = T,
                     missing = "Unclassified")

## 3. Colar os KEGG
head(prokka_X_KEGG_table)
annot_table <- merge(annot_table, prokka_X_KEGG_table,
                     by = "prk_id",
                     all = T)

## 4. Colar os GOs
# Do eggnog estão very, very zona!
# Pegar os do IPRScan quando terminarem.

## 5. Colar os COGs
annot.eggnog.all.non.electronic <- read.delim(file = './analises/2019-04-01_eggnog_not_electronic/out.tsv',
                                              stringsAsFactors = F)

colnames(annot.eggnog.all.non.electronic)
# [1] "query_name"           "seed_eggNOG_ortholog" "seed_ortholog_evalue" "seed_ortholog_score"  "predicted_gene_name"  "GO_terms"            
# [7] "KEGG_KOs"             "BiGG_reactions"       "Annotation_tax_scope" "OGs"                  "bestOG.evalue.score"  "COG.cat"             
# [13] "eggNOG.annot"
colnames(annot.eggnog.all.non.electronic)[1] <- "prk_id"

annot_table <- merge(annot_table, annot.eggnog.all.non.electronic[, c(1, 12)],
                     by = "prk_id",
                     all = T)

## 6. Colar os produtos do Prokka
prokka_products <- read.delim(file = './analises/2019-02-18_prokka/00_parsed/products/products.tsv',
                              stringsAsFactors = F)

prokka_products <- prokka_products[,-c(1)]

annot_table <- merge(annot_table, prokka_products,
                     by = "prk_id",
                     all = T)

colnames(annot_table)[5] <- "prokka_product"

## 7. Colar os produtos do eggnog
## Anotação
annot_table <- merge(annot_table,
                     annot.eggnog.all.non.electronic[,c(1,13)],
                     by = "prk_id",
                     all = T)

colnames(annot_table)

## 8. Colar os nomes de genes do eggnog
eggnog.gene.name <- read.delim(file = './analises/2019-02-18_prokka/00_parsed/gene_name/gene_parsed.tsv',
                               stringsAsFactors = F)

annot_table <- merge(annot_table,
                     eggnog.gene.name,
                     by = "contig",
                     all = T)

## Remover linhas bichadas
annot.bkp <- annot_table
annot_table <- annot_table[-c(1,2,3),]

## 9. Adicionar l2fc
# objeto: results.sig.p.tryhard mod LFC > 1, p < 0.05
wgs.dab.tryhard <- results.sig.p.tryhard[,c(2,5)]
wgs.dab.tryhard$contig <- row.names(wgs.dab.tryhard)

## Adiciona os significativos
annot_table <- merge(annot_table, as.data.frame(wgs.dab.tryhard),
                     by = "contig",
                     all = T)

## 10. Adicionar anotação (proteínas) do Kaiju

tax.annot.wgs.kaiju <- read.delim(file = './analises/2019-03-12_kaiju-megahit-prodigal/output_assembly_proteins/parsed.kaiju.out',
                                  stringsAsFactors = F)

# Remove linhas vazias
tax.annot.wgs.kaiju <- subset(tax.annot.wgs.kaiju, kaiju_tax != "")

annot_table <- merge(annot_table,
                     tax.annot.wgs.kaiju,
                     by = "prk_id",
                     all = T)

# Separar taxonomia por linhas
annot_table <- separate(data = annot_table,
                        col = "taxon1",
                        into = sc("taxon1", "taxon2", "taxon3"),
                        sep = ",")

#### saídas finais ####
# matriz completa
annot_table

# Somente p-value
annot_sig <- subset(annot_table, ! is.na(annot_table$pvalue))
# Adicionar contagem
annot_sig_new <- merge(annot_sig, cpm_table, by = "contig")

# Adicionar matriz do kraken
kraken_tax_assembly <- read.delim(file = './analises/2019-04-16_kraken-05/kraken-parsed-to-R.tsv',
                                  stringsAsFactors = F)

annot_sig_new <- merge(annot_sig_new, kraken_tax_assembly, by = "contig")

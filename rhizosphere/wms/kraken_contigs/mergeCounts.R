## Merge kraken contigs and normalized counts
setwd("/data/rcsilva/projects/saccharome-wgs/")

## Leia a classificação
kr_class = read.delim(
  file = "./analises/2019-08-28_new_classifications/kraken_out_assembly/2019-08-28_megahit.fa_conf0.51_kraken2def/classified_out_per_contig.tsv",
  header = TRUE
)

## Contagens normalizadas
counts = as.data.frame(deseq.norm.counts)


#### 1. Unir as contagens com contigs anotados ####
counts$contig = row.names(counts)
counts_taxonomy = merge(x = kr_class, y = counts, by = "contig")

write.table(counts_taxonomy,
            file="/data/AEM/dados/2019-08-28_novo_kraken_confidence51/contigs/otu_counts.txt",
            quote = F,
            sep = '\t',
            row.names = F)


#### 2. Remover contaminantes de fungos ####
fungi_taxa = read.delim(
  file = "/data/AEM/dados/2019-08-28_novo_kraken_confidence51/3_lab_sample/parse/taxonomy.txt"
)

fungi_counts = read.delim(
  file = "/data/AEM/dados/2019-08-28_novo_kraken_confidence51/3_lab_sample/parse/counts.txt"
)

fungi_contaminants = merge(
  fungi_taxa,
  fungi_counts,
  by = "OTU"
)

write.table(fungi_contaminants,
            file="/data/AEM/dados/2019-08-28_novo_kraken_confidence51/3_lab_sample/fungi_contaminants.txt",
            quote = F,
            sep = '\t',
            row.names = F)

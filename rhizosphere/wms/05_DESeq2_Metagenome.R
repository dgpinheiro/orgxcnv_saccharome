### DESeq2 para Metagenomas!!!! ###
library("BiocParallel")
library("tximport")
library("DESeq2")
library("clValid")
library("ggplot2")
library("xlsx")
library("stringr")

## Capturando a pasta de trabalho
setwd("/work/rcsilva/projects/saccharome-wgs")

## Entradas:
dir <- "/work/rcsilva/projects/saccharome-wgs/analises/2019-02-26_kallisto_count_megahit01"
setwd(dir)

# Amostras
samples <- read.table("/work/rcsilva/projects/saccharome-wgs/analises/2019-02-26_kallisto_count_megahit01/samples.txt",
                      header=T)

# Número de processadores, não usar mais que 20 no Thor
parameter.threads <- 30
register(MulticoreParam(parameter.threads))

# Importar os dados
files <- file.path(dir, "kallisto", samples$folder, "abundance.tsv")
names(files) <- (samples$folder)

## Objeto do tximport para importar no DESeq2
txi <- tximport(files, type = "kallisto", txOut=TRUE)

## Importando metagenoma
dds <- DESeqDataSetFromTximport(txi, colData = samples, design =~ treatment)

## Para gerar contagens:
deseq.counts <- estimateSizeFactors(dds)
deseq.norm.counts <- counts(deseq.counts, normalized=TRUE)

## Então, roda o DESeq
deseq.dds <- DESeq(dds, parallel=T, fitType="local")


## Colhe resultados
dds.res <- results(deseq.dds, parallel = T)

## Tira os significativos
results.sig <- subset(dds.res, (dds.res$log2FoldChange > 2 &
                         dds.res$padj < 0.05) |
                        dds.res$log2FoldChange < -2 &
                        dds.res$padj < 0.05 )

# Significativos
drownames(results.sig)
results.sig$log2FoldChange

## Junta resultados significativos e L2FC
df <- cbind(row.names(results.sig), results.sig$log2FoldChange)
            
# Escrever a saída
write.table(df, 
            file="/work/rcsilva/projects/saccharome-wgs/analises/2019-03-01_DESeq2_Megahit_out/outputs/out-lfc.tab",
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = F)


## Tira os significativos
results.sig.pval <- subset(dds.res,
                           (dds.res$log2FoldChange > 2 & dds.res$padj < 0.2 |
                              dds.res$log2FoldChange < -2 & dds.res$padj < 0.2 )
)

#### Extração das contagens ####

counts_sigmax <- subset(deseq.dds,
               names %in% rownames(results.sig))

counts_sigmax_table <- counts_sigmax@assays$data$counts

write.table(counts_sigmax_table,
            file="/work/rcsilva/projects/saccharome-wgs/analises/2019-03-08_DESeq2_GOcounts/counts/counts.txt",
            quote = F,
            sep = '\t')
        

#### Extração dos significativos sem correção de FDR (p-value) #####
sig.pvalue <- subset(dds.res,
                     (dds.res$log2FoldChange > 2 & dds.res$pvalue < 0.05 |
                        dds.res$log2FoldChange < -2 & dds.res$pvalue < 0.05 )
                     )

### Pega os nomes
names <- as.vector(row.names(sig.pvalue))

### Exporta os nomes
write.table(names,
            file="/work/rcsilva/projects/saccharome-wgs/analises/2019-03-11_GOcounts_pval_toSTAMP/names.txt",
            quote = F,
            sep = '\t',
            row.names = F)

  


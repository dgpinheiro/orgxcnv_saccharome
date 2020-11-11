# Neste script:
## 1. Obtenção das contagens normalizadas
## 2. Geração da tabela KEGG para normalizados
## 3. Extração do subconjunto P (significativo a p-valor) com EC number

## Separar p-value por EC

# Bibliotecas
library("R.utils")
library('DESeq2')

# Fixando diretório
setwd('/work/rcsilva/projects/saccharome-wgs/')

## Carregar os dados na memória do R
## Carregar -> load('./dados/RData/all-plots.RData')
## Salvar -> save.image("/data/rcsilva/projects/saccharome-wgs/dados/RData/all-plots.RData")

## Carrega a matriz do DESeq2
dds.res <- loadToEnv("./dados/RData/megahit-deseq.RData")[["dds.res"]];
dds <- loadToEnv("./dados/RData/megahit-deseq.RData")[["dds"]];

## Estimar fatores de tamanho para extração dos valores de contagem normalizados
sf.sac.wgs <- estimateSizeFactors(dds)
sf.sac.wgs.counts <- counts(sf.sac.wgs, normalized=T)
sac.wgs.normalized.counts <- as.data.frame(sf.sac.wgs.counts)
sac.wgs.normalized.counts$contig <- NA
sac.wgs.normalized.counts$contig <- rownames(sac.wgs.normalized.counts)

## Separa os significativos
results.sig.p <- subset(dds.res,
                      (dds.res$log2FoldChange >= 1 &
                         dds.res$pvalue < 0.05) |
                        dds.res$log2FoldChange <= -1 &
                        dds.res$pvalue < 0.05 )

results.sig.q <- subset(dds.res,
                      (dds.res$log2FoldChange >= 2 &
                         dds.res$padj < 0.05 ) |
                        dds.res$log2FoldChange <= -2 &
                        dds.res$padj < 0.05 )

## Significativos independente de valor?
results.sig.p.tryhard <- subset(dds.res,
                               (dds.res$log2FoldChange >= 1 &
                                  dds.res$pvalue < 0.05) |
                                 dds.res$log2FoldChange <= 1 &
                                 dds.res$pvalue < 0.05 )

## Extrai os valores
sig.values <- results.sig[,c(2,5)]

## Cola os nomes de linhas como id
sig.values$contig <- NA
sig.values$contig <- row.names(sig.values)

## Cola os IDs de ORFs
pval_annot_table <- merge(as.data.frame(sig.values), contig_X_prokka_table, by = "contig")


#### KEGG: Trocar os valores pelos normalizados ####

# Remover as contagens
colnames(cpm_KO_dict)
orf_ko_dict <- cpm_KO_dict[,-c(4,5,6,7,8,9)]

# Colar as contagens normalizadas
cpm_normalized.orf.ko.dict <- merge(sac.wgs.normalized.counts, orf_ko_dict, by = "contig")

# Exportando para parse
write.table(cpm_normalized.orf.ko.dict,
            file="./analises/2019-04-05_normalized_kegg_plots_STAMP/norm_ko.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### KEGG: Subconjunto p-significativo ####
orf.ko.pval <- subset(cpm_normalized.orf.ko.dict, contig %in% sig.values$contig)

## Separando sem cabeçalhos sujos
colnames(orf.ko.pval)
orf.ko.pval.ok <- orf.ko.pval[,c(10,11,12,13,2,3,4,5,6,7)]

## Exportando a tabela pval-sig
write.table(orf.ko.pval.ok,
            file="./analises/2019-04-05_normalized_kegg_plots_STAMP/pval.sig.norm_ko.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### KEGG counts para o Microbiome Analyst ####
# Abre a tabela do eggnog não-eletrônico
annot.eggnog.non.electro <- read.delim(file = './analises/2019-04-01_eggnog_not_electronic/parsed_KEGG/prk_kegg.tsv',
                                       stringsAsFactors = F)

# Limpa os vazios
annot.eggnog.non.electro.clean <- subset(annot.eggnog.non.electro, KO != "")

# Cola as contagens
head(annot.eggnog.non.electro.clean)

# Parsear a tabela (virgula -> newline)
write.table(annot.eggnog.non.electro.clean,
            file="/tmp/parse-ko",
            quote = F,
            sep = '\t',
            row.names = F)

# Re-importa
parsed.eggnog.ko <- read.delim(file = '/tmp/parsed-ko',
                               stringsAsFactors = F)


kegg_manalyst <- merge(parsed.eggnog.ko,
                       cpm_normalized.orf.ko.dict[,c(2,3,4,5,6,7,9)],
                       by = "prk_id")

## Pronto?
head(kegg_manalyst)

## Remove PRK_ID
kegg_manalyst.ready <- kegg_manalyst[,-c(1)]

## Soma os termos KEGG repetidos
kegg_manalyst.ready <- ddply(kegg_manalyst.ready, "KO", numcolwise(sum))

## Exportar!
write.table(kegg_manalyst.ready,
            file="./analises/2019-04-05_normalized_kegg_plots_STAMP/kegg_manalyst.ready.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### KEGG, origem: Koala ####
# Importando a matriz Koala
koala.annot.150 <- read.delim(file = './analises/2019-04-05_ghost-koala-annotation/query.ko',
                              stringsAsFactors = F)

colnames(koala.annot.150) <- c("prk_id", "KO")

# Limpa os vazios
koala.clean <- subset(koala.annot.150, KO != "")

# Recuperar as contagens brutas
koala_manalyst.step1 <- merge(koala.clean,
                              contig_X_prokka_table,
                              by = "prk_id")

koala_manalyst.step2 <- merge(koala_manalyst.step1,
                             cpm_table,
                             by = "contig")

## Somente KO e contagem
colnames(koala_manalyst.step2)
koala_ready <- ddply(koala_manalyst.step2, "KO", numcolwise(sum))
head(koala_ready)

length_header <- grep("length", colnames(koala_ready))

## Remove length
koala_ready <- koala_ready[,-c(length_header)]
head(koala_ready)
             
## Exporta
write.table(koala_ready,
            file="./analises/2019-04-05_ghost-koala-annotation/koala-manalyst.txt",
            quote = F,
            sep = '\t',
            row.names = F)         

#### Gene ontology com as contagens normalizadas para heatmap ####
head(eggnog_go_prk)
head(eggnog_dict)

## 1. Cola o contig com os GOs
# Lê os novos GOs (não-eletrônico)
sac.enog.input.nonelec <- read.delim(file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-04-01_eggnog_not_electronic/parsed_GO/prk_go.tsv',
                                     stringsAsFactors = F)

sac.enog.input.nonelec.clean <- subset(sac.enog.input.nonelec, GO_terms != "")

## Exporta de volta (parse)
write.table(sac.enog.input.nonelec.clean,
            file="./analises/2019-04-01_eggnog_not_electronic/parsed_GO/to_clean.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

## Lê de volta
sac.enog.clean <- read.delim(file = './analises/2019-04-01_eggnog_not_electronic/parsed_GO/cleaned.tsv',
                             stringsAsFactors = F)

## Cola os contigs
sac.go.contigs <- merge(sac.enog.clean, contig_X_prokka_table, by = "prk_id")

## Cola as contagens
sac.go.contigs.counts <- merge(sac.go.contigs, sac.wgs.normalized.counts, by = "contig")
colnames(sac.go.contigs.counts)[3] <- "GO_terms"

## Cola o dicionário
sac.go.counts.dict <- merge(sac.go.contigs.counts, eggnog_dict, by = "GO_terms")

## Filtra o dicionário
sac.go.counts.dict <- subset(sac.go.counts.dict, description != "molecular_function")
sac.go.counts.dict <- subset(sac.go.counts.dict, description != "biological_process")
sac.go.counts.dict <- subset(sac.go.counts.dict, description != "cellular_component")

## Deixa somente GO e contagens: subconjunto p-value
sac.go.pval <- subset(sac.go.counts.dict, contig %in% rownames(results.sig))

## Remove colunas a mais
colnames(sac.go.pval)
sac.go.pval <- sac.go.pval[,c(10,11,4,5,6,7,8,9)]

## Exporta para analise
write.table(sac.go.pval,
            file="./analises/2019-03-29_STAMP_GO_eggnog_subset_P/go.pval.non-electro.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### Gene ontology para todos por categoria ####
sac.go.counts.dict.bp <- subset(sac.go.counts.dict, ontology == "biological_process")
sac.go.counts.dict.mf <- subset(sac.go.counts.dict, ontology == "molecular_function")
sac.go.counts.dict.cc <- subset(sac.go.counts.dict, ontology == "cellular_component")

### Soma por descrição
sac.go.counts.dict.bp <- ddply(sac.go.counts.dict.bp, "description", numcolwise(sum))
sac.go.counts.dict.mf <- ddply(sac.go.counts.dict.mf, "description", numcolwise(sum))
sac.go.counts.dict.cc <- ddply(sac.go.counts.dict.cc, "description", numcolwise(sum))

### Exporta as três
write.table(sac.go.counts.dict.cc,
            file="./analises/2019-03-29_STAMP_GO_eggnog_subset_P/go.cc.pval.non-electro.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

  
#### EC number com contagens/LFC para o grupo p-val ####
ec_table <- read.delim(file = './analises/2019-02-18_prokka/00_parsed/ec/prkid_ecno.tsv',
                       stringsAsFactors = F)

### Cruza somente aqueles que são p-sig
dim(results.sig.p)

## Cola header de contig
ec_prk_k <- merge(ec_table, contig_X_prokka_table, by="prk_id")

## Somente aqueles p-sig
ec_p_sig <- subset(ec_prk_k, contig %in% row.names(results.sig.p))
print(dim(ec_p_sig)[1])

## Cola os nomes
results.sig.p$contig <- rownames(results.sig.p)
df.sig.p <- as.data.frame(results.sig.p[,c(2,7)])

## Adiciona os L2FC
ec_l2fc <- merge(ec_p_sig, df.sig.p, by = "contig")

## Embaralhar
ec_l2fc <- ec_l2fc[with(ec_l2fc, order(contig, log2FoldChange)), ]

## Somar por EC
ec_l2fc_sum <- ddply(ec_l2fc, "ec_no", numcolwise(sum))

#### EC number com contagens/LFC para geral? ####
# WIP
ec_all <- merge(ec_table, df.sig.p, by = "contig")
ec_all_l2fc_sum <- ddply(ec_all, "ec_no", numcolwise(sum))

# Interrompido.
# Os |LFC| são menos que mil a mais.
write.table(ec_l2fc_sum,
            file="./dados/R_tables/ec_l2fc_sum.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

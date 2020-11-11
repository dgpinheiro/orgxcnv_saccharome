#### DESeq counts for paper

#### 1. Pré-processamento ####
## Bibliotecas
library('xtable')
library('pheatmap')
library('plyr')
library('dendextend')

#### Função: teste de z-score ####

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

### Contagem da taxonomia - saccharome
setwd('/work/rcsilva/projects/saccharome-wgs/')
qval_table <- read.delim(file='analises/2019-03-21_useful_R_observations/out-L2FC.test1.tab',
                         stringsAsFactors = F)


## Quantas receberam anotação funcional
length(qval_table$product) - 453

## Sem proteínas hipotéticas
qval_not_hypo <- subset(qval_table, product != "hypothetical protein")

# xtable(qval_not_hypo)
row.names(qval_not_hypo) <- qval_not_hypo$prokka_id
qval_not_hypo <- qval_not_hypo[,c(1,3,4,5)]
qval_not_hypo <- na.omit(qval_not_hypo)
xtable(qval_not_hypo)

#### 2. Seleção dos sub-conjuntos top 10 up e down ####

## Extrair os menores e os maiores valores
lesser10_l2fc <- sort(qval_not_hypo[,c(4)])[1:10]
top10_l2fc <- tail(sort(qval_not_hypo[,c(4)]), 10)

## Subsetting
qval_down <- subset(qval_not_hypo, L2FC %in% lesser10_l2fc)
qval_up <- subset(qval_not_hypo, L2FC %in% top10_l2fc)

## qval_paper
qval_plot_all <- rbind(qval_down, qval_up)

## Tirar CSV
write.table(qval_plot_all,
            file="./dados/plots-tables/qval_plot.tsv",
            quote = F,
            sep = '\t')

#### 3. Categorias GO para o conjunto Q ####
eggnog_go <- subset(eggnog_out, query_name %in% qval_table$prokka_id)

## 412 com termos
dim(eggnog_go)

## Somente termos e go
eggnog_go <- eggnog_go[,c(1,6)]

## Remove os vazios
eggnog_go <- subset(eggnog_go, GO_terms != "")

## Escreve a tabela
write.table(eggnog_go,
            file="/tmp/saccharome-eggnog-go.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

## Lê novamente a tabela parseada e o dicionário
eggnog_go_prk <- read.delim(file = '/tmp/enog/ready/enog.prk.go.ready.tsv',
                            stringsAsFactors = F)

eggnog_dict <- read.delim(file = '/tmp/enog/ready/out_ready_dict.tsv',
                          stringsAsFactors = F)

colnames(dict_go_merged)[2] <- "prokka_id"

## Junta as tabelas
dict_go_merged <- merge(eggnog_go_prk, eggnog_dict, by = "GO_terms")
go_qval_ready <- merge(qval_table, dict_go_merged, by = "prokka_id")
cpm_go_qval <- merge(go_qval_ready, cpm_table, by = "contig")
cpm_go_qval <- cpm_go_qval[,-c(9)]

## Veja o resultado
View(head(cpm_go_qval))

#### 4. Heatmap do conjunto Q para GO ####
## Separar por processo
cpm_go_bp_q <- subset(cpm_go_qval, ontology == "biological_process")
cpm_go_mf_q <- subset(cpm_go_qval, ontology == "molecular_function")
cpm_go_cc_q <- subset(cpm_go_qval, ontology == "cellular_component")

#### A. Processo biológico ####
bp.ddply <- ddply(cpm_go_bp_q,"description",numcolwise(sum))

### Somente valores numéricos com nome de linha
bp.rows <- bp.ddply[,c(3:8)]
rownames(bp.rows) <- bp.ddply$description
bp.go.q.zscore <- t(apply(bp.rows, 1, cal_z_score))

### Primeiro plot
pheatmap(bp.go.q.zscore,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

### Separando os 10 top | z-score | por tratamento ###
crz.go.q.vec <- as.double(bp.go.q.zscore[,c(1)] + bp.go.q.zscore[,c(2)] + bp.go.q.zscore[,c(3)])
orz.go.q.vec <- as.double(bp.go.q.zscore[,c(4)] + bp.go.q.zscore[,c(5)] + bp.go.q.zscore[,c(6)])

go_diffs <- crz.go.q.vec - orz.go.q.vec
go_outliers <- (-4 > go_diffs | go_diffs > 4)

pheatmap(subset(bp.go.q.zscore, go_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### B. Função molecular ####
mf.ddply <- ddply(cpm_go_mf_q,"description",numcolwise(sum))

### Somente valores numéricos com nome de linha
mf.rows <- mf.ddply[,c(3:8)]
rownames(mf.rows) <- mf.ddply$description
mf.go.q.zscore <- t(apply(mf.rows, 1, cal_z_score))

### Primeiro plot
pheatmap(mf.go.q.zscore,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

### Separando os 10 top | z-score | por tratamento ###
crz.go.mf.q.vec <- as.double(mf.go.q.zscore[,c(1)] + mf.go.q.zscore[,c(2)] + mf.go.q.zscore[,c(3)])
orz.go.mf.q.vec <- as.double(mf.go.q.zscore[,c(4)] + mf.go.q.zscore[,c(5)] + mf.go.q.zscore[,c(6)])

go.mf_diffs <- crz.go.mf.q.vec - orz.go.mf.q.vec
go.mf_outliers <- (-4 > go.mf_diffs | go.mf_diffs > 4)

pheatmap(subset(mf.go.q.zscore, go.mf_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))



#### Análise para o subconjunto P: Pré-processamento ####
### Motivação: Nem BP nem MF deram interessantes! ###
## Separar para o subconjunto P
# Leitura dos cabeçalhos
pval_headers <- read.delim(file='analises/2019-03-21_useful_R_observations/pvalue_protein_headers.1c.txt',
                           stringsAsFactors = F)


eggnog_pval_go <- subset(eggnog_out, query_name %in% pval_headers$prokka_id)

## 412 com termos
dim(eggnog_pval_go)

## Somente termos e go
eggnog_pval_go <- eggnog_pval_go[,c(1,6)]

## Remove os vazios
eggnog_pval_go <- subset(eggnog_pval_go, GO_terms != "")

## Escreve a tabela
write.table(eggnog_pval_go,
            file="/tmp/pval/saccharome-eggnog-pval-go.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

## Lê novamente a tabela parseada e o dicionário
eggnog_go_pval_prk <- read.delim(file = '/tmp/pval/en.2c.tsv',
                            stringsAsFactors = F)

eggnog_pval_dict <- read.delim(file = '/tmp/pval/dict_ready.tsv',
                          stringsAsFactors = F)



## Junta as tabelas
dict_go_pval_merged <- merge(eggnog_go_pval_prk, eggnog_pval_dict, by = "GO_terms")
colnames(dict_go_pval_merged)[2] <- "prk_id"

## Adiciona o nome do contig
go_pval_ready <- merge(dict_go_pval_merged, contig_X_prokka_table, by = "prk_id")

## Adiciona as contagens
cpm_go_pval <- merge(go_pval_ready, cpm_table, by = "contig")

## Remove a coluna length
cpm_go_pval <- cpm_go_pval[,-c(6)]

## Veja o resultado
View(head(cpm_go_pval))

#### Análise para o subconjunto P: Heatmap ####
cpp_go_bp_q <- subset(cpm_go_pval, ontology == "biological_process")
cpp_go_mf_q <- subset(cpm_go_pval, ontology == "molecular_function")
cpp_go_cc_q <- subset(cpm_go_pval, ontology == "cellular_component")

#### A. Processo biológico ####
bp.p.ddply <- ddply(cpp_go_bp_q,"description",numcolwise(sum))

### Somente valores numéricos com nome de linha
bp.p.rows <- bp.p.ddply[,c(2:7)]
rownames(bp.p.rows) <- bp.p.ddply$description
bp.go.p.zscore <- t(apply(bp.p.rows, 1, cal_z_score))

### Primeiro plot
pheatmap(bp.go.p.zscore,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

### Separando os 10 top | z-score | por tratamento ###
crz.go.p.vec <- as.double(bp.go.p.zscore[,c(1)] + bp.go.p.zscore[,c(2)] + bp.go.p.zscore[,c(3)])
orz.go.p.vec <- as.double(bp.go.p.zscore[,c(4)] + bp.go.p.zscore[,c(5)] + bp.go.p.zscore[,c(6)])

go_p_diffs <- crz.go.p.vec - orz.go.p.vec
go_p_outliers <- (-5 > go_p_diffs | go_p_diffs > 5)

pheatmap(subset(bp.go.p.zscore, go_p_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

### Figura pronta!

#### Exportando heatmap (processo biológico) para STAMP ####
head(cpm_go_pval)

### Salva processo biológico ###
write.table(cpp_go_bp_q,
            file="./analises/2019-03-29_STAMP_GO_eggnog_subset_P/cpm_go_bp_pval.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### Análise para o subconjunto T: Pré-processamento ####
### Motivação: Nem BP nem MF deram interessantes! ###
## Separar para o subconjunto P
# Leitura dos cabeçalhos

## 412 com termos
dim(eggnog_out)

## Somente termos e go
eggnog_tval_go <- eggnog_out[,c(1,6)]

## Remove os vazios
eggnog_tval_go <- subset(eggnog_tval_go, GO_terms != "")

## Escreve a tabela
write.table(eggnog_tval_go,
            file="/tmp/tval/saccharome-eggnog-tval-go.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

### TBC

## Lê novamente a tabela parseada e o dicionário
eggnog_go_tval_prk <- read.delim(file = '/tmp/tval/2c.txt',
                                 stringsAsFactors = F)

eggnog_tval_dict <- read.delim(file = '/tmp/tval/dict_out.tsv',
                               stringsAsFactors = F)

## Junta as tabelas
dict_go_tval_merged <- merge(eggnog_go_tval_prk, eggnog_tval_dict, by = "GO_terms")

## Renomeia para consistência
colnames(dict_go_tval_merged)[2] <- "prk_id"

## Adiciona o nome do contig
go_tval_ready <- merge(dict_go_tval_merged, contig_X_prokka_table, by = "prk_id")

## Adiciona as contagens
cpm_go_tval <- merge(go_tval_ready, cpm_table, by = "contig")

## Remove a coluna length
cpm_go_tval <- cpm_go_tval[,-c(6)]

## Veja o resultado
View(head(cpm_go_tval))

#### Análise para o subconjunto T: Heatmap ####
cpp_go_bp_t <- subset(cpm_go_tval, ontology == "biological_process")
cpp_go_mf_t <- subset(cpm_go_tval, ontology == "molecular_function")
cpp_go_cc_t <- subset(cpm_go_tval, ontology == "cellular_component")

#### A. Processo biológico ####
bp.t.ddply <- ddply(cpp_go_bp_t_sig,"description",numcolwise(sum))

### Somente valores numéricos com nome de linha
bp.t.rows <- bp.t.ddply[,c(2:7)]
rownames(bp.t.rows) <- bp.t.ddply$description
bp.go.t.zscore <- t(apply(bp.t.rows, 1, cal_z_score))

### Primeiro plot
pheatmap(bp.go.t.zscore,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

### Separando os 10 top | z-score | por tratamento ###
crz.go.t.vec <- as.double(bp.go.t.zscore[,c(1)] + bp.go.t.zscore[,c(2)] + bp.go.t.zscore[,c(3)])
orz.go.t.vec <- as.double(bp.go.t.zscore[,c(4)] + bp.go.t.zscore[,c(5)] + bp.go.t.zscore[,c(6)])

go_t_diffs <- crz.go.t.vec - orz.go.t.vec
go_t_outliers <- (-4 > go_t_diffs | go_t_diffs > 4)

pheatmap(subset(bp.go.t.zscore, go_t_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

### Figura pronta!

#### Exportando heatmap (processo biológico) para STAMP ####
head(cpm_go_pval)

### Salva processo biológico ###
write.table(cpp_go_bp_t,
            file="./analises/2019-03-29_STAMP_GO_eggnog_subset_P/cpm_go_bp_tval.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

### Teste: soma de colunas ###
rowSums(head(cpp_go_bp_t[,c(6:11)])) > 0.005

### Filtrar
cpp_go_bp_t_sig <- subset(cpp_go_bp_t, rowSums(cpp_go_bp_t[,c(6:11)]) > 1)
dim(cpp_go_bp_t)
dim(cpp_go_bp_t_sig)

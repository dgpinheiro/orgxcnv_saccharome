## Plot com novos KOs
#### Pré-processamento ####

## Setup
setwd('/work/rcsilva/projects/saccharome-wgs')
load("/data/rcsilva/projects/saccharome-wgs/dados/RData/KO-plots.RData")

## Tabela de contagens
head(cpm_table)

## Carrega o dicionário do KEGG
kegg_dict_main <- read.delim('/data/db/metagenomics/db-KEGG/ready_R/kegg_dictionary.tsv',
                             stringsAsFactors = F)

## Carrega a nova análise
eggnog_out <- read.delim('/work/rcsilva/projects/saccharome-wgs/analises/2019-03-25_eggnog_stringent/annotations.tsv',
                         stringsAsFactors = F)

## Separa somente ID e kegg
eggnog_id_ko <- eggnog_out[,c(1,7)]

## Função z-score
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

## Salva a tabela com os IDs para parse no shell
write.table(x = eggnog_id_ko,
            file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-25_eggnog_stringent/to_parse.tsv',
            quote = F,
            row.names = F,
            sep = '\t'
            )

## Re-escreve a tabela separada
# Colocar caminhos, rodar fora do contrário:
# system('~/bin/parse-eggnog-two-columns ./to_parse.tsv > ./out_parsed.tsv')

## Lê de novo a tabela de KOs
eggnog_id_ko_parsed <- read.delim('./analises/2019-03-25_eggnog_stringent/out_parsed_ready.tsv',
                                  stringsAsFactors = F)

colnames(eggnog_id_ko_parsed)[1] <- "prk_id"

#### Geração da tabela para STAMP ####
# Tabela do prokka
head(contig_X_prokka_table)

# Unir nomes dos contigs com ids
en_prk_id_ko <- merge(contig_X_prokka_table, eggnog_id_ko_parsed,
                      by = "prk_id")

head(en_prk_id_ko)
colnames(en_prk_id_ko)[3] <- "KO"

# Unir a tabela com dicionário
en_pre_counts <- merge(en_prk_id_ko, kegg_dict_main, by = "KO")
head(en_pre_counts)

# Remover os blanks
en_pre_counts <- subset(en_pre_counts, KO != "")
head(en_pre_counts)
dim(en_pre_counts)


# Unir o dicionário com contagens
head(cpm_table)
head(en_pre_counts)
colnames(en_pre_counts)[3]
colnames(cpm_table)[1] <- "contig"

## Adiciona as contagens
en_out <- merge(en_pre_counts, cpm_table, by = "contig")
head(en_out)

## Bug: definição (shell)
en_out[3494,]$definition
grep("\n", en_out$definition)

## Limpa as definições quebradas
en_out_clean <- en_out[-c(grep("\n", en_out$definition)), ]
grep("\tK", en_out_clean$definition)

## Salva para o STAMP
write.table(
  x=en_out_clean,
  file='/work/rcsilva/projects/saccharome-wgs/analises/2019-03-27_STAMP_KO_eggnog_stringent/out_clean.tsv',
  quote = F,
  row.names = F,
  sep = '\t'
)

## Remove o comprimento
head(en_out_clean)
en_out_clean <- en_out_clean[,-c(8)]

## Redução da tabela para análise no STAMP
en_def_unique <- ddply(en_out_clean, "definition", numcolwise(sum))
en_hd_unique <- ddply(en_out_clean, "header", numcolwise(sum))
en_hd2_unique <- ddply(en_out_clean, "header2", numcolwise(sum))
en_ptw_unique <- ddply(en_out_clean, "pathway", numcolwise(sum))

## Salva separadamente
write.table(
  x=en_def_unique,
  file='./analises/2019-03-27_STAMP_KO_eggnog_stringent/stamp_ready/def.tsv',
  quote = F,
  row.names = F,
  sep = '\t'
)

write.table(
  x=en_hd_unique,
  file='./analises/2019-03-27_STAMP_KO_eggnog_stringent/stamp_ready/hd.tsv',
  quote = F,
  row.names = F,
  sep = '\t'
)

write.table(
  x=en_hd2_unique,
  file='./analises/2019-03-27_STAMP_KO_eggnog_stringent/stamp_ready/hd2.tsv',
  quote = F,
  row.names = F,
  sep = '\t'
)

write.table(
  x=en_ptw_unique,
  file='./analises/2019-03-27_STAMP_KO_eggnog_stringent/stamp_ready/ptw.tsv',
  quote = F,
  row.names = F,
  sep = '\t'
)



#### Geração do heatmap #####
## Gera matriz sem repetições
en_definitions_unique <- ddply(en_out_clean, "definition", numcolwise(sum))
en_definitions_unique <- en_definitions_unique[,-c(2)]

## Remove coluna de nome, adiciona externamente
rownames(en_definitions_unique) <- en_definitions_unique$definition
en_definitions_unique <- en_definitions_unique[,-c(1)]

## Aplica z-score
en_norm <- t(apply(en_definitions_unique, 1, cal_z_score))

pheatmap(en_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "KO all",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

#### ko_description -5 > Z > 5
crz_en_vec <- as.double(en_norm[,c(1)] + en_norm[,c(2)] + en_norm[,c(3)])
orz_en_vec <- as.double(en_norm[,c(4)] + en_norm[,c(5)] + en_norm[,c(6)])

en_ko_diffs <- crz_en_vec - orz_en_vec
en_outliers <- (-5 > en_ko_diffs | en_ko_diffs > 5)
#View(subset(bp_heatmap_norm, bp_outliers))

pheatmap(subset(en_norm, en_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         #main = "Etiquetas GO de maior desvio para processo biológico",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

### Para rotas: en_ptw_unique

# Arruma os nomes de linhas
rownames(en_ptw_unique) <- en_ptw_unique$pathway
en_ptw_unique <- en_ptw_unique[,-c(1)]
en_ptw_norm <- t(apply(en_ptw_unique, 1, cal_z_score))

crz_en_ptw_vec <- as.double(en_ptw_norm[,c(1)] + en_ptw_norm[,c(2)] + en_ptw_norm[,c(3)])
orz_en_ptw_vec <- as.double(en_ptw_norm[,c(4)] + en_ptw_norm[,c(5)] + en_ptw_norm[,c(6)])

en_ko_ptw_diffs <- crz_en_ptw_vec - orz_en_ptw_vec
en_ptw_outliers <- (-4 > en_ko_ptw_diffs | en_ko_ptw_diffs > 4)

pheatmap(subset(en_ptw_norm, en_ptw_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         #main = "Etiquetas GO de maior desvio para processo biológico",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

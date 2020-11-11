# Parse KO tables into heatmaps/cluster pics
library('pheatmap')
library('plyr')
library('dendextend')

#### Teste de z-score ####

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

#### Dados de entrada ####

## 1. Matriz de contagens
cpm_table <- read.delim("/work/rcsilva/projects/saccharome-wgs/dados/subsets/cpm_all_samples.tsv",
                        stringsAsFactors = F)

## 2. IDs (prokka) com anotações KEGG
prokka_X_KEGG_table <- read.delim("/work/rcsilva/projects/saccharome-wgs/analises/2019-03-22_parse_kegg_table/kegg_headers_correct4.tsv",
                                  stringsAsFactors = F)

## 3. Dicionário de termos KEGG
kegg_dict_table <- read.delim("/work/rcsilva/projects/saccharome-wgs/analises/2019-03-22_parse_kegg_table/kegg_dict_uniq_v2.txt",
                              stringsAsFactors = F)

## 4. Relação entre contigs e IDs do prokka
contig_X_prokka_table <- read.delim('/work/rcsilva/projects/saccharome-wgs/analises/2019-02-18_prokka/gene-to-protein_ids.txt',
                                    stringsAsFactors = F)

#### Subsetting ####
subset_cxp_table <- subset(contig_X_prokka_table,
                           contig_X_prokka_table$prk_id %in% prokka_X_KEGG_table$prk_id)

subset_cpm_table <- subset(cpm_table,
                           cpm_table$target_id %in% subset_cxp_table$contig)

subset_cpm_table <- subset_cpm_table[,c(1,3,4,5,6,7,8)]
colnames(subset_cpm_table)[1] <- c("contig")
# ead(subset_cpm_table)

#### Merging tables ####
cpm_prokka <- merge(subset_cxp_table, subset_cpm_table, by="contig")
cpm_KO <- merge(cpm_prokka, prokka_X_KEGG_table, by="prk_id")
cpm_KO_dict <- merge(cpm_KO, kegg_dict_table, by="ko_number")

# Final output (pré-heatmap): cpm_KO_dict

#### Heatmap por rotas ####
cpm_pathways <- ddply(cpm_KO_dict,"pathway",numcolwise(sum))

### Somente valores numéricos com nome de linha
cpm_pathways_names <- cpm_pathways[,c(2:7)]
rownames(cpm_pathways_names) <- cpm_pathways$pathway
cpm_pathways_norm <- t(apply(cpm_pathways_names, 1, cal_z_score))

### Primeiro plot
pheatmap(cpm_pathways_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

#### Separando os 10 top | z-score | por tratamento ####
crz_ko_vec <- as.double(cpm_pathways_norm[,c(1)] + cpm_pathways_norm[,c(2)] + cpm_pathways_norm[,c(3)])
orz_ko_vec <- as.double(cpm_pathways_norm[,c(4)] + cpm_pathways_norm[,c(5)] + cpm_pathways_norm[,c(6)])

ko_diffs <- crz_ko_vec - orz_ko_vec
ko_outliers <- (-4 > ko_diffs | ko_diffs > 4)

pheatmap(subset(cpm_pathways_norm, ko_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Heatmap por nomes ####
cpm_gnames <- ddply(cpm_KO_dict,"name",numcolwise(sum))

### Somente valores numéricos com nome de linha
cpm_gnames_names <- cpm_gnames[,c(2:7)]
rownames(cpm_gnames_names) <- cpm_gnames$name
cpm_gnames_norm <- t(apply(cpm_gnames_names, 1, cal_z_score))

### Primeiro plot
pheatmap(cpm_gnames_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Test",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

#### Separando os 10 top | z-score | por tratamento ####
crz_kog_vec <- as.double(cpm_gnames_norm[,c(1)] + cpm_gnames_norm[,c(2)] + cpm_gnames_norm[,c(3)])
orz_kog_vec <- as.double(cpm_gnames_norm[,c(4)] + cpm_gnames_norm[,c(5)] + cpm_gnames_norm[,c(6)])

kog_diffs <- crz_kog_vec - orz_kog_vec
kog_outliers <- (-5 > kog_diffs | kog_diffs > 5)

pheatmap(subset(cpm_gnames_norm, kog_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Ortólogos do KEGG com |Z-score| maior que 5",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Exportando para testes no STAMP ####
write.table(x = cpm_pathways,
            file="/work/rcsilva/projects/saccharome-wgs/analises/2019-03-24_stamp_KO_pathways/kotable.tsv",
            quote = F,
            row.names = F,
            sep = "\t")

#### Para o microbiome analyst ####
ko_to_MA <- cpm_KO_dict[,c(1,4,5,6,7,8,9)]
ko_to_MA_uniq <- unique(ko_to_MA)

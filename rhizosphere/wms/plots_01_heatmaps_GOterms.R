## Manual heatmap: GO terms ##
library('pheatmap')
library(plyr)
library(dendextend)
# ddply(df,"x",numcolwise(sum))

#### Tutoriais ####
# Dave Tang: https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

## Leitura das matrizes GO
bp_mat <- read.delim(
  file = '/data/rcsilva/projects/saccharome-wgs/analises/2019-03-20_heatmap/inputs/bp.out.txt',
  header = T,
  stringsAsFactors = F)

cc_mat <- read.delim(
  file = '/data/rcsilva/projects/saccharome-wgs/analises/2019-03-20_heatmap/inputs/cc.out.txt',
  header = T)

mf_mat <- read.delim(
  file = '/data/rcsilva/projects/saccharome-wgs/analises/2019-03-20_heatmap/inputs/mf.out.txt',
  header = T)

#### Processo biológico ####

## Agregar os valores
### Agrega pela anotação
bp_uniques <- ddply(bp_mat,"GO_Annotation",numcolwise(sum))

### Separa as colunas de contagem, adicionando os nomes como nomes de linha
bp_heatmap_mat <- bp_uniques[,c(2:7)]
rownames(bp_heatmap_mat) <- bp_uniques[,1]

## Gera o heatmap
## Gera o heatmap
pheatmap(bp_heatmap_mat,
         clustering_method = "average",
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Biological Process",
         cellheight = 10,
         cellwidth = 20,
         #kmeans_k = 70,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Função molecular ####

## Agregar os valores
### Agrega pela anotação
mf_uniques <- ddply(mf_mat,"GO_Annotation",numcolwise(sum))

### Separa as colunas de contagem, adicionando os nomes como nomes de linha
mf_heatmap_mat <- mf_uniques[,c(2:7)]
rownames(mf_heatmap_mat) <- mf_uniques[,1]

## Gera o heatmap
pheatmap(mf_heatmap_mat,
         clustering_method = "average",
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Molecular Function",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

## Agregar os valores
### Agrega pela anotação
cc_uniques <- ddply(cc_mat,"GO_Annotation",numcolwise(sum))

### Separa as colunas de contagem, adicionando os nomes como nomes de linha
cc_heatmap_mat <- cc_uniques[,c(2:7)]
rownames(cc_heatmap_mat) <- cc_uniques[,1]
cc_go_tags <- cc_mat[,4]
names(cc_go_tags) <- cc_mat[,2]

## Gera o heatmap
pheatmap(head(cc_heatmap_mat, 80),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Molecular Function",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

  
#### Teste de z-score ####

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

## Componente celular
cc_heatmap_norm <- t(apply(cc_heatmap_mat, 1, cal_z_score))
pheatmap(cc_heatmap_norm)

pheatmap(cc_heatmap_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Cellular Component",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "tomato"))(100))

## Função molecular
mf_heatmap_norm <- t(apply(mf_heatmap_mat, 1, cal_z_score))
pheatmap(mf_heatmap_norm)

pheatmap(mf_heatmap_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Molecular Function",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

## PB
bp_heatmap_norm <- t(apply(bp_heatmap_mat, 1, cal_z_score))
pheatmap(bp_heatmap_norm)

pheatmap(bp_heatmap_norm,
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "All GO tags by Biological Process",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Separando os 10 top | z-score | por função ####

## Matrizes
bp_heatmap_norm
cc_heatmap_norm
mf_heatmap_norm

#### Processo biológico - módulo de Z maior que 5 ####
crz_bp_vec <- as.double(bp_heatmap_norm[,c(1)] + bp_heatmap_norm[,c(2)] + bp_heatmap_norm[,c(3)])
orz_bp_vec <- as.double(bp_heatmap_norm[,c(4)] + bp_heatmap_norm[,c(5)] + bp_heatmap_norm[,c(6)])

bp_diffs <- crz_bp_vec - orz_bp_vec
bp_outliers <- (-5 > bp_diffs | bp_diffs > 5)
#View(subset(bp_heatmap_norm, bp_outliers))

pheatmap(subset(bp_heatmap_norm, bp_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Etiquetas GO de maior desvio para processo biológico",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Função molecular - módulo de Z maior que 5 ####
crz_mf_vec <- as.double(mf_heatmap_norm[,c(1)] + mf_heatmap_norm[,c(2)] + mf_heatmap_norm[,c(3)])
orz_mf_vec <- as.double(mf_heatmap_norm[,c(4)] + mf_heatmap_norm[,c(5)] + mf_heatmap_norm[,c(6)])

mf_diffs <- crz_mf_vec - orz_mf_vec
mf_outliers <- (-5 > mf_diffs | mf_diffs > 5)
#View(subset(mf_heatmap_norm, mf_outliers))

pheatmap(subset(mf_heatmap_norm, mf_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Etiquetas GO de maior desvio para função molecular",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#### Componente - módulo de Z maior que 5 ####
crz_cc_vec <- as.double(cc_heatmap_norm[,c(1)] + cc_heatmap_norm[,c(2)] + cc_heatmap_norm[,c(3)])
orz_cc_vec <- as.double(cc_heatmap_norm[,c(4)] + cc_heatmap_norm[,c(5)] + cc_heatmap_norm[,c(6)])

cc_diffs <- crz_cc_vec - orz_cc_vec
cc_outliers <- (-4 > cc_diffs | cc_diffs > 4)
#View(subset(cc_heatmap_norm, cc_outliers))

pheatmap(subset(cc_heatmap_norm, cc_outliers),
         clustering_method = "average",
         #annotation_row = cc_go_tags,
         scale = "row",
         legend = T,
         display_numbers = T,
         main = "Etiquetas GO de maior desvio para componente celular, |Z| > 4",
         cellheight = 10,
         cellwidth = 20,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))


## Diagrame o Venn!
library("plyr")

##### 1. Soma total #####

# Leia a matriz de gênero
migf_genus = read.delim(
  file = "/data/AEM/dados/2019-08-21_kraken_classification_rcsilva/migf/wms_taxonomy_genus_final.tsv",
  stringsAsFactors = F
)

## Agregar por taxonomia
collapsed_genus = ddply(migf_genus,
                        "taxonomy",
                        numcolwise(sum))

## Inserir soma das colunas
collapsed_genus$CRZ_sum <- rowSums(x = collapsed_genus[,c(2,3,4)])
collapsed_genus$ORZ_sum <- rowSums(x = collapsed_genus[,c(5,6,7)])

## Exportar essa tabela
write.table(collapsed_genus,
            file = "../saccharome-wgs/dados/R_tables/2019_08_28_collapsed_genus_for_venn.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

##### 2. Três contagens #####
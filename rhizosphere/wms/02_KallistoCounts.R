## Capturando a pasta de trabalho
setwd("/work/rcsilva/projects/saccharome-wgs")

## CPM matrix, diretório:
# "analises/2019-02-26_kallisto_count_megahit01/count/cpm_all_samples.tsv"

## Importando a matriz
cpm_mat <- read.table(
  "analises/2019-02-26_kallisto_count_megahit01/count/cpm_all_samples.tsv",
  header = T,
  stringsAsFactors = F)

## Retirando o tamanho da sequência, deixando só valores numéricos
cpm_clean <- cpm_mat[,c(3,4,5,6,7,8)]

## Colando os nomes das proteínas
row.names(cpm_clean) <- cpm_mat$target_id

## Obter a soma das linhas? (maybe too large)
cpm_clean[,7] <- rowSums(cpm_clean)
names(cpm_clean$V7) <- "Soma"

## Filtrar: Somente com mais de 1 CPM
cpm_mto <- subset(cpm_clean, V7 > 1)

### Permanova: dados de solo em relação aos termos
setwd('/work/rcsilva/projects/saccharome-wgs/dados/R_permanova')

## Objetivo:
# - Cruzar os dados físico-químicos com a análise de termos (KO)

## Bibliotecas
library('vegan')

## Leitura das matrizes
soil_params <- read.delim('./dados/R_tables/R_permanova/params.tsv')

### Filtrar a matriz
soil_tab <- soil_params[,c(1,5,6,7,8,9,11,12,13,14,15)]

### Manter o nome
rownames(soil_tab) <- soil_tab$param

### Remover a primeira coluna
soil_tab <- soil_tab[,c(4,5,6,9,10,11)]
soil_tab <- t(soil_tab)

### Nomear as colunas como na amostragem
colnames(soil_tab) <- c("CRZ2",
                        "CRZ3",
                        "CRZ4",
                        "ORZ2",
                        "ORZ3",
                        "ORZ4")


### Test: simpleRDA
plot(rda(t(cpm_gnames[,c(2:7)]), soil_tab),
     display = c("sites", "cn"),
     scaling = 1
)

### Test all factors
plot(rda(t(cpm_gnames[,c(8:12)]) ~ . , data = as.data.frame(soil_tab)),
     display = c("sites", "cn"),
     scaling = 2
     )

        
### RDA
rda.soil <- rda(t(cpm_gnames[,c(2:7)]), soil_tab)


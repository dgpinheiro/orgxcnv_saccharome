## Mantel Test: Saccharome

# Tutoriais indicados
# Vegan: https://www.rdocumentation.org/packages/vegan/versions/2.4-2/topics/mantel
# Legendre & Legendre (classiqueira), tanto o numérico quanto o in R
# G. Zahn: http://geoffreyzahn.com/getting-your-otu-table-into-r/

# Diretórios e bibliotecas
setwd('/work/rcsilva/projects/saccharome-wgs/')
library('ade4')
library('vegan')
library('adespatial')
library('cluster')
library('FD')
library('gclus')
library('biomformat')

# Importando BIOM - instalação de pacotes
# devtools::install_github("biomformat", "joey711")


  #### Importando e massageando os dados ####

# Leitura das matrizes
# Matriz de taxonomia: Kraken
biom_path <- '/work/rcsilva/projects/saccharome-wgs/analises/2019-04-29_kraken_all_official/all.biom'
biom_matrix <- read_biom(biom_path)
colnames(biom_matrix) <- c("CRZ2", "CRZ3", "CRZ4", "ORZ2", "ORZ3", "ORZ4")

# Ver os dados:
biom_data(biom_matrix)
observation_metadata(biom_matrix)


  
# Matriz físico-química
soil.df.wgs <- soil.df[,c(1,2,6,7,8,12,13,14)]
colnames(soil.df.wgs) <- c("parameters", "description",
                           "CRZ2", "CRZ3", "CRZ4",
                           "ORZ2", "ORZ3", "ORZ4")

## Teste de mantel ##
mantel.partial()
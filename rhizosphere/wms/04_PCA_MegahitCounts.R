## Script 4
## Realizar o PCA nas contagens Kallisto <- Megahit

## Separar somente as contagens, sem as somas
cpm <- cpm_mto[,c(1,2,3,4,5,6)]

# t
# princomp(~ ., data = cpm, cor = TRUE)
# princomp(cpm)

### Somente para 06 observações
# Separa somente 10 observações
test10 <- head(cpm, 20)

#### 1. Experimento: PCA com dados toys ####
setwd("/work/rcsilva/projects/saccharome-wgs/analises/2019-02-28_aprenda_PCA/data")

## Instalação das faltantes
install.packages("missMDA")
install.packages("FactoMineR")

## Bibliotecas do livro
library(ade4)
library(vegan)
library(gclus)
library(ape)
library(missMDA)
library(FactoMineR)

# Importando as funções do livro
source("./Functions/cleanplot.pca.R")
source("./Functions/PCA.newr.R")
source("./Functions/CA.newr.R")

# Importando os dados do livro
load("./Data/Doubs.RData")
load("./Data/mite.RData")

# Removendo coluna vazia
spe <- spe[-8, ]
env <- env[-8, ]
spa <- spa[-8, ]

# PCA da matriz de correlação
## 
env.pca <- rda(env, scale=T)
summary(env.pca)

## Extrair eigenvalores
ev <- env.pca$CA$eig
screeplot(env.pca,
          bstick = T,
          npcs = length(env.pca$CA$eig))

par(mfrow = c(1,2))
biplot(env.pca, scaling = 1, main = "PCA - escala 1")
biplot(env.pca, main = "PCA - escala 2")

#### 2. Nosso PCA: Somente 10 observações ####
# Objeto de análise: test10

## Transpor a matriz para ficar idêntica a env
t.test10 <- transpose(test10)

### Tentar recuperar os nomes
colnames(t.test10) = rownames(test10)
rownames(t.test10) = colnames(test10)

## Realiza a análise de redundância
test10.pca <- rda(t.test10, scale=T)

## Sumário dos resultados obtidos
summary(test10.pca)

## Extrai os eigenvalores
t10.ev <- test10.pca$CA$eig

## Plota a importância dos eigenvalores
screeplot(test10.pca,
          bstick = T,
          npcs = length(test10.pca$CA$eig))

biplot(test10.pca, scaling=3)

#### 3. PCOA: Todas as observações ####

## Aplica a transformação de Hellinger
spe.h <- decostand(spe, "hellinger")



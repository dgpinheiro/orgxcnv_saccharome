# Top level

## Bibliotecas
library('xtable')

### Contagem da taxonomia - saccharome
setwd('/work/rcsilva/projects/saccharome-wgs/')
qval_table <- read.delim(file='analises/2019-03-21_useful_R_observations/out-L2FC.test1.tab',
                         stringsAsFactors = F)

count(qval_table$L2FC > 0)

## Quantas receberam anotação funcional
length(qval_table$product) - 453

## Sem proteínas hipotéticas
qval_not_hypo <- subset(qval_table, product != "hypothetical protein")

# xtable(qval_not_hypo)
row.names(qval_not_hypo) <- qval_not_hypo$prokka_id
qval_not_hypo <- qval_not_hypo[,c(1,3,4,5)]
qval_not_hypo <- na.omit(qval_not_hypo)

xtable(qval_not_hypo)

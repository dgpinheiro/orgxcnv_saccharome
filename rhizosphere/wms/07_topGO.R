## TopGO para metagenomas

## Verifica se precisa instalar
# source("https://bioconductor.org/biocLite.R")
# biocLite("topGO")
# biocLite("ALL")
# install.packages("GOplot")

## Tutorial:
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/topGO/inst/doc/topGO.pdf
# https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

## Evoca a biblioteca
# library('topGO')
# library('ALL')
library('GOplot')

# Tutorial
data(EC)
head(EC$david)
head(EC$genelist)


circ <- circle_dat(EC$david, EC$genelist)


# Gráfico de bolhas
GOBubble(circ,
         title = 'Bubble plot with background colour',
         display = 'multiple',
         bg.col = T,
         labels = 3)

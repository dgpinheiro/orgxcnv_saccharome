## Merge Taxa - Kaiju

# Bibliotecas
library(data.table)

# Diretório
setwd('/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/output_per_sample/taxa-levels/1_phy/prok/out-test')

#### Phylum ####
taxlevel = "Phylum"

# Tabelas
pCRZ2 <- read.delim(file = "./CRZ2.tsv.out",
                    stringsAsFactors = F)

pCRZ3 <- read.delim(file = "./CRZ3.tsv.out",
                    stringsAsFactors = F)

pCRZ4 <- read.delim(file = "./CRZ4.tsv.out",
                    stringsAsFactors = F)

pORZ2 <- read.delim(file = "./ORZ2.tsv.out",
                    stringsAsFactors = F)

pORZ3 <- read.delim(file = "./ORZ3.tsv.out",
                    stringsAsFactors = F)

pORZ4 <- read.delim(file = "./ORZ4.tsv.out",
                    stringsAsFactors = F)

#### Juntando as tabelas ####
tmp <- merge(pCRZ2, pCRZ3, by = taxlevel, all = T)
tmp <- merge(tmp, pCRZ4, by = taxlevel, all = T)
tmp <- merge(tmp, pORZ2, by = taxlevel, all = T)
tmp <- merge(tmp, pORZ3, by = taxlevel, all = T)
tmp <- merge(tmp, pORZ4, by = taxlevel, all = T)

## Remove todas as colunas extras de reino
# tmp <- tmp[,-c(grep("Kingdom", colnames(tmp))[1:length(grep("Kingdom", colnames(tmp)))-1])]
tmp <- tmp[,-c(grep("Kingdom", colnames(tmp)))]

# Cria o caminho
path <- paste(taxlevel, "ready", "tsv", sep = ".")

tmp[is.na(tmp)] <- 0

# Salva a matriz
write.table(
  x=tmp,
  file=path,
  quote = F,
  row.names = F,
  sep = '\t'
)

#### Gênero ####
setwd('/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/output_per_sample/taxa-levels/5_gen/prok/ready')

#### Phylum ####
taxlevel = "Genus"

# Tabelas
gCRZ2 <- read.delim(file = "./CRZ2.tsv.out",
                    stringsAsFactors = F)

gCRZ3 <- read.delim(file = "./CRZ3.tsv.out",
                    stringsAsFactors = F)

gCRZ4 <- read.delim(file = "./CRZ4.tsv.out",
                    stringsAsFactors = F)

gORZ2 <- read.delim(file = "./ORZ2.tsv.out",
                    stringsAsFactors = F)

gORZ3 <- read.delim(file = "./ORZ3.tsv.out",
                    stringsAsFactors = F)

gORZ4 <- read.delim(file = "./ORZ4.tsv.out",
                    stringsAsFactors = F)

# Somente os gêneros
gCRZ2 <- gCRZ2[,c(1,dim(gCRZ2)[2])]
gCRZ3 <- gCRZ3[,c(1,dim(gCRZ3)[2])]
gCRZ4 <- gCRZ4[,c(1,dim(gCRZ4)[2])]
gORZ2 <- gORZ2[,c(1,dim(gORZ2)[2])]
gORZ3 <- gORZ3[,c(1,dim(gORZ3)[2])]
gORZ4 <- gORZ4[,c(1,dim(gORZ4)[2])]

# Limpa os gêneros vazios
gCRZ2 <- subset(gCRZ2, Genus != "" )
gCRZ3 <- subset(gCRZ3, Genus != "" )
gCRZ4 <- subset(gCRZ4, Genus != "" )
gORZ2 <- subset(gORZ2, Genus != "" )
gORZ3 <- subset(gORZ3, Genus != "" )
gORZ4 <- subset(gORZ4, Genus != "" )

#### Juntando as tabelas ####
tmp <- merge(gCRZ2, gCRZ3, by = "Genus", all = T)
tmp <- merge(tmp, gCRZ4, by = taxlevel, all = T)
tmp <- merge(tmp, gORZ2, by = taxlevel, all = T)
tmp <- merge(tmp, gORZ3, by = taxlevel, all = T)
tmp <- merge(tmp, gORZ4, by = taxlevel, all = T)

### Método DT ###
dt.c2 <- data.table(gCRZ2, key="Genus")
dt.c3 <- data.table(gCRZ3, key="Genus")
dt.c4 <- data.table(gCRZ4, key="Genus")
dt.o2 <- data.table(gORZ2, key="Genus")
dt.o3 <- data.table(gORZ3, key="Genus")
dt.o4 <- data.table(gORZ4, key="Genus")

df.tmp <- merge(dt.c2, dt.c3, all = T)
df.tmp <- merge(df.tmp, dt.c4, all = T)
df.tmp <- merge(df.tmp, dt.o2, all = T)
df.tmp <- merge(df.tmp, dt.o3, all = T)
df.tmp <- merge(df.tmp, dt.o4, all = T)

# Cria o caminho
path <- paste(taxlevel, "ready", "tsv", sep = ".")

df.tmp[is.na(df.tmp)] <- 0

# Salva a matriz
write.table(
  x=df.tmp,
  file=path,
  quote = F,
  row.names = F,
  sep = '\t'
)



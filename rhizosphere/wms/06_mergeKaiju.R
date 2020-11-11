## Kaiju tab to matrix

#### Abrindo os arquivos ####

ORZ2 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/ORZ2.ok'
)

ORZ3 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/ORZ3.ok'
)

ORZ4 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/ORZ4.ok'
)

CRZ2 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/CRZ2.ok'
)

CRZ3 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/CRZ3.ok'
)

CRZ4 <- read.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/taxout/genus/CRZ4.ok'
)


#### Juntando os gêneros ####
all.genus <- as.vector(unique(CRZ2$V2, CRZ3$V2, CRZ4$V2, ORZ2$V2, ORZ3$V2, ORZ4$V2))

## Teste do subconjunto
subset(CRZ2, CRZ2$V2 %in% all.genus)

## Cria a tabela vazia
genre.table <- data.frame(
  all.genus,
  "CRZ2" = 0,
  "CRZ3" = 0,
  "CRZ4" = 0,
  "ORZ2" = 0,
  "ORZ3" = 0,
  "ORZ4" = 0,
  stringsAsFactors = F)

### Adiciona os valores à tabela

counter <- 1
for (genus in genre.table$all.genus) {
  ## CRZ2
  if (genus %in% CRZ2$V2) {
    genre.table$CRZ2[counter] <- subset(CRZ2$V1, CRZ2$V2 == genus)[1]
  } else {
    genre.table$CRZ2[counter] <- 0
  }
  counter <- counter + 1
}


counter <- 1
for (genus in genre.table$all.genus) {
  ## CRZ3
  if (genus %in% CRZ3$V2) {
    genre.table$CRZ3[counter] <- subset(CRZ3$V1, CRZ3$V2 == genus)[1]
  } else {
    genre.table$CRZ3[counter] <- 0
  }
  counter <- counter + 1
  print(paste("Counter is", counter, "and genus is:", genus))
}

counter <- 1
for (genus in genre.table$all.genus) {
  ## CRZ4
  if (genus %in% CRZ4$V2) {
    genre.table$CRZ4[counter] <- subset(CRZ4$V1, CRZ4$V2 == genus)[1]
  } else {
    genre.table$CRZ4[counter] <- 0
  }
  counter <- counter + 1
}

counter <- 1
for (genus in genre.table$all.genus) {
  ## ORZ4
  if (genus %in% ORZ4$V2) {
    genre.table$ORZ4[counter] <- subset(ORZ4$V1, ORZ4$V2 == genus)[1]
  } else {
    genre.table$ORZ4[counter] <- 0
  }
  counter <- counter + 1
}

counter <- 1
for (genus in genre.table$all.genus) {
  ## ORZ3
  if (genus %in% ORZ3$V2) {
    genre.table$ORZ3[counter] <- subset(ORZ3$V1, ORZ3$V2 == genus)[1]
  } else {
    genre.table$ORZ3[counter] <- 0
  }
  counter <- counter + 1
}

counter <- 1
for (genus in genre.table$all.genus) {
  ## ORZ2
  if (genus %in% ORZ2$V2) {
    genre.table$ORZ2[counter] <- subset(ORZ2$V1, ORZ2$V2 == genus)[1]
  } else {
    genre.table$ORZ2[counter] <- 0
  }
  counter <- counter + 1
}

#### Escrever tabela para o STAMP ####

write.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/final_out.tsv',
  x = genre.table,
  col.names = T,
  quote = F,
  row.names = F,
  sep='\t'
)

#### Subconjunto ao menos 1000 sequências ####
seqs1000 <- subset(genre.table, rowSums(genre.table[,c(2,3,4,5,6,7)]) > 1000)

write.table(
  file = '/work/rcsilva/projects/saccharome-wgs/analises/2019-03-12_kaiju-megahit-prodigal/final_out_min100.tsv',
  x = seqs1000,
  col.names = T,
  quote = F,
  row.names = F,
  sep='\t'
)


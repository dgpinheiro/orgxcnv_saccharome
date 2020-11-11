### KOs e testes ###
library(ggplot2)        # plotting & data
library(dplyr)          # data manipulation
library(tidyr)          # data re-shaping
library(magrittr)       # pipe operator
library(gridExtra)      # provides side-by-side plotting


## Matriz:
head(cpm_KO_dict)

#### Todos com infos. uteis ####
cpm_KO_dict_ok <- cpm_KO_dict[,c(1,11,12,13,4,5,6,7,8,9)]

## Exportar a rota do enxofre
write.table(cpm_KO_dict_ok,
            file="/tmp/sakaromi/ko_all_lowmem.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### KO: rota do Enxofre ####

## Separa só o metabolismo de enxofre
cpm_ko_sulfur <- subset(cpm_KO_dict,
                        pathway == " Sulfur metabolism [PATH:ko00920]")


## Agregar por nome
sulfur_ko <- ddply(cpm_ko_sulfur,"name",numcolwise(sum))
path_ko <- sulfur_ko <- ddply(cpm_KO_dict,"pathway",numcolwise(sum))

## Inverter a tabela
cpm_ko_sulfur_ok <- cpm_ko_sulfur[,c(3,2,1,10,11,12,13,4,5,6,7,8,9)]

## Exportar a rota do enxofre
write.table(cpm_ko_sulfur_ok,
            file="/tmp/sakaromi/sulfur.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### KO: rota do Tolueno ####

## Separa só o metabolismo de enxofre
cpm_ko_toluene <- subset(cpm_KO_dict,
                        pathway == " Toluene degradation [PATH:ko00623]")

## Inverter a tabela
cpm_ko_toluene_ok <- cpm_ko_toluene[,c(3,2,1,10,11,12,13,4,5,6,7,8,9)]

## Exportar a rota do enxofre
write.table(cpm_ko_toluene_ok,
            file="/tmp/sakaromi/ko-toluene.tsv",
            quote = F,
            sep = '\t',
            row.names = F)


#### KO: somente ID (teste) ####

only_ko <- ddply(cpm_KO_dict,
                 "ko_number",
                 numcolwise(sum))

## tirar CPM
only_ko[,c(2,3,4,5,6,7)] <- only_ko[,c(2,3,4,5,6,7)] * 1000000

# Check
head(only_ko)

# Export
write.table(only_ko,
            file="/tmp/sakaromi/only_ko.tsv",
            quote = F,
            sep = '\t',
            row.names = F)

#### Teste T na mão: small size ####
# Tabela filtrada: cpm_KO_dict_ok
small <- head(cpm_KO_dict_ok, 30)

# Extrai as médias
small$CRZ_mean <- apply(small[,c(5,6,7)], 1, mean)
small$ORZ_mean <- apply(small[,c(8,9,10)], 1, mean)

# Teste?
obj <- t.test(x = small$CRZ_mean,
       y = small$ORZ_mean,
       alternative = c("two.sided", "less", "greater"),
       mu = 0, 
       paired = FALSE,
       var.equal = FALSE,
       conf.level = 0.95)




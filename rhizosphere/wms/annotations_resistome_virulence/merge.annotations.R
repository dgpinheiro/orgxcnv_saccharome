## Virulence & Resistome

setwd(dir = "/work/rcsilva/projects/saccharome-wgs/analises/2019-07-27_vfdb/")
fc_table <- read.delim(file = "./pvalue/contig_l2fc.tsv")

#### Virulencia ####
## Ler os dados
virulence_table <- read.delim(file = "./ardb/parsed.results.tsv")
colnames(virulence_table)[1] <- "contig"

## Junta as tabelas, somente significativos
fc_virulence <- merge(fc_table, virulence_table, by = "contig")

## Somente os unicos!
unique_fc_virulence <- subset(fc_virulence,
                              !duplicated(fc_virulence$contig))

#### Resistoma ####
resistome_table <- read.delim(file = "./megares/resistome.parsed.tsv")

## Troca o nome
colnames(resistome_table)[1] <- "contig"

## Junta com L2FC
fc_resistome <- merge(fc_table, resistome_table, by = "contig")

## Somente os unicos!
unique_fc_resistome <- subset(fc_resistome,
                              !duplicated(fc_resistome$contig))

#### Salvar as tabelas ####

## Virulência!
write.table(x = unique_fc_virulence,
            file = "./R_output/virulence.tsv",
            quote = F,
            sep = "\t",
            row.names = F)


## Resistoma
write.table(x = unique_fc_resistome,
            file = "./R_output/resistome.tsv",
            quote = F,
            sep = "\t",
            row.names = F)







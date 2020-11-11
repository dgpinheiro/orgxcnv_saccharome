## Matriz BIOM de bactérias - Kraken
library('vegan')
library('ade')
library('biom')
library('biomformat')
# library('BiodiversityR')

#### Testar as matrizes físico-químicas ####

# Carregar os dados, ver nomes de colunas
load('./dados/RData/soil.dgpinheiro.RData')
colnames(soil.df)

# Separar só as amostras de WGS #
soil.df.wgs <- soil.df[,-c(3,4,5,9,10,11)]

# Renomear as variáveis #
colnames(soil.df.wgs) <- c("Parameters",
                           "Description",
                           "CRZ2",
                           "CRZ3",
                           "CRZ4",
                           "ORZ2",
                           "ORZ3",
                           "ORZ4")

#### Processar o arquivo BIOM ####

## Importando a matriz completa ##
# kraken_biom <- biomformat::read_biom('./analises/2019-04-29_kraken_all_official/all.biom')
kraken_tsv <- read.delim(file = './analises/2019-04-29_kraken_all_official/tabletax.noplants.txt',
                         stringsAsFactors = F)

## Le o arquivo de taxonomia como taxmaps ##
saccharome_kraken <- parse_tax_data(kraken_tsv,
                                    class_cols = "taxonomy",
                                    class_sep = ";",
                                    class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                                    class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

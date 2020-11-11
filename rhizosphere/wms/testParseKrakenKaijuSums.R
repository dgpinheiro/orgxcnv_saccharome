colnames(masterlist)
masterlist.p.sum <- masterlist[,c(9,12)]
masterlist.p.sum.kaiju <- masterlist[,c(9,11)]
masterlist.kraken.sum <- ddply(masterlist.p.sum, "kraken_taxonomy",numcolwise(sum))

masterlist.kraken.sum.kaiju <- ddply(masterlist.p.sum.kaiju, 
                                     "kaiju_tax",
                                     numcolwise(sum))

masterlist.p.sum.kaiju <- subset(masterlist.p.sum.kaiju, ! is.na(log2FoldChange))



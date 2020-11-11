require(rJava)
require(xlsxjars)
library(xlsx)
library(ggfortify)
library(RColorBrewer)
library(NMF)
library(ggplot2)
library(reshape)
library(amap)
library(plotrix)
library(scales)
library(ade4)
library(ggord)
library(ggcorrplot)
library(factoextra)
library("FactoMineR")

in.dir <- "~/Dropbox/LB/Saccharome/Analises_Solos"
setwd(in.dir)
out.dir <- "/home/michelligf/Dropbox/LB/Saccharome/Analises_Solos/resultados_hipotest"

soil.df <- read.xlsx("./parametros_soil_df.xlsx",1)

# width x height (in inches)
# A4:  8.27in x 11.69in
param.width <- 8.27
param.height <- 11.69
#the default pointsize of plotted text (in big points)
param.pointsize <- 8

# "br" brazilian portuguese or "en" - english
lang="br"

########### END CONFIGURATION #####################

head(soil.df)
summary(soil.df)

param.col <- colnames(soil.df)[1]
#param.desc <- colnames(soil.df)[2]

if (lang == "br") {
  param.desc <- 'description'
  tr <- c("O"="Orgânico","C"="Convencional","OX"="Orgânico (Externo)","CX"="Convencional (Externo)")
  src.tr <- "Origem"
  bar.x <- "Parâmetros do Solo"
  bar.y <- "Valores"
  contrib.desc <- "Contribuição"
} else {
  if (lang == "en") {
    param.desc <- 'description.en'
    tr <- c("O"="Organic","C"="Conventional","OX"="Organic (External)", "CX"="Conventional (External)")
    src.tr <- "Source"
    bar.x <- "Soil parameters"
    bar.y <- "Values"
    contrib.desc <- "Contribution"
  }
}

if (! param.desc %in% colnames(soil.df)) {
  soil.df[[param.desc]] <- soil.df[[param.col]]
}

rownames(soil.df) <- soil.df[,param.col]

soil.df[[param.desc]] <- as.factor(
  sapply(as.character(soil.df[,param.desc]), 
         function(x) { max=40; 
         return(
           ifelse(nchar(x)>max,paste(strtrim(x,max-3),"...",sep=""),x) 
         )
         }
  ) 
)

head(soil.df)

samps <- setdiff( colnames(soil.df), c(param.col, 'description','description.en'))

samps

samps.o <- samps[ grep("^OBS[0-9]+", samps) ]
samps.c <- samps[ grep("^CBS[0-9]+", samps) ]

samps.ox <- samps[ grep("^OBSX",samps) ]
samps.cx <- samps[ grep("^CBSX",samps) ]

samps.o
samps.c
samps.ox
samps.cx

# x - parâmetro obrigatória para a função que será passada para
# o método apply receberá a linha neste caso
# os demais parâmetros servem para indicar os nomes das amostras 
# (orgânicas [a] e convencionais [b]), respectivamente, neste caso

# mantendo backup de soil.df e permitindo modificações em soil.df.test
soil.df.test <- soil.df

# 95% confidence intervals of the mean
soil.df.test$mean.o <- apply(soil.df.test[,samps.o], 1, mean)
soil.df.test$mean.c <- apply(soil.df.test[,samps.c],1,mean)

soil.df.test$sd.o <- apply(soil.df.test[,samps.o],1,sd)
soil.df.test$sd.c <- apply(soil.df.test[,samps.c],1,sd)
soil.df.test$se.o <- apply(soil.df.test[,samps.o],1,sem)
soil.df.test$se.c <- apply(soil.df.test[,samps.c],1,sem)

# 95% confidence intervals of the mean
soil.df.test$ci.low.o <-  apply(soil.df.test[,c("mean.o","se.o")],1, function(x) {return(x[[1]]-qnorm(0.975)*x[[2]])} )
soil.df.test$ci.high.o <- apply(soil.df.test[,c("mean.o","se.o")],1, function(x) {return(x[[1]]+qnorm(0.975)*x[[2]])} )
soil.df.test$ci.low.c <-  apply(soil.df.test[,c("mean.c","se.c")],1, function(x) {return(x[[1]]-qnorm(0.975)*x[[2]])} )
soil.df.test$ci.high.c <- apply(soil.df.test[,c("mean.c","se.c")],1, function(x) {return(x[[1]]+qnorm(0.975)*x[[2]])} )

soil.df.test$t.test <- apply(soil.df.test[,samps], 1, t.test.pvalue, samps.o, samps.c )
soil.df.test$w.test <- apply(soil.df.test[,samps], 1, w.test.pvalue, samps.o, samps.c )
soil.df.test$s.test.o <- apply(soil.df.test[,samps], 1, s.test.pvalue, samps.o )
soil.df.test$s.test.c <- apply(soil.df.test[,samps], 1, s.test.pvalue, samps.c )
soil.df.test$k.test.o <- apply(soil.df.test[,samps], 1, k.test.pvalue, samps.o )
soil.df.test$k.test.c <- apply(soil.df.test[,samps], 1, k.test.pvalue, samps.c )

#http://wiki.icmc.usp.br/images/7/73/Testenorm2013.pdf
soil.df.test$test.selected <- apply(soil.df.test[,c('s.test.o','s.test.c','k.test.o','k.test.c')], 1, test.selection )
soil.df.test[['p.value.selected']] <- apply(soil.df.test[,c('w.test','t.test','test.selected')], 1, function(x) { return(x[x['test.selected']]) } )


head(soil.df.test)

# copiando a transposta de soil.df para soil.df.manova para multivariate anova (manova)
# contendo somente as variáveis dependentes
soil.df.manova <- t(as.matrix(soil.df[,c(samps.c,samps.o)]))

# criando um array com a variável resposta (condition)
condition <- as.factor( gsub("[0-9]","",rownames(soil.df.manova)) )

# http://geog.uoregon.edu/bartlein/old_courses/geog414s05/topics/manovaex1.htm
# Testing the hypothesis that the mean values of some of 
# the soil variables differ from soil sources: O (Organic) or C (Conventional)

# Univariate analysis of variance with Iron (Ferro Fe)
summary(aov(soil.df.manova[,"Fe"] ~ condition))

# examine the hypothesis that the vector of means of the soil variables 
# is similar among soil sources using MANOVA:
#soil.manova <- summary(manova(soil.df.manova ~ condition),tol=1e-19)
#soil.manova

#summary(aov(cbind(condition) ~ Sand,data=as.data.frame(soil.df.manova)))

#summary(aov(as.formula(paste('cbind(condition)', 
#                             paste(colnames(soil.df.manova)[1:3],collapse=" * "), sep=" ~ ")), 
#            data=as.data.frame(soil.df.manova)))



# ANOVA DE CADA VARIÁVEL
summary(aov(soil.df.manova ~ condition))


# REFERÊNCIA PARA ESTUDOS POSTERIORES
#https://www.r-bloggers.com/multiple-analysis-of-variance-manova/

soil.df.test$anova <- sapply(rownames(soil.df.test),
                             function(x, manova.summary ) { 
                               return( 
                                 manova.summary[[ paste(" Response",x,sep=" ") ]][["Pr(>F)"]][1] 
                               ) }, summary(aov(soil.df.manova ~ condition)) )

# KRUSKAL WALLIS - não paramétrico
soil.df.test$kruskal <- sapply(rownames(soil.df.test),
                               function(x, c) { 
                                 return( 
                                   kruskal.test(cbind(as.data.frame(soil.df.manova)[[x]]) ~ c)$p.value
                                 ) }, condition )


write.xlsx(soil.df.test,file=paste(out.dir,"/hipotest_results.xlsx",sep=""),sheetName="Testes de Hipótese")
write.xlsx(as.data.frame(soil.manova$stats),file=paste(out.dir,"/hipotest_results.xlsx",sep=""),sheetName="MANOVA",append=TRUE)

soil.matrix.samps.t <- as.matrix(t(soil.df[,samps]))


svg(filename = paste(out.dir,"/hipotest_heatmap_",lang,".svg",sep=""), 
    width = param.width, 
    height = param.height,
    pointsize = param.pointsize
)

#png(filename=paste(out.dir,"hipotest_heatmap_",lang,".png",sep=""),
#    width = param.width,
#    height = param.height,
#    res = 300,            # 300 pixels per inch
#    pointsize = param.pointsize,
#     units="in"
#     )        

aheatmap(as.matrix(soil.df[,samps]),
         labRow = gsub('[\\{\\}\\[\\]\\^]',"", as.character(soil.df[[param.desc]]),perl=TRUE),
         color=colorRampPalette(c("blue", "white", "red"))(n = 47),
         border_color=NA,
         Colv = TRUE,
         Rowv = NA,
         #distfun=function(x) { return(dist(x, method = "euclidean")) }, # euclidean
         distfun=function(x) { return(Dist(x,method="pearson")) },
         hclustfun=function(x) { return(hclust(x, method = "ward.D")) }, 
         annColors=setNames(list(setNames(c('green2','orange'),tr[c('O','C')])),
                            src.tr
         ), 
         main="", 
         annRow=NA, 
         annCol=setNames(data.frame(matrix(as.character(tr[gsub("BS[0-9X]+","",samps)])),row.names=samps),src.tr), 
         cellwidth=20, 
         cellheight=10,
         breaks=NMF:::generate_breaks(c(-5,5) , 49, center=0 ),
         scale="row"
)

dev.off()



#png(filename=paste(out.dir,"hipotest_PCA_",lang,".png",sep=""),
#     width = param.width,
#     height = param.height,
#     res = 300,            # 300 pixels per inch
#     pointsize = param.pointsize,
#     units="in"
#    )        

sel <- c("pH","P","S","Ca","Mg","K","Al","B","Cu","Fe","Mn","Zn","SOM","Clay","Silt","Sand")
#pca.x <- data.frame(t(soil.df[sel,samps]))
#pca.x.scaled <- as.data.frame(scale(pca.x[,sel],scale=TRUE,center=TRUE))
#pca.res <- prcomp(pca.x.scaled, scale.=FALSE, center=FALSE)
#sum(((pca.res$sdev/sum(pca.res$sdev))*100)[1:11])
#sel <- unique(rownames(which(abs(pca.res$rotation)>=max(abs(pca.res$rotation[,1:11]))-(0.3*max(abs(pca.res$rotation[,1:6]))), arr.ind = TRUE)))
#sel <- c("B","S","pH","Silt")

pca.x <- data.frame(t(soil.df[sel,samps]))
#pca.sel<-as.data.frame(t(pca.x))
#pca.sel$var <- apply(pca.sel,1,var)
#sel<-rownames(pca.sel[order(pca.sel$var,decreasing=TRUE),])[1:8]
pca.x.scaled <- as.data.frame(scale(pca.x[,sel],scale=TRUE,center=TRUE))

pca.res <- prcomp(pca.x.scaled, scale.=FALSE, center=FALSE)


# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/
pca.f <- PCA(pca.x, scale.unit = TRUE, ncp = 5, graph = TRUE)
fviz_eig(pca.f, addlabels = TRUE, ylim = c(0, 50))
fviz_pca_var(pca.f, col.var = "black")
corrplot(pca.f$var$cos2, is.corr=FALSE)
fviz_cos2(pca.f, choice = "var", axes = 1:2)

svg(filename = paste(out.dir,"/hipotest_PCA_corrplot_",lang,".svg",sep=""), 
    width = param.width, 
    height = param.height,
    pointsize = param.pointsize
)
corrplot.x <- corrplot(pca.f$var$cos2, is.corr=FALSE)
dev.off()



svg(filename = paste(out.dir,"/hipotest_PCA_contrib_",lang,".svg",sep=""), 
    width = param.width, 
    height = param.height,
    pointsize = param.pointsize
)

gg.cp <- ggcorrplot(pca.f$var$cos2, method = "square", outline.col = "white", type="full", pch.cex = 1)
gg.cp
# Contributions of variables to PC1
fviz.cp.1 <- fviz_contrib(pca.f, choice = "var", axes = 1, top = 10)+ggtitle("Contribution",subtitle="Dim1")
# Contributions of variables to PC2
fviz.cp.2 <- fviz_contrib(pca.f, choice = "var", axes = 2, top = 10)+ggtitle("Contribution",subtitle="Dim2")
fviz.cp.12 <- fviz_contrib(pca.f, choice = "var", axes = 1:2, top = 10)+ggtitle("Contribution",subtitle="Dim1-2")
fviz.cp.12


fviz.biplot <- fviz_pca_biplot(pca.f,geom.ind = c("point","text"),pointshape = 21, pointsize = "cos2",
                               habillage = as.factor(as.character(tr[gsub("BS[0-9]*","",samps)])),
                               fill.ind = as.factor(as.character(tr[gsub("BS[0-9]*","",samps)])),
                               addEllipses = TRUE,
                               ellipse.type="t",
                               ellipse.level=0.95,
                               label="all",
                               col.var="black",
                               fill.var="white",
                               #gradient.cols =colorRampPalette(brewer.pal(9,"Greys"))(3),
                               alpha.var ="contrib",
                               col.circle="black",
                               repel=TRUE,
                               legend.title = list(fill = src.tr, color = src.tr,
                                                   alpha = contrib.desc, size= contrib.desc)
)+coord_fixed(ratio = 1)+
  scale_color_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
                                                                     tr[c("O","C","OX","CX")]))+
  scale_fill_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
                                                                    tr[c("O","C","OX","CX")]))


svg(filename = paste(out.dir,"/hipotest_PCA_",lang,".svg",sep=""), 
    width = param.width, 
    height = param.width,
    pointsize = param.pointsize
)

ggarrange(
  ggarrange(fviz.biplot, ggarrange( fviz.cp.1, fviz.cp.2 ,nrow=2, labels=c("(c)","(d)")) ,ncol=2, widths = c(.7, .3), labels=c("(a)","")),
  ggarrange(gg.cp, fviz.cp.12, ncol=2, widths=c(.7,.3), labels=c("(b)","(e)")),
  nrow=2,
  heights=c(.65,.25),
  widths=1
)

dev.off()






#require(vegan)
#plot(pca(pca.x.scaled))
#biplot(pca.res, pc.biplot=FALSE,expand=3)
#pca.x.scaled[,c("Fe","P","Mn","S")]

#pcax.limit <- max(c(abs(min(pca.res$x[,c("PC1")])),abs(max(pca.res$x[,c("PC1")]))))
#pcay.limit <- max(c(abs(min(pca.res$x[,c("PC2")])),abs(max(pca.res$x[,c("PC2")]))))

#x11()
#plot(pca.x.scaled[,c("Fe","Zn")])
#text(pca.x.scaled[,c("Fe","Zn")], row.names(pca.x), cex=0.6, pos=4, col="red")

#ggord(dudi.pca(pca.x,scannf=FALSE,center = TRUE, scale = FALSE ), as.factor(tr[gsub('BS[0-9]*','',rownames(pca.x))]),
#      txt=2,obslab = TRUE,veccol = 'darkgray' )+
#  scale_color_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
#                                                                     tr[c("O","C","OX","CX")]))+
#  scale_fill_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
#                                                                          tr[c("O","C","OX","CX")]))+
#  theme_classic()+coord_fixed(ratio = 1)

# autoplot(pca.res, 
#          scale=0,
#          data=setNames(data.frame(as.matrix(pca.x), 
#                                   as.factor(as.character(tr[gsub("BS[0-9]*","",samps)]))),
#                        c(colnames(pca.x),src.tr)
#            ),
#          frame=TRUE,
#          colour=src.tr,
#          shape=FALSE,
#          label.size=4,
#          loadings=TRUE,
#          loadings.colour="lightgray",
#          loadings.label.size=3,
#          loadings.label.colour="black",
#          loadings.label=TRUE,
#          variance_percentage = TRUE
# )+
#   scale_color_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
#                                                                      tr[c("O","C","OX","CX")]))+
#   scale_fill_manual(labels=tr[c("O","C","OX","CX")],values=setNames(c('green2','orange','greenyellow','orange4'),
#                                                                      tr[c("O","C","OX","CX")]))+
#   theme_classic()+
#   coord_fixed(ratio = 1)+
#   geom_point(data=rbind(pca.res$x[grep("^OBS",rownames(pca.res$x)),c("PC1","PC2")]),
#                   aes(x=PC1,y=PC2),colour = "darkgreen")+
#   geom_point(data=rbind(pca.res$x[grep("^CBS",rownames(pca.res$x)),c("PC1","PC2")]),
#              aes(x=PC1,y=PC2),colour = "darkorange")

#scale_color_manual(labels=tr[c("O","C")],values=setNames(c('green2','orange'),tr[c("O","C")]))
#scale_x_continuous(limit=c(-0.6,1))+
#scale_y_continuous(limit=c(-0.6,1))+



soil.df.test.t <- as.data.frame(t(soil.df.test[c("mean.o","mean.c")]))
soil.df.test.t$Id <- gsub("[^.]+.","",rownames(soil.df.test.t))
soil.df.test.t$Source <- as.factor(tr[toupper(gsub("[^.]+.","",rownames(soil.df.test.t)))])

soil.df.test.m <- melt(soil.df.test.t,id=c("Id","Source"))

soil.df.test.m$ci.low <- apply(soil.df.test.m, 1, 
                               function(x,df) {
                                 return( df[x[['variable']] , paste("ci.low",x[['Id']],sep=".") ])   
                               }, 
                               soil.df.test )
soil.df.test.m$ci.high <- apply(soil.df.test.m, 1, 
                                function(x,df) {
                                  return( df[x[['variable']] , paste("ci.high",x[['Id']],sep=".") ])   
                                }, 
                                soil.df.test )
soil.df.test.m[,param.desc] <- soil.df[soil.df.test.m$variable,param.desc]

#T-test letters with p<=0.05 
p.value.cutoff <- 0.05


for (param in levels(soil.df.test.m$variable)) {
  print(param)
  
  iop <- which(soil.df.test.m$Id=="o" & soil.df.test.m$variable==param)
  icp <- which(soil.df.test.m$Id=="c" & soil.df.test.m$variable==param)
  
  if (soil.df.test[param,'p.value.selected'] <= p.value.cutoff) {
    soil.df.test.m[iop,"Letter"] <- "a"
    soil.df.test.m[icp,"Letter"] <- "b"
  } else {
    soil.df.test.m[iop,"Letter"] <- "a"
    soil.df.test.m[icp,"Letter"] <- "a"
  }
}

if ((length(samps.cx)==1)&(length(samps.ox)==1)) {
  for (param in levels(soil.df.test.m$variable)) {
    soil.df.test.m[which(soil.df.test.m$Id=="o" & soil.df.test.m$variable==param),"X"] <- soil.df[param,samps.ox]
    soil.df.test.m[which(soil.df.test.m$Id=="c" & soil.df.test.m$variable==param),"X"] <- soil.df[param,samps.cx]
  }
}



svg(filename = paste(out.dir,"/hipotest_bars_",lang,".svg",sep=""), 
    width = param.width, 
    height = param.height,
    pointsize = param.pointsize
)

# png(filename=paste(out.dir,"hipotest_bars_",lang,".png",sep=""),
#     width = param.width,
#     height = param.height,
#     res = 300,            # 300 pixels per inch
#     pointsize = param.pointsize,
#     units="in"
#     )

mylog2_trans <- function (base = 2) 
{
  trans <- function(x) log(x + 1, base)
  inv <- function(x) base^x
  trans_new(paste0("log-", format(base)), trans, inv, log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

soil.df.test.m$desc <- factor(gsub("([^ %])%","\\1 %", as.character(soil.df.test.m[[param.desc]]),perl=TRUE))
soil.df.test.m$desc <- factor(gsub(" ","~",  as.character(soil.df.test.m$desc)))
soil.df.test.m$desc <- factor(gsub("%","paste('%')", as.character(soil.df.test.m$desc)))
soil.df.test.m$desc <- factor(gsub("\\.",".", as.character(soil.df.test.m$desc)))

# Default bar plot
p<- ggplot(soil.df.test.m, aes(x=desc, y=value, fill=Source)) + 
  scale_y_continuous(trans=mylog2_trans(),breaks=trans_breaks('log2', function(x) 2^(x+1)))+
  coord_flip() +
  geom_bar(stat="identity", color="black", position=position_dodge(), width=.7,size=.5) +
  geom_errorbar(aes(ymin=ci.low, ymax=ci.high), width=.2, size=.5,
                position=position_dodge(.7))+
  labs(title="", y = bar.y, x= bar.x)+
  theme_classic()+
  theme(
    axis.text = element_text(angle=0, vjust=0.5,size=12,face="bold"),
    legend.title=element_blank(),legend.position="bottom",
    legend.text=element_text(face="italic")
  )+
  scale_x_discrete(limits=unique(as.character(soil.df.test.m$desc)),labels=function(l) parse(text=sprintf('%s',l)) )+
  geom_text(aes(label=Letter,colour=Source),position=position_dodge(1),vjust=.9,hjust=-1,size=4)+
  scale_fill_manual(labels=tr[c("O","C")],values=setNames(c('green2','orange'),tr[c("O","C")]))+
  scale_colour_manual(labels=tr[c("O","C")],values=setNames(c('green2','orange'),tr[c("O","C")]))

if ((length(samps.cx)==1)&(length(samps.ox)==1)) {
  p <- p+geom_point(aes(x=desc,y=X,fill=Source),stat="identity", color="red", position=position_dodge(.7), size=.5)
}

print(p)

dev.off()

if ( sum(c("Sand","Silt","Clay") %in% rownames(soil.df))==3 ) {
  
  print("Performing soil texture analysis")
  
  my.soil.texture <- function (soiltexture = NULL, main = "", at = seq(0.1, 0.9, by = 0.1), 
                               axis.labels = c("percent sand", "percent silt", "percent clay"), 
                               tick.labels = list(l = seq(10, 90, by = 10), r = seq(10, 
                                                                                    90, by = 10), b = seq(10, 90, by = 10)), show.names = TRUE, 
                               show.lines = TRUE, col.names = "gray", bg.names = par("bg"), 
                               show.grid = FALSE, col.axis = "black", col.lines = "gray", 
                               col.grid = "gray", lty.grid = 3, show.legend = FALSE, label.points = FALSE, 
                               point.labels = NULL, col.symbols = "black", pch = par("pch"), 
                               ...) 
  {
    par(xpd = TRUE)
    plot(0.5, type = "n", axes = FALSE, xlim = c(0, 1), ylim = c(0, 
                                                                 1), main = NA, xlab = NA, ylab = NA)
    triax.plot(x = NULL, main = main, at = at, axis.labels = axis.labels, cex.axis=2,
               tick.labels = tick.labels, col.axis = col.axis, show.grid = show.grid, 
               col.grid = col.grid, lty.grid = lty.grid)
    arrows(0.12, 0.41, 0.22, 0.57, length = 0.15)
    arrows(0.78, 0.57, 0.88, 0.41, length = 0.15)
    arrows(0.6, -0.1, 0.38, -0.1, length = 0.15)
    if (show.lines) {
      triax.segments <- function(h1, h3, t1, t3, col) {
        segments(1 - h1 - h3/2, h3 * sin(pi/3), 1 - t1 - 
                   t3/2, t3 * sin(pi/3), col = col)
      }
      h1 <- c(85, 70, 80, 52, 52, 50, 20, 8, 52, 45, 45, 65, 
              45, 20, 20)/100
      h3 <- c(0, 0, 20, 20, 7, 0, 0, 12, 20, 27, 27, 35, 40, 
              27, 40)/100
      t1 <- c(90, 85, 52, 52, 43, 23, 8, 0, 45, 0, 45, 45, 
              0, 20, 0)/100
      t3 <- c(10, 15, 20, 7, 7, 27, 12, 12, 27, 27, 55, 35, 
              40, 40, 60)/100
      triax.segments(h1, h3, t1, t3, col.lines)
    }
    if (show.names) {
      xpos <- c(0.5, 0.7, 0.7, 0.73, 0.73, 0.5, 0.275, 0.275, 
                0.27, 0.27, 0.25, 0.135, 0.18, 0.055, 0.49, 0.72, 
                0.9)
      ypos <- c(0.66, 0.49, 0.44, 0.36, 0.32, 0.35, 0.43, 0.39, 
                0.3, 0.26, 0.13, 0.072, 0.032, 0.024, 0.18, 0.15, 
                0.06) * sin(pi/3)
      snames <- c("Argilosa", "Argilo-", "siltosa", "Franco-argilo-", "siltosa", 
                  "Franco-argilosa", "Argilo-", "arenosa", "Franco-argilo-", "arenosa", 
                  "Franco-arenosa", "Areia", "franca", "Arenosa", "Franca", "Franco-siltosa", 
                  "Siltosa")
      boxed.labels(xpos, ypos, snames, border = FALSE, xpad = 0.5, cex=1.6)
    }
    par(xpd = FALSE)
    if (is.null(soiltexture)) 
      return(NULL)
    soilpoints <- triax.points(soiltexture, show.legend = show.legend, 
                               label.points = label.points, point.labels = point.labels, 
                               col.symbols = col.symbols, pch = pch, ...)
    invisible(soilpoints)
  }
  
  
  
  svg(filename = paste(out.dir,"/hipotest_soils_",lang,".svg",sep=""), 
      width = param.width, 
      height = param.height,
      pointsize = param.pointsize
  )
  
  # png(filename=paste(out.dir,"hipotest_soils_",lang,".png",sep=""),
  #     width = param.width,
  #     height = param.height,
  #     res = 300,            # 300 pixels per inch
  #     pointsize = param.pointsize,
  #     units="in"
  #     )
  
  
  # simple soil texture
  if (lang == "br") {
    # Em português
    my.soil.texture(show.lines=T, show.names=T, col.lines='black', col.names='black', main='Classificação Textural de Solos (USDA)',
                    axis.labels=c("Areia (%)", "Silte (%)", "Argila (%)"))
  } else {
    if (lang == "en") {
      soil.texture(show.lines=T, show.names=T, col.lines='black', col.names='black', main='USDA Textural Soil Classification',
                   axis.labels=c("Sand (%)", "Silt (%)", "Clay (%)"))
    }
  }
  
  
  # ORGÂNICO
  x <- as.data.frame(t(soil.df[c("Sand","Silt","Clay"),samps.o]))
  triax.points(x, cex=2, pch=16, col.symbol=rgb(98,189, 24, alpha=100, max=255))
  
  triax.points(cbind(median(x$Sand), median(x$Silt), median(x$Clay)),
               col.symbols='green2', pch=16, cex=2)
  
  # CONVENCIONAL
  x <- as.data.frame(t(soil.df[c("Sand","Silt","Clay"),samps.c]))
  triax.points(x, cex=2, pch=16, col.symbol=rgb(255,159, 0, alpha=100, max=255))
  
  triax.points(cbind(median(x$Sand), median(x$Silt), median(x$Clay)),
               col.symbols='orange', pch=16, cex=2)
  
  
  if ((length(samps.cx)==1)&(length(samps.ox)==1)) {
    triax.points(cbind(soil.df["Sand",samps.ox], soil.df["Silt",samps.ox], soil.df["Clay",samps.ox]),
                 col.symbols='palegreen4', pch=16, cex=2)
    
    triax.points(cbind(soil.df["Sand",samps.cx], soil.df["Silt",samps.cx], soil.df["Clay",samps.cx]),
                 col.symbols='orange4', pch=16, cex=2)
    
    legend("topright", c(tr[c("O","C","OX","CX")]), cex=2.5,fill=c("green2","orange","palegreen4","orange4"),
           border=NA,box.col=NA)
  } else {
    legend("topright", c(tr[c("O","C")]), cex=2.5,fill=c("green2","orange"),
           border=NA,box.col=NA)
  }
  dev.off()
}


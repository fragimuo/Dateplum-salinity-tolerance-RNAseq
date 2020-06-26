#' ---
#' title: "Annotation analyses"
#' author: "Nicolas Delhomme & Francisco Gil-Muñoz"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Libraries
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(wordcloud))


#' Helper function
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' set the wd
setwd("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/persimmon")
#' ```

#' define a palette
pal <- brewer.pal(8,"Dark2")

#' open the db connection
con <- dbConnect(dbDriver("SQLite"),dbname="/mnt/picea/storage/reference/Taxonomy/20181107/taxonomy.sqlite")

#' read the metadata
annot <- read.delim("annotation/B2G/Persimon_20181113_blast2go_partial-export.txt",
                    stringsAsFactors = FALSE)

#' SOME RESULTS CAN BE LOADED FROM 20190319.rda
load("20190319.rda")

#' # Annotation 
#' ## Taxonomy
#' subset to only those that have description
desc.sel <- grepl("TaxID",annot$Sequence.Description)
#' and those that have a blast hit
blast.sel <- grepl("gi\\|",annot$Blast.Top.Hit.Description..HSP.)

# keep only the desc.sel if both selectors match
blast.sel[blast.sel & desc.sel] <- FALSE

#' ### TaxID based
#' collect division using the description
taxID <- as.integer(gsub(".*TaxID=| RepID=.*", "", annot$Sequence.Description[desc.sel]))

#' query the database
taxid_divnames <- dbGetQuery(con, paste0("SELECT d.div_nam, n.tax_id FROM division AS d ",
                                         "LEFT JOIN node AS n ON d.div_id = n.div_id ",
                                         "WHERE n.tax_id IN (",
                                         paste(unique(taxID), collapse = ","),
                                         ")"))

#' update the annot
annot$Division[desc.sel] <- taxid_divnames[match(taxID,taxid_divnames$tax_id),1]


#' ### GI based
gi <- sapply(strsplit(annot$Blast.Top.Hit.Description..HSP.[blast.sel],"\\|"),"[",2)

#' query the database
gi_divnames <- dbGetQuery(con, paste0("SELECT DISTINCT(g.gid),d.div_nam FROM division AS d ",
                                      "LEFT JOIN node AS n ON d.div_id = n.div_id ",
                                      "LEFT JOIN gi_taxid AS g ON n.tax_id = g.tid ",
                                      "WHERE g.gid IN (",
                                      paste(unique(gi), collapse = ","),
                                      ")"))

annot$Division[blast.sel] <- gi_divnames[match(gi,gi_divnames$gid),2]

annot$Division[is.na(annot$Division)] <- "NA"

head(annot[blast.sel,])

#' disconnect
dbDisconnect(con)

#' ## GC content

#' read the sequences
fa <- readDNAStringSet("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/trinity/Trinity.format.fasta")

gc <- rowSums(alphabetFrequency(fa,as.prob = TRUE)[,c("C","G")])

names(gc) <- sub(" .*","",names(fa))

annot$GC <- gc[annot$Sequence.Name]

gcList <- split(annot$GC,annot$Division)

plot.multidensity(gcList,col=c(1,pal),legend.x="topright",
                  xlim = c(0,1),lwd=2,legend.cex = .8,xlab = "GC proportion",
                  legend.lwd = 2)

rm(fa)

#' ## Expression values
samples <- read.delim("trinity/samples.txt",header=FALSE)
mat <- read.delim("trinity/Trinity_trans.isoform.counts.matrix",row.names = 1)
rn <- rownames(mat)
mat <- apply(mat,2,as.integer)
rownames(mat) <- rn
cd <- DataFrame(Tissue=factor(substr(samples$V1,start = 1,stop = 1)),
                        Condition=factor(substr(samples$V1,start = 2,stop = 2)))

dds <- DESeqDataSetFromMatrix(mat,
                              cd,
                              ~Condition*Tissue)

#' variance stabilising transformation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
rownames(vst) <- rn

annot <- cbind(annot,mat[annot$Sequence.Name,])
annot$meanExp <- rowMeans2(vst[annot$Sequence.Name,])

dev.null <- sapply(unique(annot$Division),function(div){
  sel <- annot$Division == div
  heatscatter(annot[sel,"GC"],
    annot[sel,"meanExp"],xlim=c(0,1),
    xlab="GC proportion",
    ylim=range(annot$meanExp),
    ylab="expression values (VST)",
    main=div,sub=paste0("n=",sum(sel)))
})

#' Save annotation file

write.table(annot, file="annotation.txt")

#' ## Variance Stabilising Transformation
#' 
#' Since there is almost no difference in library size, a VST is
#' perfectly applicable.

#' This is the situation prior to normalisation
meanSdPlot(log2(counts(dds)[rowSums(counts(dds)) > 0,]))

#' Normalisation - blind, we give no prior as we want to assess quality
PlPCA <- varianceStabilizingTransformation(dds ,blind=TRUE)
#' Extract the normalised counts
vsd <- assay(PlPCA)
#' Look at the VST fit. It looks good, almost flat, around 1 sd on average
meanSdPlot(vsd[rowSums(counts(ddsV)) > 0,])
#' The VST introduces an offset
range(vsd)
#' Which we remove, so that 0 means no expression
vsd <- vsd - min(vsd)

#' # Quality Assessment
#' ## Principal Component Analysis
pc <- prcomp(t(vsd))
percent <- round(summary(pc)$importance[2,]*100)

pal <- brewer.pal(8,"Dark2")

#Extract the names of the samples
samples <- colnames(vsd)
samples <- gsub('[0-9]+', '', samples)


plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     pch=19,
     col=pal[as.factor(samples)])
    legend("top",pch=19,
       col=pal[1:nlevels(as.factor(samples))],
       legend=levels(as.factor(samples)))

#' Try ML or DL for identifying putative plant new genes



#' 
#' Write the IDs / expression of the Viridiplatae only
head(mat[annot[annot$Division=="Viridiplantae",1],])

rownames(dds) <- rn
ddsV <- dds[annot[annot$Division=="Viridiplantae",1],]
ddsV

ddsV <- DESeq(ddsV)
resultsNames(ddsV)

resSR <- results(ddsV,name = "ConditionS.TissueR")
sum(resSR$padj <= 0.01 & !is.na(resSR$padj))

resRT <- results(ddsV,name = "ConditionT.TissueR")
sum(resRT$padj <= 0.01 & !is.na(resRT$padj))

resTS <- results(ddsV,contrast = c("Condition","T","S"))
sum(resTS$padj <= 0.01 & !is.na(resTS$padj))

resTSdf <- rownames(resTS)
resTSdf <- cbind(resTSdf, resTS$log2FoldChange)
resTSdf <- cbind(resTSdf, resTS$pvalue)
resTSdf <- cbind(resTSdf, resTS$padj)
colnames(resTSdf) <- c("Sequence.Name", "log2FoldChange", "pvalue", "padj")

annot2 <- merge(annot, resTSdf, by="Sequence.Name")




grid.newpage()

grid.draw(venn.diagram(list(rownames(resSR[which(resSR$padj <= 0.01),]),
                            rownames(resRT[which(resRT$padj <= 0.01),]),
                            rownames(resTS[which(resTS$padj <= 0.01),])),
                       filename=NULL, category.names = c("SR","TR","TS"),
                       fill=pal[1:3], resolution=100))

grid.newpage()

grid.draw(venn.diagram(list(rownames(resSR[which(resSR$padj <= 0.01),]),
                            rownames(resRT[which(resRT$padj <= 0.01),])),
                       filename=NULL, category.names = c("SR","RT"),
                       fill=pal[1:2], resolution=500))

write.table(res[which(res$padj <= 0.01),], file="TvS.txt")

sum(res$padj <= 0.01 & !is.na(res$padj))

res <- results(ddsV,contrast = list("Condition_T_vs_C","Condition_S_vs_C"))

res

res <- results(ddsV,contrast = c("Condition","T","S"))

write.table(res[which(res$padj <= 0.01),], file="TvS.txt")

res

TvS <- merge(annot, as(res,"data.frame"), by.x="Sequence.Name", by.y="row.names")

write.table(TvS, file="TvS_annot.txt")

write.table(TvS[which(TvS$padj <= 0.01),], file="TvS_annotSig.txt")


dds <- DESeq(dds)

resT <- results(dds,contrast = c("Condition","T","S"))

TvST <- merge(annot, as(resT,"data.frame"), by.x="Sequence.Name", by.y="row.names")

write.table(TvST, file="TvST_annot.txt")

write.table(TvST[which(TvST$padj <= 0.01),], file="TvST_annotSig.txt")



#' # Validation of contrasts
#'
#setdiff(resSR[which(resSR$padj <= 0.01),], union(resSR[which(resSR$padj <= 0.01),],resSR[which(resSR$padj <= 0.01),]))

#TolROOT <- setdiff(rownames(resRT[which(resRT$padj <= 0.01),]), 
#                    union(rownames(resTS[which(resTS$padj <= 0.01),]),
#                          rownames(resSR[which(resSR$padj <= 0.01),])))

# Significant genes for tolerant roots

#s.vst <- t(scale(t(vsd[as.character(TolROOT), grep("R",colnames(vsd))])))
#sum(is.na(s.vst))
#which(rowSums(is.na(s.vst))>0)
#s.vst <- s.vst[-466,]


hpal <- colorRampPalette(c("blue","white","red"))(100)
pear.D <- function(x) dist(x, method=pearson.dist)
warD <- function(x) hclust(x, method=ward.D)

B2GO_annot <- read.csv("/mnt/picea/projects/dateplum/rgarcia/persimmon-rootstocks-salinity-tolerance/annotation/B2GO_annot.txt", row.names=1)
rownames(B2GO_annot) <- B2GO_annot$Row.names

B2GO_annot <- merge(B2GO_annot, resTSdf, by="Sequence.Name")

#heatResRT <- heatmap.2(s.vst,
#                     distfun=pearson.dist,
#                     hclustfun=function(X){hclust(X,method="ward.D2")},
#                     trace="none",
#                     col=hpal,
#                     labCol = colnames(s.vst),
#                     labRow = NA)
                     # key=FALSE)

#heatResRT

#' Cut the tree into branches and look at the genes inside each branch
#' 
#cutree(as.hclust(heatResRT$rowDendrogram),5)

#' Significant genes for tolerant roots 

#sa.vst <- t(scale(t(vsd[as.character(TolROOT), colnames(vsd)])))
#sum(is.na(s.vst))

#heatResRTA <- heatmap.2(sa.vst,
#                     distfun=pearson.dist,
#                     hclustfun=function(X){hclust(X,method="ward.D2")},
#                     trace="none",
#                     col=hpal,
#                     labCol = colnames(sa.vst),
#                     labRow = NULL)

#' Significant genes for tolerant roots

TolROOT <- setdiff(rownames(resRT[which(resRT$padj <= 0.01),]), 
                    union(rownames(resTS[which(resTS$padj <= 0.01),]),
                          rownames(resSR[which(resSR$padj <= 0.01),])))


tr.vst <- t(scale(t(vsd[as.character(TolROOT), ])))
sum(is.na(tr.vst))

heatTolR <- heatmap.2(tr.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(tr.vst),
                      labRow = NA)

cltr <- cutree(as.hclust(heatTolR$rowDendrogram),20)


TolROOT <- cbind(TolROOT,clusterID=cltr)

colnames(TolROOT)[1] <- colnames(B2GO_annot)[1]


TolROOT <- merge(TolROOT, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(TolROOT, file = "TolROOT.csv")
write.csv(TolROOT[which(TolROOT$padj.x<=0.05),], file = "TolROOTpadj.csv")

#' Significant genes for sensitive roots 
 
 SensROOT <- setdiff(rownames(resSR[which(resSR$padj <= 0.01),]), 
                   union(rownames(resTS[which(resTS$padj <= 0.01),]),
                         rownames(resRT[which(resRT$padj <= 0.01),])))


#ss.vst <- t(scale(t(vsd[as.character(SensROOT), grep("R",colnames(vsd))])))
#sum(is.na(s.vst))

#heatResRS <- heatmap.2(ss.vst,
#                     distfun=pearson.dist,
#                     hclustfun=function(X){hclust(X,method="ward.D2")},
#                     trace="none",
#                     col=hpal,
#                     labCol = colnames(ss.vst),
#                     labRow = NULL)

sr.vst <- t(scale(t(vsd[as.character(SensROOT), colnames(vsd)])))
sum(is.na(sr.vst))

heatResRS <- heatmap.2(sr.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(sr.vst),
                      labRow = NA)

clrs <- cutree(as.hclust(heatResRS$rowDendrogram),20)


SensROOT <- cbind(SensROOT,clusterID=clrs)

colnames(SensROOT)[1] <- colnames(B2GO_annot)[1]


SensROOT <- merge(SensROOT, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(SensROOT, file = "SensROOT.csv")
write.csv(SensROOT[which(SensROOT$padj.x<=0.05),], file = "SensROOTpadj.csv")

#' Significant genes for T vs S

TolvSens <- rownames(resTS[which(resTS$padj <= 0.01),]) 

sts.vst <- t(scale(t(vsd[as.character(TolvSens), ])))
sum(is.na(s.vst))

heatResTS <- heatmap.2(sts.vst,
                        distfun=pearson.dist,
                        hclustfun=function(X){hclust(X,method="ward.D2")},
                        trace="none",
                        col=hpal,
                        labCol = colnames(sts.vst),
                        labRow = NA)

clts <- cutree(as.hclust(heatResTS$rowDendrogram),20)


TolvSens <- cbind(TolvSens,clusterID=clts)

colnames(TolvSens)[1] <- colnames(B2GO_annot)[1]


TolvSens <- merge(TolvSens, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(TolvSens, file = "TolvSens.csv")

grep("Auxin", TolvSens$Sequence.Description.x, ignore.case = TRUE)
?grep
#' Significant genes for T vs S only in leaves

TolvSensL <- setdiff(rownames(resTS[which(resTS$padj <= 0.01),]), 
                    union(rownames(resSR[which(resSR$padj <= 0.01),]),
                          rownames(resRT[which(resRT$padj <= 0.01),])))

stsl.vst <- t(scale(t(vsd[as.character(TolvSensL), ])))
sum(is.na(stsl.vst))

heatResTSL <- heatmap.2(stsl.vst,
                       distfun=pearson.dist,
                       hclustfun=function(X){hclust(X,method="ward.D2")},
                       trace="none",
                       col=hpal,
                       labCol = colnames(stsl.vst),
                       labRow = NA)

cltsl <- cutree(as.hclust(heatResTSL$rowDendrogram),20)


TolvSensL <- cbind(TolvSensL,clusterID=cltsl)

colnames(TolvSensL)[1] <- colnames(B2GO_annot)[1]


TolvSensL <- merge(TolvSensL, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(TolvSensL, file = "TolvSensL.csv")
write.csv(TolvSensL[which(TolvSensL$padj.x<=0.05),], file = "TolvSensLpadj.csv")


#' Salt induced only in roots

Rind <- setdiff(intersect(rownames(resSR[which(resSR$padj <= 0.01),]),
                          rownames(resRT[which(resRT$padj <= 0.01),])),
                rownames(resTS[which(resTS$padj <= 0.01),]))

sri.vst <- t(scale(t(vsd[as.character(Rind), ])))
sum(is.na(sri.vst))

heatRi <- heatmap.2(sri.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(sri.vst),
                      labRow = NA)

clir <- cutree(as.hclust(heatRi$rowDendrogram),20)


Rind <- cbind(Rind,clusterID=clir)

colnames(Rind)[1] <- colnames(B2GO_annot)[1]


Rind <- merge(Rind, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(Rind, file = "Rind.csv")
write.csv(Rind[which(Rind$padj.x<=0.05),], file = "Rindpadj.csv")

#' Tolerants > effect

TolE <- setdiff(intersect(rownames(resTS[which(resTS$padj <= 0.01),]),
                          rownames(resRT[which(resRT$padj <= 0.01),])),
                rownames(resSR[which(resSR$padj <= 0.01),]))

ste.vst <- t(scale(t(vsd[as.character(TolE), ])))
sum(is.na(ste.vst))

heatTolE <- heatmap.2(ste.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(ste.vst),
                      labRow = NA)

cltole <- cutree(as.hclust(heatTolE$rowDendrogram),20)


TolE <- cbind(TolE,clusterID=cltole)

colnames(TolE)[1] <- colnames(B2GO_annot)[1]


TolE <- merge(TolE, B2GO_annot, by="row.names", sort=FALSE)
write.csv(TolE, file = "TolE.csv")


#' Sensitive > effect

SensE <- setdiff(intersect(rownames(resTS[which(resTS$padj <= 0.01),]),
                          rownames(resSR[which(resSR$padj <= 0.01),])),
                rownames(resRT[which(resRT$padj <= 0.01),]))

sse.vst <- t(scale(t(vsd[as.character(SensE), ])))
sum(is.na(sse.vst))

heatSensE <- heatmap.2(sse.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(sse.vst),
                      labRow = NA)

clsense <- cutree(as.hclust(heatSensE$rowDendrogram),20)


SensE <- cbind(SensE,clusterID=clsense)

colnames(SensE)[1] <- colnames(annot)[1]


SensE <- merge(SensE, B2GO_annot, by="row.names", sort=FALSE)
write.csv(SensE, file = "SensE.csv")

source("~/Git/UPSCb/src/R/gopher.R")
suppressPackageStartupMessages(source("~delhomme/Git/UPSCb/src/R/gopher.R"))

# dat = change for your list of genes of interest (a subset of a population)
# dat <- scan("~/Git/UPSCb/projects/facility/doc/tiggy-gene-example.txt",what="character")
# sub("^EC:","",TolvSens$Enzyme.Code)
# sub("|EC:"," ",TolvSens$Enzyme.Code)

NamTolL <- TolLeave$Row.names
NamSensL <- SensLeave$Row.names
NamTvSensl <- TolvS$Row.names
NamTvSR <- TolvSensR$Row.names
NamLind <- Lind$Row.names
NamTolEl <- TolEl$Row.names
NamSensEl <- SensEl$Row.names
NamSL <- rownames(resSL[which(resSL$padj <= 0.01),])
NamTL <- rownames(resTL[which(resTL$padj <= 0.01),])
NamTvS <- TolvSens$Row.names

#Select all KEGGS                 
#TestEnr <- B2GO_annot[!(B2GO_annot$Enzyme.Code == "" | is.na(B2GO_annot$Enzyme.Code)), ]
#KEGGstats <- table(B2GO_annot$Enzyme.Code)
#TestEnr <- B2GO_annot[B2GO_annot$Enzyme.Code == "[EC:2.4.2.31]", ]
#TestEnr <- B2GO_annot[grep("EC:2.4.", B2GO_annot$Enzyme.Code), ]

# bg = your population (think what defines it.)
# bg <- scan("~/Git/UPSCb/projects/facility/doc/Tiggy.spruce.background.txt",what="character",skip=1)

#bg <- annot[annot$Division=="Viridiplantae",1]
#TestE <- merge(TestEnr, annot, by="Sequence.Name", sort=FALSE)



#TestEv <- TestE[which(TestE$Division.x=="Viridiplantae"),1]

#length(TestEv)

#TestEv2 <- I(intersect(TestEv, bg))
# bg$Enzyme.code <- sub("^EC:","",TolvSens$Enzyme.Code)
# we just quantify the run time
# enrichment <- gopher(dat,task = list("go","kegg","pfam"),background = bg,url="pabies")
ignore.case=TRUE

system.time(enrichmentNamTolL <- gopher(NamTolL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSensL <- gopher(NamSensL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvSensl <- gopher(NamTvSensl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvSR <- gopher(NamTvSR ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamLind <- gopher(NamLind ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTolEl <- gopher(NamTolEl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSensEl <- gopher(NamSensEl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSL <- gopher(NamSL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTL <- gopher(NamTL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvS <- gopher(NamTvS ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))

#system.time(enrichmentNamTest <- gopher(TestEv ,task = "kegg" ,background = bg,url="persimmon"))
#enrichmentNamTest$kegg[,c("id","padj")]

#TestE

# result is a list with the selected enrichment

enrichmentNamTolR$go[,c("id","padj")]
enrichmentNamSensR$go[,c("id","padj")]
enrichmentNamTvSens$go[,c("id","padj")]
enrichmentNamTvSL$go[,c("id","padj")]
enrichmentNamRind$go[,c("id","padj")]
enrichmentNamTolE$go[,c("id","padj")]
enrichmentNamSensE$go[,c("id","padj")]
enrichmentNamSR$go[,c("id","padj")]
enrichmentNamTR$go[,c("id","padj")]


ltest <- enrichmentNamTolR$pfam
alpha=.05
sel <- ltest$padj <= alpha
ltest$name <- sub(".*eucine *ich *epeat.*", "Leucine Rich Repeat", ltest$name)
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamSensR$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTvSensl$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTvSL$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamRind$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTolE$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamSensE$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()


#' We rearrange the analysis for creating a dataset LT vs LC, LS vs LC, RT vs RC and RS vs RC
#resSL <- results(ddsV,name = "ConditionT.TissueH")
#sum(resSL$padj <= 0.01 & !is.na(resSL$padj))

#resTL <- results(ddsV,name = "ConditionT.TissueH")
#sum(resTL$padj <= 0.01 & !is.na(resTL$padj))

ddsVL <- dds[annot[annot$Division=="Viridiplantae",1],]

colData(ddsVL)$Tissue <- relevel(colData(ddsVL)$Tissue,ref = "R")

ddsVL <- DESeq(ddsVL)


resultsNames(ddsVL)

resSL <- results(ddsVL,name = "ConditionS.TissueH")
sum(resSL$padj <= 0.01 & !is.na(resSL$padj))

resTL <- results(ddsVL,name = "ConditionT.TissueH")
sum(resTL$padj <= 0.01 & !is.na(resTL$padj))

resTSl <- results(ddsVL,contrast = c("Condition","T","S"))
sum(resTSl$padj <= 0.01 & !is.na(resTSl$padj))


plot.new()

grid.draw(venn.diagram(list(rownames(resSL[which(resSL$padj <= 0.01),]),
                            rownames(resTL[which(resTL$padj <= 0.01),]),
                            rownames(resTSl[which(resTSl$padj <= 0.01),])),
                       filename=NULL, category.names = c("SL","TL","TS"),
                       fill=pal[1:3], resolution=100))




#' Significant genes for tolerant roots

TolLeave <- setdiff(rownames(resTL[which(resTL$padj <= 0.01),]), 
                   union(rownames(resTSl[which(resTSl$padj <= 0.01),]),
                         rownames(resSL[which(resSL$padj <= 0.01),])))


tl.vst <- t(scale(t(vsd[as.character(TolLeave), ])))
sum(is.na(tl.vst))

heatTolL <- heatmap.2(tl.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(tl.vst),
                      labRow = NA)

cltl <- cutree(as.hclust(heatTolL$rowDendrogram),20)


TolLeave <- cbind(TolLeave,clusterID=cltl)

colnames(TolLeave)[1] <- colnames(B2GO_annot)[1]

ddsVR <- ddsV
colData(ddsVR)$Tissue <- relevel(colData(ddsVR)$Tissue, "R")
#colData(ddsVR)$Tissue <- ordered(colData(ddsVR)$Tissue,levels = c("R","H"))

ddsVR <- DESeq(ddsVR)
resultsNames(ddsVR)

resSRR <- results(ddsVR,name = "ConditionS.TissueH")
sum(resSRR$padj <= 0.01 & !is.na(resSRR$padj))

resRTR <- results(ddsVR,name = "ConditionT.TissueH")
sum(resRTR$padj <= 0.01 & !is.na(resRTR$padj))

SensLeave <- setdiff(rownames(resSL[which(resSL$padj <= 0.01),]), 
                    union(rownames(resTSl[which(resTSl$padj <= 0.01),]),
                          rownames(resTL[which(resTL$padj <= 0.01),])))



sl.vst <- t(scale(t(vsd[as.character(SensLeave), colnames(vsd)])))
sum(is.na(sl.vst))

heatResLS <- heatmap.2(sl.vst,
                       distfun=pearson.dist,
                       hclustfun=function(X){hclust(X,method="ward.D2")},
                       trace="none",
                       col=hpal,
                       labCol = colnames(sl.vst),
                       labRow = NA)

clls <- cutree(as.hclust(heatResLS$rowDendrogram),20)


SensLeave <- cbind(SensLeave,clusterID=clls)

colnames(SensLeave)[1] <- colnames(B2GO_annot)[1]


SensLeave <- merge(SensLeave, B2GO_annot, by="row.names", sort=FALSE)
write.csv(SensLeave, file = "SensLeave.csv")

#' Significant genes for T vs S

TolvS <- rownames(resTSl[which(resTSl$padj <= 0.01),]) 

stvs.vst <- t(scale(t(vsd[as.character(TolvS), ])))
sum(is.na(stvs.vst))

heatResTvS <- heatmap.2(stvs.vst,
                       distfun=pearson.dist,
                       hclustfun=function(X){hclust(X,method="ward.D2")},
                       trace="none",
                       col=hpal,
                       labCol = colnames(stvs.vst),
                       labRow = NA)

cltvs <- cutree(as.hclust(heatResTvS$rowDendrogram),20)


TolvS <- cbind(TolvS,clusterID=cltvs)

colnames(TolvS)[1] <- colnames(B2GO_annot)[1]


TolvS <- merge(TolvS, B2GO_annot, by="Sequence.Name", sort=FALSE)
write.csv(TolvS, file = "TolvS.csv")

#' Significant genes for T vs S only in leaves

TolvSensR <- setdiff(rownames(resTSl[which(resTSl$padj <= 0.01),]), 
                     union(rownames(resSL[which(resSL$padj <= 0.01),]),
                           rownames(resTL[which(resTL$padj <= 0.01),])))

stsr.vst <- t(scale(t(vsd[as.character(TolvSensR), ])))
sum(is.na(stsr.vst))

heatResTSR <- heatmap.2(stsr.vst,
                        distfun=pearson.dist,
                        hclustfun=function(X){hclust(X,method="ward.D2")},
                        trace="none",
                        col=hpal,
                        labCol = colnames(stsr.vst),
                        labRow = NA)

cltsr <- cutree(as.hclust(heatResTSL$rowDendrogram),20)


TolvSensR <- cbind(TolvSensR,clusterID=cltsr)

colnames(TolvSensR)[1] <- colnames(B2GO_annot)[1]


TolvSensR <- merge(TolvSensR, B2GO_annot, by="row.names", sort=FALSE)
write.csv(TolvSensR, file = "TolvSensR.csv")



#' Salt induced only in roots

Lind <- setdiff(intersect(rownames(resSL[which(resSL$padj <= 0.01),]),
                          rownames(resTL[which(resTL$padj <= 0.01),])),
                rownames(resTSl[which(resTSl$padj <= 0.01),]))

sli.vst <- t(scale(t(vsd[as.character(Lind), ])))
sum(is.na(sli.vst))

heatLi <- heatmap.2(sli.vst,
                    distfun=pearson.dist,
                    hclustfun=function(X){hclust(X,method="ward.D2")},
                    trace="none",
                    col=hpal,
                    labCol = colnames(sli.vst),
                    labRow = NA)

clir <- cutree(as.hclust(heatLi$rowDendrogram),20)


Lind <- cbind(Lind,clusterID=clir)

colnames(Lind)[1] <- colnames(B2GO_annot)[1]


Lind <- merge(Lind, B2GO_annot, by="row.names", sort=FALSE)
write.csv(Lind, file = "Lind.csv")

#' Tolerants > effect

TolEl <- setdiff(intersect(rownames(resTSl[which(resTSl$padj <= 0.01),]),
                          rownames(resTL[which(resTL$padj <= 0.01),])),
                rownames(resSL[which(resSL$padj <= 0.01),]))

stel.vst <- t(scale(t(vsd[as.character(TolEl), ])))
sum(is.na(stel.vst))

heatTolE <- heatmap.2(stel.vst,
                      distfun=pearson.dist,
                      hclustfun=function(X){hclust(X,method="ward.D2")},
                      trace="none",
                      col=hpal,
                      labCol = colnames(stel.vst),
                      labRow = NA)

cltole <- cutree(as.hclust(heatTolE$rowDendrogram),20)


TolEl <- cbind(TolEl,clusterID=cltole)

colnames(TolEl)[1] <- colnames(B2GO_annot)[1]


TolEl <- merge(TolEl, B2GO_annot, by="row.names", sort=FALSE)
write.csv(TolEl, file = "TolEl.csv")


#' Sensitive > effect

SensEl <- setdiff(intersect(rownames(resTSl[which(resTSl$padj <= 0.01),]),
                           rownames(resSL[which(resSL$padj <= 0.01),])),
                 rownames(resTL[which(resTL$padj <= 0.01),]))

ssel.vst <- t(scale(t(vsd[as.character(SensEl), ])))
sum(is.na(ssel.vst))

heatSensEl <- heatmap.2(ssel.vst,
                       distfun=pearson.dist,
                       hclustfun=function(X){hclust(X,method="ward.D2")},
                       trace="none",
                       col=hpal,
                       labCol = colnames(ssel.vst),
                       labRow = NA)

clsensel <- cutree(as.hclust(heatSensEl$rowDendrogram),20)


SensEl <- cbind(SensEl,clusterID=clsensel)

colnames(SensEl)[1] <- colnames(annot)[1]


SensEl <- merge(SensEl, B2GO_annot, by="row.names", sort=FALSE)
write.csv(SensEl, file = "SensEl.csv")

#' go enrichment example - to be implemented after christmas
#' (early january once the annotation are complete)
source("~/Git/UPSCb/src/R/gopher.R")
suppressPackageStartupMessages(source("~delhomme/Git/UPSCb/src/R/gopher.R"))

# dat = change for your list of genes of interest (a subset of a population)
# dat <- scan("~/Git/UPSCb/projects/facility/doc/tiggy-gene-example.txt",what="character")
# sub("^EC:","",TolvSens$Enzyme.Code)
# sub("|EC:"," ",TolvSens$Enzyme.Code)

NamTolL <- TolLeave$Sequence.Name.x
NamSensL <- SensLeave$Sequence.Name.x
NamTvSensl <- TolvS$Sequence.Name.x
NamTvSR <- TolvSensR
intersect(NamTvSR)
NamLind <- Lind$Sequence.Name.x
NamTolEl <- TolEl$Sequence.Name.x
NamSensEl <- SensEl$Sequence.Name.x
NamSL <- rownames(resSL[which(resSL$padj <= 0.01),])
NamTL <- rownames(resTL[which(resTL$padj <= 0.01),])
NamTvS <- TolvSens$Sequence.Name.x

#Select all KEGGS                 
#TestEnr <- B2GO_annot[!(B2GO_annot$Enzyme.Code == "" | is.na(B2GO_annot$Enzyme.Code)), ]
 #KEGGstats <- table(B2GO_annot$Enzyme.Code)
#TestEnr <- B2GO_annot[B2GO_annot$Enzyme.Code == "[EC:2.4.2.31]", ]
#TestEnr <- B2GO_annot[grep("EC:2.4.", B2GO_annot$Enzyme.Code), ]

# bg = your population (think what defines it.)
# bg <- scan("~/Git/UPSCb/projects/facility/doc/Tiggy.spruce.background.txt",what="character",skip=1)

#bg <- annot[annot$Division=="Viridiplantae",1]
#TestE <- merge(TestEnr, annot, by="Sequence.Name", sort=FALSE)



#TestEv <- TestE[which(TestE$Division.x=="Viridiplantae"),1]

#length(TestEv)

#TestEv2 <- I(intersect(TestEv, bg))
# bg$Enzyme.code <- sub("^EC:","",TolvSens$Enzyme.Code)
# we just quantify the run time
# enrichment <- gopher(dat,task = list("go","kegg","pfam"),background = bg,url="pabies")
ignore.case=TRUE

system.time(enrichmentNamTolL <- gopher(NamTolL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSensL <- gopher(NamSensL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvSensl <- gopher(NamTvSensl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvSR <- gopher(NamTvSR ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamLind <- gopher(NamLind ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTolEl <- gopher(NamTolEl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSensEl <- gopher(NamSensEl ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamSL <- gopher(NamSL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTL <- gopher(NamTL ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))
system.time(enrichmentNamTvS <- gopher(NamTvS ,task = list("go","kegg","pfam"),background = bg,url="persimmon"))

#system.time(enrichmentNamTest <- gopher(TestEv ,task = "kegg" ,background = bg,url="persimmon"))
#enrichmentNamTest$kegg[,c("id","padj")]

#TestE

# result is a list with the selected enrichment

enrichmentNamTolL$go[,c("id","padj")]
enrichmentNamSensL$go[,c("id","padj")]
enrichmentNamTvSensl$go[,c("id","padj")]
enrichmentNamTvSR$go[,c("id","padj")]
enrichmentNamLind$go[,c("id","padj")]
enrichmentNamTolEl$go[,c("id","padj")]
enrichmentNamSensEl$go[,c("id","padj")]
enrichmentNamSL$go[,c("id","padj")]
enrichmentNamTL$go[,c("id","padj")]
enrichmentNamTvS$go[,c("id","padj")]

ltest <- enrichmentNamTvS$pfam
ltest$name <- sub(".*eucine *ich *epeat.*", "Leucine Rich Repeat", ltest$name)
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTolL$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamSensL$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTvSensl$pfam
alpha=.05
sel <- ltest$padj <= alpha
ltest$name <- sub("*eucine *ich *epeat", "Leucine Rich Repeat", ltest$name)
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTvSR$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamLind$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.0,1),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamTolEl$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

ltest <- enrichmentNamSensEl$pfam
alpha=.05
sel <- ltest$padj <= alpha
wordcloud(ltest[sel,"name",drop=TRUE],
          unlist(ltest[sel,"nt"] / sum(ltest[sel,"nt"])),
          colors = brewer.pal(8,"Dark2"),scale = c(1.5,.5),rot.per = 0, res=300)
dev.null()

#' take that table to revigo (http://revigo.irb.hr/)

#' cluster based approach
library(WGCNA)

# info from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html

library(Mfuzz)
vignette("Mfuzz")

#' Get the data from KEGG
# to work with enzyme codes or once you got kegg enrichment
# find in files KEGGREST ; look at e.g. Git/UPSCb/projects/spruce-needles/src/R/seasonalDE.R



#'#'#' ================================================================
#'#' check that we get all annotation
#'#'#' ================================================================
stopifnot(sum(!unique(gsub(" \\(.*| var\\..*","",tax[,"tax"])) %in% nam.div$nam)==0)

#'#'#' ================================================================
#'#' extend the metadata
#'#'#' ================================================================
sprintf("These are %s different taxons represented in the data",nrow(nam.div))

#'#' these are distributed as follow:
par(mar=c(7.1,4.1,4.1,0.1))
barplot(table(annot$Division),
        main="Occurence of the different (super)kingdoms",las=2)
par(mar=c(5.1,4.1,4.1,2.1))

#' We have lost lot of plants in comparison with the previous version




sprintf("These originate from %s different UniRef90 IDs",nrow(tax))

#'#' and are distributed as follow
barplot(table(divison_name_df$div_nam[match(gsub(" \\(.*| var\\..*","",tax[,"tax"]),nam.div$nam)]))

sprintf("These UniRef IDs are matched by %s Wood FLcDNAs",sum(!is.na(mat$subject.id)))

#'#' and the taxon are distributed as follow
barplot(table(nam.div$div_nam[match(gsub(" \\(.*| var\\..*","",tax[match(mat$subject.id,tax[,"id"]),"tax"]),nam.div$nam)]))

#'#' and here are the 20 most common taxons
par(mar=c(4.1,15.1,0.1,0.1))
barplot(rev(rev(sort(table(tax[match(mat$subject.id,tax[,"id"]),"tax"],useNA="always")))[1:20]),horiz=TRUE,las=2,log="x")

mat$Taxon <- tax[match(mat$subject.id,tax[,"id"]),"tax"]
mat$Division <- nam.div$div_nam[match(gsub(" \\(.*| var\\..*","",tax[match(mat$subject.id,tax[,"id"]),"tax"]),nam.div$nam)]
head(mat)
#'#' and the most 10 common taxons for the different (super)kingdoms
lapply(lapply(lapply(lapply(split(mat$Taxon,mat$Division),table),sort),rev),"[",1:10)

#'#' and their occurence represented in a bar graph
barplot(do.call(cbind,lapply(lapply(lapply(lapply(split(mat$Taxon,mat$Division),table),sort),rev),"[",1:10)),beside=TRUE,log="y")

#'#' Interestingly most "Metazoan" are butterfly related (9 out of the 10). The 9th
#'#' is more interesting, Odobenus rosmarus divergens is: http://en.wikipedia.org/wiki/Walrus

#'#'#' =======================================================
#'#' write the updated meta-information
#'#'#' =======================================================

str(annot)

str(divison_name_df)


total <- merge(taxid_df, divison_name_df,by="Sequence.Description")

division <- left_join(annot, divison_name_df, by = "Sequence.Description", all=T)

write.csv(division,file="Filtering.csv")

#' # Session Info
#' ```{r, session info, echo=FALSE}
#' sessionInfo()
#' ```

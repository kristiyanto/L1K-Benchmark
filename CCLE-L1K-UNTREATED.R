####################### ABOUT THE SCRIPT #######################

# THIS SCRIPT COMPARES GENE EXPRESSION FOR UNPERTURB (UNTREATED)
# CELL LINES BETWEEN BROAD CCLE DATA AND 
# LINCS L1000 DATA (LEVEL 2.5 ESEMBLE, 2.5 INDIVIDUAL, LEVEL 3 AND LEVEL 4)
# 
# COMPARASION DONE BY MEASURING CORRELATION, SPEARMAN AND PEARSON
# 
# BY DANIEL KRISTIYANTO (DANIELKR@UW.EDU)
# AUTUMN 2015

################################################################


breakit <- function(x)
{
  y <- unlist(strsplit(x,"_"))[1]
  return(y)
}
removex <- function(x)
{
  y <- substring(x,2)
  return(y)
}

trimCCLE <- function(x)
{
  y                   <- as.data.frame(CCLE[CCLE$Description %in% names(x),c("Description",paste(CELL.LINES))])
  row.names(y)         <- y$Description
  y$Description        <- NULL
  return(y)
}

####################### CELL LINES TO COMPARE #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/")
CELL.LINES            <- read.csv(file="INTERSECT-CLINES-UNTRT.txt")
CELL.LINES            <- unlist(CELL.LINES)

####################### READ CCLE TABLE #######################

CCLE                  <- read.table("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/DATA/FROM KEBRA/V2_CCLE_Expression_Entrez_2012-09-29.gct", sep = "\t", header = T)    # File for CCLE
ccle.clines           <- names(CCLE)
ccle.clines           <- lapply(ccle.clines, breakit)
ccle.clines           <- unlist(ccle.clines)
names(CCLE)           <- ccle.clines

# 


####################### LEVEL 2.5 INDIVIDUAL #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-BENCHMARK/UNTRT/Lv25i/")

l2i             <- NULL
l2i.samples     <- NULL
l2i.spearman    <- NULL
l2i.pearson     <- NULL
l2i.average     <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l2i.samples   <- c(l2i.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l2i.average   <- rbind(l2i.average, colMeans(x))
  l2i           <- rbind(l2i,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="spearman"))
  l2i.spearman  <- c(l2i.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="pearson"))
  l2i.pearson   <- c(l2i.pearson,cor.pr)
}
l2i             <- cbind(l2i.samples,l2i)
l2i.cor         <- cbind("Lv.2.5i",l2i.samples,as.data.frame(l2i.spearman),as.data.frame(l2i.pearson))
row.names(l2i.average) <- CELL.LINES

####################### LEVEL 2.5  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/Lv25/")

l2          <- NULL
l2.cor      <- NULL
l2.spearman <- NULL
l2.pearson  <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  CCLE.trimmed  <- trimCCLE(x)
  x             <- x[1,row.names(CCLE.trimmed)]
  l2.spearman   <- c(l2.spearman,cor(as.numeric(x[1,]), CCLE.trimmed[,curr.clines], use="complete.obs", method="spearman"))
  l2.pearson    <- c(l2.pearson,cor(as.numeric(x[1,]), CCLE.trimmed[,curr.clines], use="complete.obs", method="pearson"))
  l2            <- rbind(l2,x)
}

l2              <- cbind(CELL.LINES,l2)
row.names(l2)   <- CELL.LINES
l2.cor          <- cbind("Lv.2.5",CELL.LINES,as.data.frame(l2.pearson), as.data.frame(l2.spearman))


####################### LEVEL 3  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/L3/")

l3             <- NULL
l3.samples     <- NULL
l3.spearman    <- NULL
l3.pearson     <- NULL
l3.average     <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l3.samples    <- c(l3.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l3.average   <- rbind(l3.average, colMeans(x))
  l3            <- rbind(l3,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="spearman"))
  l3.spearman   <- c(l3.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="pearson"))
  l3.pearson    <- c(l3.pearson,cor.pr)
}
l3             <- cbind(l3.samples,l3)
l3.cor         <- cbind("Lv.3",l3.samples,as.data.frame(l3.spearman),as.data.frame(l3.pearson))
row.names(l3.average) <- CELL.LINES


####################### LEVEL 4  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/L4/")

l4             <- NULL
l4.samples     <- NULL
l4.spearman    <- NULL
l4.pearson     <- NULL
l4.average     <- NULL 
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  CCLE.trimmed[CCLE.trimmed < 0]  <- NA
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l4.samples    <- c(l4.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l4.average    <- rbind(l4.average, colMeans(x))
  l4            <- rbind(l4,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="spearman"))
  l4.spearman   <- c(l4.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[names(a),curr.clines], use="complete.obs", method="pearson"))
  l4.pearson    <- c(l4.pearson,cor.pr)
}
l4              <- cbind(l4.samples,l4)
l4.cor          <- cbind("Lv.4",l4.samples,as.data.frame(l4.spearman),as.data.frame(l4.pearson))
row.names(l4.average) <- CELL.LINES

####################### EXPORT TO TABLEU  #######################

## AVERAGED VALUE FOR SCATTER PLOT
lv.2  <- cbind("Lv 2.5",l2)
lv.2i <- cbind("Lv 2.5i",row.names(l2i.average),as.data.frame(l2i.average))
lv.3  <- cbind("Lv 3",row.names(l3.average),as.data.frame(l3.average))
lv.4  <- cbind("Lv 4",row.names(l4.average),as.data.frame(l4.average))


names(lv.2)[1] <- "LEVEL"
names(lv.2)[2] <- "CELL-LINE"

names(lv.2i)[1] <- "LEVEL"
names(lv.2i)[2] <- "CELL-LINE"

names(lv.3)[1] <- "LEVEL"
names(lv.3)[2] <- "CELL-LINE"

names(lv.4)[1] <- "LEVEL"
names(lv.4)[2] <- "CELL-LINE"

raw.tableu    <- rbind(lv.2,lv.2i,lv.3,lv.4)
library(reshape2)
raw.tableu    <- melt(raw.tableu, id=c("LEVEL", "CELL-LINE"))
names(raw.tableu)[3] <- "GENE"
names(raw.tableu)[4] <- "L1K"
write.csv(raw.tableu, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/TABLEU/AllWells-L1K-un.csv", row.names = F)

#### CCLE DATA
CCLE.tableu <- CCLE[CCLE$Description %in% levels(raw.tableu$GENE),c("Description",levels(CELL.LINES))]
CCLE.tableu <- melt(CCLE.tableu, id=c("Description"))
names(CCLE.tableu)  <- c("GENE","CELL-LINE","CCLE")
write.csv(CCLE.tableu, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/TABLEU/CCLE.csv", row.names = F,na="")


## CORRELATION DATA

names(l2.cor)   <- c("Level","CELL-LINE","SPEARMAN","PEARSON")
names(l2i.cor)  <- c("Level","CELL-LINE","SPEARMAN","PEARSON")
names(l3.cor)   <- c("Level","CELL-LINE","SPEARMAN","PEARSON")
names(l4.cor)   <- c("Level","CELL-LINE","SPEARMAN","PEARSON")
cor.tableu      <- rbind(l2.cor,l2i.cor,l3.cor,l4.cor)
write.csv(cor.tableu, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/TABLEU/UNTREAT-CORR.csv", row.names = F)

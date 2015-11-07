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
  CCLE.trimmed         <- CCLE[,paste(CELL.LINES)]
  y                    <- CCLE.trimmed[CCLE$Description %in% names(x),]
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
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l2i.samples   <- c(l2i.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l2i           <- rbind(l2i,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="spearman"))
  l2i.spearman  <- c(l2i.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="pearson"))
  l2i.pearson   <- c(l2i.pearson,cor.pr)
}
l2i             <- cbind(l2i.samples,l2i)
l2i.cor         <- cbind(l2i.samples,l2i.spearman,l2i.pearson)

####################### LEVEL 2.5  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/Lv25/")

l2          <- NULL
l2.cor      <- NULL
l2.spearman <- NULL
l2.pearson  <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[1,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  l2.spearman   <- c(l2.spearman,cor(unlist(x), CCLE.trimmed[,curr.clines], use="complete.obs", method="spearman"))
  l2.pearson    <- c(l2.pearson,cor(unlist(x), CCLE.trimmed[,curr.clines], use="complete.obs", method="pearson"))
  l2            <- rbind(l2,x)
}

l2              <- cbind(CELL.LINES,l2)
row.names(l2)   <- CELL.LINES
l2.cor          <- cbind(CELL.LINES,l2.pearson, l2.spearman)


####################### LEVEL 3  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/L3/")

l3             <- NULL
l3.samples     <- NULL
l3.spearman    <- NULL
l3.pearson     <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l3.samples   <- c(l3.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l3           <- rbind(l3,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="spearman"))
  l3.spearman  <- c(l3.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="pearson"))
  l3.pearson   <- c(l3.pearson,cor.pr)
}
l3             <- cbind(l3.samples,l3)
l3.cor         <- cbind(l3.samples,l3.spearman,l3.pearson)


####################### LEVEL 4  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/UNTRT/L4/")

l4             <- NULL
l4.samples     <- NULL
l4.spearman    <- NULL
l4.pearson     <- NULL
for(curr.clines in CELL.LINES)
{
  x             <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  x             <- x[,intersect(names(x),CCLE$Description)]
  CCLE.trimmed  <- trimCCLE(x)
  well.no       <- seq(1,nrow(x),1)
  well.no       <- paste(curr.clines,well.no,sep="_")
  l4.samples   <- c(l4.samples,rep(curr.clines,nrow(x)))
  row.names(x)  <- well.no
  l4           <- rbind(l4,x)
  cor.sp        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="spearman"))
  l4.spearman  <- c(l4.spearman,cor.sp)
  cor.pr        <- apply(x,1, function(a) cor(a, CCLE.trimmed[,curr.clines], use="complete.obs", method="pearson"))
  l4.pearson   <- c(l4.pearson,cor.pr)
}
l4             <- cbind(l4.samples,l4)
l4.cor         <- cbind(l4.samples,l3.spearman,l4.pearson)


####################### EXPORT TO TABLEU  #######################

lv.2  <- cbind("Lv.2.5","NA",l2)
lv.2i <- cbind("Lv.2.5i",row.names(l2i),l2i)
lv.3  <- cbind("Lv.3",row.names(l3),l3)
lv.4  <- cbind("Lv.4",row.names(l4),l4)

names(lv.2)[1] <- "LEVEL"
names(lv.2)[2] <- "WELL-NO"
names(lv.2)[3] <- "CELL-LINE"


names(lv.2i)[1] <- "LEVEL"
names(lv.2i)[2] <- "WELL-NO"
names(lv.2i)[3] <- "CELL-LINE"


names(lv.3)[1] <- "LEVEL"
names(lv.3)[2] <- "WELL-NO"
names(lv.3)[3] <- "CELL-LINE"

names(lv.4)[1] <- "LEVEL"
names(lv.4)[2] <- "WELL-NO"
names(lv.4)[3] <- "CELL-LINE"

raw.tableu    <- rbind(lv.2,lv.2i,lv.3,lv.4)
library(reshape2)
raw.tableu    <- melt(raw.tableu, id=c("LEVEL","WELL-NO", "CELL-LINE"))
names(raw.tableu)[4] <- "GENE"
names(raw.tableu)[5] <- "EXP"

write.csv(raw.tableu, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/TABLEU/AllWells-L1K-un.csv", row.names = F)


CCLE.tableu <- CCLE[CCLE$Description %in% levels(raw.tableu$GENE),c("Description",levels(CELL.LINES))]
CCLE.tableu <- melt(CCLE.tableu, id=c("Description"))
names(CCLE.tableu) <- c("GENE","CELL-LINE","CCLE")
write.csv(raw.tableu, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/TABLEU/CCLE.csv", row.names = F)

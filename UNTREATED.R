
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


do.corr <- function(x)
{
  s2i <- cor(l2i[,x], CCLE.trimmed[,x], use="complete.obs", method="spearman") 
  p2i <- cor(l2i[,x], CCLE.trimmed[,x], use="complete.obs", method="pearson") 
  
  s2 <- cor(l2[,x], CCLE.trimmed[,x], use="complete.obs", method="spearman") 
  p2 <- cor(l2[,x], CCLE.trimmed[,x], use="complete.obs", method="pearson") 
  
  s3 <- cor(as.numeric(l3[,x]), CCLE.trimmed[,x], use="complete.obs", method="spearman") 
  p3 <- cor(as.numeric(l3[,x]), CCLE.trimmed[,x], use="complete.obs", method="pearson") 
  
  s4 <- cor(as.numeric(l4[,x]), CCLE.trimmed[,x], use="complete.obs", method="spearman") 
  p4 <- cor(as.numeric(l4[,x]), CCLE.trimmed[,x], use="complete.obs", method="pearson") 
  
  y <- rbind(s2i,p2i,s2,p2,s3,p3,s4,p4)
  return(y)
}
####################### READ CCLE TABLE #######################

CCLE                  <- read.table("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/DATA/FROM KEBRA/V2_CCLE_Expression_Entrez_2012-09-29.gct", sep = "\t", header = T)    # File for CCLE
ccle.clines           <- names(CCLE)
ccle.clines           <- lapply(ccle.clines, breakit)
ccle.clines           <- unlist(ccle.clines)
names(CCLE)           <- ccle.clines


CCLE.trimmed                   <- CCLE[CCLE$Description %in% row.names(l2i),]
row.names(CCLE.trimmed)        <- CCLE.trimmed$Description
CCLE.trimmed                   <- CCLE.trimmed[,paste(CELL.LINES)]
CCLE.trimmed                   <- CCLE.trimmed[order(row.names(CCLE.trimmed)),]


####################### CELL LINES TO COMPARE #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/")
CELL.LINES            <- read.csv(file="INTERSECT-CLINES-UNTRT.txt")
CELL.LINES            <- unlist(CELL.LINES)
####################### LEVEL 2.5 INDIVIDUAL #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/UNTRT/Lv25i/")

raw.l2i <- NULL
for(curr.clines in CELL.LINES)
{
  x       <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  raw.l2i <- rbind(raw.l2i,x[1,])  
}
l2i             <- raw.l2i
l2i.samples     <- l2i[,1]
l2i$Sample      <- NULL
l2i             <- l2i[,intersect(names(l2i),CCLE$Description)]
l2i             <- as.data.frame(t(l2i))
colnames(l2i)   <- CELL.LINES
l2i             <- l2i[order(row.names(l2i)),]

####################### LEVEL 2.5  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/UNTRT/Lv25/")

raw.l2 <- NULL
for(curr.clines in CELL.LINES)
{
  x       <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  raw.l2 <- rbind(raw.l2,x[1,])  
}

l2             <- raw.l2
l2.samples     <- l2[,1]
l2$Sample      <- NULL
l2             <- l2[,intersect(names(l2),CCLE$Description)]
l2             <- as.data.frame(t(l2))
colnames(l2)   <- CELL.LINES
l2             <- l2[order(row.names(l2)),]


####################### LEVEL 3  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/UNTRT/L3/")

raw.l3 <- NULL
for(curr.clines in CELL.LINES)
{
  x       <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  y       <- colMeans(x[2:ncol(x)])
  raw.l3  <- rbind(raw.l3,c(curr.clines,y))
}

l3             <- as.data.frame(raw.l3)
l3             <- l3[,intersect(names(l3),CCLE$Description)]
l3             <- as.data.frame(t(l3))
colnames(l3)   <- CELL.LINES
l3             <- l3[order(row.names(l3)),]


####################### LEVEL 4  #######################
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/UNTRT/L4/")

raw.l4 <- NULL
for(curr.clines in CELL.LINES)
{
  x       <- read.csv(file = paste0(curr.clines,".txt"),sep="\t")
  y       <- colMeans(x[2:ncol(x)])
  raw.l4  <- rbind(raw.l4,c(curr.clines,y))
}

l4             <- as.data.frame(raw.l4)
l4             <- l4[,intersect(names(l4),CCLE$Description)]
l4             <- as.data.frame(t(l4))
colnames(l4)   <- CELL.LINES
l4             <- l4[order(row.names(l4)),]


####################### CORRELATION  #######################

RESULT <- NULL
for(curr.clines in CELL.LINES)
{
  RESULT      <- cbind(RESULT,do.corr(curr.clines))
}
RESULT        <- as.data.frame(RESULT)
names(RESULT) <- CELL.LINES


####################### EXPORT TO TABLEU  #######################

res.tab <- cbind(c("Lv2.5-Ind","Lv2.5-Ind","Lv2.5","Lv2.5","Lv3","Lv3","Lv4","Lv4"), RESULT)
res.tab <- cbind(rep(c("Spearman","Pearson"),4),res.tab)
names(res.tab)[1] <- "Method"
names(res.tab)[2] <- "Level"
library(reshape2)
res.tab         <- melt(res.tab, id=c("Method","Level"))
names(res.tab)  <- c("Method","Level","Cell.Line","Correlation")
write.csv(res.tab, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/TABLEU/CCLEvsL1K-Untreated.csv", row.names = F)


l2i.tab           <- cbind("Level2.5.Ind", row.names(l2i.tab),l2i)
l2.tab            <- cbind("Level2.5", row.names(l2i.tab), l2)
l3.tab            <- cbind("Level3", row.names(l2i.tab),l3)
l4.tab            <- cbind("Level4", row.names(l2i.tab),l4)
CCLE.trimmed.tab  <- cbind("CCLE",row.names(l2i.tab),CCLE.trimmed)
names(l2i.tab)[1] <- "Level"
names(l2.tab)[1]  <- "Level"
names(l3.tab)[1]  <- "Level"
names(l4.tab)[1]  <- "Level"
names(CCLE.trimmed.tab)[1]  <- "Level"

names(l2i.tab)[2] <- "Gene"
names(l2.tab)[2]  <- "Gene"
names(l3.tab)[2]  <- "Gene"
names(l4.tab)[2]  <- "Gene"
names(CCLE.trimmed.tab)[2]  <- "Gene"
raw.tab           <- rbind(l2i.tab,l2.tab,l3.tab,l4.tab)
raw.tab           <- melt(raw.tab, id=c("Level","Gene"))
raw.tab.ccle      <- melt(CCLE.trimmed.tab, id=c("Level","Gene"))
names(raw.tab)    <- c("Level","Gene","CellLine","Expression")
names(raw.tab.ccle)    <- c("Level","Gene","CellLine","Expression")
write.csv(raw.tab, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/TABLEU/CCLEvsL1K-Untreated-Raw.csv", row.names = F)
write.csv(raw.tab, file="/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/TABLEU/CCLEvsL1K-Untreated-ccle.csv", row.names = F)


#########################################################################
# THIS SCRIPT IS TO FIND THE CELL LINES TO COMPARE BETWEEN LINCS L1K 
# AND CCLE. 
# THE OUTPUT IS THE CELL LINES THAT ARE FOUND BOTH IN CCLE AND L1K.
# BY DANIEL KRISTIYANTO (DANIELKR@UW.EDU)
# AUTUMN 2015
#########################################################################
breakit <- function(x)
{
  y <- unlist(strsplit(x,"_"))[1]
  return(y)
}


run.this=0
if(run.this==1)
  
#########################################################################
# THIS CELL LINES LIST IS ACCUIRED BY RUNNING A QUERY IN LINCS TOOL:
# SELECTING ALL CELL LINES GIVEN PERTUBATION DESCRIPTION UNTREATED.
#########################################################################

l1k.untrt <- factor(c("A375",
        "A375",
        "A375",
        "A375",
        "A375",
        "A549",
        "A549",
        "ASC",
        "FIBRNPC",
        "FIBRNPC",
        "HA1E",
        "HA1E",
        "HA1E",
        "HCC515",
        "HEK293T",
        "HEKTE",
        "HEPG2",
        "HEPG2",
        "HT29",
        "JURKAT",
        "MCF7",
        "MCF7",
        "MCF7",
        "NEU",
        "NEU",
        "NEU.KCL",
        "NEU.KCL",
        "NKDBA",
        "NPC",
        "NPC",
        "NPC",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "PC3",
        "SHSY5Y",
        "SKL",
        "SW480",
        "VCAP",
        "VCAP",
        "VCAP"))

#########################################################################
# DO INTERSECTION OF CELL LINES FOUND IN L1K WITH CCLE
#########################################################################

l1k.clines            <- levels(l1k.untrt)
CCLE.data             <- read.table("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K/DATA/FROM KEBRA/V2_CCLE_Expression_Entrez_2012-09-29.gct", sep = "\t", header = T)    # File for CCLE
ccle.clines           <- names(CCLE.data)
ccle.clines           <- lapply(ccle.clines, breakit)       # ONLY USE THE FIRST WORD BEFORE '_' IN CCLE
ccle.clines           <- unlist(ccle.clines)                
intersect.clines      <- intersect(ccle.clines,l1k.clines)  # THIS IS THE RESULT

setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/")     
write.csv(intersect.clines,file = "INTERSECT-CLINES-UNTRT.txt", row.names = F)

#########################################################################
# DO INTERSECTION OF CELL LINES FOUND IN L1K WITH ARRAY EXPRESS
#########################################################################
library(stringr)
ax.meta               <- read.table("DATA/E_MTAB_2706_sdrf.txt", header = T, sep="\t",stringsAsFactors = F, quote = NULL, comment = "")
ax.celines            <- levels(factor(ax.meta$Characteristics.cell.line.))
ax.celines            <- str_replace_all(ax.celines, "-", "")
ax.celines            <- str_replace_all(ax.celines, " ", "")
ax.celines            <- str_replace_all(ax.celines, "'.'", "")
ax.celines            <- str_replace_all(ax.celines, "/", "")

ax.intersect          <- intersect(l1k.clines,ax.celines)
setwd("/Users/Daniel/Google Drive/BIOINFORMATICS/L1K-Benchmark/")     
write.csv(ax.intersect,file = "INTERSECT-CLINES-ARRAYX-UNTRT.txt", row.names = F)

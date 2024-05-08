## ========================================== ##
##             Psi Conservation               ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 11/28/2023
#@ content: Figure 3
#@ version: 2.0

# housekeeping ------------------------------------------------------------
list.of.packages <- c("ggplot2", "ggExtra","dplyr", "tidyr", "plyr", "readr", "stringr", "tibble", "stringr", "reshape2", "limma", "plotly", "hrbrthemes", "caret", "venneuler")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) # from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them. Thanks!
library(ggplot2)
library(ggExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(reshape2)
library(limma)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
library(caret)
library(janitor)
library(gridExtra)
library(viridis)
library(epos)
library(pracma)
library(tidyverse)
library(plotly)
library(rtracklayer)

###--- Clear Space ---###
rm(list = ls()) #clear global environment except for packages & libraries 
cat("\f") #clears console

###--- General Path & Working Directory ---###
# path <- scan(what = "character", n=1) #actual directory needs to be on the following line
# "/home/cam/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

path <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

dirPath <- print(path)
setwd(dirPath) #set working directory

dataPath <- print(paste(dirPath, "DATA", sep = "/"))
pvalPath <- print(paste(dataPath, "pValues", sep = "/"))
figPath <- print(paste(dirPath, "Figures", "Figure-3","FIG3-elements", sep = "/"))
bamPath <- print(paste(dirPath, "BAM", "hg38v10", sep = "/"))
tpmPath <- print(paste(bamPath, "stringtie",sep = "/"))
pysamPath <- print(paste(dataPath, "pysamstats", sep = "/"))
pvalSOURCE <- print(paste(dirPath,"pValues","DATA", sep = "/"))
kpPath <- print(paste(dataPath, "seqLOGO", sep = "/"))
kdePath <- print(paste(dataPath, "KDE", sep = "/"))
rdsPath <- print(paste(dataPath, "RData", sep = "/"))
orthoPATH <- print(paste(dataPath, "orthogonal", sep = "/"))


###--- Colors ---###
cp1 <- c("#FF595E", "#FF924C", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93")

CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")

cat("\f") #clears console

# DATA FRAME: Ranking -----------------------------------------------------

if (T) {
  master.df <- read_csv(paste(pvalPath,"sixCellLine.pval.filtered.master.csv", sep = "/"))
  ranking <- master.df
  
  ranking$A549.mm.Direct[ranking$A549.psi == 0] <- 0
  ranking$HeLa.mm.Direct[ranking$HeLa.psi == 0] <- 0
  ranking$HepG2.mm.Direct[ranking$HepG2.psi == 0] <- 0
  ranking$Jurkat.mm.Direct[ranking$Jurkat.psi == 0] <- 0
  ranking$NTERA.mm.Direct[ranking$NTERA.psi == 0] <- 0
  ranking$SHSY5Y.mm.Direct[ranking$SHSY5Y.psi == 0] <- 0
  
  #panIVT.df <- read_csv(paste(pvalPath, "panIVT.mm.csv", sep="/"))
  #panIVT.df <- panIVT.df[,c("ACP", "mm.IVT")]
  #ranking$mm.IVT <- panIVT.df$mm.IVT[which(panIVT.df$ACP == ranking$ACP)]
  
  ranking.1 <- ranking[which(ranking$rank == 1),]
  ranking.2 <- ranking[which(ranking$rank == 2),]
  ranking.3 <- ranking[which(ranking$rank == 3),]
  ranking.4 <- ranking[which(ranking$rank == 4),]
  ranking.5 <- ranking[which(ranking$rank == 5),]
  ranking.6 <- ranking[which(ranking$rank == 6),]
  
  r1 <- nrow(ranking.1)
  r2 <- nrow(ranking.2)
  r3 <- nrow(ranking.3)
  r4 <- nrow(ranking.4)
  r5 <- nrow(ranking.5)
  r6 <- nrow(ranking.6)
  
  print(paste("1 Cell Line",r1, sep = " -- "))
  print(paste("2 Cell Line",r2, sep = " -- "))
  print(paste("3 Cell Line",r3, sep = " -- "))
  print(paste("4 Cell Line",r4, sep = " -- "))
  print(paste("5 Cell Line",r5, sep = " -- "))
  print(paste("6 Cell Line",r6, sep = " -- "))
  
}

# Individual Ranking Plots ------------------------------------------------

##--BASELINE CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct","psi", "N_reads_Direct", "p.value.Direct"), paste, sep="."))
CL.mm.sum <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- ranking.1[,c("ACP","Annotation", "chr", "position", "kmer",CL.mm)]
mm$filter <- 0

min.mm <- 30

#-ONE-#
mm$filter[rowSums(mm[, CL.mm.sum] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]

if (T) {
  mm.one$cell <- 0
  mm.one$CP <- paste(mm.one$chr,mm.one$position, sep =":")
  for (i in 1:nrow(mm.one)) {
    if (mm.one$A549.psi[i] > 0){mm.one$cell[i] = 1} 
    if (mm.one$HeLa.psi[i] > 0) {mm.one$cell[i] = 2} 
    if (mm.one$HepG2.psi[i] > 0) {mm.one$cell[i] = 3}
    if (mm.one$Jurkat.psi[i] > 0) {mm.one$cell[i] = 4}
    if (mm.one$NTERA.psi[i] > 0) {mm.one$cell[i] = 5} 
    if (mm.one$SHSY5Y.psi[i] > 0) {mm.one$cell[i] = 6}
  }
}


A549.1 <- mm.one[which(mm.one$A549.psi == 1 & mm.one$A549.mm.Direct >= min.mm),]
HeLa.1 <- mm.one[which(mm.one$HeLa.psi == 1 & mm.one$HeLa.mm.Direct >= min.mm),]
HepG2.1 <- mm.one[which(mm.one$HepG2.psi == 1 & mm.one$HepG2.mm.Direct >= min.mm),]
Jurkat.1 <- mm.one[which(mm.one$Jurkat.psi == 1 & mm.one$Jurkat.mm.Direct >= min.mm),]
NTERA.1 <- mm.one[which(mm.one$NTERA.psi == 1 & mm.one$NTERA.mm.Direct >= min.mm),]
SHSY5Y.1 <- mm.one[which(mm.one$SHSY5Y.psi == 1 & mm.one$SHSY5Y.mm.Direct >= min.mm),]

# A549.1 <- mm.one[which(mm.one$A549.psi == 1),]
# HeLa.1 <- mm.one[which(mm.one$HeLa.psi == 1),]
# HepG2.1 <- mm.one[which(mm.one$HepG2.psi == 1),]
# Jurkat.1 <- mm.one[which(mm.one$Jurkat.psi == 1),]
# NTERA.1 <- mm.one[which(mm.one$NTERA.psi == 1),]
# SHSY5Y.1 <- mm.one[which(mm.one$SHSY5Y.psi == 1),]
# 



# write.csv(A549.1,print(paste(pvalPath,"A549.1.csv", sep = "/")), row.names = F)
# write.csv(HeLa.1,print(paste(pvalPath,"HeLa.1.csv", sep = "/")), row.names = F)
# write.csv(HepG2.1,print(paste(pvalPath,"HepG2.1.csv", sep = "/")), row.names = F)
# write.csv(Jurkat.1,print(paste(pvalPath,"Jurkat.1.csv", sep = "/")), row.names = F)
# write.csv(NTERA.1,print(paste(pvalPath,"NTERA.1.csv", sep = "/")), row.names = F)
# write.csv(SHSY5Y.1,print(paste(pvalPath,"SHSY5Y.1.csv", sep = "/")), row.names = F)


print(nrow(A549.1))
print(nrow(HeLa.1))
print(nrow(HepG2.1))
print(nrow(Jurkat.1))
print(nrow(NTERA.1))
print(nrow(SHSY5Y.1))


A549.1$SSrank <- rank(A549.1$A549.mm.Direct, ties.method="random")
HeLa.1$SSrank <- rank(HeLa.1$HeLa.mm.Direct, ties.method="random")
HepG2.1$SSrank <- rank(HepG2.1$HepG2.mm.Direct, ties.method="random")
Jurkat.1$SSrank <- rank(Jurkat.1$Jurkat.mm.Direct, ties.method="random")
NTERA.1$SSrank <- rank(NTERA.1$NTERA.mm.Direct, ties.method="random")
SHSY5Y.1$SSrank <- rank(SHSY5Y.1$SHSY5Y.mm.Direct, ties.method="random")

r1.df <- rbind(data.frame(A549.1), 
               data.frame(HeLa.1), 
               data.frame(HepG2.1),
               data.frame(Jurkat.1),
               data.frame(NTERA.1),
               data.frame(SHSY5Y.1))
#r1.df$SSrank <- sprintf("%04d", r1.df$SSrank)

write.csv(r1.df,print(paste(pvalPath,"r1.df.csv", sep = "/")), row.names = F)

r1.df$ID <- paste(r1.df$cell, r1.df$SSrank, sep ="-" )


r1.out <- r1.df
r1.out$CP <- paste(r1.out$chr, r1.out$position, sep = ":")
write.csv(r1.out,print(paste(pvalPath,"rank1-CP.csv", sep = "/")), row.names = F)

r1.df$A549.mm.Direct[r1.df$A549.psi == 0] <- 0
r1.df$HeLa.mm.Direct[r1.df$HeLa.psi == 0] <- 0
r1.df$HepG2.mm.Direct[r1.df$HepG2.psi == 0] <- 0
r1.df$Jurkat.mm.Direct[r1.df$Jurkat.psi == 0] <- 0
r1.df$NTERA.mm.Direct[r1.df$NTERA.psi == 0] <- 0
r1.df$SHSY5Y.mm.Direct[r1.df$SHSY5Y.psi == 0] <- 0


r1.df$CP <- paste(r1.df$chr, r1.df$position, sep = ":")
r1.df$kmer.start <- r1.df$position - 5
r1.df$kmer.end <- r1.df$position + 5
r1.df$range <- paste(r1.df$kmer.start, r1.df$kmer.end, sep = "-")
r1.df$kmer.11 <- paste(r1.df$chr,r1.df$range , sep = ":")
write.csv(r1.df,print(paste(pvalPath,"rank1-CP.csv", sep = "/")), row.names = F)






r1.df$color <- "other"
r1.df$color[r1.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
r1.df$color[r1.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"




r1.2 <- r1.df[which(r1.df$A549.mm.Direct >= 30 |
                     r1.df$HeLa.mm.Direct >= 30 |
                     r1.df$HepG2.mm.Direct >= 30 |
                     r1.df$Jurkat.mm.Direct >= 30 |
                     r1.df$NTERA.mm.Direct >= 30 |
                     r1.df$SHSY5Y.mm.Direct >= 30),]
kmer <- r1.2$kmer.11
#kmer <- toupper(kmer)
write.table(kmer,print(paste(pvalPath,"rank1-KDErange.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)

annot <- r1.2$Annotation
#kmer <- toupper(kmer)
write.table(annot,print(paste(pvalPath,"rank1-KDEannot.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)




library(pracma)


cp2 <- c("#FFFFFF","#0450B4", "#046DC8", "#1184A7", "#15A2A2", "#6FB1A0", "#B4418E", "#D94A8C", "#EA515F",
         "#FE7434", "#FEA802", "#ffd93d")

breaks <- c(0,9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99, 99.99, 101)
labels <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100")

heatMAP.df <- r1.df[,c("ID", "A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")]
reps <- c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")
HM_melt <- melt(heatMAP.df, measure.vars = reps)
HM_melt$value[HM_melt$value < 0] <- 0
HM_melt$category <- cut(HM_melt$value, breaks = breaks, labels = labels, right = FALSE)
HM_melt$category <- factor(HM_melt$category, levels = labels)
HM_melt$ID <- factor(HM_melt$ID, levels = rev(sort(unique(HM_melt$ID))))
all.HeatMap <- HM_melt %>% 
  ggplot(aes(x = variable , y = ID, fill = category)) +
  geom_tile()+
  ylab("40+ %Mismatch for significant pvalue in 1 cell line") +
  #ggtitle("binary.TPM10.p001.mm10")+
  # scale_fill_viridis(discrete=TRUE)+
  scale_fill_manual(values = setNames(cp2, labels)) +
  labs(fill = "Mismatch (%)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype=1),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")
  )
all.HeatMap
name <- print(paste("rank1.BinaryHeatMapbinary.eps", sep = ""))
ggsave(print(paste(figPath, name, sep = "/")))


## Illustrator Annotation
g.1 <- r1.df[,c("Annotation", "chr", "position", "ID")]
g.1 %>% arrange(ID) %>% write.csv(print(paste(figPath,"rank1.csv", sep = "/")), row.names = F)

g.1.kmer <- r1.df[,c("kmer","ID")]
g.1.kmer$kmer <- gsub("T", "U", g.1.kmer$kmer)
g.1.kmer %>% arrange(ID) %>% write.csv(print(paste(figPath,"rank1.kmer.csv", sep = "/")), row.names = F)

g.1.kmer2 <- g.1.kmer %>% arrange(ID)
g.1.kmer2 <- g.1.kmer2[,c("kmer")]
write.csv(g.1.kmer2,print(paste(figPath,"rank1.kmer2.csv", sep = "/")), row.names = F)

#-KMER MOTIF COLORS-#
g.1.kmer$color <- "other"
g.1.kmer$color[g.1.kmer$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
g.1.kmer$color[g.1.kmer$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"




# DATA FRAME: Orthogonal --------------------------------------------------

if (T) {
  r.ortho <- read_csv(paste(orthoPATH, "PSI.orthogonal.MASTER.csv", sep = "/"))
  
  one <- r.ortho[which(r.ortho$rank == 1 & r.ortho$ortho > 0), ]
  two <- r.ortho[which(r.ortho$rank == 2 & r.ortho$ortho > 0), ]
  three <- r.ortho[which(r.ortho$rank == 3 & r.ortho$ortho > 0), ]
  four <- r.ortho[which(r.ortho$rank == 4 & r.ortho$ortho > 0), ]
  five <- r.ortho[which(r.ortho$rank == 5 & r.ortho$ortho > 0), ]
  
  print(paste("rank1 Ψ othrogonally confirmed: ", nrow(one), sep = ""))
  print(paste("rank2 Ψ othrogonally confirmed: ", nrow(two), sep = ""))
  print(paste("rank3 Ψ othrogonally confirmed: ", nrow(three), sep = ""))
  print(paste("rank4 Ψ othrogonally confirmed: ", nrow(four), sep = ""))
  print(paste("rank5 Ψ othrogonally confirmed: ", nrow(five), sep = ""))
  
  
}

# FIGURE: Rank5 -----------------------------------------------------------

##--BASELINE CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- r.ortho[which(r.ortho$ortho > 0),]
mm <- five

mm$filter <- 0

min.mm <- 20

#-ONE-#
mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]

cat("A549: ", length(which(mm.one$A549.psi == 1)), "\n",
    "HeLa: ", length(which(mm.one$HeLa.psi == 1)), "\n",
    "HepG2: ", length(which(mm.one$HepG2.psi == 1)), "\n",
    "Jurkat: ", length(which(mm.one$Jurkat.psi == 1)), "\n",
    "NTERA: ", length(which(mm.one$NTERA.psi == 1)), "\n",
    "SH-SY5Y: ", length(which(mm.one$SHSY5Y.psi == 1)), "\n", sep = "")

mm.one$A549.mm.Direct[mm.one$A549.psi == 0] <- 0
mm.one$HeLa.mm.Direct[mm.one$HeLa.psi == 0] <- 0
mm.one$HepG2.mm.Direct[mm.one$HepG2.psi == 0] <- 0
mm.one$Jurkat.mm.Direct[mm.one$Jurkat.psi == 0] <- 0
mm.one$NTERA.mm.Direct[mm.one$NTERA.psi == 0] <- 0
mm.one$SHSY5Y.mm.Direct[mm.one$SHSY5Y.psi == 0] <- 0

a <- mm.one[,c("Annotation", "chr", "position", "BID", "CeU", "psi", "PRAISE", "RBS", "ortho")]


mm.one$SS <- round(apply(mm.one[,c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")], 1, sd), digits = 2)
#mm.one$SSrank <- sprintf("%04f", mm.one$SS)
#mm.one <- mm.one[order(mm.one$SSrank),]
mm.one <- mm.one[order(mm.one$SS),]

SSrank <- paste(mm.one$SS, mm.one$Annotation, sep = "-")

##--PLOT--##

cp2 <- c("#D5D5D5","#0450B4", "#046DC8", "#1184A7", "#15A2A2", "#6FB1A0", "#B4418E", "#D94A8C", "#EA515F",
         "#FE7434", "#FEA802", "#ffd93d")

breaks <- c(0,9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99, 99.99, 101)
labels <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100")

heatMAP.df <- mm.one[,c("Annotation","SS", "A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")]
reps <- c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")
HM_melt <- melt(heatMAP.df, measure.vars = reps)
HM_melt$value[HM_melt$value < 0] <- 0
HM_melt$category <- cut(HM_melt$value, breaks = breaks, labels = labels, right = FALSE)
HM_melt$category <- factor(HM_melt$category, levels = labels)
HM_melt$SS <- factor(HM_melt$SS, levels = rev(sort(unique(HM_melt$SS))))
all.HeatMap <- HM_melt %>% 
  ggplot(aes(x = variable , y = Annotation, fill = category)) +
  geom_tile()+
  #ylab("40+ %Mismatch for significant pvalue in 1 cell line") +
  #ggtitle("binary.TPM10.p001.mm10")+
  # scale_fill_viridis(discrete=TRUE)+
  scale_fill_manual(values = setNames(cp2, labels)) +
  labs(fill = "Mismatch (%)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype=1),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")
  )
all.HeatMap
name <- print(paste("rank5.BinaryHeatMapbinary.eps", sep = ""))
ggsave(print(paste(figPath, name, sep = "/")))




g.1 <- r5[,c("Annotation", "chr", "position", "SS")]

g.1 %>% arrange(SS) %>% write.csv(print(paste(figPath,"rank5.csv", sep = "/")), row.names = F)

#write.csv(g.1,print(paste(figPath,"rank5.csv", sep = "/")), row.names = F)
g.1.kmer <- r5[,c("kmer","SS")]
g.1.kmer %>% arrange(SS) %>% write.csv(print(paste(figPath,"rank5.kmer.csv", sep = "/")), row.names = F)

g.1.kmer2 <- g.1.kmer %>% arrange(SS)
g.1.kmer2 <- g.1.kmer2[,c("kmer")]
write.csv(g.1.kmer2,print(paste(figPath,"rank5.kmer2.csv", sep = "/")), row.names = F)

#write.csv(g.1.kmer,print(paste(figPath,"rank5.kmer.csv", sep = "/")), row.names = F)



# FIGURE: Rank3 -----------------------------------------------------------
##--BASELINE CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- r.ortho[which(r.ortho$ortho > 0),]
mm <- three

mm$filter <- 0

min.mm <- 20

#-ONE-#
mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]

cat("A549: ", length(which(mm.one$A549.psi == 1)), "\n",
    "HeLa: ", length(which(mm.one$HeLa.psi == 1)), "\n",
    "HepG2: ", length(which(mm.one$HepG2.psi == 1)), "\n",
    "Jurkat: ", length(which(mm.one$Jurkat.psi == 1)), "\n",
    "NTERA: ", length(which(mm.one$NTERA.psi == 1)), "\n",
    "SH-SY5Y: ", length(which(mm.one$SHSY5Y.psi == 1)), "\n", sep = "")

mm.one$A549.mm.Direct[mm.one$A549.psi == 0] <- 0
mm.one$HeLa.mm.Direct[mm.one$HeLa.psi == 0] <- 0
mm.one$HepG2.mm.Direct[mm.one$HepG2.psi == 0] <- 0
mm.one$Jurkat.mm.Direct[mm.one$Jurkat.psi == 0] <- 0
mm.one$NTERA.mm.Direct[mm.one$NTERA.psi == 0] <- 0
mm.one$SHSY5Y.mm.Direct[mm.one$SHSY5Y.psi == 0] <- 0

a <- mm.one[,c("Annotation", "chr", "position", "BID", "CeU", "psi", "PRAISE", "RBS", "ortho")]


mm.one$SS <- round(apply(mm.one[,c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")], 1, sd), digits = 2)
#mm.one$SSrank <- sprintf("%04f", mm.one$SS)
#mm.one <- mm.one[order(mm.one$SSrank),]
mm.one <- mm.one[order(mm.one$SS),]

SSrank <- paste(mm.one$SS, mm.one$Annotation, sep = "-")

##--PLOT--##

cp2 <- c("#D5D5D5","#0450B4", "#046DC8", "#1184A7", "#15A2A2", "#6FB1A0", "#B4418E", "#D94A8C", "#EA515F",
         "#FE7434", "#FEA802", "#ffd93d")

breaks <- c(0,9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99, 99.99, 101)
labels <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100")

heatMAP.df <- mm.one[,c("Annotation","SS", "A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")]
reps <- c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")
HM_melt <- melt(heatMAP.df, measure.vars = reps)
HM_melt$value[HM_melt$value < 0] <- 0
HM_melt$category <- cut(HM_melt$value, breaks = breaks, labels = labels, right = FALSE)
HM_melt$category <- factor(HM_melt$category, levels = labels)
HM_melt$SS <- factor(HM_melt$SS, levels = rev(sort(unique(HM_melt$SS))))
all.HeatMap <- HM_melt %>% 
  ggplot(aes(x = variable , y = Annotation, fill = category)) +
  geom_tile()+
  #ylab("40+ %Mismatch for significant pvalue in 1 cell line") +
  #ggtitle("binary.TPM10.p001.mm10")+
  # scale_fill_viridis(discrete=TRUE)+
  scale_fill_manual(values = setNames(cp2, labels)) +
  labs(fill = "Mismatch (%)") +
  theme(axis.line.x = element_line(size = 0.5, colour = "black", linetype=1),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(), 
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")
  )
all.HeatMap
name <- print(paste("rank3.BinaryHeatMapbinary.eps", sep = ""))
ggsave(print(paste(figPath, name, sep = "/")))



g.1 <- mm.one[,c("Annotation", "chr", "position")]
g.1 %>% arrange(SS) %>% write.csv(print(paste(figPath,"rank5.csv", sep = "/")), row.names = F)
#write.csv(g.1,print(paste(figPath,"rank5.csv", sep = "/")), row.names = F)
g.1.kmer <- mm.one[,c("kmer","SS")]
g.1.kmer %>% arrange(SS) %>% write.csv(print(paste(figPath,"rank5.kmer.csv", sep = "/")), row.names = F)
g.1.kmer2 <- g.1.kmer %>% arrange(SS)
g.1.kmer2 <- g.1.kmer2[,c("kmer")]
write.csv(g.1.kmer2,print(paste(figPath,"rank5.kmer2.csv", sep = "/")), row.names = F)
#write.csv(g.1.kmer,print(paste(figPath,"rank5.kmer.csv", sep = "/")), row.names = F)










# DATA FRAME: Ionic Currents ----------------------------------------------

##--RANK 1--##

if (T) {
  
  ##--BASELINE CONDITIONS--##
  CL.mm <- as.vector(outer(CellLine, c("mm.Direct","psi"), paste, sep="."))
  CL.mm.sum <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
  mm <- ranking.1[,c("ACP","Annotation", "chr", "position", "kmer",CL.mm)]
  mm$filter <- 0
  
  min.mm <- 30
  
  #-ONE-#
  mm$filter[rowSums(mm[, CL.mm.sum] >= min.mm) > 0] <- 1
  mm.one <- mm[which(mm$filter == 1),]
  
  if (T) {
    mm.one$cell <- 0
    mm.one$CP <- paste(mm.one$chr,mm.one$position, sep =":")
    for (i in 1:nrow(mm.one)) {
      if (mm.one$A549.psi[i] > 0){mm.one$cell[i] = 1} 
      if (mm.one$HeLa.psi[i] > 0) {mm.one$cell[i] = 2} 
      if (mm.one$HepG2.psi[i] > 0) {mm.one$cell[i] = 3}
      if (mm.one$Jurkat.psi[i] > 0) {mm.one$cell[i] = 4}
      if (mm.one$NTERA.psi[i] > 0) {mm.one$cell[i] = 5} 
      if (mm.one$SHSY5Y.psi[i] > 0) {mm.one$cell[i] = 6}
    }
    
    
    A549.1 <- mm.one[which(mm.one$A549.psi == 1 & mm.one$A549.mm.Direct >= min.mm),]
    HeLa.1 <- mm.one[which(mm.one$HeLa.psi == 1 & mm.one$HeLa.mm.Direct >= min.mm),]
    HepG2.1 <- mm.one[which(mm.one$HepG2.psi == 1 & mm.one$HepG2.mm.Direct >= min.mm),]
    Jurkat.1 <- mm.one[which(mm.one$Jurkat.psi == 1 & mm.one$Jurkat.mm.Direct >= min.mm),]
    NTERA.1 <- mm.one[which(mm.one$NTERA.psi == 1 & mm.one$NTERA.mm.Direct >= min.mm),]
    SHSY5Y.1 <- mm.one[which(mm.one$SHSY5Y.psi == 1 & mm.one$SHSY5Y.mm.Direct >= min.mm),]
    
    # A549.1 <- mm.one[which(mm.one$A549.psi == 1),]
    # HeLa.1 <- mm.one[which(mm.one$HeLa.psi == 1),]
    # HepG2.1 <- mm.one[which(mm.one$HepG2.psi == 1),]
    # Jurkat.1 <- mm.one[which(mm.one$Jurkat.psi == 1),]
    # NTERA.1 <- mm.one[which(mm.one$NTERA.psi == 1),]
    # SHSY5Y.1 <- mm.one[which(mm.one$SHSY5Y.psi == 1),]
    
    
    print(nrow(A549.1))
    print(nrow(HeLa.1))
    print(nrow(HepG2.1))
    print(nrow(Jurkat.1))
    print(nrow(NTERA.1))
    print(nrow(SHSY5Y.1))
    
    
    r1.df <- rbind(data.frame(A549.1), 
                   data.frame(HeLa.1), 
                   data.frame(HepG2.1),
                   data.frame(Jurkat.1),
                   data.frame(NTERA.1),
                   data.frame(SHSY5Y.1))
    
    r1.df$A549.mm.Direct[r1.df$A549.psi == 0] <- 0
    r1.df$HeLa.mm.Direct[r1.df$HeLa.psi == 0] <- 0
    r1.df$HepG2.mm.Direct[r1.df$HepG2.psi == 0] <- 0
    r1.df$Jurkat.mm.Direct[r1.df$Jurkat.psi == 0] <- 0
    r1.df$NTERA.mm.Direct[r1.df$NTERA.psi == 0] <- 0
    r1.df$SHSY5Y.mm.Direct[r1.df$SHSY5Y.psi == 0] <- 0
    
    
    r1.df$CP <- paste(r1.df$chr, r1.df$position, sep = ":")
    r1.df$kmer.start <- r1.df$position - 5
    r1.df$kmer.end <- r1.df$position + 5
    r1.df$range <- paste(r1.df$kmer.start, r1.df$kmer.end, sep = "-")
    r1.df$kmer.11 <- paste(r1.df$chr,r1.df$range , sep = ":")
    
    r1.df$color <- "other"
    r1.df$color[r1.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
    r1.df$color[r1.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
    
  }
}
##--CREATE GENE LISTS--##

r1.annot.df <- read_table(paste(kdePath, "rank1-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r1.annot.df) <- c("Annotation")
r1.annot <- r1.annot.df$Annotation


##--GATHER IONIC CURRENTS FOR POI--##

DRS.type <- c("Direct", "IVT")
CellLine.samp <- c("SHSY5Y") 
r1.annot.samp <- c("SNAP23")


drs.list.Direct <- list()

for (drs in DRS.type) {
  
  CL.DRS.list <- list()
  
  for (CL in CellLine) {
    rep.name <- paste(CL, drs, "merged", sep = "_")
    
    annot.list <- list()
    
    for (annot in r1.annot.samp) {
      print(paste("Start:", annot, drs))
      
      in.folder <- paste(rep.name, annot, "hg38v10.txt", sep = ".")
      in.tsv <- paste(kdePath, "KDE-EA",in.folder, sep="/")
      
      EAC <-  read.table(in.tsv, header = TRUE)
      
      POI.start <- r1.df$position[which(r1.df$Annotation == annot)] - 12
      POI.end <- r1.df$position[which(r1.df$Annotation == annot)] + 2
      
      chr <- r1.df$chr[which(r1.df$Annotation == annot)]
      
      range <- POI.start:POI.end
      POI.range <- as.vector(outer(chr, range, paste, sep=":"))

      poi.list <- list()

      for (poiKMER in POI.range) {

        pos <- which(POI.range == poiKMER)
        EAC.tmp <- EAC[which(EAC$position == range[pos]),]

        poi.df <- EAC.tmp[,c("contig", "position", "reference_kmer", "event_level_mean", "event_stdv",
                             "model_kmer", "model_mean", "model_stdv")]

        poi.list[[poiKMER]] <- poi.df

      }

      annot.list[[annot]] <- poi.list
      
    }
    
    CL.DRS.list[[CL]] <- annot.list
    
  }
  
  drs.list.Direct[[drs]] <- CL.DRS.list
  
}

master.rds.file <- paste("rank1.IonicCurrent.RData", sep = "")
saveRDS(drs.list.Direct, print(paste(rdsPath, master.rds.file, sep = "/")))




##--RANK 5--##

if (T) {
  
  ranking.5$color <- "other"
  ranking.5$color[ranking.5$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  ranking.5$color[ranking.5$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  ##--BASELINE CONDITIONS--##
  CL.mm <- as.vector(outer(CellLine, c("mm.Direct","psi"), paste, sep="."))
  CL.mm.sum <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
  mm <- ranking.5[,c("ACP","Annotation", "chr", "position", "kmer", "color",CL.mm)]
  mm$filter <- 0
  
  min.mm <- 30
  
  #-ONE-#
  mm$filter[rowSums(mm[, CL.mm.sum] >= min.mm) > 0] <- 1
  mm.one <- mm[which(mm$filter == 1),]
  
  
  r5.df <- mm.one
  r5.df$A549.mm.Direct[r5.df$A549.psi == 0] <- 0
  r5.df$HeLa.mm.Direct[r5.df$HeLa.psi == 0] <- 0
  r5.df$HepG2.mm.Direct[r5.df$HepG2.psi == 0] <- 0
  r5.df$Jurkat.mm.Direct[r5.df$Jurkat.psi == 0] <- 0
  r5.df$NTERA.mm.Direct[r5.df$NTERA.psi == 0] <- 0
  r5.df$SHSY5Y.mm.Direct[r5.df$SHSY5Y.psi == 0] <- 0
  
  
  
  ##--IGV--##
  mm.five <- mm
  
  r5.out <- mm.five
  r5.out$CP <- paste(r5.out$chr, r5.out$position, sep = ":")
  r5.out$kmer.start <- r5.out$position - 5
  r5.out$kmer.end <- r5.out$position + 5
  r5.out$range <- paste(r5.out$kmer.start, r5.out$kmer.end, sep = "-")
  r5.out$kmer.11 <- paste(r5.out$chr,r5.out$range , sep = ":")
  
  r5.2 <- r5.out[which(r5.out$Annotation != "RP11-286N22.8" &
                         r5.out$Annotation != "RP11-512M8.12" &
                         r5.out$Annotation != "RP11-512M8.13" &
                         r5.out$Annotation != "RP11-641J8.4" &
                         r5.out$Annotation != "RP11-697E2.12"& 
                         r5.out$Annotation != "RP11-697E2.6" ), ]
  
  r5.2 <- r5.2[which(r5.2$A549.mm.Direct >= 30 |
                       r5.2$HeLa.mm.Direct >= 30 |
                       r5.2$HepG2.mm.Direct >= 30 |
                       r5.2$Jurkat.mm.Direct >= 30 |
                       r5.2$NTERA.mm.Direct >= 30 |
                       r5.2$SHSY5Y.mm.Direct >= 30),]
  r5.df <- r5.2
}



##--CREATE GENE LISTS--##

r5.annot.df <- read_table(paste(kdePath, "rank5-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r5.annot.df) <- c("Annotation")
r5.annot <- r5.annot.df$Annotation


##--GATHER IONIC CURRENTS FOR POI--##

DRS.type <- c("Direct", "IVT")
CellLine.samp <- c("NTERA") 
r5.annot.samp <- c("CDC42")


drs.list.Direct <- list()

for (drs in DRS.type) {
  
  CL.DRS.list <- list()
  
  for (CL in CellLine) {
    rep.name <- paste(CL, drs, "merged", sep = "_")
    
    annot.list <- list()
    
    for (annot in r5.annot) {
      print(paste("Start:", annot, drs))
      
      in.folder <- paste(rep.name, annot, "hg38v10.txt", sep = ".")
      in.tsv <- paste(kdePath, "KDE-EA",in.folder, sep="/")
      
      EAC <-  read.table(in.tsv, header = TRUE)
      
      POI.start <- r5.df$position[which(r5.df$Annotation == annot)] - 6
      POI.end <- r5.df$position[which(r5.df$Annotation == annot)]
      
      chr <- r5.df$chr[which(r5.df$Annotation == annot)]
      
      range <- POI.start:POI.end
      POI.range <- as.vector(outer(chr, range, paste, sep=":"))
      
      poi.list <- list()
      
      for (poiKMER in POI.range) {
        
        pos <- which(POI.range == poiKMER)
        EAC.tmp <- EAC[which(EAC$position == range[pos]),]
        
        poi.df <- EAC.tmp[,c("contig", "position", "reference_kmer", "event_level_mean", "event_stdv",
                             "model_kmer", "model_mean", "model_stdv")]
        
        poi.list[[poiKMER]] <- poi.df
        
      }
      
      annot.list[[annot]] <- poi.list
      
    }
    
    CL.DRS.list[[CL]] <- annot.list
    
  }
  
  drs.list.Direct[[drs]] <- CL.DRS.list
  
}

master.rds.file <- paste("rank5.IonicCurrent.RData", sep = "")
saveRDS(drs.list.Direct, print(paste(rdsPath, master.rds.file, sep = "/")))





# DATA FRAME: rank3 Ionic Current -----------------------------------------
if (T) {
  ##--BASELINE CONDITIONS--##
  CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
  mm <- r.ortho[which(r.ortho$ortho > 0),]
  mm <- three
  
  mm$filter <- 0
  
  min.mm <- 20
  
  #-ONE-#
  mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
  mm.one <- mm[which(mm$filter == 1),]
  
  cat("A549: ", length(which(mm.one$A549.psi == 1)), "\n",
      "HeLa: ", length(which(mm.one$HeLa.psi == 1)), "\n",
      "HepG2: ", length(which(mm.one$HepG2.psi == 1)), "\n",
      "Jurkat: ", length(which(mm.one$Jurkat.psi == 1)), "\n",
      "NTERA: ", length(which(mm.one$NTERA.psi == 1)), "\n",
      "SH-SY5Y: ", length(which(mm.one$SHSY5Y.psi == 1)), "\n", sep = "")
  
  mm.one$A549.mm.Direct[mm.one$A549.psi == 0] <- 0
  mm.one$HeLa.mm.Direct[mm.one$HeLa.psi == 0] <- 0
  mm.one$HepG2.mm.Direct[mm.one$HepG2.psi == 0] <- 0
  mm.one$Jurkat.mm.Direct[mm.one$Jurkat.psi == 0] <- 0
  mm.one$NTERA.mm.Direct[mm.one$NTERA.psi == 0] <- 0
  mm.one$SHSY5Y.mm.Direct[mm.one$SHSY5Y.psi == 0] <- 0

  ##--IGV--##
  mm.three <- mm.one
  
  r3.out <- mm.three
  r3.out$CP <- paste(r3.out$chr, r3.out$position, sep = ":")
  r3.out$kmer.start <- r3.out$position - 5
  r3.out$kmer.end <- r3.out$position + 5
  r3.out$range <- paste(r3.out$kmer.start, r3.out$kmer.end, sep = "-")
  r3.out$kmer.11 <- paste(r3.out$chr,r3.out$range , sep = ":")

  r3.df <- r3.out
  # kmer <- r3.out$kmer.11
  # #kmer <- toupper(kmer)
  # write.table(kmer,print(paste(kdePath,"rank3-KDErange.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)
  # 
  # annot <- r3.out$Annotation
  # #kmer <- toupper(kmer)
  # write.table(annot,print(paste(kdePath,"rank3-KDEannot.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)

}

##--CREATE GENE LISTS--##

r3.annot.df <- read_table(paste(kdePath, "rank3-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r3.annot.df) <- c("Annotation")
r3.annot <- r3.annot.df$Annotation


##--GATHER IONIC CURRENTS FOR POI--##

DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA")
# r3.annot.samp <- c("CDC42")


drs.list.Direct <- list()

for (drs in DRS.type) {
  
  CL.DRS.list <- list()
  
  for (CL in CellLine) {
    rep.name <- paste(CL, drs, "merged", sep = "_")
    
    annot.list <- list()
    
    for (annot in r3.annot) {
      print(paste("Start:", annot, drs))
      
      in.folder <- paste(rep.name, annot, "hg38v10.txt", sep = ".")
      in.tsv <- paste(kdePath, "KDE-EA","rank3", in.folder, sep="/")
      
      EAC <-  read.table(in.tsv, header = TRUE)
      
      POI.start <- r3.df$position[which(r3.df$Annotation == annot)] - 6
      POI.end <- r3.df$position[which(r3.df$Annotation == annot)]
      
      chr <- r3.df$chr[which(r3.df$Annotation == annot)]
      
      range <- POI.start:POI.end
      POI.range <- as.vector(outer(chr, range, paste, sep=":"))
      
      poi.list <- list()
      
      for (poiKMER in POI.range) {
        
        pos <- which(POI.range == poiKMER)
        EAC.tmp <- EAC[which(EAC$position == range[pos]),]
        
        poi.df <- EAC.tmp[,c("contig", "position", "reference_kmer", "event_level_mean", "event_stdv",
                             "model_kmer", "model_mean", "model_stdv")]
        
        poi.list[[poiKMER]] <- poi.df
        
      }
      
      annot.list[[annot]] <- poi.list
      
    }
    
    CL.DRS.list[[CL]] <- annot.list
    
  }
  
  drs.list.Direct[[drs]] <- CL.DRS.list
  
}

master.rds.file <- paste("rank3.IonicCurrent.RData", sep = "")
saveRDS(drs.list.Direct, print(paste(rdsPath, master.rds.file, sep = "/")))




# DATA FRAME: rank1 Ionic Currents - ortho --------------------------------


if (T) {
  ##--BASELINE CONDITIONS--##
  CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
  mm <- r.ortho[which(r.ortho$ortho > 0),]
  mm <- one
  
  mm$filter <- 0
  
  min.mm <- 20
  
  #-ONE-#
  mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
  mm.one <- mm[which(mm$filter == 1),]
  
  ##--IGV--##
  mm.three <- mm.one
  
  r5.out <- mm.three
  r5.out$CP <- paste(r5.out$chr, r5.out$position, sep = ":")
  r5.out$kmer.start <- r5.out$position - 5
  r5.out$kmer.end <- r5.out$position + 5
  r5.out$range <- paste(r5.out$kmer.start, r5.out$kmer.end, sep = "-")
  r5.out$kmer.11 <- paste(r5.out$chr,r5.out$range , sep = ":")
  
  # kmer <- r5.out$kmer.11
  # #kmer <- toupper(kmer)
  # write.table(kmer,print(paste(kdePath,"rank1-ortho-KDErange.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)
  # 
  # annot <- r5.out$Annotation
  # #kmer <- toupper(kmer)
  # write.table(annot,print(paste(kdePath,"rank1-ortho-KDEannot.txt", sep = "/")), row.names = F,quote = FALSE,col.names=FALSE)
  # 
}

##--CREATE GENE LISTS--##

r1.ortho.annot.df <- read_table(paste(kdePath, "rank1-ortho-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r1.ortho.annot.df) <- c("Annotation")
r1.ortho.annot <- r1.ortho.annot.df$Annotation


##--BASELINE CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- r.ortho[which(r.ortho$ortho > 0),]
mm <- one

mm$filter <- 0
min.mm <- 20

#-ONE-#
mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]
r1.ortho.df <- mm.one


##--GATHER IONIC CURRENTS FOR POI--##

DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA")
# r1.ortho.annot.samp <- c("CDC42")


drs.list.Direct <- list()

for (drs in DRS.type) {
  
  CL.DRS.list <- list()
  
  for (CL in CellLine) {
    rep.name <- paste(CL, drs, "merged", sep = "_")
    
    annot.list <- list()
    
    for (annot in r1.ortho.annot) {
      print(paste("Start:", annot, drs))
      
      in.folder <- paste(rep.name, annot, "hg38v10.txt", sep = ".")
      in.tsv <- paste(kdePath, "KDE-EA","rank1-ortho", in.folder, sep="/")
      
      EAC <-  read.table(in.tsv, header = TRUE)
      
      POI.start <- r1.ortho.df$position[which(r1.ortho.df$Annotation == annot)] - 6
      POI.end <- r1.ortho.df$position[which(r1.ortho.df$Annotation == annot)]
      
      chr <- r1.ortho.df$chr[which(r1.ortho.df$Annotation == annot)]
      
      range <- POI.start:POI.end
      POI.range <- as.vector(outer(chr, range, paste, sep=":"))
      
      poi.list <- list()
      
      for (poiKMER in POI.range) {
        
        pos <- which(POI.range == poiKMER)
        EAC.tmp <- EAC[which(EAC$position == range[pos]),]
        
        poi.df <- EAC.tmp[,c("contig", "position", "reference_kmer", "event_level_mean", "event_stdv",
                             "model_kmer", "model_mean", "model_stdv")]
        
        poi.list[[poiKMER]] <- poi.df
        
      }
      
      annot.list[[annot]] <- poi.list
      
    }
    
    CL.DRS.list[[CL]] <- annot.list
    
  }
  
  drs.list.Direct[[drs]] <- CL.DRS.list
  
}

master.rds.file <- paste("rank1-ortho.IonicCurrent.RData", sep = "")
saveRDS(drs.list.Direct, print(paste(rdsPath, master.rds.file, sep = "/")))




# FIGURE: rank1 - Ionic Currents ------------------------------------------

##--INDIVIDUAL CELL LINES--##

r1.annot.df <- read_table(paste(kdePath, "rank1-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r1.annot.df) <- c("Annotation")
r1.annot <- r1.annot.df$Annotation

rank1.rds.file <- paste("rank1.IonicCurrent.RData", sep = "")
drs.list.Direct <- readRDS(print(paste(rdsPath, rank1.rds.file, sep = "/")))


DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA") 
# r5.annot.samp <- c("CDC42")


for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r1.annot) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          ylim(0, 1) +
          xlim(0,200) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.75, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
        
    }
  }
}

##--SANGER CONFIRMED--##

##--SET--##

if (T) {
  
  coY.min <- 0
  coY.max <- .35
  
  coX.min <- 65
  coX.max <- 90
  
  annot <- "SET"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SET)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SET)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SET)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SET)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SET)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SET)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SET), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SET), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SET), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SET), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SET), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SET), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "#FF924C", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--MINPP1--##

if (T) {
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 65
  coX.max <- 110
  
  annot <- "MINPP1"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$MINPP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$MINPP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$MINPP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$MINPP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$MINPP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$MINPP1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$MINPP1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$MINPP1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$MINPP1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$MINPP1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$MINPP1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$MINPP1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "#FF595E", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--NUTF2--##

if (T) {
  
  coY.min <- 0
  coY.max <- .25
  
  coX.min <- 65
  coX.max <- 95
  
  annot <- "NUTF2"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$NUTF2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$NUTF2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$NUTF2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$NUTF2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$NUTF2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$NUTF2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$NUTF2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$NUTF2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$NUTF2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$NUTF2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$NUTF2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$NUTF2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "#FFCA3A", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--CENPA--##

if (T) {
  
  coY.min <- 0
  coY.max <- .25
  
  coX.min <- 105
  coX.max <- 130
  
  annot <- "CENPA"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$CENPA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$CENPA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$CENPA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$CENPA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$CENPA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$CENPA)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$CENPA), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$CENPA), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$CENPA), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$CENPA), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$CENPA), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$CENPA), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "#FFCA3A", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--TLE5--##

if (T) {
  
  coY.min <- 0
  coY.max <- .25
  
  coX.min <- 60
  coX.max <- 90
  
  annot <- "TLE5"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$TLE5)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$TLE5)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$TLE5)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$TLE5)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$TLE5)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$TLE5)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$TLE5), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$TLE5), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$TLE5), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$TLE5), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$TLE5), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$TLE5), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "#8AC926", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--BCCIP--##

if (T) {
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 65
  coX.max <- 125
  
  annot <- "BCCIP"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$BCCIP)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$BCCIP), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$BCCIP), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$BCCIP), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$BCCIP), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$BCCIP), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$BCCIP), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "#8AC926", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--STRN4--##

if (T) {
  
  coY.min <- 0
  coY.max <- .25
  
  coX.min <- 60
  coX.max <- 90
  
  annot <- "STRN4"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$STRN4)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$STRN4), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$STRN4), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$STRN4), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$STRN4), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$STRN4), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$STRN4), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "#8AC926", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--SNAP23--##


if (T) {
  
  coY.min <- 0
  coY.max <- .15
  
  coX.min <- 0
  coX.max <- 200
  
  annot <- "SNAP23"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SNAP23)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SNAP23), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SNAP23), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SNAP23), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SNAP23), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SNAP23), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SNAP23), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(annot, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}





##--YET TO BE SANGERED--##

##--BCCIP--##
if (T) {
  
  annot <- "BCCIP"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$BCCIP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$BCCIP)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$BCCIP), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$BCCIP), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$BCCIP), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$BCCIP), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$BCCIP), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$BCCIP), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      geom_density(fill = "transparent") +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .5) +
      xlim(50,150) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--STRN4--##
if (T) {
  
  annot <- "STRN4"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$STRN4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$STRN4)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$STRN4), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$STRN4), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$STRN4), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$STRN4), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$STRN4), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$STRN4), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      geom_density(fill = "transparent", linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .5) +
      xlim(50,150) +
      coord_cartesian(xlim = c(50,80)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--SEC23IP--##

if (T) {
  
  annot <- "SEC23IP"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SEC23IP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SEC23IP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SEC23IP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SEC23IP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SEC23IP)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SEC23IP)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SEC23IP), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SEC23IP), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SEC23IP), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SEC23IP), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SEC23IP), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SEC23IP), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "#FF595E", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .15) +
      xlim(0,200) +
      coord_cartesian(xlim = c(80,150)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--ATAD2--##

if (T) {
  
  annot <- "ATAD2"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$ATAD2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$ATAD2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$ATAD2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$ATAD2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$ATAD2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$ATAD2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$ATAD2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$ATAD2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$ATAD2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$ATAD2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$ATAD2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$ATAD2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "#FF595E", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .2) +
      xlim(0,200) +
      coord_cartesian(xlim = c(65,110)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--KLF16--##

if (T) {
  
  annot <- "KLF16"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$KLF16)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$KLF16)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$KLF16)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$KLF16)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$KLF16)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$KLF16)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$KLF16), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$KLF16), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$KLF16), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$KLF16), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$KLF16), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$KLF16), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "#FF924C", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .25) +
      xlim(0,200) +
      coord_cartesian(xlim = c(70,120)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--SMAD4--##

if (T) {
  
  annot <- "SMAD4"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SMAD4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SMAD4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SMAD4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SMAD4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SMAD4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SMAD4)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SMAD4), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SMAD4), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SMAD4), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SMAD4), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SMAD4), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SMAD4), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "#FFCA3A", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .15) +
      xlim(0,200) +
      coord_cartesian(xlim = c(60,130)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--PTMS--##

if (T) {
  
  annot <- "PTMS"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$PTMS)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$PTMS)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$PTMS)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$PTMS)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$PTMS)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$PTMS)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$PTMS), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$PTMS), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$PTMS), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$PTMS), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$PTMS), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$PTMS), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "#8AC926", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .25) +
      xlim(0,200) +
      coord_cartesian(xlim = c(60,90)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--CALCOCO2--##

if (T) {
  
  annot <- "CALCOCO2"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$CALCOCO2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$CALCOCO2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$CALCOCO2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$CALCOCO2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$CALCOCO2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$CALCOCO2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$CALCOCO2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$CALCOCO2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$CALCOCO2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$CALCOCO2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$CALCOCO2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$CALCOCO2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "#1982C4", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .15) +
      xlim(0,200) +
      coord_cartesian(xlim = c(75,140)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--PDLIM1--##

if (T) {
  
  annot <- "PDLIM1"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$PDLIM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$PDLIM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$PDLIM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$PDLIM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$PDLIM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$PDLIM1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$PDLIM1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$PDLIM1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$PDLIM1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$PDLIM1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$PDLIM1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$PDLIM1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .15) +
      xlim(0,200) +
      coord_cartesian(xlim = c(70,110)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--KPNA4--##

if (T) {
  
  annot <- "KPNA4"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$KPNA4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$KPNA4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$KPNA4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$KPNA4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$KPNA4)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$KPNA4)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$KPNA4), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$KPNA4), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$KPNA4), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$KPNA4), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$KPNA4), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$KPNA4), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .3) +
      xlim(0,200) +
      coord_cartesian(xlim = c(70,100)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--PAPOLA--##

if (T) {
  
  annot <- "PAPOLA"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$PAPOLA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$PAPOLA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$PAPOLA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$PAPOLA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$PAPOLA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$PAPOLA)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$PAPOLA), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$PAPOLA), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$PAPOLA), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$PAPOLA), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$PAPOLA), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$PAPOLA), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .25) +
      xlim(0,200) +
      coord_cartesian(xlim = c(70,100)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--SNAP23--##

if (T) {
  
  annot <- "SNAP23"
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SNAP23)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SNAP23)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SNAP23), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SNAP23), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SNAP23), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SNAP23), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SNAP23), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SNAP23), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, .15) +
      xlim(0,200) +
      coord_cartesian(xlim = c(65,115)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}













##--MANUAL--##

ggplot() +
  geom_density(aes(x = drs.list.Direct$Direct$A549$MINPP1),fill = "blue", color = "black", alpha = 0.7) +
  geom_density(aes(x = drs.list.Direct$IVT$A549$MINPP1),fill = "red", color = "black", alpha = 0.7) +
  ggtitle("Kernel Density Estimate: MINPP1") +
  xlab("Current (pA)") +
  geom_vline(xintercept=89.28, color="black", linewidth=0.75, linetype=2)+
  theme_minimal()


# FIGURE: rank5 - Ionic Currents ------------------------------------------


##--INDIVIDUAL CELL LINES--##

r5.annot.df <- read_table(paste(kdePath, "rank5-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r5.annot.df) <- c("Annotation")
r5.annot <- r5.annot.df$Annotation

rank5.rds.file <- paste("rank5.IonicCurrent.RData", sep = "")
drs.list.Direct <- readRDS(print(paste(rdsPath, rank5.rds.file, sep = "/")))


DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA") 
# r5.annot.samp <- c("CDC42")


for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r5.annot) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          ylim(0, 1) +
          xlim(0,200) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.75, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
      
    }
  }
}


DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA") 
r5.annot.samp <- c("SEC23IP")



for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r5.annot.samp) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          ylim(0, .15) +
          xlim(0,200) +
          coord_cartesian(xlim = c(80,150)) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
      
    }
  }
}

###




##--WDR45--##

if (T) {
  
  annot <- "WDR45"
  
  coY.min <- 0
  coY.max <- .075
  
  coX.min <- 70
  coX.max <- 140
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$WDR45)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$WDR45)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$WDR45)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$WDR45)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$WDR45)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$WDR45)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$WDR45), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$WDR45), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$WDR45), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$WDR45), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$WDR45), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$WDR45), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "#FFCA3A", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--YWHAQ--##

if (T) {
  
  annot <- "YWHAQ"
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 60
  coX.max <- 90
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$YWHAQ)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$YWHAQ)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$YWHAQ)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$YWHAQ)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$YWHAQ)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$YWHAQ)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$YWHAQ), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$YWHAQ), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$YWHAQ), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$YWHAQ), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$YWHAQ), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$YWHAQ), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "#FFCA3A", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--CDC42--##

if (T) {
  
  annot <- "CDC42"
  
  coY.min <- 0
  coY.max <- .075
  
  coX.min <- 70
  coX.max <- 145
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$CDC42)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$CDC42)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$CDC42)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$CDC42)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$CDC42)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$CDC42)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$CDC42), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$CDC42), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$CDC42), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$CDC42), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$CDC42), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$CDC42), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "#8AC926", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--LAPTM4B--##

if (T) {
  
  annot <- "LAPTM4B"
  
  coY.min <- 0
  coY.max <- .25
  
  coX.min <- 55
  coX.max <- 90
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$LAPTM4B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$LAPTM4B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$LAPTM4B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$LAPTM4B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$LAPTM4B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$LAPTM4B)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$LAPTM4B), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$LAPTM4B), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$LAPTM4B), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$LAPTM4B), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$LAPTM4B), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$LAPTM4B), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "#1982C4", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--GRB2--##

if (T) {
  
  annot <- "GRB2"
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 55
  coX.max <- 100
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$GRB2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$GRB2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$GRB2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$GRB2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$GRB2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$GRB2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$GRB2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$GRB2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$GRB2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$GRB2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$GRB2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$GRB2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--M6PR--##

if (T) {
  
  annot <- "M6PR"
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 60
  coX.max <- 105
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$M6PR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$M6PR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$M6PR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$M6PR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$M6PR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$M6PR)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$M6PR), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$M6PR), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$M6PR), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$M6PR), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$M6PR), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$M6PR), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--DAZAP1--##

if (T) {
  
  annot <- "DAZAP1"
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 60
  coX.max <- 105
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$DAZAP1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$DAZAP1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$DAZAP1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$DAZAP1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$DAZAP1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$DAZAP1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$DAZAP1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--H2AZ2--##

if (T) {
  
  annot <- "H2AZ2"
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 65
  coX.max <- 110
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$H2AZ2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$H2AZ2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$H2AZ2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$H2AZ2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$H2AZ2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$H2AZ2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$H2AZ2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$H2AZ2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$H2AZ2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$H2AZ2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$H2AZ2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$H2AZ2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "#6A4C93",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank5",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}











###








##--MANUAL--##

ggplot() +
  geom_density(aes(x = drs.list.Direct$Direct$A549$MINPP1),fill = "blue", color = "black", alpha = 0.7) +
  geom_density(aes(x = drs.list.Direct$IVT$A549$MINPP1),fill = "red", color = "black", alpha = 0.7) +
  ggtitle("Kernel Density Estimate: MINPP1") +
  xlab("Current (pA)") +
  geom_vline(xintercept=89.28, color="black", linewidth=0.75, linetype=2)+
  theme_minimal()







# FIGURE: rank3 - Ionic Currents ------------------------------------------

##--INDIVIDUAL CELL LINES--##

r3.annot.df <- read_table(paste(kdePath, "rank3-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r3.annot.df) <- c("Annotation")
r3.annot <- r3.annot.df$Annotation

rank3.rds.file <- paste("rank3.IonicCurrent.RData", sep = "")
drs.list.Direct <- readRDS(print(paste(rdsPath, rank3.rds.file, sep = "/")))

DRS.type <- c("Direct", "IVT")

for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r3.annot) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          ylim(0, 1) +
          xlim(0,200) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.75, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
      
    }
  }
}


DRS.type <- c("Direct", "IVT")
# CellLine.samp <- c("NTERA") 
r3.annot.samp <- c("SEC23IP")



for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r3.annot.samp) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          ylim(0, .15) +
          xlim(0,200) +
          coord_cartesian(xlim = c(80,150)) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
      
    }
  }
}

###




##--ZC3H7A--##

if (T) {
  
  annot <- "ZC3H7A"
  
  coY.min <- 0
  coY.max <- 0.25
  
  coX.min <- 60
  coX.max <- 110
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$ZC3H7A)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$ZC3H7A)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$ZC3H7A)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$ZC3H7A)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$ZC3H7A)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$ZC3H7A)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$ZC3H7A), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$ZC3H7A), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$ZC3H7A), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$ZC3H7A), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$ZC3H7A), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$ZC3H7A), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--DUS1L--##

if (T) {
  
  annot <- "DUS1L"
  
  coY.min <- 0
  coY.max <- 0.30
  
  coX.min <- 60
  coX.max <- 110
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$DUS1L)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$DUS1L)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$DUS1L)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$DUS1L)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$DUS1L)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$DUS1L)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$DUS1L), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$DUS1L), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$DUS1L), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$DUS1L), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$DUS1L), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$DUS1L), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--DAZAP1--##

if (T) {
  
  annot <- "DAZAP1"
  
  coY.min <- 0
  coY.max <- 0.20
  
  coX.min <- 55
  coX.max <- 110
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$DAZAP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$DAZAP1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$DAZAP1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$DAZAP1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$DAZAP1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$DAZAP1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$DAZAP1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$DAZAP1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--TMEM14B--##

if (T) {
  
  annot <- "TMEM14B"
  
  
  coY.min <- 0
  coY.max <- 0.25
  
  coX.min <- 60
  coX.max <- 110
  
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$TMEM14B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$TMEM14B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$TMEM14B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$TMEM14B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$TMEM14B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$TMEM14B)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$TMEM14B), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$TMEM14B), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$TMEM14B), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$TMEM14B), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$TMEM14B), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$TMEM14B), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--TMEM126B--##

if (T) {
  
  annot <- "TMEM126B"
  
  
  coY.min <- 0
  coY.max <- 0.12
  
  coX.min <- 73
  coX.max <- 112
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$TMEM126B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$TMEM126B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$TMEM126B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$TMEM126B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$TMEM126B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$TMEM126B)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$TMEM126B), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$TMEM126B), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$TMEM126B), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$TMEM126B), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$TMEM126B), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$TMEM126B), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--SKP1--##

if (T) {
  
  annot <- "SKP1"
  
  coY.min <- 0
  coY.max <- 0.23
  
  coX.min <- 70
  coX.max <- 95
  
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SKP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SKP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SKP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SKP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SKP1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SKP1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SKP1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SKP1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SKP1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SKP1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SKP1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SKP1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--RPS29--##

if (T) {
  
  annot <- "RPS29"
  
  coY.min <- 0
  coY.max <- 0.30
  
  coX.min <- 50
  coX.max <- 110
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$RPS29)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$RPS29)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$RPS29)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$RPS29)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$RPS29)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$RPS29)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$RPS29), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$RPS29), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$RPS29), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$RPS29), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$RPS29), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$RPS29), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--RPL15--##

if (T) {
  
  annot <- "RPL15"
  
  
  coY.min <- 0
  coY.max <- 0.075
  
  coX.min <- 60
  coX.max <- 155
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$RPL15)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$RPL15)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$RPL15)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$RPL15)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$RPL15)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$RPL15)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$RPL15), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$RPL15), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$RPL15), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$RPL15), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$RPL15), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$RPL15), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--RBCK1--##

if (T) {
  
  annot <- "RBCK1"
  
  
  coY.min <- 0
  coY.max <- 0.15
  
  coX.min <- 60
  coX.max <- 100
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$RBCK1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$RBCK1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$RBCK1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$RBCK1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$RBCK1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$RBCK1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$RBCK1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$RBCK1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$RBCK1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$RBCK1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$RBCK1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$RBCK1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}



##--PGRMC2--##

if (T) {
  
  annot <- "PGRMC2"
  
  
  coY.min <- 0
  coY.max <- 0.30
  
  coX.min <- 60
  coX.max <- 140
  
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$PGRMC2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$PGRMC2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$PGRMC2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$PGRMC2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$PGRMC2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$PGRMC2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$PGRMC2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$PGRMC2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$PGRMC2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$PGRMC2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$PGRMC2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$PGRMC2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}



##--MYO1B--##

if (T) {
  
  annot <- "MYO1B"
  
  coY.min <- 0
  coY.max <- 0.28
  
  coX.min <- 67
  coX.max <- 97
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$MYO1B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$MYO1B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$MYO1B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$MYO1B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$MYO1B)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$MYO1B)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$MYO1B), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$MYO1B), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$MYO1B), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$MYO1B), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$MYO1B), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$MYO1B), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}



##--ITFG1--##

if (T) {
  
  annot <- "ITFG1"
  
  
  coY.min <- 0
  coY.max <- 0.25
  
  coX.min <- 40
  coX.max <- 130
  
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$ITFG1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$ITFG1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$ITFG1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$ITFG1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$ITFG1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$ITFG1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$ITFG1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$ITFG1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$ITFG1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$ITFG1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$ITFG1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$ITFG1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}



##--ENSA--##

if (T) {
  
  annot <- "ENSA"
  
  
  coY.min <- 0
  coY.max <- 0.1
  
  coX.min <- 70
  coX.max <- 145
  
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$ENSA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$ENSA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$ENSA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$ENSA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$ENSA)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$ENSA)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$ENSA), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$ENSA), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$ENSA), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$ENSA), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$ENSA), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$ENSA), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--EIF3F--##

if (T) {
  
  annot <- "EIF3F"
  
  
  coY.min <- 0
  coY.max <- 0.2
  
  coX.min <- 65
  coX.max <- 110
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$EIF3F)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$EIF3F)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$EIF3F)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$EIF3F)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$EIF3F)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$EIF3F)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$EIF3F), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$EIF3F), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$EIF3F), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$EIF3F), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$EIF3F), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$EIF3F), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--EIF3D--##

if (T) {
  
  annot <- "EIF3D"
  
  coY.min <- 0
  coY.max <- 0.25
  
  coX.min <- 55
  coX.max <- 100
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$EIF3D)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$EIF3D)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$EIF3D)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$EIF3D)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$EIF3D)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$EIF3D)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$EIF3D), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$EIF3D), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$EIF3D), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$EIF3D), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$EIF3D), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$EIF3D), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--CALR--##

if (T) {
  
  annot <- "CALR"
  
  coY.min <- 0
  coY.max <- 0.30
  
  coX.min <- 60
  coX.max <- 110
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$CALR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$CALR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$CALR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$CALR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$CALR)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$CALR)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$CALR), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$CALR), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$CALR), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$CALR), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$CALR), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$CALR), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--CALM1--##

if (T) {
  
  annot <- "CALM1"
  
  coY.min <- 0
  coY.max <- .20
  
  coX.min <- 60
  coX.max <- 105
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$CALM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$CALM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$CALM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$CALM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$CALM1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$CALM1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$CALM1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$CALM1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$CALM1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$CALM1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$CALM1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$CALM1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


##--BCAP31--##

if (T) {
  
  annot <- "BCAP31"
  
  coY.min <- 0
  coY.max <- 0.10
  
  coX.min <- 60
  coX.max <- 92
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$BCAP31)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$BCAP31)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$BCAP31)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$BCAP31)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$BCAP31)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$BCAP31)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$BCAP31), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$BCAP31), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$BCAP31), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$BCAP31), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$BCAP31), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$BCAP31), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank3",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}



# FIGURE: rank1 - Ionic Current - ortho -----------------------------------


##--CREATE GENE LISTS--##

r1.ortho.annot.df <- read_table(paste(kdePath, "rank1-ortho-KDEannot.txt", sep = "/"), col_names = FALSE)
colnames(r1.ortho.annot.df) <- c("Annotation")
r1.ortho.annot <- r1.ortho.annot.df$Annotation

rank1.rds.file <- paste("rank1-ortho.IonicCurrent.RData", sep = "")
drs.list.Direct <- readRDS(print(paste(rdsPath, rank1.rds.file, sep = "/")))



DRS.type <- c("Direct", "IVT")

for (CL in CellLine) {
  
  direct.list <- drs.list.Direct$Direct[[CL]]
  ivt.list <- drs.list.Direct$IVT[[CL]]
  
  for (annot in r1.ortho.annot) {
    
    POI.range <- names(direct.list[[annot]])
    
    direct.annot <- direct.list[[annot]]
    ivt.annot <- ivt.list[[annot]]
    
    
    for (poiKMER in POI.range) {
      
      direct.pos <- direct.annot[[poiKMER]]
      ivt.pos <- ivt.annot[[poiKMER]]
      
      nDirect <- nrow(direct.pos)
      nIVT <- nrow(ivt.pos)
      
      if (nDirect > 2 & nIVT > 2) {
        pos <- direct.pos$position[1]
        
        
        drs.df <- rbind(data.frame("EM" = direct.pos$event_level_mean, "ID"= "Direct"),
                        data.frame("EM" = ivt.pos$event_level_mean, "ID" ="IVT")) 
        PM.mean <- direct.pos$model_mean[1]
        PM.kmer <- direct.pos$model_kmer[1]
        
        out.file <- paste(CL, annot, pos, "KDE-IonicCurrent.eps", sep = ".")
        
        CL.plot.reads <- drs.df %>%
          ggplot(aes(x = EM, color = ID)) +
          geom_density(fill = "transparent", linewidth=0.5) +
          ggtitle(paste(CL, poiKMER, PM.kmer, sep = " ")) +
          scale_color_manual(values = c("Direct" = "red",
                                        "IVT" = "black")) +
          ylim(0, 1) +
          xlim(0,200) +
          xlab("Current (pA)") +
          geom_vline(xintercept=PM.mean, color="black", linewidth=0.75, linetype=2)+
          theme(axis.line=element_line(),
                #axis.text.x=element_blank(),
                #axis.text.y=element_blank(),
                #axis.ticks=element_blank(),
                #axis.title.x=element_blank(),
                #axis.title.y=element_blank(),
                legend.position="none",
                panel.background=element_blank(),
                panel.border=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                plot.background=element_blank())
        CL.plot.reads
        ggsave(print(paste(kdePath,"KDE-plot", "rank1-ortho",annot,out.file, sep = "/")), width = 4, height = 4)
      }
      
      else {print(paste(CL, poiKMER, "doesn't exist", sep = " "))}
      
    }
  }
}



##--IPO7--##

if (T) {
  
  annot <- "IPO7"
  
  coY.min <- 0
  coY.max <- .075
  
  coX.min <- 70
  coX.max <- 150
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$IPO7)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$IPO7)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$IPO7)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$IPO7)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$IPO7)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$IPO7)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$IPO7), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$IPO7), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$IPO7), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$IPO7), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$IPO7), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$IPO7), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1-ortho",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--SLC25A1--##

if (T) {
  
  annot <- "SLC25A1"
  
  coY.min <- 0
  coY.max <- .20
  
  coX.min <- 50
  coX.max <- 100
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$SLC25A1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$SLC25A1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$SLC25A1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$SLC25A1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$SLC25A1)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$SLC25A1)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$SLC25A1), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$SLC25A1), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$SLC25A1), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$SLC25A1), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$SLC25A1), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$SLC25A1), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1-ortho",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--MRPL46--##

if (T) {
  
  annot <- "MRPL46"
  
  coY.min <- 0
  coY.max <- 0.30
  
  coX.min <- 65
  coX.max <- 100
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$MRPL46)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$MRPL46)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$MRPL46)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$MRPL46)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$MRPL46)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$MRPL46)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$MRPL46), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$MRPL46), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$MRPL46), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$MRPL46), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$MRPL46), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$MRPL46), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1-ortho",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}

##--DRG2--##

if (T) {
  
  annot <- "DRG2"
  
  coY.min <- 0
  coY.max <- 0.065
  
  coX.min <- 75
  coX.max <- 140
  
  
  ivt.pan <- rbind(data.frame(do.call(rbind, drs.list.Direct$IVT$A549$DRG2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HeLa$DRG2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$HepG2$DRG2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$Jurkat$DRG2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$NTERA$DRG2)),
                   data.frame(do.call(rbind, drs.list.Direct$IVT$SHSY5Y$DRG2)))
  ivt.pan$ID <- "panIVT"
  
  
  direct.df <-rbind(data.frame(do.call(rbind, drs.list.Direct$Direct$A549$DRG2), "ID" = "A549"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HeLa$DRG2), "ID" = "HeLa"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$HepG2$DRG2), "ID" = "HepG2"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$Jurkat$DRG2), "ID" = "Jurkat"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$NTERA$DRG2), "ID" = "NTERA"),
                    data.frame(do.call(rbind, drs.list.Direct$Direct$SHSY5Y$DRG2), "ID" = "SHSY5Y"))
  
  
  drs.df <- rbind(ivt.pan, direct.df)
  drs.df <- rownames_to_column(drs.df, var = "RowID")
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    out.file <- paste(annot, drs.pos, "overlay.KDE-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill = ID)) +
      scale_color_manual(values = c("A549" = "#FF595E", 
                                    "HeLa" = "#FF924C", 
                                    "HepG2" = "#FFCA3A", 
                                    "Jurkat" = "#8AC926", 
                                    "NTERA" = "#1982C4", 
                                    "SHSY5Y" = "#6A4C93",
                                    "panIVT" = "black")) +
      scale_fill_manual(values = c("A549" = "transparent", 
                                   "HeLa" = "transparent", 
                                   "HepG2" = "transparent", 
                                   "Jurkat" = "transparent", 
                                   "NTERA" = "transparent", 
                                   "SHSY5Y" = "transparent",
                                   "panIVT" = "transparent")) +
      geom_density(linewidth=0.5) +
      labs(title = paste(annot), subtitle = paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
            #axis.text.x=element_blank(),
            #axis.text.y=element_blank(),
            #axis.ticks=element_blank(),
            #axis.title.x=element_blank(),
            #axis.title.y=element_blank(),
            legend.position="none",
            plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8),
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    CL.plot.reads
    ggsave(print(paste(kdePath,"KDE-plot", "rank1-ortho",annot,out.file, sep = "/")), width = 3, height = 3)
    
    
  }
}


# FIGURE: SANGER ----------------------------------------------------------

library(sangerseqR)
sangerPath <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU/DATA/sanger"

dirPath <- print(path)
setwd(dirPath) #set working directory

if (T) {
  sangerIN = paste(sangerPath, "CENPA_HepG2_F_CENPA_FWD__2024-01-29_C10.ab1", sep = "/")
  sangerOUT = paste(sangerPath, "rank1-CENPA.pdf", sep = "/")
  
  ITS<-read.abif(paste(sangerIN))
  ITSseq <- sangerseq(ITS)
  str(ITSseq)
  chromatogram(ITSseq, width = 20, height = 2, trim5=155 , trim3=289 ,showcalls = "primary", filename = paste(sangerOUT))
}


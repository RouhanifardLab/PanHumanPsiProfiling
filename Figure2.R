## ========================================== ##
##             Psi Conservation               ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 11/28/2023
#@ content: Figure 2
#@ version: 4.0

# housekeeping ------------------------------------------------------------
list.of.packages <- c("ggplot2", "ggExtra","dplyr", "tidyr", "plyr", "readr", "ape", "stringr", "tibble", "stringr", "reshape2", "limma", "plotly", "rtracklayer", "VennDiagram", "Rsamtools", "hrbrthemes", "caret", "venneuler")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) # from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them. Thanks!
library(ggplot2)
library(ggExtra)
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(tibble)
library(stringr)
library(reshape2)
library(limma)
library(plotly)
library(rtracklayer)
library(VennDiagram)
library(Rsamtools)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()
library(caret)
library(venneuler)
library(vtable)
library(janitor)
library(gridExtra)
library(viridis)
library(epos)
library(pracma)
library(tidyverse)

###--- Clear Space ---###
rm(list = ls()) #clear global environment except for packages & libraries 
cat("\f") #clears console

###--- General Path & Working Directory ---###
path <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/home/cam/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

path <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

dirPath <- print(path)
setwd(dirPath) #set working directory

dataPath <- print(paste(dirPath, "DATA", sep = "/"))
pvalPath <- print(paste(dataPath, "pValues", sep = "/"))
figPath <- print(paste(dirPath, "Figures", "Figure-2","FIG2-elements", sep = "/"))
bamPath <- print(paste(dirPath, "BAM", "hg38v10", sep = "/"))
tpmPath <- print(paste(bamPath, "stringtie",sep = "/"))
pysamPath <- print(paste(dataPath, "pysamstats", sep = "/"))
pvalSOURCE <- print(paste(dirPath,"pValues","DATA", sep = "/"))
kpPath <- print(paste(dataPath, "seqLOGO", sep = "/"))
bootPath <- print(paste(dirPath,"CODE","bootSTRAP", sep = "/"))
bsPath <- print(paste(bootPath, "pVALUES", sep="/"))
TPMbsPath <- print(paste(bootPath, "stringtie/DATA", sep = "/"))
orthoPATH <- print(paste(dataPath, "orthogonal", sep = "/"))
qpPath <- print(paste(dataPath, "quantProteomics", sep = "/"))

###--- Colors ---###
cp1 <- c("#FF595E", "#FF924C", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93")

CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")

cat("\f") #clears console

# DATA FRAME: pVal MASTER -------------------------------------------------


m.df <- master.df[,c(stable, CL.ACP)]

CL.ACP <- as.vector(outer(CellLine, c("N_reads_Direct","mm.Direct", "p.value.Direct", "psi", "which.IVT"), paste, sep="."))
CL.psi <- as.vector(outer(CellLine, c("psi"), paste, sep="."))

if (T){
  master.df <- read_csv(paste(pvalPath,"sixCellLine.pval-IVT.mmDirect.PSI.master.csv", sep = "/"))

  #-Column Vectors-#
  stable <- c("ACP", "Annotation", "chr", "position", "kmer")
  CL.ACP <- as.vector(outer(CellLine, c("N_reads_Direct", "N_reads_IVT","mm.Direct", "mm.IVT", "p.value.Direct","pval.NA", "psi", "true", "which.IVT"), paste, sep="."))
  CL.psi <- as.vector(outer(CellLine, c("psi"), paste, sep="."))
  
  master.list <- list()
  
  for (CL in CellLine) {
    
    CL.col <- as.vector(outer(CL, c("N_reads_Direct", "N_reads_IVT","mm.Direct", "mm.IVT", "p.value.Direct",
                                    "pval.NA", "psi", "true", "which.IVT"), paste, sep="."))
    CL.tmp <- master.df[,c(stable,CL.col)]
    
    prefix_to_remove <- paste0(CL,".")
    columns_to_exclude <- c(stable)
    names(CL.tmp) <- ifelse(names(CL.tmp) %in% columns_to_exclude, names(CL.tmp), 
                            gsub(paste0("^", prefix_to_remove), "", names(CL.tmp)))
    
    
    CL.tmp$mm.Direct[CL.tmp$mm.Direct < 10] <- 0
    CL.tmp <- CL.tmp[which(CL.tmp$N_reads_Direct >= 30 &
                             CL.tmp$N_reads_IVT >= 10 & 
                             CL.tmp$mm.IVT <= 10),]
  
    CL.tmp$psi <- 0 
    CL.tmp$psi[which(CL.tmp$p.value.Direct < 0.001 & CL.tmp$mm.Direct >= 10)] <- 1
    CL.tmp$pIVT <- 0
    CL.tmp$pIVT[which(CL.tmp$which.IVT == "pIVT")] <- 1
    
    CL.tmp <- subset(CL.tmp, select = -c(which.IVT))
    
    CL.tmp <- CL.tmp %>% rename_with(~ paste0(CL,".", .), -c(ACP, Annotation, chr, position, kmer))
    master.list[[CL]] <- CL.tmp
  }
  

  ACP.master.UF <- master.list %>% reduce(full_join, by=stable)
  ACP.priority <- ACP.master.UF[,c("ACP", CL.psi)]
  ACP.priority[is.na(ACP.priority)] <- 0
  ACP.priority$rank <- rowSums(ACP.priority[CL.psi])
  ACP.priority <- ACP.priority[which(ACP.priority$rank > 0),]
  
  ACP.priority.map <- map(master.list, ~ .x %>% 
                            filter(ACP %in% ACP.priority$ACP))
  
  ACP.master <- ACP.priority.map %>% reduce(full_join, by=stable)
  ACP.master$rank <- rowSums(ACP.master[CL.psi])
  
  
  ACP.master <- ACP.master[complete.cases(ACP.master), ]  
  
  pval.file <- paste("sixCellLine",".pval.filtered.master.csv", sep = "")
  write.csv(ACP.master,print(paste(pvalPath, pval.file, sep = "/")), row.names = F)
}

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
  
  # 
  #   write.csv(ranking.1,print(paste(pvalPath,"ranking1.csv", sep = "/")), row.names = F)
  #   write.csv(ranking.2,print(paste(pvalPath,"ranking2.csv", sep = "/")), row.names = F)
  #   write.csv(ranking.3,print(paste(pvalPath,"ranking3.csv", sep = "/")), row.names = F)
  #   write.csv(ranking.4,print(paste(pvalPath,"ranking4.csv", sep = "/")), row.names = F)
  #   write.csv(ranking.5,print(paste(pvalPath,"ranking5.csv", sep = "/")), row.names = F)
  #   write.csv(ranking.6,print(paste(pvalPath,"ranking6.csv", sep = "/")), row.names = F)
  
  
  
  
  ranking.6$color <- "other"
  ranking.6$color[ranking.6$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  ranking.6$color[ranking.6$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  
  ##--pIVT only--##
  ranking.pIVT <- ranking[which(ranking$A549.pIVT == 1 & 
                                  ranking$HeLa.pIVT == 1 & 
                                  ranking$HepG2.pIVT == 1 & 
                                  ranking$Jurkat.pIVT == 1 & 
                                  ranking$NTERA.pIVT == 1 & 
                                  ranking$SHSY5Y.pIVT == 1),]
  ranking.pIVT.1 <- ranking.pIVT[which(ranking.pIVT$rank == 1),]
  ranking.pIVT.2 <- ranking.pIVT[which(ranking.pIVT$rank == 2),]
  ranking.pIVT.3 <- ranking.pIVT[which(ranking.pIVT$rank == 3),]
  ranking.pIVT.4 <- ranking.pIVT[which(ranking.pIVT$rank == 4),]
  ranking.pIVT.5 <- ranking.pIVT[which(ranking.pIVT$rank == 5),]
  ranking.pIVT.6 <- ranking.pIVT[which(ranking.pIVT$rank == 6),]
  
  r1 <- nrow(ranking.pIVT.1)
  r2 <- nrow(ranking.pIVT.2)
  r3 <- nrow(ranking.pIVT.3)
  r4 <- nrow(ranking.pIVT.4)
  r5 <- nrow(ranking.pIVT.5)
  r6 <- nrow(ranking.pIVT.6)
  
  print(paste("1 Cell Line",r1, sep = " -- "))
  print(paste("2 Cell Line",r2, sep = " -- "))
  print(paste("3 Cell Line",r3, sep = " -- "))
  print(paste("4 Cell Line",r4, sep = " -- "))
  print(paste("5 Cell Line",r5, sep = " -- "))
  print(paste("6 Cell Line",r6, sep = " -- "))
  
  ##--BASELINE CONDITIONS--##
  CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
  mm <- ranking.pIVT.6[,c("ACP","Annotation", "chr", "position", "kmer", CL.mm)]
  mm$filter <- 0
  
  min.mm <- 30
  
  #-ONE-#
  mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
  mm.one <- mm[which(mm$filter == 1),]
  
  #-ALL-#
  mm$filter <- 0
  mm$filter[rowSums(mm[, CL.mm] > min.mm) == length(CL.mm)] <- 1
  mm.all <- mm[which(mm$filter == 1),]
  
  
  ranking$color <- "other"
  ranking$color[ranking$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  ranking$color[ranking$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  mm <- ranking[,c("rank", "ACP","Annotation", "chr", "position", "kmer", "color", CL.mm)]
  
}

# FIGURE: Conserved Ψ HeatMAP ---------------------------------------------

##--EVERYTHING--##
CL.col <- as.vector(outer(CellLine, c("mm.Direct", "N_reads_Direct", "psi"), paste, sep="."))
# mm <- ranking[,c("ACP","rank",CL.mm)]


##--BASELINE CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- ranking.6[,c("ACP","Annotation", "chr", "position", "kmer", "color", CL.col)]
mm$filter <- 0

min.mm <- 30

#-ONE-#
mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]

#-ALL-#
mm$filter <- 0
mm$filter[rowSums(mm[, CL.mm] > min.mm) == length(CL.mm)] <- 1
mm.all <- mm[which(mm$filter == 1),]


##--IGV--##
r6.out <- mm.all
r6.out$CP <- paste(r6.out$chr, r6.out$position, sep = ":")
r6.out$kmer.start <- r6.out$position - 5
r6.out$kmer.end <- r6.out$position + 5
r6.out$range <- paste(r6.out$kmer.start, r6.out$kmer.end, sep = "-")
r6.out$kmer.11 <- paste(r6.out$chr,r6.out$range , sep = ":")
write.csv(r6.out,print(paste(pvalPath,"rank6-CP.csv", sep = "/")), row.names = F)



#------------------------------------------------------------------------------------------#
###---PLOT---###

r6 <- mm.all
r6$SS <- round(apply(r6[,c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")], 1, sd), digits = 2)
#r6$SSrank <- sprintf("%04f", r6$SS)
#r6 <- r6[order(r6$SSrank),]
r6 <- r6[order(r6$SS),]

SSrank <- paste(r6$SS, r6$Annotation, sep = "-")


cp2 <- c("#FFFFFF","#0450B4", "#046DC8", "#1184A7", "#15A2A2", "#6FB1A0", "#B4418E", "#D94A8C", "#EA515F",
         "#FE7434", "#FEA802", "#ffd93d")

breaks <- c(0,9.99, 19.99, 29.99, 39.99, 49.99, 59.99, 69.99, 79.99, 89.99, 99.99, 101)
labels <- c("0-9","10-19","20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90-99", "100")

heatMAP.df <- r6[,c("SS", "A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")]
reps <- c("A549.mm.Direct", "HeLa.mm.Direct", "HepG2.mm.Direct", "Jurkat.mm.Direct", "NTERA.mm.Direct", "SHSY5Y.mm.Direct")
HM_melt <- melt(heatMAP.df, measure.vars = reps)
HM_melt$value[HM_melt$value < 0] <- 0
HM_melt$category <- cut(HM_melt$value, breaks = breaks, labels = labels, right = FALSE)
HM_melt$category <- factor(HM_melt$category, levels = labels)
HM_melt$SS <- factor(HM_melt$SS, levels = rev(sort(unique(HM_melt$SS))))
all.HeatMap <- HM_melt %>% 
  ggplot(aes(x = variable , y = SS, fill = category)) +
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
name <- print(paste("rank6.BinaryHeatMapbinary.eps", sep = ""))
ggsave(print(paste(figPath, name, sep = "/")))






g.1 <- r1.df[,c("annot", "chr", "position", "ID")]

g.1 %>% arrange(ID) %>% write.csv(print(paste(figPath,"rank1.csv", sep = "/")), row.names = F)

write.csv(g.1,print(paste(figPath,"rank1.csv", sep = "/")), row.names = F)
g.1.kmer <- r1.df[,c("kmer","ID")]
g.1.kmer %>% arrange(ID) %>% write.csv(print(paste(figPath,"rank1.kmer.csv", sep = "/")), row.names = F)

g.1.kmer2 <- g.1.kmer %>% arrange(ID)
g.1.kmer2 <- g.1.kmer2[,c("kmer")]
write.csv(g.1.kmer2,print(paste(figPath,"rank1.kmer2.csv", sep = "/")), row.names = F)


write.csv(g.1.kmer,print(paste(figPath,"rank1.kmer.csv", sep = "/")), row.names = F)

g.6 <- r6[,c("annot", "chr", "position", "SS", "SSrank")]
write.csv(g.6,print(paste(figPath,"rank6.csv", sep = "/")), row.names = F)
g.6.kmer <- r6[,c("kmer")]
write.csv(g.6.kmer,print(paste(figPath,"rank6.kmer.csv", sep = "/")), row.names = F)




g.6 <- r6[,c("annot", "chr", "position", "SS", "SSrank")]
write.csv(g.6,print(paste(figPath,"rank6.csv", sep = "/")), row.names = F)
g.6.kmer <- r6[,c("kmer")]
write.csv(g.6.kmer,print(paste(figPath,"rank6.kmer.csv", sep = "/")), row.names = F,)



# DATA FRAME: LOCATION ----------------------------------------------------
loco <- read_csv(paste(dataPath, "ST_loco.csv", sep = "/"))

g = readGFF("/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU/CODE/bootSTRAP/REF/gencode.v38.annotation.gff3")
gff3=as.data.frame(g)

gff3 <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]


genes.df <- read_csv(paste(dataPath, "genesDF.csv", sep = "/"))
result <- r6.out %>%
  left_join(genes.df[,c("Annotation", "strand")], by = "Annotation")



##--ALL CONDITIONS--##
CL.mm <- as.vector(outer(CellLine, c("mm.Direct"), paste, sep="."))
mm <- ranking.6[,c("ACP","Annotation", "chr", "position", "kmer", "color",CL.mm)]

##--IGV--##
r6.out <- mm
# r6.out$CP <- paste(r6.out$chr, r6.out$position, sep = ":")
# r6.out$kmer.start <- r6.out$position - 5
# r6.out$kmer.end <- r6.out$position + 5
# r6.out$range <- paste(r6.out$kmer.start, r6.out$kmer.end, sep = "-")
# r6.out$kmer.11 <- paste(r6.out$chr,r6.out$range , sep = ":")


r6.out$CDS = F
r6.out$Three_prime_UTR = F
r6.out$Five_Prime_UTR = F
r6.out$start_codon = F
r6.out$stop_codon = F
r6.out$stop_codon_redefined_as_selenocysteine = F

for (row in 1:nrow(r6.out)){
  # type = gff3$type[which(gff3$start <= r6.out$position[row] & 
  #                          gff3$end >= r6.out$position[row] &
  #                          gff3$gene_name == r6.out$Annotation[row])]
  
  
  type = gff3$type[which(gff3$start <= r6.out$position[row] & 
                           gff3$end >= r6.out$position[row] &
                           gff3$seqid == r6.out$chr[row])]
  
  if (length(type)>0){
    if ("CDS" %in% type) {
      r6.out$CDS[row] = T
    }
    if ("three_prime_UTR" %in% type) {
      r6.out$Three_prime_UTR[row] = T
    }
    if ("five_prime_UTR" %in% type) {
      r6.out$Five_Prime_UTR[row] = T
    }
    if ("start_codon" %in% type) {
      r6.out$start_codon[row] = T
    }
    if ("stop_codon" %in% type) {
      r6.out$stop_codon[row] = T
    }
    if ("stop_codon_redefined_as_selenocysteine" %in% type) {
      r6.out$stop_codon_redefined_as_selenocysteine[row] = T
    }
  }
}


start.codon <- r6.out[which(r6.out$start_codon == 1),]

all.false <- r6.out[which(r6.out$CDS == 0 &
                            r6.out$Three_prime_UTR == 0 &
                            r6.out$Five_Prime_UTR == 0 &
                            r6.out$start_codon == 0 &
                            r6.out$stop_codon == 0),]



write.csv(pus7, print(paste(pvalPath, "pus7-conserved.10reads-10mm.csv", sep = "/")), row.names = F)






pos <- 197486726
a <- hspd1[which(hspd1$start <= pos & hspd1$end>= pos),]






hspd1 <- gff3[which(gff3$gene_name == "HSPD1" & gff3$gene_type == "protein_coding" & ((gff3$type == "exon") | (gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR"))), c("type", "start", "end", "gene_name")]
pos <- 197486726
a <- hspd1[which(hspd1$start <= pos & hspd1$end>= pos),]


KIAA0930 <- gff3[which(gff3$gene_name == "KIAA0930" & gff3$gene_type == "protein_coding" & ((gff3$type == "exon") | (gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR"))), c("type", "start", "end", "gene_name")]
pos <- 45196436
a <- KIAA0930[which(KIAA0930$start <= pos & KIAA0930$end>= pos),]





# BOOTSTRAP: Variable Conserved -------------------------------------------
if (T) {
  PUS.TPM.BS.list <- list()
  
  for (CL in CellLine) {
    print(paste("start:CL",CL, sep=" "))
    
    BS.list <- list()
    
    for (i in 1:10) {
      in.file.name <- paste0(CL,"_Direct_merged.Sampled-",i,".hg38v10.abundance.tsv")
      TPM <- read_tsv(print(paste(TPMbsPath, in.file.name, sep = "/")))
      
      colnames(TPM) <- c("Gene_ID","gene_name","Reference","Strand",
                         "Start","End", "Coverage", "FPKM", "TPM")
      BS.temp <- data.frame(CL = CL,
                            BS = i,
                            PLD3 = TPM$TPM[which(TPM$gene_name == "PLD3")],
                            RNF167 = TPM$TPM[which(TPM$gene_name == "RNF167")],
                            FKBP4 = TPM$TPM[which(TPM$gene_name == "FKBP4")],
                            PMPCB = TPM$TPM[which(TPM$gene_name == "PMPCB")],
                            SLC2A1 = TPM$TPM[which(TPM$gene_name == "SLC2A1")])
      
      BS.list[[i]] <- BS.temp
    }
    
    tpm.temp <- data.table::rbindlist(BS.list)
    PUS.TPM.BS.list[[CL]] <- tpm.temp
  }
  
  PUS.TPM.BS <- data.table::rbindlist(PUS.TPM.BS.list)
  #write.csv(PUS.TPM.BS, print(paste(pvalPath,"BS-PUStpm.csv", sep = "/")), row.names = F)
  
}

PUS.TPM.BS$BS <- paste0("BS",PUS.TPM.BS$BS)

PUS.TPM.BS.melt <- melt(PUS.TPM.BS, id = c("CL", "BS"))

##--BOX PLOT--##

psi.syn.box <- PUS.TPM.BS.melt %>% 
  ggplot(aes(x = factor(variable), y = log10(value), fill = CL)) + 
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(position = position_dodge(width = 0.8), size = .25) +
  scale_fill_manual(values = cp1) +
  ylab("log10(TPM)")+
  theme(axis.line=element_line(linewidth = 0.25, colour = "black", linetype=1),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        #axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
psi.syn.box
ggsave(print(paste(figPath, "BS-conserved.TPM-box.eps", sep = "/")))


##--BAR PLOT--##
psi.syn.box <- PUS.TPM.BS.melt %>% 
  ggplot(aes(x = factor(variable), y = value, fill = CL)) + 
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), width = 0.5, color="grey") +
  geom_point(position = position_dodge(width = 0.8), size = 0.5, color="black") +
  scale_fill_manual(values = cp1) +
  theme(axis.line=element_line(linewidth = 0.25, colour = "grey", linetype=1),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
psi.syn.box
ggsave(print(paste(figPath, "BS-conserved.TPM-bar.eps", sep = "/")), width = 6, height = 4)






##--BAR-ABSOLUTE PLOT--##
summary_df <- PUS.TPM.BS.melt %>%
  group_by(variable, CL) %>%
  summarize(mean_value = mean(value), sd_value = sd(value))

psi.syn.bar <- summary_df %>% 
  ggplot(aes(x = variable, y = mean_value, fill = CL)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cp1)+
  geom_errorbar(
    aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
    position = position_dodge(width = 0.7),
    width = 0.25
  ) +
  theme(axis.line=element_line(linewidth = 0.25, colour = "black", linetype=1),
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        #legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
psi.syn.bar
ggsave(print(paste(figPath, "BS-conserved.TPM-barABS.eps", sep = "/")))




fkbp4.melt <- PUS.TPM.BS.melt[which(PUS.TPM.BS.melt$variable == "FKBP4"),]
fkbp4.df <- pivot_wider(fkbp4.melt, names_from = BS, values_from = value)
fkbp4.stats <- fkbp4.melt %>%
  group_by(CL) %>%
  summarize(
    Mean_TPM = mean(value),
    Median_TPM = median(value),
    SD_TPM = sd(value),
    Min_TPM = min(value),
    Max_TPM = max(value)
  )


slc.melt <- PUS.TPM.BS.melt[which(PUS.TPM.BS.melt$variable == "SLC2A1"),]
slc.df <- pivot_wider(slc.melt, names_from = BS, values_from = value)
slc.stats <- slc.melt %>%
  group_by(CL) %>%
  summarize(
    Mean_TPM = mean(value),
    Median_TPM = median(value),
    SD_TPM = sd(value),
    Min_TPM = min(value),
    Max_TPM = max(value)
  )



##--Conserved Occupancy--##

if (T) {
  PUS.TPM.BS.list <- list()
  
  for (CL in CellLine) {
    print(paste("start:CL",CL, sep=" "))
    
    BS.list <- list()
    
    for (i in 1:10) {
      in.file.name <- paste0(CL,"_Direct_merged.Sampled-",i,".hg38v10.abundance.tsv")
      TPM <- read_tsv(print(paste(TPMbsPath, in.file.name, sep = "/")))
      
      colnames(TPM) <- c("Gene_ID","gene_name","Reference","Strand",
                         "Start","End", "Coverage", "FPKM", "TPM")
      BS.temp <- data.frame(CL = CL,
                            BS = i,
                            NIP7 = TPM$TPM[which(TPM$gene_name == "NIP7")],
                            IDI1 = TPM$TPM[which(TPM$gene_name == "IDI1")],
                            GTF3C2 = TPM$TPM[which(TPM$gene_name == "GTF3C2")],
                            AK2 = TPM$TPM[which(TPM$gene_name == "AK2")],
                            POGK = TPM$TPM[which(TPM$gene_name == "POGK")])
      
      BS.list[[i]] <- BS.temp
    }
    
    tpm.temp <- data.table::rbindlist(BS.list)
    PUS.TPM.BS.list[[CL]] <- tpm.temp
  }
  
  PUS.TPM.BS <- data.table::rbindlist(PUS.TPM.BS.list)
  #write.csv(PUS.TPM.BS, print(paste(pvalPath,"BS-PUStpm.csv", sep = "/")), row.names = F)
  
}

PUS.TPM.BS$BS <- paste0("BS",PUS.TPM.BS$BS)

PUS.TPM.BS.melt <- melt(PUS.TPM.BS, id = c("CL", "BS"))

nip7.melt <- PUS.TPM.BS.melt[which(PUS.TPM.BS.melt$variable == "NIP7"),]
nip7.df <- pivot_wider(nip7.melt, names_from = BS, values_from = value)
nip7.stats <- nip7.melt %>%
  group_by(CL) %>%
  summarize(
    Mean_TPM = mean(value),
    Median_TPM = median(value),
    SD_TPM = sd(value),
    Min_TPM = min(value),
    Max_TPM = max(value)
  )


IDI1.melt <- PUS.TPM.BS.melt[which(PUS.TPM.BS.melt$variable == "IDI1"),]
IDI1.df <- pivot_wider(IDI1.melt, names_from = BS, values_from = value)
IDI1.stats <- IDI1.melt %>%
  group_by(CL) %>%
  summarize(
    Mean_TPM = mean(value),
    Median_TPM = median(value),
    SD_TPM = sd(value),
    Min_TPM = min(value),
    Max_TPM = max(value)
  )

IDI1.df <- pivot_wider(PUS.TPM.BS.melt, names_from = BS, values_from = value)
IDI1.stats <- PUS.TPM.BS.melt %>%
  group_by(CL, variable) %>%
  summarize(
    Mean_TPM = mean(value),
    Median_TPM = median(value),
    SD_TPM = sd(value),
    Min_TPM = min(value),
    Max_TPM = max(value)
  )






# DATA FRAME: Orthogonal Parsing ------------------------------------------

ortho.df <- read_csv(paste(orthoPATH,"PSI.orthogonal.MASTER.csv", sep = "/"))

df <- ortho.df[which(ortho.df$ortho > 0),]

##--BASELINE CONDITIONS--##
mm <- ortho.df
mm$filter <- 0

min.mm <- 20

#-ONE-#
mm$filter[rowSums(mm[, CL.mm] >= min.mm) > 0] <- 1
mm.one <- mm[which(mm$filter == 1),]
mm.one <- mm[which(mm$filter == 1 | mm$ortho > 0),]

#-ALL-#
mm$filter <- 0
mm$filter[rowSums(mm[, CL.mm] > min.mm) == length(CL.mm)] <- 1
mm.all <- mm[which(mm$filter == 1),]

mm.ortho <- mm.one[which(mm.one$ortho > 0),]

mm.one <- mm[which(mm$filter == 1 | mm$ortho > 0),]



df <- ortho.df[which(ortho.df$ortho > 0), c("Annotation", CL.mm)]
CL.sim.TRUB1 <- mm.all[which(mm.all$color == "trub1"),c("Annotation",CL.mm)]
CL.sim.TRUB1 <- subset(CL.sim.TRUB1, !(Annotation == "FLAD1" | Annotation == "GRWD1"))

df <- ortho.df[which(ortho.df$ortho > 0), c("Annotation", CL.mm)]
df[, 2:7] <- df[, 2:7] / 100

# Numeric encoding of Transcripts
transcript_mapping <- df %>%
  dplyr::distinct(Annotation) %>%
  dplyr::mutate(Transcript_ID = row_number()) %>%
  dplyr::select(Annotation, Transcript_ID)


data_encoded <- df %>%
  left_join(transcript_mapping, by = "Annotation") %>%
  dplyr::select(-Annotation)

cor_matrix <-  data_encoded%>%
  dplyr::select(-Transcript_ID, ) %>%
  cor(method = "pearson")
hc <- hclust(as.dist(1 - cor_matrix), method = "complete")

plot(hc, main = "Hierarchical Clustering Dendrogram")


row_names <- col_names <- c(CL.mm)

# Convert the correlation matrix into a tidy data frame
cor_df$row_names <- row_names
cor_df <- melt(cor_df, id.vars = "row_names")
colnames(cor_df) <- c("Row", "Column", "Correlation")

# Create a heatmap using ggplot2
ggplot(cor_df, aes(x = Row, y = Column, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Correlation Heatmap") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))













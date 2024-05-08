## ========================================== ##
##             Psi Conservation               ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 11/28/2023
#@ content: Figure 4
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
# path <- scan(what = "character", n=1) #actual directory needs to be on the following line
# "/home/cam/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

path <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"

dirPath <- print(path)
setwd(dirPath) #set working directory

dataPath <- print(paste(dirPath, "DATA", sep = "/"))
pvalPath <- print(paste(dataPath, "pValues", sep = "/"))
figPath <- print(paste(dirPath, "Figures", "Figure-4","FIG4-elements", sep = "/"))
bamPath <- print(paste(dirPath, "BAM", "hg38v10", sep = "/"))
tpmPath <- print(paste(bamPath, "stringtie",sep = "/"))
pysamPath <- print(paste(dataPath, "pysamstats", sep = "/"))
pvalSOURCE <- print(paste(dirPath,"pValues","DATA", sep = "/"))
kpPath <- print(paste(dataPath, "seqLOGO", sep = "/"))
kdePath <- print(paste(dataPath, "KDE", sep = "/"))
rdsPath <- print(paste(dataPath, "RData", sep = "/"))
dmodPath <- print(paste(dataPath, "strandANALYSIS", sep = "/"))
bootPath <- print(paste(dirPath,"CODE","bootSTRAP", sep = "/"))
bsPath <- print(paste(bootPath, "pVALUES", sep="/"))
TPMbsPath <- print(paste(bootPath, "stringtie/DATA", sep = "/"))
orthoPATH <- print(paste(dataPath, "orthogonal", sep = "/"))

###--- Colors ---###
cp1 <- c("#FF595E", "#FF924C", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93")

CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")

cat("\f") #clears console

# DATA FRAME: pVal MASTER -------------------------------------------------

if (T){
  master.df <- read_csv(paste(pvalPath,"sixCellLine.pval-IVT.mmDirect.PSI.master.csv", sep = "/"))
  # master.df <- master.bed
  
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
    CL.tmp <- CL.tmp[which(CL.tmp$N_reads_Direct >= 10 &
                             CL.tmp$N_reads_IVT >= 10 & 
                             CL.tmp$mm.IVT <= 10),]
    
    CL.tmp$psi <- 0 
    CL.tmp$psi[which(CL.tmp$p.value.Direct < 0.001 & CL.tmp$mm.Direct >= 10 & CL.tmp$N_reads_Direct >= 30)] <- 1
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
  
  
  ##--REQUIRE 30 READS ACROSS--##
  #ACP.master <- ACP.master[complete.cases(ACP.master), ]  
  
  
  ##--ALLOW INCOMPLETE CASES--##
  ACP.master[is.na(ACP.master)] <- 0
  ACP.master$A549.p.value.Direct[which(ACP.master$A549.N_reads_Direct == 0)] <- 1
  ACP.master$HeLa.p.value.Direct[which(ACP.master$HeLa.N_reads_Direct == 0)] <- 1
  ACP.master$HepG2.p.value.Direct[which(ACP.master$HepG2.N_reads_Direct == 0)] <- 1
  ACP.master$Jurkat.p.value.Direct[which(ACP.master$Jurkat.N_reads_Direct == 0)] <- 1
  
  
  #pval.file <- paste("sixCellLine",".pval.filtered.master.csv", sep = "")
  pval.file <- paste("sixCellLine.pval.IC.filtered.master.csv", sep = "")
  
  write.csv(ACP.master,print(paste(pvalPath, pval.file, sep = "/")), row.names = F)
}


# DATA FRAME: Conserved REwork --------------------------------------------

if (T) {
  stable <- c("ACP", "Annotation", "chr", "position", "kmer")
  CL.ACP <- as.vector(outer(CellLine, c("N_reads_Direct", "N_reads_IVT","mm.Direct", "mm.IVT", "p.value.Direct"), paste, sep="."))
  CL.psi <- as.vector(outer(CellLine, c("psi"), paste, sep="."))
  
  master <- read_csv(paste(pvalPath,"sixCellLine.pval-IVT.mmDirect.PSI.master.csv", sep = "/"))
  master.df <- master[,c(stable, CL.ACP)]
  
  master.df[is.na(master.df)] <- 0
  master.df$A549.p.value.Direct[which(master.df$A549.N_reads_Direct == 0)] <- 1
  master.df$HeLa.p.value.Direct[which(master.df$HeLa.N_reads_Direct == 0)] <- 1
  master.df$HepG2.p.value.Direct[which(master.df$HepG2.N_reads_Direct == 0)] <- 1
  master.df$Jurkat.p.value.Direct[which(master.df$Jurkat.N_reads_Direct == 0)] <- 1
  
  
  master.list <- list()
  
  for (CL in CellLine) {
    
    CL.col <- as.vector(outer(CL, c("N_reads_Direct", "N_reads_IVT","mm.Direct", "mm.IVT", "p.value.Direct"), paste, sep="."))
    CL.tmp <- master.df[,c(stable,CL.col)]
    
    prefix_to_remove <- paste0(CL,".")
    columns_to_exclude <- c(stable)
    names(CL.tmp) <- ifelse(names(CL.tmp) %in% columns_to_exclude, names(CL.tmp), 
                            gsub(paste0("^", prefix_to_remove), "", names(CL.tmp)))
    
    CL.tmp$mm.Direct[CL.tmp$mm.Direct < 10] <- 0
    CL.tmp <- CL.tmp[which(CL.tmp$mm.IVT <= 10),]
    
    CL.tmp$psi <- 0 
    CL.tmp$psi[which(CL.tmp$p.value.Direct < 0.001 & CL.tmp$mm.Direct >= 20 & CL.tmp$N_reads_Direct >= 20)] <- 1
    
    CL.tmp$nReadMIN <- 0
    CL.tmp$nReadMIN[which(CL.tmp$N_reads_Direct >= 20)] <- 1
    
    CL.tmp$psi.nReadMIN <- 0
    CL.tmp$psi.nReadMIN[which(CL.tmp$psi == 1 & CL.tmp$nReadMIN == 1)] <- 1
    
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
  
  
  ACP.master[is.na(ACP.master)] <- 0
  ACP.master$A549.p.value.Direct[which(ACP.master$A549.N_reads_Direct == 0)] <- 1
  ACP.master$HeLa.p.value.Direct[which(ACP.master$HeLa.N_reads_Direct == 0)] <- 1
  ACP.master$HepG2.p.value.Direct[which(ACP.master$HepG2.N_reads_Direct == 0)] <- 1
  ACP.master$Jurkat.p.value.Direct[which(ACP.master$Jurkat.N_reads_Direct == 0)] <- 1
  
  # pval.file <- paste("sixCellLine.pval.IC.filtered.master.csv", sep = "")
  pval.file <- paste("sixCellLine.pval.IC-20.filtered.master.csv", sep = "")
  
  write.csv(ACP.master,print(paste(pvalPath, pval.file, sep = "/")), row.names = F)
}

##


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
  
  
  # write.csv(ranking.1,print(paste(pvalPath,"ranking1.csv", sep = "/")), row.names = F)
  # write.csv(ranking.2,print(paste(pvalPath,"ranking2.csv", sep = "/")), row.names = F)
  # write.csv(ranking.3,print(paste(pvalPath,"ranking3.csv", sep = "/")), row.names = F)
  # write.csv(ranking.4,print(paste(pvalPath,"ranking4.csv", sep = "/")), row.names = F)
  # write.csv(ranking.5,print(paste(pvalPath,"ranking5.csv", sep = "/")), row.names = F)
  # write.csv(ranking.6,print(paste(pvalPath,"ranking6.csv", sep = "/")), row.names = F)
  
}




# FIGURE: TYPE II -- Conserved --------------------------------------------
##--Column Vectors--##

if (T) {
  stable <- c("ACP", "Annotation", "chr", "position", "kmer")
  
  master.df <- read_csv(paste(pvalPath, "sixCellLine.pval.filtered.master.csv", sep = "/"))
  # master.df <- read_csv(paste(pvalPath, "sixCellLine.pval.IC.filtered.master.csv", sep = "/"))
  #master.df <- read_csv(paste(pvalPath, "sixCellLine.pval.IC-20.filtered.master.csv", sep = "/"))
  
  ranking <- master.df
  
  annot.list <- list()
  
  for (CL in CellLine) {
    CL.filter <- as.vector(outer(CL, c("N_reads_Direct", "mm.Direct", "p.value.Direct", "psi"), paste, sep="."))
    col.filter <- c(stable, CL.filter, "rank")
    CL.tmp <- ranking[, which(names(ranking) %in% col.filter)]
    
    prefix_to_remove <- paste0(CL,".")
    columns_to_exclude <- c(stable, "rank")
    names(CL.tmp) <- ifelse(names(CL.tmp) %in% columns_to_exclude, names(CL.tmp), 
                            gsub(paste0("^", prefix_to_remove), "", names(CL.tmp)))
    
    CL.tmp$id <- CL
    CL.tmp <- CL.tmp[which(CL.tmp$psi == 1),]
    
    count_data <- table(CL.tmp$id, CL.tmp$Annotation)
    multiPOS <- as.data.frame(count_data)
    multiPOS <- multiPOS %>%
      rename("Var1" = "Cell",
             "Var2" = "Gene",
             "Freq" = "Count")
    multiPOS <- multiPOS[multiPOS$Count != 0, c("Gene", "Count")]
    
    multiPOS <- multiPOS %>% rename_with(~ paste0(CL,".", .), -c(Gene))
    
    
    annot.list[[CL]] <- multiPOS
  }
  
  
  
  multiPOS <- annot.list %>% reduce(full_join, by=c("Gene"))
  multiPOS[is.na(multiPOS)] <- 0
  
  
  multiPOS$hyper <- 0
  multiPOS$hyper[which(multiPOS$A549.Count > 1 |
                         multiPOS$HeLa.Count > 1 |
                         multiPOS$HepG2.Count > 1 |
                         multiPOS$Jurkat.Count > 1 |
                         multiPOS$NTERA.Count > 1 |
                         multiPOS$SHSY5Y.Count > 1)] <- 1
  multiPOS <- multiPOS[which(multiPOS$hyper == 1),]
  
  CL.count <- as.vector(outer(CellLine, c("Count"), paste, sep="."))
  multiPOS$SS <- round(apply(multiPOS[,c(CL.count)], 1, sd), digits = 2)
  
  
  
  hyper.list <- list()
  
  for (CL in CellLine) {
    CL.filter <- as.vector(outer(CL, c("N_reads_Direct", "mm.Direct", "p.value.Direct", "psi"), paste, sep="."))
    col.filter <- c(stable, "rank",CL.filter)
    CL.tmp <- ranking[, which(names(ranking) %in% col.filter)]
    
    prefix_to_remove <- paste0(CL,".")
    columns_to_exclude <- c(stable, "rank")
    names(CL.tmp) <- ifelse(names(CL.tmp) %in% columns_to_exclude, names(CL.tmp), 
                            gsub(paste0("^", prefix_to_remove), "", names(CL.tmp)))
    CL.tmp <- CL.tmp[which(CL.tmp$psi == 1),]
    CL.tmp$id <- CL 
    
    CL.tmp <- CL.tmp[CL.tmp$Annotation %in% multiPOS$Gene,]
    
    CL.tmp <- CL.tmp[,c("ACP", "Annotation", "chr", "position", "kmer", "rank", "psi", "mm.Direct")]
    
    CL.tmp <- CL.tmp %>% rename_with(~ paste0(CL,".", .), -c(ACP, Annotation, chr, position, kmer, rank))
    
    hyper.list[[CL]] <- CL.tmp
  }
  
  
  hyper.df.all <- hyper.list %>% reduce(full_join, by=c("ACP", "Annotation", "chr", "position", "kmer", "rank"))
  hyper.df.all[is.na(hyper.df.all)] <- 0
  
  # CL.nReadMIN <- as.vector(outer(CellLine, c("nReadMIN"), paste, sep="."))
  # nReadMIN <- ranking[,c(stable, CL.nReadMIN)]
  # nReadMIN$rank <- rowSums(nReadMIN[CL.nReadMIN])
  # nReadMIN <- nReadMIN[which(nReadMIN$rank > 1),]
  # 
  # hyper.df.all <- hyper.df.all %>%
  #   filter(ACP %in% nReadMIN$ACP)
  
  
  
  CL.psi <- as.vector(outer(CellLine, c("psi"), paste, sep="."))
  hyper.df.all$rank <- rowSums(hyper.df.all[CL.psi])
  hyper.df.all <- hyper.df.all[which(hyper.df.all$rank > 0),]
  
  hyper.df.all$color <- "other"
  hyper.df.all$color[hyper.df.all$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  hyper.df.all$color[hyper.df.all$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
}


HDF <- hyper.df.all

count_data <- table(HDF$Annotation)
HDF1 <- as.data.frame(count_data)
HDF1 <- HDF1 %>%
  rename("Var1" = "Annotation",
         "Freq" = "Count")
HDF1 <- HDF1[HDF1$Count > 0, c("Annotation", "Count")]

hyper2 <- hyper.df.all %>%
  filter(Annotation %in% HDF1$Annotation)

hyper.trub1 <- hyper2[which(hyper2$color == "trub1"),]
hyper.pus7 <- hyper2[which(hyper2$color == "pus7"),]



##--LOCO ANNOTATION GFF3--##
library(rtracklayer)
if (T) {
  g = readGFF(paste(refPath,"gencode.v38.annotation.gff3",sep = "/"))
  gff3=as.data.frame(g)

  gff3 <- gff3[which(gff3$gene_type == "protein_coding" & ((gff3$type == "CDS") | (gff3$type == "three_prime_UTR") | (gff3$type == "five_prime_UTR") | (gff3$type == "start_codon") | (gff3$type == "stop_codon") | (gff3$type == "stop_codon_redefined_as_selenocysteine"))), c("type", "seqid","start", "end", "gene_name")]
}




if (T) {
  
  hyper.df.all$CDS = F
  hyper.df.all$Three_prime_UTR = F
  hyper.df.all$Five_Prime_UTR = F
  hyper.df.all$start_codon = F
  hyper.df.all$stop_codon = F
  hyper.df.all$stop_codon_redefined_as_selenocysteine = F
  
  for (row in 1:nrow(hyper.df.all)){
    type = gff3$type[which(gff3$start <= hyper.df.all$position[row] & 
                             gff3$end >= hyper.df.all$position[row] &
                             gff3$seqid == hyper.df.all$chr[row])]
    
    if (length(type)>0){
      if ("CDS" %in% type) {
        hyper.df.all$CDS[row] = T
      }
      if ("three_prime_UTR" %in% type) {
        hyper.df.all$Three_prime_UTR[row] = T
      }
      if ("five_prime_UTR" %in% type) {
        hyper.df.all$Five_Prime_UTR[row] = T
      }
      if ("start_codon" %in% type) {
        hyper.df.all$start_codon[row] = T
      }
      if ("stop_codon" %in% type) {
        hyper.df.all$stop_codon[row] = T
      }
      if ("stop_codon_redefined_as_selenocysteine" %in% type) {
        hyper.df.all$stop_codon_redefined_as_selenocysteine[row] = T
      }
    }
  }
  
  
  start.codon <- hyper.df.all[which(hyper.df.all$start_codon == 1),]
  cds <- hyper.df.all[which(hyper.df.all$CDS == 1),]
  stop.codon <- hyper.df.all[which(hyper.df.all$stop_codon == 1),]
  five.prime <- hyper.df.all[which(hyper.df.all$Five_Prime_UTR == 1),]
  
  all.false <- hyper.df.all[which(hyper.df.all$CDS == 0 &
                                hyper.df.all$Three_prime_UTR == 0 &
                                hyper.df.all$Five_Prime_UTR == 0 &
                                hyper.df.all$start_codon == 0 &
                                hyper.df.all$stop_codon == 0),]
}



pval.file <- paste("sixCellLine.IC-20.hyperII.csv", sep = "")


write.csv(hyper.df.all,print(paste(pvalPath, pval.file, sep = "/")), row.names = F)

# FIGURE: strand analysis -------------------------------------------------

##--DOUBLE MOD--##
if (T) {
  
  LAPTM4B.Direct.list <- list()
  
  for (CL in CellLine) {
    
    file.name <- paste0(CL,"_Direct_LAPTM4B-StrandAnalysis.txt")
    file.path <- paste(dmodPath, file.name, sep = "/")
    
    dMOD.table <- read.table(file.path)
    
    colnames(dMOD.table) <- c("position", "read", "base")
    dMOD.poi <- dMOD.table[which(dMOD.table$position == 97776037 | 
                                   dMOD.table$position == 97776039),]
    dMOD.poi$position <- dMOD.poi$position + 1
    
    
    dMOD.df <- pivot_wider(dMOD.poi, 
                           names_from = position, 
                           values_from = base)
    #dMOD.df <- dMOD.df[complete.cases(dMOD.df), ]
    
    LAPTM4B.Direct.list[[CL]] <- dMOD.df
    
  } 
  
  
  LAPTM4B.IVT.list <- list()
  
  for (CL in CellLine) {
    
    file.name <- paste0(CL,"_IVT_LAPTM4B-StrandAnalysis.txt")
    file.path <- paste(dmodPath, file.name, sep = "/")
    
    dMOD.table <- read.table(file.path)
    
    colnames(dMOD.table) <- c("position", "read", "base")
    dMOD.poi <- dMOD.table[which(dMOD.table$position == 97776037 | 
                                   dMOD.table$position == 97776039),]
    dMOD.poi$position <- dMOD.poi$position + 1
    
    
    dMOD.df <- pivot_wider(dMOD.poi, 
                           names_from = position, 
                           values_from = base)
    dMOD.df <- dMOD.df[complete.cases(dMOD.df), ]
    
    LAPTM4B.IVT.list[[CL]] <- dMOD.df
    
  } 
}

lap.a549 <- LAPTM4B.Direct.list$A549
lap.hela <- LAPTM4B.Direct.list$HeLa
lap.hepg2 <- LAPTM4B.Direct.list$HepG2
lap.jurkat <- LAPTM4B.Direct.list$Jurkat
lap.ntera <- LAPTM4B.Direct.list$NTERA
lap.shsy5y <- LAPTM4B.Direct.list$SHSY5Y


lap <- rbind(lap.a549, lap.hela, lap.hepg2, lap.jurkat, lap.ntera, lap.shsy5y)
lap.cov <- lap[complete.cases(lap), ]


lap.filter <- lap.cov[which(lap.cov$`97776038` != "A"),]
lap.filter <- lap.filter[which(lap.filter$`97776040` != "A"),]
lap.filter$combo <- paste0(lap.filter$`97776038`, lap.filter$`97776040`)

nrow(lap.filter[which(lap.filter$combo == "CC"),])
nrow(lap.filter[which(lap.filter$combo == "TC"),])
nrow(lap.filter[which(lap.filter$combo == "CT"),])
nrow(lap.filter[which(lap.filter$combo == "TT"),])


# 
# 


# DATA FRAME: Ionic Currents ----------------------------------------------
##--GATHER IONIC CURRENTS FOR POI--##

DRS.type <- c("Direct", "IVT")
CellLine.samp <- c("HepG2") 
dm.annot <- c("LAPTM4B")

if (T) {
  drs.list.Direct <- list()
  
  for (drs in DRS.type) {
    
    CL.DRS.list <- list()
    
    for (CL in CellLine) {
      
      print(paste("Start:", CL, dm.annot, drs))
      
      rep.name <- paste(CL, drs, dm.annot, sep = "_")
      
      
      in.folder <- paste(rep.name, "txt", sep = ".")
      in.tsv <- paste(dmodPath, "EA",in.folder, sep="/")
      
      EAC <-  read.table(in.tsv, header = TRUE)
      EAC.POI <- EAC[which(EAC$position == (97776037-2) | EAC$position == (97776039-2)),
                     c("contig", "position", "read_name","reference_kmer", "event_level_mean", "event_stdv", 
                       "model_kmer", "model_mean", "model_stdv")]
      
      
      CL.DRS.list[[CL]] <- EAC.POI
      
    }
    
    drs.list.Direct[[drs]] <- CL.DRS.list
    
  }
  
  master.rds.file <- paste("doublemod.IonicCurrent.RData", sep = "")
  saveRDS(drs.list.Direct, print(paste(rdsPath, master.rds.file, sep = "/")))
}




# DATA FRAME: Ortho multi -------------------------------------------------

r.ortho <- read_csv(paste(orthoPATH, "PSI.orthogonal.MASTER.csv", sep = "/"))
ortho.df <- r.ortho[which(r.ortho$ortho > 0),]

##--multiPOS stats--##
multiPOS_summary <- multiPOS %>%
  group_by(CL, variable) %>%
  summarise(mean_value = mean(value),
            sd_value = sd(value))

count_data <- table(ortho.df$Annotation)
ortho.multiPOS <- as.data.frame(count_data)

multiPOS <- multiPOS %>%
  rename("Var1" = "Cell",
         "Var2" = "Gene",
         "Freq" = "Count")

multiPOS <- multiPOS %>%
  rename("Cell" = "Var1",
         "Gene" = "Var2",
         "Count" = "Freq")
multiPOS <- multiPOS[multiPOS$Count >= 2,]

summary_data <- aggregate(Gene ~ Cell + Count, data = multiPOS, FUN = function(x) length(unique(x)))
summary_data <- summary_data[summary_data$Count != 0,]


# FIGURE: Double Mod Strand Analysis --------------------------------------


##--DOUBLE MOD STRAND ANALYSIS--##
if (T) {
  
  LAPTM4B.Direct.list <- list()
  
  for (CL in CellLine) {
    
    file.name <- paste0(CL,"_Direct_LAPTM4B-StrandAnalysis.txt")
    file.path <- paste(dmodPath, file.name, sep = "/")
    
    dMOD.table <- read.table(file.path)
    
    colnames(dMOD.table) <- c("position", "read", "base")
    dMOD.poi <- dMOD.table[which(dMOD.table$position == 97776037 | 
                                   dMOD.table$position == 97776039),]
    dMOD.poi$position <- dMOD.poi$position + 1
    
    
    dMOD.df <- pivot_wider(dMOD.poi, 
                           names_from = position, 
                           values_from = base)
    dMOD.df <- dMOD.df[complete.cases(dMOD.df), ]
    
    LAPTM4B.Direct.list[[CL]] <- dMOD.df
    
  } 
  
  
  LAPTM4B.IVT.list <- list()
  
  for (CL in CellLine) {
    
    file.name <- paste0(CL,"_IVT_LAPTM4B-StrandAnalysis.txt")
    file.path <- paste(dmodPath, file.name, sep = "/")
    
    dMOD.table <- read.table(file.path)
    
    colnames(dMOD.table) <- c("position", "read", "base")
    dMOD.poi <- dMOD.table[which(dMOD.table$position == 97776037 | 
                                   dMOD.table$position == 97776039),]
    dMOD.poi$position <- dMOD.poi$position + 1
    
    
    dMOD.df <- pivot_wider(dMOD.poi, 
                           names_from = position, 
                           values_from = base)
    dMOD.df <- dMOD.df[complete.cases(dMOD.df), ]
    
    LAPTM4B.IVT.list[[CL]] <- dMOD.df
    
  } 
}

##--DOUBLE MOD IONIC CURRENT--##
if (T){
  dm.rds.file <- paste("doublemod.IonicCurrent.RData", sep = "")
  dm.list <- readRDS(print(paste(rdsPath, dm.rds.file, sep = "/")))
  
  pos.nuc <- list(
    Direct = data.frame(do.call(rbind, LAPTM4B.Direct.list)),
    IVT = data.frame(do.call(rbind, LAPTM4B.IVT.list))
  )
  
  ic.list <- list(
    Direct = data.frame(do.call(rbind, dm.list$Direct)),
    IVT = data.frame(do.call(rbind, dm.list$IVT))
    
  )
  
  
  
  
  
  combos <- c("CT", "TC", "TT", "CC", "CA", "AC")
  DRS.type <- c("Direct", "IVT")
  
  POI.reads.list <- list()
  for (drs in DRS.type) {
    drs.tmp <- pos.nuc[[drs]]
    ic.tmp <- ic.list[[drs]]
    
    c.list <- list()
    for (c in combos){
      
      c.split <- strsplit(c,"")[[1]]
      reads <- drs.tmp$read[which(drs.tmp$X97776038 == c.split[1] & drs.tmp$X97776040 == c.split[2])]
      
      reads.df <- ic.tmp[ic.tmp$read_name %in% reads,]
      
      c.list[[c]] <- reads.df
    }
    POI.reads.list[[drs]] <- c.list
  }
}
##--INDIVIDUAL STRAND COMBOS--##
ivt.pos <- rbind(data.frame(do.call(rbind, POI.reads.list$IVT), "ID" = "IVT"))
for (c in combos){
  
  coY.min <- 0
  coY.max <- .2
  
  coX.min <- 55
  coX.max <- 100
  
  direct.38 <- POI.reads.list$Direct[[c]]
  direct.38 <- direct.38[which(direct.38$position == 97776035),]
  direct.38$ID <- "Direct"
  df.38 <- rbind(data.frame(direct.38),
                 data.frame(ivt.pos[which(ivt.pos$position == 97776035),]))
  
  PM.mean <- df.38$model_mean[1]
  PM.kmer <- df.38$model_kmer[1]
  drs.pos <- (df.38$position[1]+3)
  
  out.file <- paste("LAPTM4B", drs.pos, c, "merged.DoubleMod-IonicCurrent.eps", sep = ".")
  
  CL.plot.reads <- df.38 %>%
    ggplot(aes(x = event_level_mean, color = ID)) +
    scale_color_manual(values = c("Direct" = "red", "IVT" = "black")) +
    geom_density(fill = "transparent", linewidth=.5) +
    ggtitle(paste(PM.kmer, drs.pos, c ,sep = " ")) +
    ylim(0, 1) +
    xlim(0,200) +
    coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
    xlab("Current (pA)") +
    geom_vline(xintercept=PM.mean, color="black", linewidth=0.35, linetype=2)+
    theme(axis.line=element_line(),
          #axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
          #axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  CL.plot.reads
  ggsave(print(paste(dmodPath,"plot",out.file, sep = "/")), width = 3, height = 3)
  
  
  
  direct.40 <- POI.reads.list$Direct[[c]]
  direct.40 <- direct.40[which(direct.40$position == 97776037),]
  direct.40$ID <- "Direct"
  df.40 <- rbind(data.frame(direct.40),
                 data.frame(ivt.pos[which(ivt.pos$position == 97776037),]))
  
  
  
  PM.mean <- df.40$model_mean[1]
  PM.kmer <- df.40$model_kmer[1]
  drs.pos <- (df.40$position[1]+3)
  
  out.file <- paste("LAPTM4B", drs.pos, c, "merged.DoubleMod-IonicCurrent.eps", sep = ".")
  
  CL.plot.reads <- df.40 %>%
    ggplot(aes(x = event_level_mean, color = ID)) +
    scale_color_manual(values = c("Direct" = "red", "IVT" = "black")) +
    geom_density(fill = "transparent", linewidth=.5) +
    ggtitle(paste(PM.kmer, drs.pos, c ,sep = " ")) +
    ylim(0, 1) +
    xlim(0,200) +
    coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
    xlab("Current (pA)") +
    geom_vline(xintercept=PM.mean, color="black", linewidth=0.35, linetype=2)+
    theme(axis.line=element_line(),
          #axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
          #axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  CL.plot.reads
  ggsave(print(paste(dmodPath,"plot",out.file, sep = "/")), width = 3, height = 3)
  
  
}


##--OVERLAY--##
if (T) {
  
  coY.min <- 0
  coY.max <- .22
  
  coX.min <- 65
  coX.max <- 95
  
  ivt.pos <- rbind(data.frame(do.call(rbind, POI.reads.list$IVT), "ID" = "IVT"))
  
  direct.df <-rbind(data.frame(POI.reads.list$Direct$TT, "ID" = "TT"),
                    data.frame(POI.reads.list$Direct$CT, "ID" = "CT"),
                    data.frame(POI.reads.list$Direct$TC, "ID" = "TC"),
                    data.frame(POI.reads.list$Direct$CC, "ID" = "CC"))
  
  
  drs.df <- rbind(ivt.pos, direct.df)
  drs.df$POI <- paste0(drs.df$contig, ":",drs.df$position)
  
  chr <- drs.df$contig[1]
  pos <- unique(drs.df$position)
  POI.range <- as.vector(outer(chr, pos, paste, sep=":")) 
  
  for (poiKMER in POI.range) {
    
    tmp.df <- drs.df[which(drs.df$POI == poiKMER),]
    
    PM.mean <- tmp.df$model_mean[1]
    PM.kmer <- tmp.df$model_kmer[1]
    drs.pos <- tmp.df$position[1]
    
    pos.kmer <- which(POI.range == poiKMER)
    
    out.file <- paste("LAPTM4B", pos.kmer, "overlay-merged.DoubleMod-IonicCurrent.eps", sep = ".")
    
    CL.plot.reads <- tmp.df %>%
      ggplot(aes(x = event_level_mean, color = ID, fill=ID)) +
      scale_color_manual(values = c("TT" = "red", 
                                    "TC" = "purple", 
                                    "CT" = "green", 
                                    "CC" = "blue",
                                    "IVT" = "black")) +
      scale_fill_manual(values = c("TT" = "transparent", 
                                    "TC" = "transparent", 
                                    "CT" = "transparent", 
                                    "CC" = "transparent",
                                    "IVT" = "lightgrey")) +
      geom_density(linewidth=0.5) +
      ggtitle(paste(poiKMER, PM.kmer, sep = " ")) +
      ylim(0, 1) +
      xlim(0,200) +
      coord_cartesian(xlim = c(coX.min,coX.max), ylim = c(coY.min,coY.max)) +
      xlab("Current (pA)") +
      geom_vline(xintercept=PM.mean, color="black", linewidth=0.25, linetype=2)+
      theme(axis.line=element_line(),
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
    CL.plot.reads
    ggsave(print(paste(dmodPath,"plot",out.file, sep = "/")), width = 6, height = 3)
    
    
  }
}





























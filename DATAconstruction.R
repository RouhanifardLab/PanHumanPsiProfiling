## ========================================== ##
##             Psi Conservation               ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 06/21/2023
#@ version: 10.0

# housekeeping ------------------------------------------------------------
list.of.packages <- c("pval.rds.file", "dplyr", "tidyr", "plyr", "readr", "stringr", "tibble", "stringr", "reshape2", "limma", "caret")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) # from http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them. Thanks!
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(stringr)
library(reshape2)
library(limma)
library(caret)
library(janitor)
library(pracma)
library(tidyverse)
library(purrr)

###--- Clear Space ---###
rm(list = ls()) #clear global environment except for packages & libraries 
cat("\f") #clears console

###--- General Path & Working Directory ---###
path <- scan(what = "character", n=1) #actual directory needs to be on the following line
"/Users/carolinemccormick/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"
# 
# path <- scan(what = "character", n=1) #actual directory needs to be on the following line
# "/home/cam/Rouhanifard lab Dropbox/Papers/Conservation_PseudoU"


dirPath <- print(path)
setwd(dirPath) #set working directory

dataPath <- print(paste(dirPath, "DATA", sep = "/"))
pvalPath <- print(paste(dataPath, "pValues", sep = "/"))
figPath <- print(paste(dirPath, "Figures", "FIG1", sep = "/"))
bamPath <- print(paste(dirPath, "hg38v10", "BAM", sep = "/"))
pysamPath <- print(paste(dataPath, "pysamstats", sep = "/"))
pvalSOURCE <- print(paste(dirPath,"pValues","DATA", sep = "/"))
tpmPath <- print(paste(dataPath, "stringtie", sep = "/"))
rdsPath <- print(paste(dataPath, "RData", sep = "/"))


cat("\f") #clears console


colNAMES <- c("Annotation", "chr", "position", "kmer", 
              "T_Direct", "C_Direct", "T_IVT", "C_IVT", 
              "N_reads_Direct", "N_reads_IVT",
              "expected.mm", "mm.Direct", "mm.IVT", "p.value.Direct")



CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")
stable <- c("ACP", "Annotation", "chr", "position", "kmer")



# PVAL filter: Direct !0 --------------------------------------------------

##--CELL LINE LOOP--##

#CellLine = c("SHSY5Y")


for (CL in CellLine) {
  print(paste("start:CL",CL, sep=" "))
  
  #-Files-#
  inFile.list <- c("panIVT", "pIVT")
  panIVT <- paste(CL,".pVALs.csv", sep = "")
  pIVT <-  paste(CL,"-paired.pVALs.csv", sep = "")
  
  for (f in inFile.list) {
    print(paste("start:",CL,f, sep=" "))
    
    in.file <- get(f)
    
    pval.UF <- read_csv(paste(pvalPath, in.file, sep = "/"))
    pval.UF <- pval.UF[,c(colNAMES)]
    pval <- pval.UF[which(pval.UF$N_reads_Direct >= 10 & pval.UF$N_reads_IVT >= 10),]
    
    out.pval.file <- paste(CL, f, "pVALs.csv", sep = ".")
    write.csv(pval,print(paste(pvalPath, out.pval.file, sep = "/")), row.names = F)
  }
}



# DATA FRAME: CL-pval Loop ------------------------------------------------
#CellLine = c("HeLa.sample","Jurkat.sample","SHSY5Y.sample")
#stable <- c("ACP", "Annotation", "chr", "position", "kmer")
#CellLine = c("SHSY5Y")


#------------------------------------------------------------------------------------------#
##--MASTER FILES--##

master.list <- list()

master.rds.file <- paste("sixCellLine",".MASTER.RData", sep = "")
pval.rds.file <- paste("sixCellLine",".pval-IVT.master.RData", sep = "")
pval.file <- paste("sixCellLine",".pval-IVT.PSI.master.csv", sep = "")

#------------------------------------------------------------------------------------------#
##--CELL LINE LOOP--##


for (CL in CellLine) {
  print(paste("start:CL",CL, sep=" "))
  
  #-Files-#
  inFile.list <- c("panIVT", "pIVT")
  panIVT <- paste(CL,".panIVT.pVALs.csv", sep = "")
  pIVT <-  paste(CL,".pIVT.pVALs.csv", sep = "")
  
  #------------------------------------------------------------#
  ##--Cell Line panIVT & pIVT MASTER--##
  
  CL.pvals.list <- list()
  CL.pvals.list.prefix <- list()
  
  for (f in inFile.list) {
    print(paste("start:",CL,f, sep=" "))
    
    in.file <- get(f)
    
    pval.UF <- read_csv(paste(pvalPath, in.file, sep = "/"))
    pval.UF <- pval.UF[,c(colNAMES)]
    pval.UF$pval.NA <- 0
    pval.UF$pval.NA[which(is.na(pval.UF$p.value.Direct))] = 1
    pval.UF$p.value.Direct[which(is.na(pval.UF$p.value.Direct))] = 0  
    
    pval <- pval.UF
    
    pval$mm.DirectMINUSmm.IVT <- pval$mm.Direct - pval$mm.IVT
    pval$total_Direct <- pval$T_Direct + pval$C_Direct
    pval$total_IVT <- pval$T_IVT + pval$C_IVT
    pval<- pval[which(pval$total_IVT != 0),]
    
    pval$psi <- 0
    pval$psi[which(((-log10(pval$p.value.Direct)) > 2) & 
                     (pval$p.value.Direct < 0.001) & 
                     (pval$N_reads_Direct >= 10) & 
                     (pval$mm.IVT < 10))] <- 1
    
    pval$ACP <- paste(pval$Annotation, pval$chr, pval$position, sep = "")
    pval$true <- 1
    pval$which.IVT <- f
    
    
    ##--write to cell line list--##
    CL.pvals.list[[f]] <- pval
    
    pval.prefix <- pval %>% rename_with(~ paste0(f,".", .), -c(ACP, Annotation, chr, position, kmer))
    CL.pvals.list.prefix[[f]] <- pval.prefix
  } 
  
  #------------------------------------------------------------#
  ##--WRITE TO RData: Cell Line panIVT & pIVT MASTER--##
  print(paste("FINISH:",CL, "pval-IVT.master.RData", sep=" "))
  
  #out.rds.file <- paste(CL,".pval-IVT.master.RData", sep = "")
  #saveRDS(CL.pvals.list, print(paste(rdsPath, out.rds.file, sep = "/")))
  
  #------------------------------------------------------------#
  ##--unfiltered MASTER--##
  
  CL.ACP <- CL.pvals.list.prefix %>% reduce(full_join, by=stable)
  ACP.priority <- CL.ACP[,c("ACP", "pIVT.true", "panIVT.true")]
  
  ACP.priority[is.na(ACP.priority)] <- 0
  ACP.priority$panIVT.true[which(ACP.priority$pIVT.true == 1)] <- 0
  
  ACP.pIVT <- ACP.priority$ACP[which(ACP.priority$pIVT.true == 1)]
  ACP.panIVT <- ACP.priority$ACP[which(ACP.priority$panIVT.true == 1)]
  
  pIVT.priority <- filter(CL.pvals.list$pIVT, ACP %in% ACP.pIVT)      
  panIVT.priority <- filter(CL.pvals.list$panIVT, ACP %in% ACP.panIVT)      
  
  CL.df <- rbind(pIVT.priority, panIVT.priority)
  CL.pvals.df <- CL.df #pval filter
  
  CL.pvals.df$CP <- paste(CL.pvals.df$chr, CL.pvals.df$position, sep = "")
 
  CL.df <- CL.pvals.df %>% rename_with(~ paste0(CL,".", .), -c(ACP, Annotation, chr, position, kmer, CP))
  master.list[[CL]] <- CL.df
  
  
  #------------------------------------------------------------#
  ##--WRITE TO RData: Cell Line panIVT & pIVT MASTER--##
  
  out.pval.file <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")
  write.csv(CL.pvals.df,print(paste(pvalPath, out.pval.file, sep = "/")), row.names = F)

  #-LOOPING-#
  print(paste("FINISH:",CL, "pval-IVT.master.RData", sep=" "))
}


#------------------------------------------------------------------------------------------#
##--WRITE TO RData: MASTER--##

master.rds.file <- paste("sixCellLine",".MASTER.RData", sep = "")
saveRDS(master.list, print(paste(rdsPath, master.rds.file, sep = "/")))

#------------------------------------------------------------------------------------------#
print(" !!!!!!!-bhahahahahaha IT WORKED-!!!!!!!")



# DATA FRAME: MASTER  -----------------------------------------------------
##--READ RData: MASTER--##

master.rds.file <- paste("sixCellLine",".MASTER.RData", sep = "")

master.list <- readRDS(print(paste(rdsPath, master.rds.file, sep = "/")))


colNAMES <- c("Annotation", "chr", "position", "kmer", 
              "T_Direct", "C_Direct", "T_IVT", "C_IVT", 
              "N_reads_Direct", "N_reads_IVT",
              "expected.mm", "mm.Direct", "mm.IVT", "p.value.Direct")




#------------------------------------------------------------------------------------------#
##--MASTER PSI CSV--##
if (T) {
  
  stable <- c("ACP", "Annotation", "chr", "position", "kmer")
  CL.psi <- as.vector(outer(CellLine, c("psi"), paste, sep="."))
  
  
  ACP.master.UF <- master.list %>% reduce(full_join, by=stable)
  ACP.priority <- ACP.master.UF[,c("ACP", CL.psi)]
  ACP.priority[is.na(ACP.priority)] <- 0
  ACP.priority$rank <- rowSums(ACP.priority[CL.psi])
  ACP.priority <- ACP.priority[which(ACP.priority$rank > 0),]
  ACP.priority$true <- 1
  
  ACP.priority.map <- map(master.list, ~ .x %>% 
                            filter(ACP %in% ACP.priority$ACP))
  
  ACP.master <- ACP.priority.map %>% reduce(full_join, by=stable)
  ACP.master$rank <- rowSums(ACP.master[CL.psi])
  
  CL.ACP <- as.vector(outer(CellLine, c("T_Direct", "C_Direct", "N_reads_Direct", "N_reads_IVT","mm.Direct", "mm.IVT", "p.value.Direct","pval.NA", "psi", "true", "which.IVT"), paste, sep="."))
  
  ACP.master <- ACP.master[,c(stable, CL.ACP)]
  
  pval.file <- paste("sixCellLine",".pval-IVT.mmDirect.PSI.master-kinetic.csv", sep = "")
  write.csv(ACP.master,print(paste(pvalPath, pval.file, sep = "/")), row.names = F)
  
}



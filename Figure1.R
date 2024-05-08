## ========================================== ##
##             Psi Conservation               ##
## ========================================== ##

#@ author: C.A.McCormick
#@ date: 11/28/2023
#@ content: Figure 1
#@ version: 15.0

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
figPath <- print(paste(dirPath, "Figures", "Figure-1","FIG1-elements", sep = "/"))
bamPath <- print(paste(dirPath, "BAM", "hg38v10", sep = "/"))
tpmPath <- print(paste(bamPath, "stringtie",sep = "/"))
pysamPath <- print(paste(dataPath, "pysamstats", sep = "/"))
pvalSOURCE <- print(paste(dirPath,"pValues","DATA", sep = "/"))
kpPath <- print(paste(dataPath, "seqLOGO", sep = "/"))
bootPath <- print(paste(dirPath,"CODE","bootSTRAP", sep = "/"))
bsPath <- print(paste(bootPath, "pVALUES", sep="/"))
TPMbsPath <- print(paste(bootPath, "stringtie/DATA", sep = "/"))
orthoPATH <- print(paste(dataPath, "orthogonal", sep = "/"))


###--- Colors ---###

cp1 <- c("#FF595E", "#FF924C", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93")


cat("\f") #clears console

CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")
# FIGURE: Cell Line Global mm vs tpm --------------------------------------

##--Manual Check--##
# pval.file.name <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")

#------------------------------------------------------------------------------------------#
##--SCATTER: SPECIFIC LABELS--##
for (CL in CellLine) {
  #------------------------------------------------------------------------------------------#
  ##--SPECIFIC LABELS--##
  pval.file.name <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")
  
  out.file.tpm <- paste(CL,".pval.panIVT.log-TPMvMM-SCALE_HD.eps", sep = "")
  out.file.reads <- paste(CL,".pval.panIVT.log-nREADSvMM-SCALE_HD.eps", sep = "")
  
  pval.df <- read_csv(paste(pvalPath, pval.file.name, sep = "/"))
  # pval.df <- pval.df[which(pval.df$which.IVT == "pIVT"),]
  
  pval.df <- pval.df[which(pval.df$N_reads_Direct >= 10 &
                             pval.df$mm.IVT <= 10 &
                             pval.df$p.value.Direct < 0.001 &
                             # pval.df$which.IVT == "pIVT" &
                             pval.df$N_reads_IVT >= 10),]
  
  pval.pan <- pval.df
  

  pval.df$color <- "other"
  pval.df$color[pval.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  pval.df$color[pval.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  CL.plot.reads <- pval.df %>%
    ggplot(aes(x = log10(N_reads_Direct), y = mm.Direct, label = Annotation)) +
    geom_point(data = . %>% filter(color == "other"), color = "azure4", stroke = NA) +
    geom_point(data = . %>% filter(color == "pus7"), color = "magenta4", stroke = NA) +
    geom_point(data = . %>% filter(color == "trub1"), color = "springgreen4", stroke = NA) +
    geom_text(data=subset(pval.df,
                          CP == "chr2235424407" |
                          CP == "chr135603333" | 
                          CP =="chr101044099" | 
                          CP == "chr1469380270")) +
    #geom_text(aes(label = Annotation), hjust=0, vjust=0) +
    #xlim(1, 5) +
    #ylim(0, 110) +
    scale_x_continuous(
      #breaks = seq(0, 110, by = 5)
      breaks = c(1,3,5), 
      labels = c("1", "3", "5"), 
      limits = c(1, 5),
      #expand = c(0, 0.05)
    )+
    scale_y_continuous(
      #breaks = seq(0, 110, by = 5)
      breaks = c(0, 10, 25, 50, 75, 100), 
      labels = c("0", "10", "25", "50", "75", "100"), 
      limits = c(0, 100),
      expand = c(0, 0.05)
    ) +
    geom_hline(yintercept=10, color="red", linewidth=0.75, linetype=2) + 
    geom_hline(yintercept=c(25,50,75), color="black", linewidth=0.25) + 
    geom_hline(yintercept=100, color="black", linewidth=0.5) + 
    geom_vline(xintercept=log10(30), color="red", linewidth=0.75, linetype=2)+
    geom_vline(xintercept=5, color="black", linewidth=0.5)+
    #annotate("rect", xmin=c(1), xmax = c(log10(40)), ymin=c(0) , ymax=c(100))+
    #annotate("rect", xmin=c(log10(40)), xmax = c(5), ymin=c(0) , ymax=c(10))+
    theme(axis.line=element_line(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          #axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          #legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  CL.plot.reads
  ggsave(print(paste(figPath, out.file.reads, sep = "/")), width = 2, height = 4)
  
  x_hist <- ggplot(pval.df, aes(x = log10(N_reads_Direct))) +
    geom_histogram(binwidth = 0.2, fill = "gray", color = "black") +
    theme_minimal() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  x_hist
  
  y_hist <- ggplot(pval.df, aes(x = mm.Direct)) +
    geom_histogram(binwidth = 5, fill = "gray", color = "black") +
    coord_flip() +
    theme_minimal() +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  y_hist
}  

#------------------------------------------------------------------------------------------#
##--HISTOGRAMS--##

#-MERGED-#
for (CL in CellLine) {
  pval.file.name <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")
  
  out.file.histX <- paste(CL,".pval.panIVT.log-TPM-histX.eps", sep = "")
  out.file.histY <- paste(CL,".pval.panIVT.MM-histY.eps", sep = "")
  
  pval.df <- read_csv(paste(pvalPath, pval.file.name, sep = "/"))
  # pval.df <- pval.df[which(pval.df$which.IVT == "pIVT"),]
  
  pval.df <- pval.df[which(pval.df$N_reads_Direct >= 10 &
                             pval.df$mm.IVT <= 10 &
                             pval.df$p.value.Direct < 0.001 &
                             # pval.df$which.IVT == "pIVT" &
                             pval.df$N_reads_IVT >= 10),]
  
  pval.df$color <- "other"
  pval.df$color[pval.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  pval.df$color[pval.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"

  x_hist <- ggplot(pval.df, aes(x = log10(N_reads_Direct))) +
    geom_histogram(fill = "lightskyblue4", color = "lightskyblue2",binwidth = 0.2) +
    xlim(0,5)+
    theme(axis.line=element_line(),
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
  x_hist
  ggsave(print(paste(figPath, out.file.histX, sep = "/")), width = 2, height = 4)



  y_hist <- ggplot(pval.df,  aes(x = mm.Direct)) +
    geom_histogram(fill = "lightskyblue4", color = "lightskyblue2",binwidth = 3) +
    
    xlim(0,100)+
    #ylim(0,1000)+
    #coord_flip() +
    theme(axis.line=element_line(),
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
  y_hist
  ggsave(print(paste(figPath, out.file.histY, sep = "/")), width = 2, height = 4)
  
}  


#------------------------------------------------------------------------------------------#
#-MOTIF-#
for (CL in CellLine) {
  pval.file.name <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")
  
  out.file.histX <- paste(CL,".pval.pIVT_panIVT.log-TPM-histX.eps", sep = "")
  out.file.histY <- paste(CL,".pval.pIVT_panIVT.MM-histY.eps", sep = "")
  
  pval.df <- read_csv(paste(pvalPath, pval.file.name, sep = "/"))
  pval.df <- pval.df[which(pval.df$which.IVT == "pIVT"),]
  
  pval.df <- pval.df[which(pval.df$N_reads_Direct >= 10 &
                             pval.df$mm.IVT <= 10 &
                             pval.df$p.value.Direct < 0.001 &
                             pval.df$which.IVT == "pIVT" &
                             pval.df$N_reads_IVT >= 10),]
  
  pval.df$color <- "other"
  pval.df$color[pval.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  pval.df$color[pval.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  # x_hist <- ggplot(pval.df, aes(x = log10(N_reads_Direct))) +
  #   #geom_histogram(binwidth = 0.2, fill = "gray", color = "black") +
  #   geom_histogram(data = . %>% filter(color == "other"), 
  #                  fill = "lightskyblue4", color = "lightskyblue2",binwidth = 0.2) +
  #   geom_histogram(data = . %>% filter(color == "pus7"), fill = "magenta4", color = "magenta2", binwidth = 0.2) +
  #   geom_histogram(data = . %>% filter(color == "trub1"), fill = "springgreen4",  color = "springgreen2", binwidth = 0.2) +
  #   # scale_x_continuous(
  #   #   #breaks = seq(0, 110, by = 5)
  #   #   breaks = c(1,3,5), 
  #   #   labels = c("1", "3", "5"), 
  #   #   limits = c(0, 5),
  #   #   #expand = c(0, 0.05)
  #   # )+
  #   xlim(0,5)+
  #   theme(axis.line=element_line(),
  #         #axis.text.x=element_blank(),
  #         #axis.text.y=element_blank(),
  #         #axis.ticks=element_blank(),
  #         axis.title.x=element_blank(),
  #         axis.title.y=element_blank(),
  #         #legend.position="none",
  #         panel.background=element_blank(),
  #         panel.border=element_blank(),
  #         panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank(),
  #         plot.background=element_blank())
  # x_hist
  # ggsave(print(paste(figPath, out.file.histX, sep = "/")), width = 2, height = 4)
  # 
  # 
  
  y_hist <- ggplot(pval.df,  aes(x = mm.Direct)) +
    #geom_histogram(binwidth = 1, fill = "gray", color = "black") +
    geom_histogram(data = . %>% filter(color == "other"), 
                   fill = "lightskyblue4", color = "lightskyblue2",binwidth = 2.5) +
    geom_histogram(data = . %>% filter(color == "pus7"), fill = "magenta4", color = "magenta2", binwidth = 2.5) +
    geom_histogram(data = . %>% filter(color == "trub1"), fill = "springgreen4",  color = "springgreen2", binwidth = 2.5) +
    xlim(0,100)+
    #ylim(0,1000)+
    #coord_flip() +
    theme(axis.line=element_line(),
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
  y_hist
  ggsave(print(paste(figPath, out.file.histY, sep = "/")), width = 2, height = 4)
  
}  


  #------------------------------------------------------------------------------------------#
  ##--SPECIFIC SITES ONLY--##
for (CL in CellLine) {
  
  out.file.reads2 <- paste(CL,".LABEL.pval.pIVT_panIVT.log-nREADSvMM-SCALE_HD.eps", sep = "")
  
  pval.file.name <- paste(CL,".pval-IVT.PSI.master.csv", sep = "")

  pval.df <- read_csv(paste(pvalPath, pval.file.name, sep = "/"))
  # pval.df <- pval.df[which(pval.df$which.IVT == "pIVT"),]
  
  pval.df <- pval.df[which(pval.df$N_reads_Direct >= 10 &
                             pval.df$mm.IVT <= 10 &
                             pval.df$p.value.Direct < 0.001 &
                             # pval.df$which.IVT == "pIVT" &
                             pval.df$N_reads_IVT >= 10),]

  pval.df$color <- "other"
  pval.df$color[pval.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  pval.df$color[pval.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  subPATH <- print(paste(dataPath, "submissionDATA", sep = "/"))
  
  
  

  pval.df <- pval.df[which(pval.df$CP == "chr2235424407" | #MCM5
                             pval.df$CP == "chr135603333" | #PSMB2
                             pval.df$CP == "chr101044099" | #IDI1
                             pval.df$CP == "chr1469380270" #ERH
                           ),]



  CL.plot.reads <- pval.df %>%
    ggplot(aes(x = log10(N_reads_Direct), y = mm.Direct, label = Annotation)) +
    geom_point(data = . %>% filter(color == "other"), color = "azure4", stroke = NA) +
    geom_point(data = . %>% filter(color == "pus7"), color = "magenta4", stroke = NA) +
    geom_point(data = . %>% filter(color == "trub1"), color = "springgreen4", stroke = NA) +
    geom_text(data=subset(pval.df,
                          CP == "chr2235424407" |
                          CP == "chr135603333" |
                          CP =="chr101044099" |
                          CP == "chr1469380270")) +
    #geom_text(aes(label = Annotation), hjust=0, vjust=0) +
    scale_x_continuous(
      #breaks = seq(0, 110, by = 5)
      breaks = c(1,3,5), 
      labels = c("1", "3", "5"), 
      limits = c(1, 5),
      #expand = c(0, 0.05)
    )+
    scale_y_continuous(
      #breaks = seq(0, 110, by = 5)
      breaks = c(0, 10, 25, 50, 75, 100), 
      labels = c("0", "10", "25", "50", "75", "100"), 
      limits = c(0, 100),
      expand = c(0, 0.05)
    ) +
    geom_hline(yintercept=10, color="red", linewidth=0.75, linetype=2) + 
    geom_hline(yintercept=c(25,50,75), color="black", linewidth=0.25) + 
    geom_hline(yintercept=100, color="black", linewidth=0.5) + 
    geom_vline(xintercept=log10(30), color="red", linewidth=0.75, linetype=2)+
    geom_vline(xintercept=5, color="black", linewidth=0.5)+
    #annotate("rect", xmin=c(1), xmax = c(log10(40)), ymin=c(0) , ymax=c(100))+
    #annotate("rect", xmin=c(log10(40)), xmax = c(5), ymin=c(0) , ymax=c(10))+
    theme(axis.line=element_line(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          #axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          #legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  CL.plot.reads
  ggsave(print(paste(subPATH, out.file.reads2, sep = "/")), width = 2, height = 4)


}








# DATA FRAME: CL stats ----------------------------------------------------
CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")

CL.list <- list()
for (cl in CellLine) {
  UF.file.name <- paste0(cl,".pval-IVT.PSI.master.csv")
  UF.df <- read_csv(paste(pvalPath, UF.file.name, sep = "/"))
  
  U.df <- UF.df[which(UF.df$N_reads_IVT >= 10 & UF.df$N_reads_Direct >= 30),]
  U.df$color <- "other"
  U.df$color[U.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
  U.df$color[U.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
  
  U <- nrow(U.df)
  U.trub1 <- length(which(U.df$color == "trub1"))
  U.pus7 <- length(which(U.df$color == "pus7"))
  
  pval.df <- U.df[which(U.df$N_reads_Direct >= 30 &
                          U.df$mm.IVT <= 10 &
                          U.df$mm.Direct >= 10 &
                          U.df$p.value.Direct < 0.001 &
                          U.df$N_reads_IVT >= 10),]
  
  
  CL.temp <- data.frame(CL = cl, 
                        U.sites = U,
                        U.trub1 = U.trub1, 
                        U.pus7 = U.pus7,
                        psi.sites = nrow(pval.df), 
                        psi.U = ((nrow(pval.df)) / (U)),
                        TRUB1.sites = length(which(pval.df$color == "trub1")), 
                        TRUB1.U = ((length(which(pval.df$color == "trub1"))) / (U.trub1)) , 
                        TRUB1.psi = ((length(which(pval.df$color == "trub1"))) / (nrow(pval.df))),
                        PUS7.sites = length(which(pval.df$color == "pus7")), 
                        PUS7.U = ((length(which(pval.df$color == "pus7"))) / (U.pus7)), 
                        PUS7.psi = ((length(which(pval.df$color == "pus7"))) / (nrow(pval.df))))
  CL.list[[cl]] <- CL.temp
}


CL.stats <- data.table::rbindlist(CL.list)




CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")

CL.list <- list()

for (CL in CellLine) {
  
  tpm.file.name <- paste(CL,"_Direct_merged.hg38v10.abundance.tsv", sep="")
  TPM <- read_tsv(print(paste(tpmPath, tpm.file.name, sep = "/")))
  colnames(TPM) <- c("Gene_ID","gene_name","chr","Strand","Start","End", "Coverage", "FPKM", "TPM")
  
  temp.df <- data.frame(CL = CL, 
                        pus7.tpm = TPM$TPM[which(TPM$gene_name == "PUS7")],
                        trub1.tpm = TPM$TPM[which(TPM$gene_name == "TRUB1")])
  CL.list[[CL]] <- temp.df
}
CL.tpm <- data.table::rbindlist(CL.list)

CL.stats$PUS7.tpm <- CL.tpm$pus7.tpm[which(CL.stats$CL == CL.tpm$CL)]
CL.stats$TRUB1.tpm <- CL.tpm$trub1.tpm[which(CL.stats$CL == CL.tpm$CL)]


write.csv(CL.stats, print(paste(pvalPath,"CL.stats.csv", sep = "/")), row.names = F)

stats <- read_csv(print(paste(pvalPath,"CL.stats.csv", sep = "/")))

# BOOTSTRAP: Weighted PUS TPM prop ----------------------------------------

##--DATA FRAME: BS STATS--## 
 # Run on LINUX
if (T) {
  
  CellLine = c("A549", "HeLa", "HepG2", "Jurkat", "NTERA", "SHSY5Y")
  U.col <- c("Annotation", "chr", "position", "kmer", "T_Direct", "C_Direct", "mm.Direct", "mm.IVT",
             "p.value.Direct", "total_Direct", "total_IVT", "ACP", "which.IVT", "color")
  
  BSit <- as.vector(outer("BS", c(1:10), paste, sep=""))
  
  BS.stats.list <- list()
  
  for (bs in BSit) {
    
    CL.list <- list()
    for (CL in CellLine) {
      UF.file.name <- paste(CL,bs,"pval-IVT.PSI.master.csv", sep = ".")
      UF.df <- read_csv(paste(bsPath, UF.file.name, sep = "/"))
      
      U.df <- UF.df[which(UF.df$N_reads_IVT >= 10 & 
                            UF.df$N_reads_Direct >= 10 & 
                            # UF.df$which.IVT == "pIVT" & 
                            UF.df$mm.IVT <= 10),]
      U.df$color <- "other"
      U.df$color[U.df$kmer %in% c("TGTAG","TCTAG","TATAG","TTTAG","TGTAA","TCTAA","TATAA","TTTAA")] <- "pus7"
      U.df$color[U.df$kmer %in% c("GTTCA", "GTTCT", "GTTCC", "GTTCG")] <- "trub1"
      U.df <- U.df[which(U.df$color == "trub1" | U.df$color == "pus7"),c(U.col)]
      
      U.df$psi <- 0
      U.df$psi[which(U.df$mm.IVT <= 10 & 
                       U.df$mm.Direct >= 10 & 
                       U.df$p.value.Direct < 0.001
                       )] <- 1
      
      trub1.df <- U.df[which(U.df$color == "trub1"),]
      trub1.U <- sum(trub1.df$total_Direct)
      trub1.C <- sum(trub1.df$C_Direct[which(trub1.df$psi == 1)])
      pus7.df <- U.df[which(U.df$color == "pus7"),]
      pus7.U <- sum(pus7.df$total_Direct)
      pus7.C <- sum(pus7.df$C_Direct[which(pus7.df$psi == 1)])
      
      CL.temp <- data.frame("CL" = CL,
                            "BS"=bs, 
                            trub1.pos = nrow(trub1.df),
                            trub1.noMOD = nrow(trub1.df[which(trub1.df$psi != 1),]),
                            trub1.psi = nrow(trub1.df[which(trub1.df$psi == 1),]),
                            trub1.kmer = trub1.U,
                            trub1.C = trub1.C,
                            trub1.frac = trub1.C / trub1.U,
                            pus7.pos = nrow(pus7.df),
                            pus7.noMOD = nrow(pus7.df[which(pus7.df$psi != 1),]),
                            pus7.psi = nrow(pus7.df[which(pus7.df$psi == 1),]),
                            pus7.kmer = pus7.U,
                            pus7.C = pus7.C,
                            pus7.frac = pus7.C / pus7.U)
      CL.list[[CL]] <- CL.temp
    }
    
    CL.stats <- data.table::rbindlist(CL.list)
    
    BS.stats.list[[bs]] <- CL.stats
  }
  
  
  BS.stats <- data.table::rbindlist(BS.stats.list)
  write.csv(BS.stats, print(paste(pvalPath,"PUS.panIVT.stats-BS.csv", sep = "/")), row.names = F)
  
}




##--DATA FRAME: PLOT--##
if (T) {
  BS.stats <- read_csv(paste(pvalPath, "PUS.panIVT.stats-BS.csv", sep = "/"))
  
  #-switch back to binary-#
  
  BS.stats$trub1.frac <- BS.stats$trub1.psi/BS.stats$trub1.pos
  BS.stats$pus7.frac <- BS.stats$pus7.psi/BS.stats$pus7.pos

  
  
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
                            pus7.TPM = TPM$TPM[which(TPM$gene_name == "PUS7")],
                            trub1.TPM = TPM$TPM[which(TPM$gene_name == "TRUB1")])
      
      BS.list[[i]] <- BS.temp
    }
    
    tpm.temp <- data.table::rbindlist(BS.list)
    PUS.TPM.BS.list[[CL]] <- tpm.temp
  }
  
  PUS.TPM.BS <- data.table::rbindlist(PUS.TPM.BS.list)
  write.csv(PUS.TPM.BS, print(paste(pvalPath,"BS-PUStpm.csv", sep = "/")), row.names = F)
  
  
  PUS.TPM.BS$BS <- paste0("BS",PUS.TPM.BS$BS)
  
  PUS.TPM.BS.melt <- melt(PUS.TPM.BS, id = c("CL", "BS"))
  
  
  summary_df <- PUS.TPM.BS.melt %>%
    group_by(variable, CL) %>%
    summarize(mean_value = mean(value), sd_value = sd(value))
  
  BS.psi <- merge(BS.stats, PUS.TPM.BS, by = c("CL", "BS"))
  
}



##--SCATTER PLOTS W V+H ERROR BARS--##
if (T){
  
  cp1 <- c("#FF595E", "#FF924C", "#FFCA3A", "#8AC926", "#1982C4", "#6A4C93")
  
  BS.psi.pus7 <- BS.psi[,c("CL", "BS", "pus7.frac", "pus7.TPM")]
  summary_psi <- BS.psi.pus7 %>%
    group_by(CL) %>%
    summarize(mean.pus7 = mean(pus7.frac), sd.pus7 = sd(pus7.frac))
  summary_tpm <- BS.psi.pus7 %>%
    group_by(CL) %>%
    summarize(mean.tpm = mean(pus7.TPM), sd.tpm = sd(pus7.TPM))
  BS.pus7 <- merge(summary_psi, summary_tpm, by = c("CL"))
  
  #-individual cell lines-#
  model <- lm(mean.tpm ~ mean.pus7, data = BS.pus7)
  r_squared <- summary(model)$r.squared
  print(r_squared)
  
  
  a <- ggplot(BS.pus7, aes(x = round((mean.pus7 * (10^4)), digits = 3), y = mean.tpm, color = as.factor(CL))) +
    geom_point(size = 2) +
    scale_color_manual(values = cp1) + 
    scale_y_continuous(breaks = c(0, 30, 65), 
                       labels = c("0", "30", "65"), 
                       limits = c(0, 65),
                       expand = c(0, 0.05)) +
    theme(axis.line = element_line(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank()) +
    geom_errorbar(aes(ymin = mean.tpm - sd.tpm, ymax = mean.tpm + sd.tpm), width = 0.025) + 
    geom_errorbarh(aes(xmin = round(((mean.pus7 - sd.pus7) * (10^4)), digits = 3),
                       xmax = round(((mean.pus7 + sd.pus7) * (10^4)), digits = 3)), height = 0.05) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.25, linetype=2)
  print(a)
  ggsave(print(paste(figPath, "BS-PUS7.binary.prop-tpm.scatter.eps", sep = "/")),width = 3, height = 3)
  

  
  
  BS.psi.trub1 <- BS.psi[,c("CL", "BS", "trub1.frac", "trub1.TPM")]
  summary_psi <- BS.psi.trub1 %>%
    group_by(CL) %>%
    summarize(mean.trub1 = mean(trub1.frac), sd.trub1 = sd(trub1.frac))
  summary_tpm <- BS.psi.trub1 %>%
    group_by(CL) %>%
    summarize(mean.tpm = mean(trub1.TPM), sd.tpm = sd(trub1.TPM))
  BS.trub1 <- merge(summary_psi, summary_tpm, by = c("CL"))

  model <- lm(mean.tpm ~ mean.trub1, data = BS.trub1)
  r_squared <- summary(model)$r.squared
  print(r_squared)
  
  
  a <- ggplot(BS.trub1, aes(x = round((mean.trub1 * (10^4)), digits = 3), y = mean.tpm, color = as.factor(CL))) +
    geom_point(size = 2) +
    scale_color_manual(values = cp1) + 
    scale_y_continuous(breaks = c(0, 30, 65), 
                       labels = c("0", "30", "65"), 
                       limits = c(0, 65),
                       expand = c(0, 0.05)) +
    theme(axis.line = element_line(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank()) +
    geom_errorbar(aes(ymin = mean.tpm - sd.tpm, ymax = mean.tpm + sd.tpm), width = 0.025) + 
    geom_errorbarh(aes(xmin = round(((mean.trub1 - sd.trub1) * (10^4)), digits = 3),
                       xmax = round(((mean.trub1 + sd.trub1) * (10^4)), digits = 3)), height = 0.05) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.25, linetype=2)
  print(a)
  ggsave(print(paste(figPath, "BS-TRUB1.binary.prop-tpm.scatter.eps", sep = "/")),width = 3, height = 3)
  
  
  
 
}


##--LINEAR MODEL--## 


pus7.df <- BS.psi[,c("CL", "BS", "pus7.frac", "pus7.TPM")]
pus7.model <- lm(y ~ x, data = data)

#-individual cell lines-#
models <- pus7.df %>%
  group_by(CL) %>%
  do(model = lm(pus7.TPM ~ pus7.frac, data = .))
r_squared <- models %>%
  summarise(r_squared = summary(model)$r.squared)


# BOOTSTRAP: PUS TPM ------------------------------------------------------



PUS <- c("PUS1", "PUSL1", "PUS3", "PUS7", "PUS7L", "PUS10", "TRUB1", "TRUB2", "DKC1", "RPUSD1", "RPUSD2")

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
                            PUS1 = TPM$TPM[which(TPM$gene_name == "PUS1")],
                            PUSL1 = TPM$TPM[which(TPM$gene_name == "PUSL1")],
                            PUS3 = TPM$TPM[which(TPM$gene_name == "PUS3")],
                            PUS7 = TPM$TPM[which(TPM$gene_name == "PUS7")],
                            PUS7L = TPM$TPM[which(TPM$gene_name == "PUS7L")],
                            PUS10 = TPM$TPM[which(TPM$gene_name == "PUS10")],
                            TRUB1 = TPM$TPM[which(TPM$gene_name == "TRUB1")],
                            TRUB2 = TPM$TPM[which(TPM$gene_name == "TRUB2")],
                            DKC1 = TPM$TPM[which(TPM$gene_name == "DKC1")],
                            RPUSD1 = TPM$TPM[which(TPM$gene_name == "RPUSD1")],
                            RPUSD2 = TPM$TPM[which(TPM$gene_name == "RPUSD2")])
      
      BS.list[[i]] <- BS.temp
    }
    
    tpm.temp <- data.table::rbindlist(BS.list)
    PUS.TPM.BS.list[[CL]] <- tpm.temp
  }
  
  PUS.TPM.BS <- data.table::rbindlist(PUS.TPM.BS.list)
  # write.csv(PUS.TPM.BS, print(paste(pvalPath,"BS-PUStpm.csv", sep = "/")), row.names = F)
  
  
  
  PUS.TPM.BS$BS <- paste0("BS",PUS.TPM.BS$BS)
  
  PUS.TPM.BS.melt <- melt(PUS.TPM.BS, id = c("CL", "BS"))
  PUS.TPM.BS.melt <- PUS.TPM.BS.melt %>%
    rename("variable" = "PUS", 
           "value" = "TPM")
  
  summary_df <- PUS.TPM.BS.melt %>%
    group_by(PUS, CL) %>%
    summarize(mean_value = mean(TPM), sd_value = sd(TPM))
}


##--ANOVA: PUS--##

if (T) {
  #-ANOVA + TUKEY POST HOC-#
  split_data <- split(PUS.TPM.BS.melt, PUS.TPM.BS.melt$PUS)
  analyze_transcript <- function(data) {
    data$CL <- as.factor(data$CL)
    anova_result <- aov(TPM ~ CL, data = data)
    tukey_result <- TukeyHSD(anova_result)
    return(list(ANOVA = summary(anova_result), Tukey = tukey_result))
  }
  results <- lapply(split_data, analyze_transcript)
  
  anova.df <- data.frame(matrix(nrow = 1, 
                                ncol = length(names(unlist(results[[1]]$ANOVA[[1]])))), 
                         stringsAsFactors = FALSE)
  colnames(anova.df) <- names(unlist(results[[1]]$ANOVA[[1]]))
  
  for (i in seq_along(results)) {
    pus_value <- names(results)[i]
    values <- unlist(results[[i]]$ANOVA[[1]])
    anova.df[i, names(values)] <- values
    anova.df$PUS[i] <- pus_value
  }
  anova.df$p_adjusted <- p.adjust(anova.df$`Pr(>F)1`, method = "BH")
  row.names(anova.df) <- anova.df$PUS
  
  Tukey.df <- data.frame(matrix(nrow = length(PUS), ncol = 15))
  colnames(Tukey.df) <- rownames(results$PUS1$Tukey$CL)
  rownames(Tukey.df) <- PUS
  for (i in seq_along(results)) {
    pus_value <- names(results)[i]
    tmp.df <- as.data.frame(t(results[[i]]$Tukey$CL))
    values <- tmp.df[4,]
    Tukey.df[i,colnames(tmp.df)] <- values
  }
  
  
  
  
  
  
  
  
  
  
  PUS.sig.df <- Tukey.df
  PUS.sig.df$PUS.sig <- anova.df[intersect(rownames(PUS.sig.df), rownames(anova.df)),
                                 c("p_adjusted")]
  PUS.sig.df.CL <- as.data.frame(t(PUS.sig.df[, !(colnames(PUS.sig.df) == "PUS.sig")]))
  
  
  PUS.sig.stars <- apply(PUS.sig.df, 2, function(x) {
    ifelse(x >= 0.05, 0,
           ifelse(x < 0.05 & x >= 0.01, 1,
                  ifelse(x < 0.01 & x >= 0.001, 2, 3)))
  })
  PUS.sig.stars <- as.data.frame(PUS.sig.stars)
  
  
  #-WRITING-#
  
  # PUS # 
  print(paste(rownames(PUS.sig.df)[which.min(PUS.sig.df$PUS.sig)], 
              ":most significant differences in expression levels among the cell lines.")) 
  print(paste(rownames(PUS.sig.df)[which.max(PUS.sig.df$PUS.sig)], 
              ":least significant differences in expression levels among the cell lines."))
  
  # CELL LINE #
  min_values <- apply(PUS.sig.df.CL, 2, min)
  min_row_names <- apply(PUS.sig.df.CL, 2, function(col) {
    rownames(PUS.sig.df.CL)[which.min(col)]
  })
  for (i in seq_along(min_values)) {
    print(paste0(names(min_values)[i],": ", min_row_names[i], ", p-value = ",  min_values[i]))
    cat("\n")
  }
}


sh <- as.vector(outer(c("SHSY5Y"),CellLine, paste, sep="-"))
sh.stat <- PUS.sig.stars[c("DKC1"), c(sh[1:5])]
sh.stat <- PUS.sig.df[c("DKC1"), c(sh[1:5])]

  
  
  
##--PUS w/ POINTS + SD--##
if (T) {
  psi.syn.box <- PUS.TPM.BS.melt %>% 
    ggplot(aes(x = factor(variable), y = value, fill = CL)) + 
    geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.8), 
                  width = .5, color="grey") +
    geom_point(position = position_dodge(width = 0.8), size = 0.05, color="black") +
    scale_fill_manual(values = cp1) +
    theme(axis.line=element_line(linewidth = 0.25, colour = "grey", linetype=1),
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
  psi.syn.box
  ggsave(print(paste(figPath, "BS-PUS.TPM-bar.eps", sep = "/")), width = 6, height = 2)
}

##--ABSOLUTE PUS--##
if (T) {
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
  ggsave(print(paste(figPath, "BS-PUS.TPM.eps", sep = "/")))
}

##--TRUB1 & PUS7--##
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
                            PUS7 = TPM$TPM[which(TPM$gene_name == "PUS7")],
                            TRUB1 = TPM$TPM[which(TPM$gene_name == "TRUB1")])
      
      BS.list[[i]] <- BS.temp
    }
    
    tpm.temp <- data.table::rbindlist(BS.list)
    PUS.TPM.BS.list[[CL]] <- tpm.temp
  }
  
  PUS.TPM.BS <- data.table::rbindlist(PUS.TPM.BS.list)
  write.csv(PUS.TPM.BS, print(paste(pvalPath,"BS-PUStpm.csv", sep = "/")), row.names = F)
  
  
  PUS.TPM.BS$BS <- paste0("BS",PUS.TPM.BS$BS)
  
  PUS.TPM.BS.melt <- melt(PUS.TPM.BS, id = c("CL", "BS"))
  
  
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
  ggsave(print(paste(figPath, "BS-PUS.TRUB1-PUS7.TPM.eps", sep = "/")))
  
}






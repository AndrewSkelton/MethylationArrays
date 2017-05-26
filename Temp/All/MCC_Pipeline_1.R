e#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Read in Infinium 450K Methylations Array IDAT Files and setup pheno table  |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
setwd("~/2015Apr_450K_Reynard/Collaboration/All/")

source('http://bioconductor.org/biocLite.R')
# biocLite("minfi")

library(minfi)
library(DMRcate)
library(lumi)
library(limma)
library(ggplot2)
library(reshape2)
library(scales)
library(IlluminaHumanMethylation450kmanifest)
# library(doParallel)
# library(gridExtra)
# library(cowplot)
# library(ggbio)
# library(biovizBase)
library(GenomicFeatures)
library(biomaRt)
library(dplyr)
library(ggdendro)
library(readr)

pheno_in <- read.table("Pheno/pheno_7_qc3.txt",
                       row.names=1,
                       header=T,
                       sep="\t",
                       stringsAsFactors=F)

pheno_ex <- pheno_in[ (pheno_in$KL_Score == "-") |
                        (pheno_in$Discrepancy_KL == 1) |
                        (pheno_in$Discrepancy_Other == 1) |
                        (pheno_in$Discrepancy_Gender == 1) |
                        (pheno_in$Gender == "ND") |
                        (pheno_in$QC_Exclude_3 == 1), ]

pheno_in <- pheno_in[ (pheno_in$KL_Score != "-") &
                        (pheno_in$Discrepancy_KL != 1) &
                        (pheno_in$Discrepancy_Other != 1) &
                        (pheno_in$Discrepancy_Gender != 1) &
                        (pheno_in$Gender != "ND") &
                        (pheno_in$QC_Exclude_3 != 1), ]

pheno_in <- pheno_in[grep("Blood", 
                          pheno_in$SampleType, 
                          invert = T),]

# pheno_in <- pheno_in[grep("Newc", 
#                           pheno_in$Source),]

pheno_in$Basename <- paste0(getwd(), 
                            "/idat/",
                            rownames(pheno_in))



# write.csv(pheno_in, file="Analysis_6/pheno_in.csv")
# write.csv(pheno_ex, file="Analysis_6/pheno_ex.csv")
##'-----------------------------------------------------------------------------------------#



##'Read in Raw IDAT Files
##'-----------------------------------------------------------------------------------------#
raw.idat      <- read.metharray.exp(targets = pheno_in)
##'-----------------------------------------------------------------------------------------#



##'Detection P Value Filter
##' At least 50% of samples must have a detection P Val < 0.01
##'-----------------------------------------------------------------------------------------#
lumi.dpval        <- detectionP(raw.idat, type = "m+u")
lumi.failed       <- lumi.dpval > 0.01
lumi.dpval.remove <- names(which(rowMeans(lumi.failed)>0.5, TRUE))
##'-----------------------------------------------------------------------------------------#


# ctrl_add          <- getControlAddress(raw.idat)
# ctlWide           <- log2(getRed(raw.idat)[ctrl_add, ])
# ctlR              <- melt(ctlWide, 
#                           varnames = c("address", 
#                                        "sample"))
# ctlWide           <- log2(getGreen(raw.idat)[ctrl_add, ])
# ctlG              <- melt(ctlWide, 
#                           varnames = c("address", 
#                                        "sample"))
# ctl               <- rbind(cbind(channel = "Red",   ctlR), 
#                            cbind(channel = "Green", ctlG))
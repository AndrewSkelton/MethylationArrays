#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Read in Infinium 450K Methylations Array IDAT Files and setup pheno table  |
#-------------------------------------------------------------------------------------------#


##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
source('http://bioconductor.org/biocLite.R')
# biocLite("MethylAid")

library(ggplot2)
library(reshape2)
library(dplyr)
library(readr)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(limma)
library(biomaRt)
library(minfi)
library(openxlsx)
library(MethylAid)

setwd("~/2016Jul_450KSanDie_Reynard/")
##'-----------------------------------------------------------------------------------------#



##'Import Pheno table and manupulate
##'-----------------------------------------------------------------------------------------#
pheno_tmp              <- read_tsv("Pheno/GSE80071.txt", col_names = F) %>% as.data.frame
pheno_in               <- read_tsv("Pheno/GSE80071_GEO.txt") %>% 
                          as.data.frame %>% 
                          dplyr::select(geo_accession,organism_ch1) %>% 
                          left_join(pheno_tmp, by = c("geo_accession"="X1"))
pheno_in$Joint         <- sapply(pheno_in$organism_ch1,function(x){strsplit(x,"_")[[1]][2]})
pheno_in               <- pheno_in %>% 
                          mutate(Disease  = ifelse(grepl("RA",geo_accession),"RA","OA"),
                                 Basename = paste0(getwd(),"/Raw_Data/",X2),
                                 JD       = paste0(Disease,"_",Joint)) %>% 
                          filter(Joint != "NA", geo_accession != "RA_29")
##'-----------------------------------------------------------------------------------------#


##'Read in Raw IDAT Files
##'-----------------------------------------------------------------------------------------#
raw.idat               <- read.metharray.exp(targets = pheno_in)
##'-----------------------------------------------------------------------------------------#


##'Detection P Value Filter
##' At least 50% of samples must have a detection P Val < 0.01
##'-----------------------------------------------------------------------------------------#
lumi.dpval.remove      <- raw.idat %>% 
                          detectionP(type = "m+u") %>% 
                          as.data.frame %>% 
                          add_rownames("CpG") %>% 
                          mutate_each(funs(. > 0.01), -CpG) %>% 
                          mutate(RowAv = rowMeans(.[,-1]),
                                 Check = ifelse(RowAv > 0.5, "True", "False")) %>% 
                          dplyr::select(CpG,RowAv,Check) %>%
                          filter(Check == "True") %>% .[["CpG"]]
##'-----------------------------------------------------------------------------------------#



##'Normalise Data - Functional Normalisation (Minfi)
##'-----------------------------------------------------------------------------------------#
# lumi.norm.base         <- preprocessFunnorm(raw.idat,
#                                             nPCs    = 2,
#                                             sex     = pData(raw.idat)$Gender,
#                                             bgCorr  = T,
#                                             dyeCorr = T,
#                                             verbose = T)
set.seed(73)
lumi.norm.base         <- preprocessSWAN(raw.idat, verbose = T) %>% mapToGenome
##'-----------------------------------------------------------------------------------------#


##'Filter SNPs and Sex Chromosome Probes
##'-----------------------------------------------------------------------------------------#
annotation             <- lumi.norm.base %>% getAnnotation %>% as.data.frame
lumi.norm.filtered     <- lumi.norm.base %>% 
                          addSnpInfo %>% 
                          dropLociWithSnps(snps = c("SBE","CpG","Probe"), maf  = 0.05)
annotation             <- lumi.norm.filtered %>% getAnnotation %>% as.data.frame
lumi.norm.filtered     <- lumi.norm.filtered[-(match(lumi.dpval.remove, 
                                                     annotation$Name) %>% na.omit),]
annotation             <- lumi.norm.filtered %>% getAnnotation %>% as.data.frame
lumi.norm.filtered     <- lumi.norm.filtered[-grep("chrX|chrY", annotation$chr),]
annotation             <- lumi.norm.filtered %>% getAnnotation %>% as.data.frame
##'-----------------------------------------------------------------------------------------#


##'PCA Plot
##'-----------------------------------------------------------------------------------------#
pca          <- lumi.norm.filtered %>% getM %>% as.data.frame %>% t %>% prcomp
d            <- pca$x %>% as.data.frame %>% add_rownames("X2") %>% 
                left_join(pData(lumi.norm.filtered) %>% as.data.frame)
pcv          <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

gg <- ggplot(d,aes(PC1,PC2,colour = JD)) +
      geom_point(size = 4) +
      xlab(label=paste0("PC1 (", pcv[1], "%)")) +
      ylab(label=paste0("PC2 (", pcv[2], "%)")) +
      theme_bw() +
      ggtitle("PCA") +
      theme(axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15))

# png(filename =  "PCA_O5Ex.png",
#     height   =  1024,
#     width    =  1024,
#     units    =  "px")
print(gg)
# dev.off()
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation - Model Fit
##'-----------------------------------------------------------------------------------------#
Treat                  <- pData(lumi.norm.filtered)$JD %>% factor
design                 <- model.matrix(~0 + Treat)
colnames(design)       <- colnames(design) %>% gsub("Treat","",.)
cm                     <- makeContrasts(OA_hip-OA_knee, RA_hip-RA_knee, levels = colnames(design))
fitB                   <- lmFit(minfi::getBeta(lumi.norm.filtered), design) %>% contrasts.fit(cm) %>% eBayes
fitM                   <- lmFit(minfi::getM(lumi.norm.filtered),    design) %>% contrasts.fit(cm) %>% eBayes
##'-----------------------------------------------------------------------------------------#


##'Differential Methylation - Extract Results
##'-----------------------------------------------------------------------------------------#
contrasts              <- cm %>% colnames
names(contrasts)       <- fitM$coefficients %>% colnames
pVal                   <- 0.01
deltaB                 <- 0.1
i                      <- 1

wb                     <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Pheno_Data")
openxlsx::writeData(wb, "Pheno_Data", pData(lumi.norm.filtered))

for(i in 1:length(contrasts)){
  gene_list_B            <- topTable(fitB, coef = contrasts[i], number = Inf, sort.by = "P")
  gene_list_M            <- topTable(fitM, coef = contrasts[i], number = Inf, sort.by = "P",
                                     p.value    = pVal)
  if(nrow(gene_list_M) > 0){
    gene_list_out          <- gene_list_M %>% 
                              add_rownames("CpG") %>% 
                              left_join(gene_list_B %>% add_rownames("CpG"), by = c("CpG" = "CpG"))
    colnames(gene_list_out)<- colnames(gene_list_out) %>% gsub(".x",".M",.) %>% gsub(".y",".B",.)
    gene_list_out          <- gene_list_out %>% 
                              filter(abs(logFC.B) > deltaB) %>% 
                              dplyr::select(CpG,logFC.B,P.Value.M,adj.P.Val.M) %>% 
                              left_join(annotation, by = c("CpG" = "Name")) %>% 
                              as.data.frame
    message(paste0("Differential Methylation in ", names(contrasts)[i], ": ", nrow(gene_list_out)))
    
    # Write xlsx sheet
    openxlsx::addWorksheet(wb, contrasts[i])
    openxlsx::writeData(wb, contrasts[i], gene_list_out)
    betas_in  <- lumi.norm.filtered %>%
                 minfi::getBeta(.) %>%
                 as.data.frame %>%
                 add_rownames("CpG") %>%
                 dplyr::slice(match(gene_list_out$CpG, CpG))
    openxlsx::addWorksheet(wb, paste0("Betas ",contrasts[i]))
    openxlsx::writeData(wb, paste0("Betas ",contrasts[i]), betas_in)
    
    # Write Bed Track
    bed_nme <- contrasts[i] %>% gsub(" ", "",.)
    bed_out <- data.frame(chr = gene_list_out$chr) %>%
               mutate(start   = gene_list_out$pos - 1,
                      end     = gene_list_out$pos %>% as.numeric,
                      name    = paste0(bed_nme, "_DeltaBeta_",
                                       gene_list_out$logFC.B %>% round(3)),
                      score   = gene_list_out$logFC.B %>% round(3),
                      strand  = gene_list_out$strand %>% as.vector,
                      box_s   = start,
                      box_e   = end,
                      colour  = ifelse(score < 0, "0,43,255", "255,0,0"))
    write.table(bed_out,
                sep       = "\t",
                row.names = F,
                col.names = F,
                quote     = F,
                file      = paste0(bed_nme,"_SanDie.bed"))
    exe     <- paste0("sed -i '1s/^/track name=\"",
                      bed_nme, "\" description=\"",
                      bed_nme, "\" visibility=3 ",
                      "itemRgb=\"On\"",
                      "\\n/' ",
                      bed_nme, "_SanDie.bed")
    system(exe)
  } else {
    message(paste0("Differential Methylation in ", names(contrasts)[i], ": 0"))
  }
}
openxlsx::saveWorkbook(wb, "DM_450K_SanDie.xlsx", overwrite = T)

##'-----------------------------------------------------------------------------------------#


##'Differential Methylation - Plot top hits
##'-----------------------------------------------------------------------------------------#
CpG_in                 <- gene_list_out$CpG %>% head
# CpG_in                 <- gene_list_out %>% filter(logFC.B > 0) %>% .[["CpG"]] %>% .[1:5]
Beta_in                <- lumi.norm.filtered %>% 
                          minfi::getBeta(.) %>% 
                          as.data.frame %>% 
                          add_rownames("CpG") %>% 
                          filter(CpG %in% CpG_in) %>% 
                          melt %>% 
                          dplyr::rename(Sentrix_ID = variable, Beta = value) %>% 
                          left_join(pData(lumi.norm.filtered) %>% as.data.frame)
gg <- ggplot(Beta_in, aes(x      = TruePair,
                          y      = Beta,
                          fill   = Type,
                          group  = Type)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_bw() + 
      facet_grid(CpG~Age_Class, scales = "free_x")
print(gg)
png(filename =  "Beta_Plot_Top_Hits_Old_Young.png",
    height   =  1024,
    width    =  1024,
    units    =  "px")
print(gg)
dev.off()
##'-----------------------------------------------------------------------------------------#

bed_nme <- names(contrasts)[i] %>% gsub(" ", "",.)
bed_out <- data.frame(chr = gene_list_out$chr) %>% 
           mutate(start   = gene_list_out$pos - 1,
                  end     = gene_list_out$pos %>% as.numeric,
                  name    = paste0(bed_nme, "_DeltaBeta_", 
                                   gene_list_out$logFC.B %>% round(3)),
                  score   = gene_list_out$logFC.B %>% round(3),
                  box_s   = start,
                  box_e   = end,
                  colour  = ifelse(score < 0, "0,43,255", "255,0,0"))
write.table(bed_out, 
            sep       = "\t", 
            row.names = F, 
            col.names = F, 
            quote     = F,   
            file      = paste0(bed_nme,".bed"))
exe     <- paste0("sed -i '1s/^/track name=\"", 
                  bed_nme, "\" description=\"",
                  bed_nme, "\" visibility=3 ",
                  "itemRgb=\"On\"",
                  "\\n/' ", 
                  bed_nme, ".bed")
system(exe)


anno_in <- annotation[match(dm[[i]]$Name,
                            rownames(annotation)),]

out        <- c()
out$chr    <- dm[[i]]$chr
out        <- as.data.frame(out)
out$start  <- (dm[[i]]$pos - 1)
out$end    <- (dm[[i]]$pos)
out$name   <- paste0(dm[[i]]$Name, "_DeltaBeta_", dm[[i]]$`B-FC` %>% round(3))
out$score  <- as.numeric(dm[[i]]$`B-FC`)
out$strand <- dm[[i]]$strand
out$box_s  <- out$start
out$box_e  <- out$end
out$colour <- as.vector("255,0,0")
out[out$score < 0,]$colour <- "0,43,255" 

write.table(out, 
            sep       = "\t", 
            row.names = F, 
            col.names = F, 
            quote     = F,   
            file      = paste0(dm_name[i],".bed"))
exe <- paste0("sed -i '1s/^/track name=\"", 
              bed_name[i], "\" description=\"",
              bed_name[i], "\" visibility=3 ",
              "itemRgb=\"On\"",
              "\\n/' ", 
              dm_name[i], ".bed")
system(exe)
# `track name="Pre Vs Les OA Hip" description="Pre Vs Les OA Hip" visibility=3 itemRgb="On"`





# Write Bed Track
bed_nme <- "CpG_In"
bed_out <- data.frame(chr = annotation$chr) %>%
           mutate(start   = annotation$pos - 1 %>% as.vector,
                  end     = annotation$pos %>% as.vector,
                  name    = annotation$Name,
                  score   = 0,
                  strand  = annotation$strand %>% as.vector,
                  box_s   = start,
                  box_e   = end,
                  colour  = "0,0,0")
write.table(bed_out,
            sep       = "\t",
            row.names = F,
            col.names = F,
            quote     = F,
            file      = paste0(bed_nme,".bed"))
exe     <- paste0("sed -i '1s/^/track name=\"",
                  bed_nme, "\" description=\"",
                  bed_nme, "\" visibility=3 ",
                  "itemRgb=\"On\"",
                  "\\n/' ",
                  bed_nme, ".bed")
system(exe)

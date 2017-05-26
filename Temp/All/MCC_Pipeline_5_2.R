#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Differential Methylation - Paired                                          |
#-------------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/DM_Paired"))
setwd(paste0(current_dir,
             "/DM_Paired"))

lumi.norm.filtered.paired <- lumi.norm.filtered[,grep("ND",
                                                      pData(lumi.norm.filtered)$SampleType_Pair,
                                                      invert = T)]
lumi.norm.filtered.paired <- lumi.norm.filtered.paired[,grep("Knee",
                                                             pData(lumi.norm.filtered.paired)$SampleType_Pair)]

##'Differential Methylation - Limma - Paired Design
##'-----------------------------------------------------------------------------------------#
treatments           <- levels(factor(pData(lumi.norm.filtered.paired)$SampleType))
treatment_arrays     <- factor(pData(lumi.norm.filtered.paired)$SampleType)
levels(treatment_arrays) <- levels(treatment_arrays)[c(2,1,3,4)]
batch_arrays         <- factor(pData(lumi.norm.filtered.paired)$Source)
pairs                <- factor(pData(lumi.norm.filtered.paired)$Sample.ID)

design               <- model.matrix(~pairs + treatment_arrays)
design               <- design[,-c(24,25)]
# colnames(design)[24:26] <- treatments[-1]
# colnames(design)[1:length(treatments)] <- treatments

# design               <- model.matrix(~pairs + treatment_arrays)
# colnames(design)[24:26] <- treatments

#Beta or M Value Selection
fitB                  <- lmFit(minfi::getBeta(lumi.norm.filtered.paired), design)
fitM                  <- lmFit(minfi::getM(lumi.norm.filtered.paired),    design)

# cont_mat              <- makeContrasts(Lesioned_OA_Knee-Preserved_OA_Knee,
#                                        Lesioned_OA_Hip-Preserved_OA_Hip,
#                                        Preserved_OA_Knee-Preserved_OA_Hip,
#                                        Lesioned_OA_Knee-Lesioned_OA_Hip,
#                                        levels = colnames(design))

fit2B                 <- eBayes(fitB)
fit2M                 <- eBayes(fitM)
contrasts             <- 24 #ncol(design)
pVal                  <- 0.05
##

##'Adjusted Betas - Adjusting for pair covariate 
##'-----------------------------------------------------------------------------------------#
# mAdj.in     <- minfi::getM(lumi.norm.filtered.paired)
# mAdj.fit    <- fitM$coefficients[,-c(1:length(treatments),ncol(fitM$coefficients))]
# mAdj        <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments), ncol(design))])
# betaAdj_p   <- ilogit2(mAdj)
##'-----------------------------------------------------------------------------------------#


##
for(i in 1:length(contrasts)) {
  
  dir.create(paste0(current_dir,
                    "/DM_Paired/",
                    gsub(" ", 
                         "",
                         contrasts[i])))
  setwd(paste0(current_dir,
               "/DM_Paired/",
               gsub(" ", 
                    "",
                    contrasts[i])))
  
  gene_list_B             <- topTable(fit2B,
                                      coef          = contrasts,
                                      number        = Inf,
                                      adjust.method = "BH",
                                      sort.by       = "P")
  
  gene_list_M             <- topTable(fit2M,
                                      coef          = contrasts,
                                      p.value       = pVal,
                                      number        = Inf,
                                      adjust.method = "BH",
                                      sort.by       = "P")
  
  gene_list_B             <- gene_list_B[match(rownames(gene_list_M), 
                                               rownames(gene_list_B)),]
  
  gene_list_out           <- cbind(gene_list_M[,c(1,4,5)], 
                                   gene_list_B[,c(1,4,5)])
  colnames(gene_list_out) <- c("M-FC", "M-Pval", "M-Adj.Pval",
                               "B-FC", "B-Pval", "B-Adj.Pval")
  gene_list_out           <- cbind(gene_list_out,
                                   annotation[match(rownames(gene_list_out),
                                                    rownames(annotation)),])
  gene_list_out           <- gene_list_out[abs(gene_list_out$`B-FC`) > 0.1,]
  assign(paste0(gsub(" ", 
                     "",
                     contrasts[i]),
                "_Pair_",
                pVal), 
         gene_list_out)
  
  write.csv(gene_list_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "Pair_DM.csv"))
  
  beta_out <- betaAdj_p[match(rownames(gene_list_out)[1:100], 
                              rownames(betaAdj_p)),]
  
  
  colnames(beta_out) <- paste0(pData(lumi.norm.filtered.paired)$Sample.ID,
                               "_",
                               pData(lumi.norm.filtered.paired)$SampleType,
                               "_",
                               pData(lumi.norm.filtered.paired)$Source) 
  write.csv(beta_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "_Paired_AdjustedBetas_top100.csv"))
  
  cpgs_in           <- eval(parse(text = paste0("`", 
                                                gsub(" - ", 
                                                     "-",
                                                     contrasts[i]), 
                                                "_Pair_", 
                                                pVal,
                                                "`")))
  if(nrow(cpgs_in) < 10) {
    cpgs_in           <- rownames(cpgs_in)[1:nrow(cpgs_in)]
  } else {
    cpgs_in           <- rownames(cpgs_in)[1:10]
  }
  
  
  df                <- melt(betaAdj_p[match(cpgs_in, 
                                            rownames(betaAdj_p)),])
  pheno             <- pData(lumi.norm.filtered.paired)[match(df$Var2, 
                                                              rownames(pData(lumi.norm.filtered.paired))),]
  df                <- as.data.frame(cbind(df, pheno))
  rownames(df)      <- 1:nrow(df)
  colnames(df)[1:3] <- c("CpG", "SentrixID", "Beta") 
  
  # df <- df[df$SampleType == "Preserved_OA_Knee",]
  gg <- ggplot(df, aes(x=CpG, y=Beta))                            + 
    geom_boxplot(aes(fill      = SampleType,
                     colour    = SampleType),
                 notch         = T,
                 outlier.shape = NA,
                 alpha         = 0.4)                         +
    geom_point(aes(shape       = SampleType, 
                   fill        = SampleType),
               position        = position_dodge(width=.75), 
               pch             = 21, 
               size            = 2,
               alpha           = 0.4)                         +
    theme_bw()                                                +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
    scale_y_continuous(limits = c(0, 1))                      +
    ggtitle(paste0("Top 10 DM Probes \n",
                   contrasts[i])) #+
  #     geom_text(label       = df$Sample.ID, 
  #               size        = 4, 
  #               vjust       = 1.2, 
  #               hjust       = -0.2)   
  
  png(paste0("DM_Vis_",
             gsub(" ", 
                  "_", 
                  contrasts[i]),
             "_Pair.png"),
      width  = 12.53, 
      height = 6.98, 
      units  = "in", 
      res    = 600)
  print(gg)
  dev.off()
  
  df <- df[grep(paste0(paste0(strsplit(gsub(" - ","-",
                                            contrasts[i]), "-")[[1]][1],
                              "|",
                              strsplit(gsub(" - ","-",
                                            contrasts[i]), "-")[[1]][2])),
                df$SampleType),]
  
  gg <- ggplot(df, aes(x=CpG, y=Beta))                        + 
    geom_boxplot(aes(fill      = SampleType,
                     colour    = SampleType),
                 notch         = T,
                 outlier.shape = NA,
                 alpha         = 0.4)                         +
    geom_point(aes(shape       = SampleType, 
                   fill        = SampleType),
               position        = position_dodge(width=.75), 
               pch             = 21, 
               size            = 2,
               alpha           = 0.4)                         +
    theme_bw()                                                +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
    scale_y_continuous(limits = c(0, 1))                      +
    ggtitle(paste0("Top 10 DM Probes \n",
                   contrasts[i]))
  
  png(paste0("DM_Vis2_",
             gsub(" ", 
                  "_", 
                  contrasts[i]),
             "_Pair.png"),
      width  = 12.53, 
      height = 6.98, 
      units  = "in", 
      res    = 600)
  print(gg)
  dev.off()
  
  setwd("../")
}
##'-----------------------------------------------------------------------------------------#
setwd(current_dir)




















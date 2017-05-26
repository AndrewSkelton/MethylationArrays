#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Differential Methylation                                                   |
#-------------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/DM"))
setwd(paste0(current_dir,
             "/DM"))

##'Differential Methylation - Limma - Grouped Design
##'-----------------------------------------------------------------------------------------#
treatments           <- levels(factor(pData(lumi.norm.filtered)$SampleType))
treatment_arrays     <- factor(pData(lumi.norm.filtered)$SampleType)
batch_arrays         <- factor(pData(lumi.norm.filtered)$Source)

design               <- model.matrix(~0 + treatment_arrays + batch_arrays)
colnames(design)     <- c(treatments, paste0("be",
                                             levels(batch_arrays)[-1]))

#Beta or M Value Selection
fitB                  <- lmFit(minfi::getBeta(lumi.norm.filtered), design)
fitM                  <- lmFit(minfi::getM(lumi.norm.filtered),    design)

cont_mat              <- makeContrasts(Preserved_OA_Hip-Preserved_OA_Knee,
                                       Preserved_OA_Hip-Preserved_Control_Hip, 
                                       Preserved_OA_Hip-Lesioned_OA_Hip,
                                       Preserved_OA_Knee-Lesioned_OA_Knee,
                                       levels = colnames(design))

fit2B                 <- eBayes(contrasts.fit(fitB,
                                              contrasts=cont_mat))
fit2M                 <- eBayes(contrasts.fit(fitM,
                                              contrasts=cont_mat))

contrasts             <- colnames(cont_mat)
pVal                  <- 0.01
##


##'Adjusted Betas - Adjusting for batch covariate 
##'-----------------------------------------------------------------------------------------#
# betaAdj.in  <- getBeta(lumi.norm.filtered)
# betaAdj.fit <- fitB$coefficients[,-c(1:3)]
# betaAdj     <- as.matrix(betaAdj.in) - betaAdj.fit %*% t(design[,-c(1:3)])
mAdj.in     <- minfi::getM(lumi.norm.filtered)
mAdj.fit    <- fitM$coefficients[,-c(1:length(treatments))]
mAdj        <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments))])
betaAdj     <- ilogit2(mAdj)

# Cpgs_in              <- c("cg27219955","cg24569447","cg10151685","cg00512726","cg25936177",
#                           "cg16179969","cg24224304","cg07726953","cg17793497","cg21525369",
#                           "cg13061237","cg27629948","cg26581449")
# beta_query           <- betaAdj[match(Cpgs_in, rownames(betaAdj)),]
# colnames(beta_query) <- paste0(pData(lumi.norm.filtered)$Sample.ID,"_",
#                              pData(lumi.norm.filtered)$SampleType,"_",
#                              pData(lumi.norm.filtered)$Source)
# write.csv(beta_query, file = "Normalised_Betas_Query_April2016.csv")

region_in   <- c("chr12", 104774416, 105232067)
anno_in     <- annotation[(annotation$chr == region_in[1]) & 
                          (annotation$pos > as.numeric(region_in[2])) & 
                          (annotation$pos < as.numeric(region_in[3])),]
Cpgs_in_tar <- match(rownames(anno_in), rownames(betaAdj))
beta_query  <- betaAdj[Cpgs_in_tar,]


# annotation[grep("cg00608661", rownames(annotation)),]
# lumi.dpval.remove[grep("cg00608661", lumi.dpval.remove)]
# remove[grep("cg00608661", names(remove))]
# anno_tmp <- as.data.frame(getAnnotation(lumi.norm.base))
# anno_tmp[grep("cg00608661", rownames(anno_tmp)),]
beta_query           <- betaAdj[match("cg00608661", rownames(betaAdj)),] %>% as.data.frame %>% t
colnames(beta_query) <- paste0(pData(lumi.norm.filtered)$Sample.ID,"_",
                               pData(lumi.norm.filtered)$SampleType,"_",
                               pData(lumi.norm.filtered)$Source)
anno_in              <- annotation[match("cg00608661", rownames(annotation)),] %>% as.data.frame
write.csv(beta_query, file = "Normalised_Betas_cg00608661_18April2016.csv")
write.csv(anno_in,    file = "Annotation_cg00608661_18April2016.csv")
##'-----------------------------------------------------------------------------------------#

wb          <- openxlsx::createWorkbook()
for(i in 1:length(contrasts)) {
  
  # dir.create(paste0(current_dir,
  #                   "/DM/",
  #                   gsub(" ", 
  #                        "",
  #                        contrasts[i])))
  # setwd(paste0(current_dir,
  #              "/DM/",
  #              gsub(" ", 
  #                   "",
  #                   contrasts[i])))
  
  gene_list_B             <- topTable(fit2B,
                                      coef          = contrasts[i],
                                      number        = Inf,
                                      adjust.method = "BH",
                                      sort.by       = "P")
  
  gene_list_M             <- topTable(fit2M,
                                      coef          = contrasts[i],
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
  wsname                  <- contrasts[i] %>% 
                             gsub("Preserved", "Pre", .) %>% 
                             gsub("Lesioned", "Les", .) %>%
                             gsub(" - ", "-", .)
  openxlsx::addWorksheet(wb, wsname)
  openxlsx::writeData(wb, wsname, gene_list_out)
  assign(paste0(gsub(" - ", 
                     "_",
                     contrasts[i])), 
         gene_list_out)
}
saveWorkbook(wb, "Differential_Expression.xlsx", overwrite = T)

  
  
  cat(paste0(contrasts[i],
             ": ",
             nrow(gene_list_out),
             "\n"))
# }
  write.csv(gene_list_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "_DM.csv"))
  
  beta_out <- betaAdj[match(rownames(gene_list_out)[1:100], 
                            rownames(betaAdj)),]
  
  
  colnames(beta_out) <- paste0(pData(lumi.norm.filtered)$Sample.ID,
                               "_",
                               pData(lumi.norm.filtered)$SampleType,
                               "_",
                               pData(lumi.norm.filtered)$Source) 
  write.csv(beta_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "_AdjustedBetas_top100.csv"))
  
  cpgs_in           <- eval(parse(text = paste0("`", 
                                                gsub(" - ", 
                                                     "_",
                                                     contrasts[i]), 
                                                # "_", 
                                                # pVal,
                                                "`")))
  cpgs_in           <- rownames(cpgs_in)[1:10]
  
  df                <- melt(betaAdj[match(cpgs_in, 
                                          rownames(betaAdj)),])
  pheno             <- pData(lumi.norm.filtered)[match(df$Var2, 
                                                       rownames(pData(lumi.norm.filtered))),]
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
             ".png"),
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
             ".png"),
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


##'Plot Beta
##'-----------------------------------------------------------------------------------------#
for(i in 1:length(contrasts)) {
  cpgs_in           <- eval(parse(text = paste0("`", 
                                                gsub(" - ", 
                                                     "-",
                                                     contrasts[i]), 
                                                "_", 
                                                pVal,
                                                "`")))
  cpgs_in           <- rownames(cpgs_in)[1:10]
  
  df                <- melt(betaAdj[match(cpgs_in, 
                                              rownames(betaAdj)),])
  pheno             <- pData(lumi.norm.filtered)[match(df$Var2, 
                                                       rownames(pData(lumi.norm.filtered))),]
  df                <- as.data.frame(cbind(df, pheno))
  rownames(df)      <- 1:nrow(df)
  colnames(df)[1:3] <- c("CpG", "SentrixID", "Beta") 
  
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
                contrasts[i]))

  png(paste0("DM_Vis_",
             gsub(" ", 
                  "_", 
                  contrasts[i]),
             ".png"),
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
             ".png"),
      width  = 12.53, 
      height = 6.98, 
      units  = "in", 
      res    = 600)
  print(gg)
  dev.off()
}
##'-----------------------------------------------------------------------------------------#

# tmp.b <- c(0.3398919356,0.3024910872,0.2807809162,0.2792369636,0.2663079128)
# tmp   <- c("cg14063191", "cg25224992", "cg11601932", "cg14901671", "cg11978634")
# `OA_Hip-Hip_Control_0.01`[grep(paste(tmp, collapse="|"), rownames(`OA_Hip-Hip_Control_0.01`)),1:8]
search_in <- c("cg00739667",
               "cg00875989",
               "cg03507326",
               "cg10287137")
foo       <- betaAdj[grep(paste(search_in,
                                collapse="|"),
                          rownames(betaAdj)),]
write.csv(foo,
          file="CpG_Query_30_07_15.csv")


foo <- apply(betaAdj, 2, mean)
anno_27 <- annotation[annotation$Methyl27_Loci == T,]
anno_27_mat <- betaAdj[match(rownames(anno_27),
                             rownames(betaAdj)),]
foo2 <- apply(anno_27_mat, 2, mean)
out <- rbind(foo, foo2)
rownames(out) <- c("Filtered probes 450K",
                   "Filtered probes 27k")
colnames(out) <- paste0(pData(lumi.norm.filtered)$Sample.ID,
                             "_",
                             pData(lumi.norm.filtered)$SampleType,
                             "_",
                             pData(lumi.norm.filtered)$Source) 
write.csv(out, file="Mean_Methy.csv") 

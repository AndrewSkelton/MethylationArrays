#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Differential Methylation - Search                                          |
#-------------------------------------------------------------------------------------------#


current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/DM_Query"))
setwd(paste0(current_dir,
             "/DM_Query"))

clusters             <- read.table("../../Pheno/clusters.txt",
                                   sep    = "\t",
                                   header = T)
nof_subgroup         <- pData(lumi.norm.filtered)[grep("Preserved_Control_Hip",
                                                       pData(lumi.norm.filtered)$SampleType),]
clusters             <- rbind(clusters,
                              data.frame(chip.postion = paste0(nof_subgroup$Sentrix_Barcode, 
                                                                "_", 
                                                                nof_subgroup$Sentrix_Position), 
                                         cluster       = 3, 
                                         type          = "Preserved_Control_Hip"))
clusters$combo       <- paste0(clusters$type, 
                               clusters$cluster)
cluster_pos          <- na.omit(match(clusters$chip.postion,
                                        rownames(pData(lumi.norm.filtered))))

lumi.norm.sub        <- lumi.norm.filtered[,cluster_pos]
clusters             <- clusters[match(colnames(lumi.norm.sub),
                                       clusters$chip.postion),]
pData(lumi.norm.sub) <- cbind(pData(lumi.norm.sub),
                              clusters[match(clusters$chip.postion,
                                             colnames(lumi.norm.sub)),])

##'Differential Methylation - Limma - Grouped Design
##'-----------------------------------------------------------------------------------------#
treatments           <- levels(factor(clusters$combo))
treatment_arrays     <- factor(clusters$combo)
batch_arrays         <- factor(pData(lumi.norm.sub)$Source)

design               <- model.matrix(~0 + treatment_arrays + batch_arrays)
colnames(design)     <- c(treatments, paste0("be",
                                             levels(batch_arrays)[-1]))

#Beta or M Value Selection
fitB                  <- lmFit(minfi::getBeta(lumi.norm.sub), design)
fitM                  <- lmFit(minfi::getM(lumi.norm.sub),    design)

cont_mat              <- makeContrasts(hip    = OA_Hip1-OA_Hip2,
                                       knee   = OA_Knee1-OA_Knee2,
                                       all    = (((OA_Knee1+OA_Hip1)/2)-((OA_Knee2+OA_Hip2)/2)),
                                       levels = colnames(design))

fit2B                 <- eBayes(contrasts.fit(fitB,
                                              contrasts = cont_mat))
fit2M                 <- eBayes(contrasts.fit(fitM,
                                              contrasts = cont_mat))

contrasts             <- colnames(cont_mat)
pVal                  <- 0.01
##

mAdj.in     <- minfi::getM(lumi.norm.sub)
mAdj.fit    <- fitM$coefficients[,-c(1:length(treatments))]
mAdj        <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments))])
betaAdj     <- ilogit2(mAdj)

##
for(i in 1:length(contrasts)) {
  
  dir.create(paste0(current_dir,
                    "/DM_Query/",
                    gsub(" ", 
                         "",
                         contrasts[i])))
  setwd(paste0(current_dir,
               "/DM_Query/",
               gsub(" ", 
                    "",
                    contrasts[i])))
  
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
  assign(paste0(gsub(" ", 
                     "",
                     contrasts[i]),
                "_",
                pVal), 
         gene_list_out)
  
  write.csv(gene_list_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "_DM.csv"))
  
  beta_out <- betaAdj[match(rownames(gene_list_out)[1:100], 
                            rownames(betaAdj)),]
  
  
  colnames(beta_out) <- paste0(pData(lumi.norm.sub)$Sample.ID,
                               "_",
                               pData(lumi.norm.sub)$combo,
                               "_",
                               pData(lumi.norm.sub)$Source) 
  write.csv(beta_out, 
            file=paste0(gsub(" ",
                             "",
                             contrasts[i]),
                        "_AdjustedBetas_top100.csv"))
  
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
  pheno             <- pData(lumi.norm.sub)[match(df$Var2, 
                                                  rownames(pData(lumi.norm.sub))),]
  df                <- as.data.frame(cbind(df, pheno))
  rownames(df)      <- 1:nrow(df)
  colnames(df)[1:3] <- c("CpG", "SentrixID", "Beta") 
  
  # df <- df[df$SampleType == "Preserved_OA_Knee",]
  gg <- ggplot(df, aes(x=CpG, y=Beta))                            + 
    geom_boxplot(aes(fill      = combo,
                     colour    = combo),
                 notch         = T,
                 outlier.shape = NA,
                 alpha         = 0.4)                         +
    geom_point(aes(shape       = combo, 
                   fill        = combo),
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
  
#   gg <- ggplot(df, aes(x=CpG, y=Beta))                        + 
#     geom_boxplot(aes(fill      = combo,
#                      colour    = combo),
#                  notch         = T,
#                  outlier.shape = NA,
#                  alpha         = 0.4)                         +
#     geom_point(aes(shape       = combo, 
#                    fill        = combo),
#                position        = position_dodge(width=.75), 
#                pch             = 21, 
#                size            = 2,
#                alpha           = 0.4)                         +
#     theme_bw()                                                +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))  +
#     scale_y_continuous(limits = c(0, 1))                      +
#     ggtitle(paste0("Top 10 DM Probes \n",
#                    contrasts[i]))
#   
#   png(paste0("DM_Vis2_",
#              gsub(" ", 
#                   "_", 
#                   contrasts[i]),
#              ".png"),
#       width  = 12.53, 
#       height = 6.98, 
#       units  = "in", 
#       res    = 600)
#   print(gg)
#   dev.off()
  
  setwd("../")
}
##'-----------------------------------------------------------------------------------------#
setwd(current_dir)

probes_of_interest      <- c("cg18603538", "cg25730670", "cg14199090",
                             rownames(annotation[grep("LEP", annotation$UCSC_RefGene_Name),]))
betaAdj_query           <- betaAdj[match(probes_of_interest,
                                         rownames(betaAdj)),]
colnames(betaAdj_query) <- paste0(pData(lumi.norm.sub)$Sample.ID,
                                  "_",
                                  pData(lumi.norm.sub)$combo,
                                  "_",
                                  pData(lumi.norm.sub)$Source) 
anno_of_interest        <- annotation[match(probes_of_interest,
                                            rownames(annotation)),]
write.csv(anno_of_interest,
          file = "Query_Probes_Anno.csv")
write.csv(betaAdj_query, 
          file = "Query_Probes_Beta.csv")

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
    geom_boxplot(aes(fill      = combo,
                     colour    = combo),
                 notch         = T,
                 outlier.shape = NA,
                 alpha         = 0.4)                         +
    geom_point(aes(shape       = combo, 
                   fill        = combo),
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

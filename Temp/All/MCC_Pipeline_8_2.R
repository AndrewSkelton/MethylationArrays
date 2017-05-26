#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Linear Variable Covariate Model - Predicted Methylation Age                |
#-------------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/mAge-cAge_Regression"))
setwd(paste0(current_dir,
             "/mAge-cAge_Regression"))

##'Data Setup - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
for(i in 1:length(treatments)) {
# lvcm.sample <- treatments[i]
# lvcm.in     <- lumi.norm.filtered[,grep(lvcm.sample, 
#                                         pData(lumi.norm.filtered)$SampleType)]
lvcm.sample <- "Hip_All"
lvcm.in     <- lumi.norm.filtered[,grep("Lesioned_OA_Hip|Preserved_Control_Hip|Preserved_OA_Hip", 
                                        pData(lumi.norm.filtered)$SampleType)]
# lvcm.sample <- "Knee_All"
# lvcm.in     <- lumi.norm.filtered[,grep("Preserved_OA_Knee|Lesioned_OA_Knee", 
#                                         pData(lumi.norm.filtered)$SampleType)]
# lvcm.sample <- "OA_All"
# lvcm.in     <- lumi.norm.filtered[,grep("Preserved_OA_Knee|Lesioned_OA_Knee|Lesioned_OA_Hip|Preserved_OA_Hip", 
#                                         pData(lumi.norm.filtered)$SampleType)]
lvcm.sample <- "All"
lvcm.in     <- lumi.norm.filtered


lvcm.design <- model.matrix(~as.numeric(pData(lvcm.in)$mcAgeDiff))
lvcm.fit    <- lmFit(minfi::getM(lvcm.in), lvcm.design)
lvcm.fit2   <- eBayes(lvcm.fit)
##'-----------------------------------------------------------------------------------------#

dir.create(paste0(current_dir,
                  "/mAge-cAge_Regression/",
                  lvcm.sample))
setwd(paste0(current_dir,
             "/mAge-cAge_Regression/",
             lvcm.sample))

##'Model Filtering
##'-----------------------------------------------------------------------------------------#
lvcm.filter    <- topTable(fit           = lvcm.fit2,
                           coef          = 2,
                           p.value       = 0.01,
                           number        = Inf,
                           adjust.method = "none",
                           sort.by       = "P")

lvcm.filter    <- cbind(lvcm.filter,
                        annotation[match(rownames(lvcm.filter),
                                         rownames(annotation)),])
write.csv(lvcm.filter,
          file=paste0(lvcm.sample,
                      "_MethAgeVsChronAge_Correlation.csv"))

beta_sub       <- betaAdj[,match(colnames(lvcm.in),
                                 colnames(betaAdj))]

beta_out       <- beta_sub[match(rownames(lvcm.filter)[1:100], 
                                 rownames(beta_sub)),]
colnames(beta_out) <- paste0(pData(lvcm.in)$Sample.ID,
                             "_",
                             pData(lvcm.in)$SampleType,
                             "_",
                             pData(lvcm.in)$Source) 
write.csv(beta_out, 
          file=paste0(lvcm.sample,
                      "_MethAgeVsChronAge_AdjustedBetas_top100.csv"))

# lvcm.filter$afc <- abs(lvcm.filter$logFC)
# lvcm.filter.fc  <- lvcm.filter[with(lvcm.filter, order(-afc)), ]
##'-----------------------------------------------------------------------------------------#



##'Plot Regression
##'-----------------------------------------------------------------------------------------#
lvcm.cpgIn  <- rownames(lvcm.filter)[1:5]
lvcm.df     <- melt(beta_sub[lvcm.cpgIn,])

if(length(lvcm.cpgIn) > 1) {
  colnames(lvcm.df)         <- c("CpG", "Sample", "Beta")
} else {
  lvcm.df$CpG               <- lvcm.cpgIn
  lvcm.df$Sample            <- rownames(lvcm.df)
  rownames(lvcm.df)         <- c(1:nrow(lvcm.df))
  lvcm.df                   <- lvcm.df[,c(2,3,1)]
  colnames(lvcm.df)[3]      <- "Beta"
}

lvcm.df     <- cbind(lvcm.df,
                     pData(lvcm.in)[match(lvcm.df$Sample,
                                          rownames(pData(lvcm.in))),])
lvcm.df     <- as.data.frame(lvcm.df)
lvcm.df$mcAgeDiff <- as.numeric(lvcm.df$mcAgeDiff)

lvcm.gg      <- ggplot(lvcm.df, 
                       aes(x                = mcAgeDiff, 
                           y                = Beta, 
                           colour           = SampleType)) + 
                  geom_point() +
                  geom_smooth(method        = lm, 
                              se            = T, 
                              aes(group     = 1)) +
                  scale_x_continuous(limits = c(min(lvcm.df$mcAgeDiff),
                                                max(lvcm.df$mcAgeDiff))) +
                  theme_bw()

if(length(lvcm.cpgIn) > 1) {
  lvcm.gg <- lvcm.gg + 
    facet_grid(SampleType ~ CpG)
}

png(paste0(lvcm.sample,
           "_MethAgeVsChronAge_Correlations.png"),
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600) 
print(lvcm.gg)
dev.off()
setwd("../")
##'-----------------------------------------------------------------------------------------#
}
setwd(current_dir)




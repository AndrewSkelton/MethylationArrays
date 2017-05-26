#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Linear Variable Covariate Model                                            |
#-------------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/mcAgeRatio_Regression_Grouped"))
setwd(paste0(current_dir,
             "/mcAgeRatio_Regression_Grouped"))

##'Data Setup - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
# for(i in 1:length(levels(factor(pData(lumi.norm.filtered)$SampleType)))) {
# lvcm.sample <- levels(factor(pData(lumi.norm.filtered)$SampleType))[i]
# lvcm.in     <- lumi.norm.filtered[,grep(lvcm.sample, 
#                                         pData(lumi.norm.filtered)$SampleType)]

pheno_tmp                 <- read_csv("../../pData_norm.csv") %>% as.data.frame
pData(lumi.norm.filtered)$mcAgeRatio <- pheno_tmp$MethAge / pheno_tmp$Age

lvcm.sample.1 <- "Hip_All"
lvcm.in.1     <- lumi.norm.filtered[,grep("Lesioned_OA_Hip|Preserved_Control_Hip|Preserved_OA_Hip", 
                                        pData(lumi.norm.filtered)$SampleType)]
lvcm.sample.2 <- "Knee_All"
lvcm.in.2     <- lumi.norm.filtered[,grep("Preserved_OA_Knee|Lesioned_OA_Knee", 
                                        pData(lumi.norm.filtered)$SampleType)]
lvcm.sample.3 <- "OA_All"
lvcm.in.3     <- lumi.norm.filtered[,grep("Preserved_OA_Knee|Lesioned_OA_Knee|Lesioned_OA_Hip|Preserved_OA_Hip", 
                                        pData(lumi.norm.filtered)$SampleType)]
lvcm.sample.4 <- "All"
lvcm.in.4     <- lumi.norm.filtered

lvcm.sample.list <- list(lvcm.sample.1, lvcm.sample.2,
                         lvcm.sample.3, lvcm.sample.4)
lvcm.in.list     <- list(lvcm.in.1,     lvcm.in.2,
                         lvcm.in.3,     lvcm.in.4)
for(i in 1:length(lvcm.sample.list)) {
  lvcm.sample <- lvcm.sample.list[[i]]
  lvcm.in     <- lvcm.in.list[[i]]


# setwd(tmp)
# tmp <- getwd()

dir.create(paste0(current_dir,
                  "/mcAgeRatio_Regression_Grouped/",
                  lvcm.sample))
setwd(paste0(current_dir,
             "/mcAgeRatio_Regression_Grouped/",
             lvcm.sample))

if(length(levels(factor(pData(lvcm.in)$Source))) > 1) {
  lvcm.design <- model.matrix(~(as.numeric(pData(lvcm.in)$mcAgeRatio)) +
                                pData(lvcm.in)$Source)
} else {
  lvcm.design <- model.matrix(~as.numeric(pData(lvcm.in)$mcAgeRatio))
}
lvcm.fit    <- lmFit(minfi::getM(lvcm.in), lvcm.design)
lvcm.fit2   <- eBayes(lvcm.fit)
##'-----------------------------------------------------------------------------------------#

if(length(levels(factor(pData(lvcm.in)$Source))) > 1) {
  mAdj.in.age     <- minfi::getM(lvcm.in)
  mAdj.fit.age    <- lvcm.fit2$coefficients[,-c(1:2)]
  mAdj.age        <- as.matrix(mAdj.in.age) - mAdj.fit.age %*% t(lvcm.design[,-c(1:2)])
  betaAdj.age     <- ilogit2(mAdj.age)
} else {
  betaAdj.age     <- minfi::getBeta(lvcm.in)
}

# betaAdj.age <- betaAdj
##'Model Filtering
##'-----------------------------------------------------------------------------------------#
lvcm.filter    <- topTable(fit           = lvcm.fit2,
                           coef          = 2,
                           p.value       = 0.001,
                           number        = Inf,
                           adjust.method = "none",
                           sort.by       = "P")
# qval           <- qvalue(p = lvcm.filter$P.Value)
# lvcm.filter.q  <- cbind(lvcm.filter,
#                         qval$qvalues,
#                         qval$pvalues,
#                         qval$lfdr)
# lvcm.filter.f  <- lvcm.filter.q[lvcm.filter.q$`qval$lfdr` < 0.01,]
# summary(qval)
cat(paste0(lvcm.sample,
           ": ",
           nrow(lvcm.filter),
           "\n")) %>% message
# }

lvcm.filter    <- cbind(lvcm.filter,
                        annotation[match(rownames(lvcm.filter),
                                         rownames(annotation)),])
write.csv(lvcm.filter,
          file=paste0(lvcm.sample,
                      "_P0.001_mcAgeRatio_Correlation.csv"))
# "_PAdj0.05_Age_Correlation.csv"))

# betaAdj_subset     <- betaAdj.age[,grep(lvcm.sample, 
#                                         pData(lumi.norm.filtered)$SampleType)]
betaAdj_subset     <- betaAdj.age

beta_out           <- betaAdj_subset[match(rownames(lvcm.filter)[1:nrow(lvcm.filter)], 
                                           rownames(betaAdj.age)),]
if(!is.null(nrow(beta_out))) {
  colnames(beta_out) <- paste0(pData(lvcm.in)$Sample.ID,
                             "_",
                             pData(lvcm.in)$SampleType,
                             "_",
                             pData(lvcm.in)$Source) 

  write.csv(beta_out, 
            file=paste0(lvcm.sample,
                        "_P0.001_mcAgeRatio_AdjustedBetas.csv"))
  # "_PAdj0.05_Age_AdjustedBetas_top100.csv"))
}
  # lvcm.filter$afc <- abs(lvcm.filter$logFC)
  # lvcm.filter.fc  <- lvcm.filter[with(lvcm.filter, order(-afc)), ]
  ##'-----------------------------------------------------------------------------------------#
  
  
  
  ##'Plot Regression
  ##'-----------------------------------------------------------------------------------------#
  if(nrow(lvcm.filter) > 10) {
    lvcm.cpgIn  <- rownames(lvcm.filter)[1:10]
  } else {
    lvcm.cpgIn  <- rownames(lvcm.filter)[1:nrow(lvcm.filter)]
  }
  lvcm.df     <- melt(betaAdj_subset[lvcm.cpgIn,])
  
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
  lvcm.df$Age <- as.numeric(lvcm.df$Age)
  
  lvcm.gg      <- ggplot(lvcm.df, 
                         aes(x                = mcAgeRatio, 
                             y                = Beta, 
                             colour           = SampleType)) + 
                    geom_point() +
                    geom_smooth(method        = lm, 
                                se            = T, 
                                aes(group     = 1)) +
                    scale_y_continuous(limits = c(0,1)) +
                    # scale_x_continuous(limits = c(0,100)) +
  #                   scale_x_continuous(limits = c(min(lvcm.df$Age),
  #                                                 max(lvcm.df$Age))) +
                    theme_bw()
  
  if(length(lvcm.cpgIn) > 1) {
    lvcm.gg <- lvcm.gg + 
               # facet_grid(SampleType ~ CpG)
               facet_grid(. ~ CpG)
  }
  
  png(paste0(lvcm.sample,
             "_mcAgeRatio_Correlations.png"),
      width  = 12.53, 
      height = 6.98, 
      units  = "in", 
      res    = 600) 
    print(lvcm.gg)
  dev.off()
}
##'-----------------------------------------------------------------------------------------#


setwd(current_dir)



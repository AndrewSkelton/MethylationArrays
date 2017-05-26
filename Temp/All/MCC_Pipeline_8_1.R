#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Linear Variable Covariate Model - A Tad more Complex -mAge Regression      |
#-------------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/mAge_Regression_Advanced_Model"))
setwd(paste0(current_dir,
             "/mAge_Regression_Advanced_Model"))

##'Data Setup - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#

lvcm.sample <- "All"
lvcm.in     <- lumi.norm.filtered

SampleType  <- factor(pData(lvcm.in)$SampleType)
mAge        <- as.numeric(pData(lvcm.in)$MethAge)
Data_Source <- factor(pData(lvcm.in)$Source)
lvcm.design <- model.matrix(~0 + SampleType:mAge + Data_Source)
colnames(lvcm.design)
lvcm.fit    <- lmFit(minfi::getM(lvcm.in), lvcm.design)
lvcm.fit2   <- eBayes(lvcm.fit)
##'-----------------------------------------------------------------------------------------#

# if(length(levels(factor(pData(lvcm.in)$Source))) > 1) {
#   mAdj.in.age     <- minfi::getM(lvcm.in)
#   mAdj.fit.age    <- lvcm.fit2$coefficients[,-c(1:2)]
#   mAdj.age        <- as.matrix(mAdj.in.age) - mAdj.fit.age %*% t(lvcm.design[,-c(1:2)])
#   betaAdj.age     <- ilogit2(mAdj.age)
# } else {
#   betaAdj.age     <- minfi::getBeta(lvcm.in)
# }
mAdj.in.age     <- minfi::getM(lvcm.in)
mAdj.fit.age    <- lvcm.fit2$coefficients[,-c(6:10)]
mAdj.age        <- as.matrix(mAdj.in.age) - mAdj.fit.age %*% t(lvcm.design[,-c(6:10)])
betaAdj.age     <- ilogit2(mAdj.age)

##'Model Filtering
##'-----------------------------------------------------------------------------------------#

for(i in c(6:10)) {
  lvcm.filter    <- topTable(fit           = lvcm.fit2,
                             coef          = i,
                             p.value       = 0.01,
                             number        = Inf,
                             adjust.method = "BH",
                             sort.by       = "P")
  
  lvcm.name      <- colnames(lvcm.design)[i]
  lvcm.name      <- strsplit(lvcm.name, ":")[[1]][1]
  lvcm.name      <- strsplit(lvcm.name, "SampleType")[[1]][2]
  
#   cat(paste0(lvcm.name,
#              ": ",
#              nrow(lvcm.filter),
#              "\n"))
# }
  
  dir.create(paste0(current_dir,
                    "/mAge_Regression_Advanced_Model/",
                    lvcm.name))
  setwd(paste0(current_dir,
               "/mAge_Regression_Advanced_Model/",
               lvcm.name))
  
  lvcm.filter    <- cbind(lvcm.filter,
                          annotation[match(rownames(lvcm.filter),
                                           rownames(annotation)),])
  write.csv(lvcm.filter,
            file=paste0(lvcm.name,
                        "_Age_Correlation.csv"))
  
  betaAdj_subset     <- betaAdj.age
  beta_out           <- betaAdj_subset[match(rownames(lvcm.filter)[1:nrow(lvcm.filter)], 
                                             rownames(betaAdj_subset)),]
  if(!is.null(nrow(beta_out))) {
    colnames(beta_out) <- paste0(pData(lvcm.in)$Sample.ID,
                                 "_",
                                 pData(lvcm.in)$SampleType,
                                 "_",
                                 pData(lvcm.in)$Source) 
    
    write.csv(beta_out, 
              file=paste0(lvcm.name,
                          "_mAge_AdjustedBetas_top100.csv"))
    ##'-----------------------------------------------------------------------------------------#
    
    
    ##'Plot Regression
    ##'-----------------------------------------------------------------------------------------#
    if(nrow(lvcm.filter) > 5) {
      lvcm.cpgIn  <- rownames(lvcm.filter)[1:5]
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
                           aes(x                = Age, 
                               y                = Beta, 
                               colour           = SampleType)) + 
      geom_point() +
      geom_smooth(method        = lm, 
                  se            = T, 
                  aes(group     = 1)) +
      scale_y_continuous(limits = c(0,1)) +
      scale_x_continuous(limits = c(0,100)) +
      #                   scale_x_continuous(limits = c(min(lvcm.df$Age),
      #                                                 max(lvcm.df$Age))) +
      theme_bw()
    
    if(length(lvcm.cpgIn) > 1) {
      lvcm.gg <- lvcm.gg + 
        facet_grid(SampleType ~ CpG)
      # facet_grid(. ~ CpG)
    }
    
    png(paste0(lvcm.name,
               "_mAge_Correlations.png"),
        width  = 12.53, 
        height = 6.98, 
        units  = "in", 
        res    = 600) 
    print(lvcm.gg)
    dev.off()
  }
  ##'-----------------------------------------------------------------------------------------#
}

setwd(current_dir)



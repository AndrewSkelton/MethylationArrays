##'Data Setup - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
lvcm.sample <- "All"
lvcm.in     <- lumi.norm.filtered

SampleType  <- factor(pData(lvcm.in)$SampleType)
Age         <- as.numeric(pData(lvcm.in)$Age)
Data_Source <- factor(pData(lvcm.in)$Source)
# lvcm.design <- model.matrix(~0 + SampleType*Age + Data_Source)
# colnames(lvcm.design)[11:14] <- levels(SampleType)[2:5]

# lvcm.design <- model.matrix(~0 + SampleType + Data_Source + SampleType:Age)
# colnames(lvcm.design)[10:14] <- levels(SampleType)

lvcm.design <- model.matrix(~0 + SampleType + SampleType:Age + Data_Source)
colnames(lvcm.design)[10:14] <- levels(SampleType)

lvcm.fit    <- lmFit(minfi::getM(lvcm.in), lvcm.design)

# cont_mat    <- makeContrasts(all       = ((Lesioned_OA_Hip+Lesioned_OA_Knee+Preserved_Control_Hip+Preserved_OA_Hip+Preserved_OA_Knee)/5),
#                              oa_only   = ((Lesioned_OA_Hip+Lesioned_OA_Knee+Preserved_OA_Hip+Preserved_OA_Knee)/4),
#                              hip_only  = ((Lesioned_OA_Hip+Preserved_Control_Hip+Preserved_OA_Hip)/3),
#                              knee_only = ((Lesioned_OA_Knee+Preserved_OA_Knee)/2),
#                              levels    = colnames(lvcm.design))
# lvcm.fit.con<- contrasts.fit(lvcm.fit, contrasts = cont_mat)
# lvcm.fit2   <- eBayes(lvcm.fit.con)

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
mAdj.fit.age    <- lvcm.fit$coefficients[,-c(11:14)]
mAdj.age        <- as.matrix(mAdj.in.age) - mAdj.fit.age %*% t(lvcm.design[,-c(11:14)])
betaAdj.age     <- ilogit2(mAdj.age)
# mAdj.in.age     <- minfi::getM(lvcm.in)
# mAdj.fit.age    <- lvcm.fit$coefficients[,-c(10:14)]
# mAdj.age        <- as.matrix(mAdj.in.age) - mAdj.fit.age %*% t(lvcm.design[,-c(10:14)])
# betaAdj.age     <- ilogit2(mAdj.age)
#   
##'Model Filtering
##'-----------------------------------------------------------------------------------------#

for(i in c(10:14)) {
  # for(i in c(10:14)) {
  # for(i in colnames(cont_mat)) {
  lvcm.filter    <- topTable(fit           = lvcm.fit2,
                             coef          = i,
                             p.value       = 0.05,
                             number        = Inf,
                             adjust.method = "BH",
                             sort.by       = "P")
  
  lvcm.filter    <- lvcm.filter[with(lvcm.filter, order(-logFC)), ]
  
  lvcm.name      <- colnames(lvcm.design)[i]
  #   lvcm.name      <- strsplit(lvcm.name, ":")[[1]][1]
  #   lvcm.name      <- strsplit(lvcm.name, "SampleType")[[1]][2]
  
  # lvcm.name <- i
  cat(paste0(lvcm.name,
             ": ",
             nrow(lvcm.filter),
             "\n"))
}

  assign(paste0(paste0("InteractionModel_cAge", 
                     "_",
                     colnames(lvcm.design)[i])), 
         lvcm.filter)
}

#   assign(lvcm.name, 
#          lvcm.filter)
# }
save(all, oa_only, hip_only, knee_only, file = "Age_Collapse_Lists.RData")



dir.create(paste0(current_dir,
                  "/cAge_Regression_Advanced_Model/",
                  lvcm.name))
setwd(paste0(current_dir,
             "/cAge_Regression_Advanced_Model/",
             lvcm.name))

lvcm.filter    <- cbind(lvcm.filter,
                        annotation[match(rownames(lvcm.filter),
                                         rownames(annotation)),])
write.csv(lvcm.filter,
          file=paste0(lvcm.name,
                      "_Age_Correlation.csv"))

betaAdj_subset     <- betaAdj
beta_out           <- betaAdj_subset[match(rownames(lvcm.filter)[1:nrow(lvcm.filter)], 
                                           rownames(betaAdj_subset)),]
if(!is.null(nrow(beta_out))) {
  colnames(beta_out) <- paste0(pData(lvcm.in)$Sample.ID,
                               "_",
                               pData(lvcm.in)$SampleType,
                               "_",
                               pData(lvcm.in)$Source) 
}    
write.csv(beta_out, 
          file=paste0(lvcm.name,
                      "_Age_AdjustedBetas_top100.csv"))
##'-----------------------------------------------------------------------------------------#


##'Plot Regression
##'-----------------------------------------------------------------------------------------#
if(nrow(lvcm.filter) > 5) {
  lvcm.cpgIn  <- rownames(lvcm.filter)[1:6]
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
  # facet_grid(CpG ~ .)
} else {
  lvcm.gg <- lvcm.gg + 
    facet_grid(SampleType ~ .)
}

png(paste0(lvcm.name,
           "_Age_Correlations.png"),
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600) 
print(lvcm.gg)
dev.off()
}
##'-----------------------------------------------------------------------------------------#
}
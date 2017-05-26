#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Fitting Splines for Age                                                    |
#-------------------------------------------------------------------------------------------#

library(splines)

##'Data Setup - Fit the Model to a Basis Matrix for Natural Cubic Splines - Using 5
##'             Degrees of Freedom
##'-----------------------------------------------------------------------------------------#
for(i in 1:length(treatments)) {
lvcm.sample <- treatments[i]
lvcm.in     <- lumi.norm.filtered[,grep(lvcm.sample, 
                                        pData(lumi.norm.filtered)$SampleType)]

lvcm.x      <- ns(as.numeric(pData(lvcm.in)$Age), df=5)

if(length(levels(factor(pData(lvcm.in)$Source))) > 1) {
  lvcm.design <- model.matrix(~lvcm.x +
                               pData(lvcm.in)$Source)
} else {
  lvcm.design <- model.matrix(~lvcm.x)
}

lvcm.fit    <- lmFit(minfi::getM(lvcm.in), lvcm.design)
lvcm.fit2   <- eBayes(lvcm.fit)

lvcm.filter <- topTable(fit           = lvcm.fit2,
                        coef          = 2:6,
                        p.value       = 0.0001,
                        number        = Inf,
                        adjust.method = "none",
                        sort.by       = "F")
##'-----------------------------------------------------------------------------------------#



##'Data Setup - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
lvcm.cpgIn  <- rownames(lvcm.filter)[21:25]
lvcm.df     <- melt(betaAdj[lvcm.cpgIn,
                            grep(lvcm.sample,
                                 pData(lumi.norm.filtered)$SampleType)])

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
                           colour           = SampleType,
                           group            = CpG)) + 
                  geom_point() +
                  theme_bw() +
                  stat_smooth(method = "lm", 
                              formula = y ~ ns(x, df=5)) +
                  scale_x_continuous(limits = c(min(lvcm.df$Age),
                                                max(lvcm.df$Age)))
if(length(lvcm.cpgIn) > 1) {
  lvcm.gg <- lvcm.gg + 
    facet_grid(SampleType ~ CpG)
}

png(paste0(lvcm.sample,
           "_Splines.png"),
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600) 
print(lvcm.gg)
dev.off()
##'-----------------------------------------------------------------------------------------#
}
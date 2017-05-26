#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Post Normalisation QC Metrics                                              |
#-------------------------------------------------------------------------------------------#

treatments           <- levels(factor(pData(lumi.norm.filtered)$SampleType))
treatment_arrays     <- factor(pData(lumi.norm.filtered)$SampleType)
batch_arrays         <- factor(pData(lumi.norm.filtered)$Source)
design               <- model.matrix(~0 + treatment_arrays + batch_arrays)
colnames(design)     <- c(treatments, paste0("be",
                                             levels(batch_arrays)[-1]))
fitA                 <- lmFit(minfi::getM(lumi.norm.filtered), design)
mAdj.in              <- minfi::getM(lumi.norm.filtered)
mAdj.fit             <- fitA$coefficients[,-c(1:length(treatments))]
mAdj                 <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments))])
betaAdj              <- ilogit2(mAdj)

##'PCA
##'-----------------------------------------------------------------------------------------#
#lumi.norm.base
pca          <- prcomp(t(minfi::getM(lumi.norm.filtered)))



pca_in       <- minfi::getM(preprocessRaw(raw.idat))
pca_in       <- minfi::getM(lumi.norm.filtered)
# pca_in       <- mAdj

remove              <- which(pca_in==-Inf, arr.ind=TRUE)[,1]
remove              <- c(remove, which(pca_in==Inf, arr.ind=TRUE)[,1])
# remove              <- c(remove, which(pca_in==NaN, arr.ind=TRUE)[,1])
remove              <- c(remove, which(is.na(pca_in), arr.ind=TRUE)[,1])
if(length(remove) > 0) {
  pca_in    <- pca_in[-remove,]
}
pca          <- prcomp(t(pca_in))

# pca          <- prcomp(t(mAdj))

d            <- as.data.frame(pca$x)
d            <- as.data.frame(cbind(d, pData(lumi.norm.base)))
pcv          <- round((pca$sdev)^2 / sum(pca$sdev^2)*100, 2)

for(i in 1:1) {
  for(k in c("Source", "SampleType")) {
    j  <- i+1
    gg <- ggplot(d, aes_string(x         = paste0("PC",i), 
                               y         = paste0("PC",j)))       + 
      geom_point(aes_string(colour = paste0(k)), 
                 size   = 3)                     + 
      #         geom_text(label       = d$SampleType, 
      #                   size        = 4, 
      #                   vjust       = 1.2, 
      #                   hjust       = -0.2)                  +
      theme_bw()                                     +
      ggtitle("PCA of Normalised Methylation Data")  +
      theme(axis.title.x    = element_text(size=15),
            axis.title.y    = element_text(size=15)) +
      xlab(label=paste0("PC", i, "(", pcv[i], "%)"))      +
      ylab(label=paste0("PC", j, "(", pcv[j], "%)"))
    
    png(paste0("Norm_PCA_all_", k, "_", i, j, ".png"),
        width  = 12.53, 
        height = 6.98, 
        units  = "in", 
        res    = 600)
    print(gg)
    dev.off()
  }
}
##'-----------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Run Copy Number Variation (CNV) Analysis                                   |
#-------------------------------------------------------------------------------------------#


##'Data In (Requires all probes)
##'-----------------------------------------------------------------------------------------#
fitA              <- lmFit(minfi::getM(lumi.norm.base), design)
mAdj.in           <- minfi::getM(lumi.norm.base)
mAdj.fit          <- fitA$coefficients[,-c(1:5)]
mAdj              <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:5)])
betaAdj           <- ilogit2(mAdj)
dat0              <- as.data.frame(betaAdj)
##'-----------------------------------------------------------------------------------------#


##'Set Up CNV Reference Annotation
##'-----------------------------------------------------------------------------------------#
cnv.anno   <- CNV.create_anno()
##'-----------------------------------------------------------------------------------------#


##'Set Up CNV Reference Annotation
##'-----------------------------------------------------------------------------------------#
cnv.load   <- CNV.load(as.data.frame(mAdj.in[na.omit(match(rownames(mAdj.in),
                                                           cnv.anno@probes@ranges@NAMES)),]))
cnv.load   <- CNV.load(as.data.frame(mAdj.in))
##'-----------------------------------------------------------------------------------------#



##'Set Up CNV fit - 1 Sample Vs Control
##'-----------------------------------------------------------------------------------------#
sample_in  <- pData(lumi.norm.filtered)[pData(lumi.norm.filtered)$SampleType == "Preserved_OA_Hip",]
sample_in  <- rownames(sample_in)[4]
disease_in <- pData(lumi.norm.filtered)$SampleType == "Preserved_Control_Hip"

cnv.x      <- CNV.fit(query = cnv.load[sample_in], 
                      ref   = cnv.load[disease_in], 
                      anno  = cnv.anno)
cnv.x      <- CNV.bin(cnv.x)
cnv.x      <- CNV.segment(cnv.x)
##'-----------------------------------------------------------------------------------------#


##'Plots
##'-----------------------------------------------------------------------------------------#
CNV.genomeplot(cnv.x)
##'-----------------------------------------------------------------------------------------#













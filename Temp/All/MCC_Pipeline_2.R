#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Basic QC - Raw Data                                                        |
#-------------------------------------------------------------------------------------------#



##'Run Minfi QC Report Tools
##'-----------------------------------------------------------------------------------------#
qcReport(raw.idat, 
         sampNames  = pData(raw.idat)$Sample.ID,
         sampGroups = pData(raw.idat)$SampleType, 
         pdf        = "qcReport.pdf")


png("Raw_Beta_Dist_source.png",
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600)
densityPlot(betaAdj, 
            # sampGroups = pData(raw.idat)$SampleType,
            sampGroups = pData(raw.idat)$Source,
            main       = "Beta", 
            xlab       = "Beta Distribution")
dev.off()
##'-----------------------------------------------------------------------------------------#
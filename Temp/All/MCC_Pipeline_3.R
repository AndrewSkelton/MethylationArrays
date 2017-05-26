  #-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Normalise and Filter                                                       |
#-------------------------------------------------------------------------------------------#

setwd("12Oct_Query/")

##'Normalise Data - Functional Normalisation (Minfi)
##'-----------------------------------------------------------------------------------------#
lumi.norm.base      <- preprocessFunnorm(raw.idat,
                                         nPCs    = 2,
                                         sex     = pData(raw.idat)$Gender,
                                         bgCorr  = T,
                                         dyeCorr = T,
                                         verbose = T)
# foo                 <- factor(pData(raw.idat)$Gender)
# levels(foo)         <- c("F","M")
# lumi.norm.base      <- preprocessQuantile(raw.idat,
#                                           sex     = as.vector(foo))
##'-----------------------------------------------------------------------------------------#



##'Filters
##' -Inf Values
##' SNPs - MAF of 5%
##' ChrX/Y Removal
##' detection PVal
##'-----------------------------------------------------------------------------------------#
Mvals               <- minfi::getM(lumi.norm.base)
remove              <- which(Mvals==-Inf, arr.ind=TRUE)[,1]
if(length(remove) > 0) {
  lumi.norm         <- lumi.norm.base[-remove,]
}else{
  lumi.norm         <- lumi.norm.base
}
lumi.norm           <- lumi.norm[-match(lumi.dpval.remove,
                                        rownames(lumi.norm)),]

manifest            <- getManifest(raw.idat)
lumi.norm.snp       <- addSnpInfo(lumi.norm)
lumi.norm.snp.drop  <- dropLociWithSnps(lumi.norm.snp, 
                                        snps = c("SBE","CpG","Probe"), 
                                        maf  = 0.05)
annotation          <- as.data.frame(getAnnotation(lumi.norm.snp.drop))
lumi.norm.filtered  <- lumi.norm.snp.drop[-grep("chrX|chrY", annotation$chr),]
# lumi.norm.filtered  <- lumi.norm.snp.drop[grep("chrX|chrY", annotation$chr),]
annotation          <- as.data.frame(getAnnotation(lumi.norm.filtered))

foo <- rownames(annotation)
save(foo, file = "SexProbes.RData")
##'-----------------------------------------------------------------------------------------#

manifest    <- getManifest(raw.idat)
snpProbesI  <- getProbeInfo(manifest, type = "SnpI")$Name
snpProbesII <- getProbeInfo(manifest, type = "SnpI")$Name
##'-----------------------------------------------------------------------------------------#

anno_in     <- as.data.frame(getAnnotation(lumi.norm.filtered))
features_in <- c("Body", "1stExon", "TSS1500", "TSS200", "5'UTR", "3'UTR")

foo <- c()
for(i in 1:length(features_in)) {
  anno_tmp <- anno_in[grep(features_in[i], 
                           anno_in$UCSC_RefGene_Group),]
  print(paste0("[Group] ",
               features_in[i],
               " - ",
               nrow(anno_tmp)))
  tmp <- betaAdj[match(rownames(anno_tmp),
                       rownames(betaAdj)),]
  foo <- rbind(foo,
               colMeans(tmp))
  rownames(foo)[i] <- paste0(features_in[i],
                             " - ",
                             nrow(anno_tmp))
}
colnames(foo)    <- paste0(pData(lumi.norm.filtered)$Sample.ID,
                           "_",
                           pData(lumi.norm.filtered)$SampleType,
                           "_",
                           pData(lumi.norm.filtered)$Source) 
write.csv(foo,
          file = "Beta_Means_RefGroup.csv")


foo <- c()
for(i in 1:length(unique(anno_in$Relation_to_Island))) {
  anno_tmp <- anno_in[grep(unique(anno_in$Relation_to_Island)[i], 
                           anno_in$Relation_to_Island),]
  print(paste0("[Group] ",
               unique(anno_in$Relation_to_Island)[i],
               " - ",
               nrow(anno_tmp)))
  tmp <- betaAdj[match(rownames(anno_tmp),
                       rownames(betaAdj)),]
  foo <- rbind(foo,
               colMeans(tmp))
  rownames(foo)[i] <- paste0(unique(anno_in$Relation_to_Island)[i],
                                 " - ",
                                 nrow(anno_tmp))
}
colnames(foo)    <- paste0(pData(lumi.norm.filtered)$Sample.ID,
                           "_",
                           pData(lumi.norm.filtered)$SampleType,
                           "_",
                           pData(lumi.norm.filtered)$Source) 
write.csv(foo,
          file = "Beta_Means_IslandRelation.csv")

##'-----------------------------------------------------------------------------------------#

foo <- unique(annotation$Relation_to_Island)

load("Probes_In.RData")
annotation <- annotation[match(array_in_probes,
                               rownames(annotation)),]

Feature.Island  <- rownames(annotation[annotation$Relation_to_Island == "Island",])
Feature.OpenSea <- rownames(annotation[annotation$Relation_to_Island == "OpenSea",])
Feature.N_Shelf <- rownames(annotation[annotation$Relation_to_Island == "N_Shelf",])
Feature.N_Shore <- rownames(annotation[annotation$Relation_to_Island == "N_Shore",])
Feature.S_Shore <- rownames(annotation[annotation$Relation_to_Island == "S_Shore",])
Feature.S_Shelf <- rownames(annotation[annotation$Relation_to_Island == "S_Shelf",])

save(Feature.Island, Feature.OpenSea, Feature.N_Shelf, 
     Feature.N_Shore, Feature.S_Shelf, Feature.S_Shore, file="FeatureIsland.Rdata")

##
Feature.Body    <- rownames(annotation[annotation$UCSC_RefGene_Group == "Body",])
Feature.1Exon   <- rownames(annotation[annotation$UCSC_RefGene_Group == "1stExon",])
Feature.3UTR    <- rownames(annotation[annotation$UCSC_RefGene_Group == "3'UTR",])
Feature.5UTR    <- rownames(annotation[annotation$UCSC_RefGene_Group == "5'UTR",])
Feature.TSS200  <- rownames(annotation[annotation$UCSC_RefGene_Group == "TSS200",])
Feature.TSS1500 <- rownames(annotation[annotation$UCSC_RefGene_Group == "TSS1500",])

save(Feature.Body, Feature.1Exon,  Feature.3UTR, 
     Feature.5UTR, Feature.TSS200, Feature.TSS1500, file="FeatureGene.Rdata")










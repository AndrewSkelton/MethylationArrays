#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : Hovarth Adaptation                                                         |
#-------------------------------------------------------------------------------------------#



##'Set the Working Directory, load essential packages, and import files
##'-----------------------------------------------------------------------------------------#
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("RPMM"))

# library(WGCNA)
# library(sqldf)
# library(impute)
# library(BMIQ)
library(dynamicTreeCut)
##'-----------------------------------------------------------------------------------------#

# current_dir <- getwd()
# dir.create(paste0(current_dir,
#                   "/Hovarth"))
# setwd(paste0(current_dir,
#              "/Hovarth"))


##'Horvath Functions and Prep
##'-----------------------------------------------------------------------------------------#
# source("~/Louise/Hovarth/NORMALIZATION.R")
source("../Scripts/Hovarth/NORMALIZATION.R")
source("../Scripts/Hovarth/Horvath_MethyAge_Normalisation_Functions.R")
probeAnnotation21kdatMethUsed <- read.csv("../Scripts/Hovarth/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k            <- read.csv("../Scripts/Hovarth/datMiniAnnotation.csv")
datClock                      <- read.csv("../Scripts/Hovarth/AdditionalFile3.csv")

trafo      <- function(x, adult.age=20) { 
  x <- (x+1)/(1+adult.age)
  y <- ifelse(x<=1, log(x),x-1)
  y 
}

anti.trafo <- function(x, adult.age=20) { 
  ifelse(x<0, 
         (1+adult.age)*exp(x)-1, 
         (1+adult.age)*x+adult.age) 
}

asnumeric1 <- function(x) {
  as.numeric(as.character(x))
}

# beta_h_in <- minfi::getM(lumi.norm.base)
fitA              <- lmFit(minfi::getM(lumi.norm.base), design)
mAdj.in           <- minfi::getM(lumi.norm.base)
mAdj.fit          <- fitA$coefficients[,-c(1:length(treatments))]
mAdj              <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments))])
betaAdj           <- ilogit2(mAdj)
dat0              <- as.data.frame(betaAdj)
dat0$CpG          <- as.character(rownames(dat0))
dat0              <- dat0[,c(ncol(dat0), 1:(ncol(dat0)-1))]
##'-----------------------------------------------------------------------------------------#


##'Horvath - Stage 1
##'-----------------------------------------------------------------------------------------#
nSamples     <- dim(dat0)[[2]]-1
nProbes      <- dim(dat0)[[1]]
dat0[,1]     <- gsub(x=dat0 [,1],pattern="\"",replacement="")

# file.remove("LogFile.txt")
# file.create("LogFile.txt")
# DoNotProceed <- F

match1       <- match(probeAnnotation21kdatMethUsed$Name, dat0[,1])
dat1         <- dat0[match1,]

set.seed(1)
##'-----------------------------------------------------------------------------------------#


##'Horvath - Stage 2
##'-----------------------------------------------------------------------------------------#
meanMethBySample      <- as.numeric(apply(as.matrix(dat1[,-1]), 2, mean, na.rm=T))
minMethBySample       <- as.numeric(apply(as.matrix(dat1[,-1]), 2, min,  na.rm=T))
maxMethBySample       <- as.numeric(apply(as.matrix(dat1[,-1]), 2, max,  na.rm=T))

datMethUsed           <- t(dat1[,-1])
noMissingPerSample    <- apply(as.matrix(is.na(datMethUsed)), 1, sum)
gs_in                 <- probeAnnotation21kdatMethUsed$goldstandard2

datMethUsedNormalized <- BMIQcalibration(datM              = datMethUsed, 
                                         goldstandard.beta = gs_in, 
                                         plots             = F)
##'-----------------------------------------------------------------------------------------#



##'Horvath - Stage 3
##'-----------------------------------------------------------------------------------------#
selectCpGsClock <- is.element(dimnames(datMethUsedNormalized)[[2]], 
                              as.character(datClock$CpGmarker[-1]))
datMethClock0   <- data.frame(datMethUsedNormalized[,selectCpGsClock])
datMethClock    <- data.frame(datMethClock0[as.character(datClock$CpGmarker[-1])])
predictedAge    <- as.numeric(anti.trafo(datClock$CoefficientTraining[1] + 
                                           as.matrix(datMethClock) %*% 
                                           as.numeric(datClock$CoefficientTraining[-1])))

Comment                      <- ifelse(predictedAge <0, 
                                       "Negative DNAm age.", 
                                       ifelse(predictedAge > 100, 
                                              "Old DNAm age.", 
                                              rep("",length(predictedAge))))
Comment[is.na(predictedAge)] <- "Age prediction was not possible. "
##'-----------------------------------------------------------------------------------------#



##'Horvath - Stage 4 - Write Output
##'-----------------------------------------------------------------------------------------#
datout <- data.frame(SampleID=colnames(dat1)[-1], 
                     DNAmAge=predictedAge, 
                     Comment, 
                     noMissingPerSample, 
                     meanMethBySample, 
                     minMethBySample,
                     maxMethBySample)

datout <- cbind(datout, as.numeric(pData(lumi.norm.base)$Age))
datout <- cbind(datout, pData(lumi.norm.base)$SampleType)
colnames(datout)[8:9] <- c("Chron_Age", "SampleType")
write.csv(datout, file="Methylation_Age.csv")

pData(lumi.norm.filtered)$MethAge   <- as.numeric(datout$DNAmAge)
pData(lumi.norm.filtered)$mcAgeDiff <- pData(lumi.norm.filtered)$MethAge - as.numeric(pData(lumi.norm.filtered)$Age)
##'-----------------------------------------------------------------------------------------#


gg <- ggplot(datout, aes(x=Chron_Age, y=DNAmAge, colour=SampleType)) +
            geom_point() + 
            theme_bw() + 
            scale_y_continuous(limits=c(0,100)) +
            scale_x_continuous(limits=c(0,100)) + 
            stat_smooth(method    = lm, 
                        fullrange = T,
                        se        = F) +
            geom_abline(intercept = 0) + 
            facet_grid(SampleType ~ .)
png("MethAge_Vs_ChronAge.png",
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600) 
print(gg)
dev.off()


datout2 <- as.data.frame(cbind(datout, pData(lumi.norm.base)))
write.csv(datout2, file="Methylation_Age_plus_pheno_data.csv")

setwd(current_dir)


clock_out           <- betaAdj[match(as.vector(datClock$CpGmarker[-1]),
                                     rownames(betaAdj)),]
colnames(clock_out) <- paste0(pData(lumi.norm.base)$Sample.ID,
                             "_",
                             pData(lumi.norm.base)$SampleType,
                             "_",
                             pData(lumi.norm.base)$Source) 
write.csv(clock_out,
          file="Hovarth_353_Adjusted_Betas.csv")




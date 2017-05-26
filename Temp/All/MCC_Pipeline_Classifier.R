#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : My First Classifier                                                        |
#-------------------------------------------------------------------------------------------#



##'Pull in Libraries
##'-----------------------------------------------------------------------------------------#
library(CORElearn)
library(ExplainPrediction)
library(FSelector)
library(evtree)
library(varSelRF)
library(party)
library(rpart)
##'-----------------------------------------------------------------------------------------#


##'
##'-----------------------------------------------------------------------------------------#
# classifier.features <- rownames(`Preserved_OA_Hip-Preserved_OA_Knee_0.01`)
# classifier.typeA    <- 'Preserved_OA_Hip'
# classifier.typeB    <- 'Preserved_OA_Knee'
classifier.features <- rownames(`Preserved_Control_Hip-Preserved_OA_Hip_0.01`)
classifier.typeA    <- 'Preserved_Control_Hip'
classifier.typeB    <- 'Preserved_OA_Hip'
classifier.outcome  <- as.vector(pData(lumi.norm.filtered)$SampleType)
classifier.samples  <- grep(paste(classifier.typeA,
                                  classifier.typeB,
                                  sep = "|"),
                            classifier.outcome)
classifier.pheno    <- pData(lumi.norm.filtered)[classifier.samples,]
classifier.df       <- betaAdj[,classifier.samples]
classifier.df       <- classifier.df[match(classifier.features,
                                           rownames(classifier.df)),]
classifier.df       <- as.data.frame(t(classifier.df))
classifier.df       <- cbind(classifier.df,
                             as.data.frame(classifier.pheno[,c(2,3,6)]))
##'-----------------------------------------------------------------------------------------#


##'First things first, what to classify....? - OA Vs Non-OA
##'-----------------------------------------------------------------------------------------#
# classifier.df       <- as.data.frame(cbind(t(betaAdj),
#                                      as.data.frame(pData(lumi.norm.filtered)[,c(2,3,6)])))
# classifier.outcome  <- as.vector(pData(lumi.norm.filtered)$SampleType)
# classifier.outcome[grep("OA",      classifier.outcome)] <- "OA"
# classifier.outcome[grep("Control", classifier.outcome)] <- "Control"
# classifier.df$SampleType <- classifier.outcome
##'-----------------------------------------------------------------------------------------#


##'Designate Test and Training Samples
##'-----------------------------------------------------------------------------------------#
set.seed(123)
trainIdxs          <- sample(x       = nrow(classifier.df), 
                             size    = 0.7 * nrow(classifier.df), 
                             replace = F)
testIdxs           <- c(1:nrow(classifier.df))[-trainIdxs]
##'-----------------------------------------------------------------------------------------#


##'Feature Selection - Optimised for large datasets
##'-----------------------------------------------------------------------------------------#
estReliefF         <- attrEval(formula          = "SampleType", 
                               data             = classifier.df[trainIdxs,], 
                               estimator        = "ReliefFexpRank", 
                               ReliefIterations = 10)
classifier.df      <- classifier.df[,estReliefF > (max(estReliefF)-0.2)]
classifier.df$SampleType <- factor(classifier.pheno$SampleType)
##'-----------------------------------------------------------------------------------------#


##'Build the Model - Random Forest
##'-----------------------------------------------------------------------------------------#
modelRF            <- CoreModel(formula            = "SampleType", 
                                data               = classifier.df[trainIdxs,], 
                                model              = "rf", 
                                selectionEstimator = "MDL",
                                minNodeWeightRF    = 5,
                                rfNoTrees          = 100, 
                                maxThreads         = 20)
# destroyModels(modelRF)

pred <- predict(object  = modelRF, 
                newdata = classifier.df[testIdxs,],
                type    = "both")

mEval <- modelEval(model          = modelRF, 
                   correctClass   = classifier.df[testIdxs,]$SampleType, 
                   predictedClass = pred$class, 
                   predictedProb  = pred$prob)

plot(modelRF)

print(mEval)

plot(pred$probabilities[,1]*100)

explainVis(model        = modelRF, 
           trainData    = classifier.df[trainIdxs,], 
           testData     = classifier.df[testIdxs,], 
           method       = "EXPLAIN",
           visLevel     = "model",
           problemName  = "classifier.df", 
           fileType     = "none", 
           classValue   = 1, 
           displayColor = "color")
##'-----------------------------------------------------------------------------------------#


##'Build the Model - Tree
##'-----------------------------------------------------------------------------------------#
modelT             <- CoreModel(formula            = "SampleType", 
                                data               = classifier.df[trainIdxs,], 
                                model              = "tree", 
                                selectionEstimator = "MDL", 
                                maxThreads         = 20)
# destroyModels(modelT)

pred               <- predict(object  = modelT, 
                              newdata = classifier.df[testIdxs,],
                              type    = "both")

mEval              <- modelEval(model          = modelT, 
                                correctClass   = classifier.df[testIdxs,]$SampleType, 
                                predictedClass = pred$class, 
                                predictedProb  = pred$prob)


print(mEval)
plot(pred$probabilities[,1]*100)

explainVis(model        = modelT, 
           trainData    = classifier.df[trainIdxs,], 
           testData     = classifier.df[testIdxs,], 
           method       = "EXPLAIN",
           visLevel     = "model",
           problemName  = "classifier.df", 
           fileType     = "none", 
           classValue   = 1, 
           displayColor = "color")

plot(modelT, classifier.df[trainIdxs,])

rpm <- getRpartModel(modelT,classifier.df[trainIdxs,])
plot(rpm, branch=0.5, minbranch=5, compress=TRUE)
text(rpm, pretty=0, digits=3)
##'-----------------------------------------------------------------------------------------#


##'Build the Model - Conditional Inference Tree
##'-----------------------------------------------------------------------------------------#
ct               <- ctree(formula = SampleType~., 
                          data    = classifier.df[trainIdxs,])
plot(ct, 
     main="Conditional Inference Tree")
table(predict(ct), 
      classifier.df[trainIdxs,]$SampleType)

tr.pred          <- predict(ct, 
                            newdata = classifier.df[testIdxs,], 
                            type    = "prob")
##'-----------------------------------------------------------------------------------------#


##'Build the Model - RPart
##'-----------------------------------------------------------------------------------------#
rpart.fit <- rpart(formula = SampleType ~ ., 
                   method  = "class", 
                   data    = classifier.df[trainIdxs,])

printcp(rpart.fit) # display the results
plotcp(rpart.fit) # visualize cross-validation results
summary(rpart.fit) # detailed summary of splits

plot(rpart.fit, 
     uniform = T, 
     main    = "Classification Tree for Sample Type")
text(rpart.fit, 
     use.n   = T,
     all     = T, 
     cex     = .8)

table(subset(classifier.df[trainIdxs,], 
             cg05877497 > 0.3748248)$SampleType)
##'-----------------------------------------------------------------------------------------#





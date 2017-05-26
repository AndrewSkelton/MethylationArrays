probes_of_interest <- c("cg23606718",
                        "cg00059225",
                        "cg23606718",
                        "cg23500537",
                        "cg11071401",
                        "cg12757011",
                        "cg12757011",
                        "cg24150153")

model_all <- betaAdj[match(probes_of_interest,
                           rownames(betaAdj)),]

pheno_all <- as.data.frame(pData(lumi.norm.filtered))

for(i in 1:length(unique(pheno_all$SampleType))) {
  tmp <- model_all[,grep(unique(pheno_all$SampleType)[i],
                         pheno_all$SampleType)]
  print(unique(pheno_all$SampleType)[i])
  print(rowMeans(tmp))
}


model_interaction_model_cAge <- betaAdj.age[match(probes_of_interest,
                                                  rownames(betaAdj.age)),]
for(i in 1:length(unique(pheno_all$SampleType))) {
  tmp <- model_interaction_model_cAge[,grep(unique(pheno_all$SampleType)[i],
                                            pheno_all$SampleType)]
  print(unique(pheno_all$SampleType)[i])
  print(rowMeans(tmp))
}











probe_in         <- "cg23606718"
raw_data_probe_m <- getM(lumi.norm.filtered[probe_in,])
raw_data_probe_b <- getBeta(lumi.norm.filtered[probe_in,])
sampleType_in    <- pData(lumi.norm.filtered)$SampleType
sampleID_in      <- rownames(pData(lumi.norm.filtered))

##'Model 1 - A vs B
coef_m1          <- mAdj.fit[probe_in,]
adjusted_m_m1    <- mAdj[probe_in,]
adjusted_b_m1    <- betaAdj[probe_in,]
  
##'Model 2 - ~SampleType:Age + Data_Source
coef_m2          <- mAdj.fit.age[probe_in,]
adjusted_m_m2    <- mAdj.age[probe_in,]
adjusted_b_m2    <- betaAdj.age[probe_in,]

##'Model 3 - ~Age + Data_Source
# coef_m2          <- mAdj.fit.age[probe_in,]
# adjusted_m_m2    <- mAdj.age[probe_in,]
# adjusted_b_m2    <- betaAdj.age[probe_in,]
  
  
df <- data.frame(Value = c(raw_data_probe_m,
                          raw_data_probe_b,
                          adjusted_m_m1,
                          adjusted_b_m1,
                          adjusted_m_m2,
                          adjusted_b_m2),
                model = c(rep("Raw",       length(sampleType_in)*2),
                          rep("Model_1",   length(sampleType_in)*2),
                          rep("Model_2",   length(sampleType_in)*2)),
                Type  = rep(sampleType_in, 6),
                data  = c(rep("M",         length(sampleType_in)),
                          rep("Beta",      length(sampleType_in)),
                          rep("M",         length(sampleType_in)),
                          rep("Beta",      length(sampleType_in)),
                          rep("M",         length(sampleType_in)),
                          rep("Beta",      length(sampleType_in))),
                ID    = rep(sampleID_in, 6))


df_in <- df[df$Type == "Preserved_Control_Hip",]

ggplot(df_in, aes(x = ID,
               y = Value,
               group = Type,
               colour = Type)) +
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  facet_grid(data ~ model, scales = "free")



SampleType <- rep(c("A", "B", "C"), 6)
Age        <- rep(45:47, each=6) # Note use of "each"
Age        <- rep(45:47, 6)
Batch      <- c(rep("X", 6),
                rep("Y", 12))

colnames(model.matrix(~0 + SampleType + Batch))
colnames(model.matrix(~0 + SampleType*Age + Batch))
colnames(model.matrix(~0 + SampleType:Age + Batch))
colnames(model.matrix(~Age + Batch))


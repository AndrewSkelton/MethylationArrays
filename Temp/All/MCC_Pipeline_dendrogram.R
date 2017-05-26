
##'Dendrogram Utilising the ggplot2 library
##'-----------------------------------------------------------------------------------------#
for(i in 1:length(levels(factor(pData(lumi.norm.filtered)$SampleType)))) {
  betaAdj_in <- betaAdj[,grep(levels(factor(pData(lumi.norm.filtered)$SampleType))[i],
                              pData(lumi.norm.filtered)$SampleType)]
  d          <- dist(t(betaAdj_in), 
                     method   = 'euclidian')
  h          <- hclust(d, 
                       method = 'complete')
  dd         <- dendro_data(h)
  labs       <- label(dd)
  labs$group <- pData(lumi.norm.filtered)[match(colnames(betaAdj_in),
                                                rownames(pData(lumi.norm.filtered))),]$SampleType
  
  g          <- ggdendrogram(h, 
                             leaf_labels = F, 
                             labels      = T, 
                             rotate      = T) +
                    scale_x_discrete(labels   = c()) +
                    scale_y_continuous(expand = c(0.1, 0)) +
                    geom_text(data    = labs, 
                              aes(label  = label, 
                                  x      = x, 
                                  y      = y,
                                  colour = as.factor(group)), 
                              hjust   = 1.0, 
                              size    = 4.0) +
                    theme(legend.text  = element_text(size=14),
                          legend.title = element_text(size=14))
  print(g)
  
  png(paste0("Dendrogram/Dendrogram_",
             levels(factor(pData(lumi.norm.filtered)$SampleType))[i],
             "_Adjusted.png"),
      width  = 12.53, 
      height = 6.98, 
      units  = "in", 
      res    = 600)
  print(g)
  dev.off()
}
##'-----------------------------------------------------------------------------------------#



##'Differential Methylation from Hierarchical Clustering
##'-----------------------------------------------------------------------------------------#
i          <- 4
betaAdj_in <- betaAdj[,grep(levels(factor(pData(lumi.norm.filtered)$SampleType))[i],
                            pData(lumi.norm.filtered)$SampleType)]
d          <- dist(t(betaAdj_in), 
                   method   = 'euclidian')
h          <- hclust(d, 
                     method = 'complete')
dd         <- dendro_data(h)

foo <- cutree(rev(h), k=4)
tmp <- data.frame(samples=names(foo), groups=foo)
tmp <- tmp[with(tmp, order(-groups)), ]
tmp[tmp$groups == 4,]$groups <- 3

tmp$groups           <- factor(tmp$groups)
levels(tmp$groups)   <- c("B","C","A")
# lumi.norm.in         <- lumi.norm.filtered[,pData(lumi.norm.filtered)$SampleType == "Preserved_OA_Hip"]
lumi.norm.in         <- lumi.norm.filtered[,grep("Hip", pData(lumi.norm.filtered)$SampleType)]

pData(lumi.norm.in)$Dendro_Group  <- "NB"
pData(lumi.norm.in)[match(as.vector(tmp$samples),
                          rownames(pData(lumi.norm.in))),]$Dendro_Group <- as.vector(tmp$groups)
# pData(lumi.norm.in)$Dendro_Group  <- factor(pData(lumi.norm.in)$Dendro_Group)

# levels(pData(lumi.norm.in)$Dendro_Group) <- c("AB","AB","C")
treatments           <- levels(factor(pData(lumi.norm.in)$Dendro_Group))
treatment_arrays     <- factor(pData(lumi.norm.in)$Dendro_Group)
batch_arrays         <- factor(pData(lumi.norm.in)$Source)

design               <- model.matrix(~0 + treatment_arrays + batch_arrays)
colnames(design)     <- c(treatments, paste0("be",
                                             levels(batch_arrays)[-1]))

#Beta or M Value Selection
fitB                  <- lmFit(minfi::getBeta(lumi.norm.in), design)
fitM                  <- lmFit(minfi::getM(lumi.norm.in),    design)

cont_mat              <- makeContrasts(#AB - C,
                                       A-B,
                                       A-C, 
                                       B-C,
                                       # Query_A.B_C = ((A + B)/2) - C,
                                       Query_A.C_B = ((A + C)/2) - B,
                                       levels = colnames(design))

fit2B                 <- eBayes(contrasts.fit(fitB,
                                              contrasts=cont_mat))
fit2M                 <- eBayes(contrasts.fit(fitM,
                                              contrasts=cont_mat))

mAdj.in               <- minfi::getM(lumi.norm.in)
mAdj.fit              <- fitM$coefficients[,-c(1:length(treatments))]
mAdj                  <- as.matrix(mAdj.in) - mAdj.fit %*% t(design[,-c(1:length(treatments))])
betaAdj.dendro        <- ilogit2(mAdj)

contrasts             <- colnames(cont_mat)
pVal                  <- 0.05
i                     <- 5


for(i in 1:length(contrasts)) {
gene_list_B             <- topTable(fit2B,
                                    coef          = contrasts[i],
                                    number        = Inf,
                                    adjust.method = "BH",
                                    sort.by       = "P")

gene_list_M             <- topTable(fit2M,
                                    coef          = contrasts[i],
                                    p.value       = pVal,
                                    number        = Inf,
                                    adjust.method = "BH",
                                    sort.by       = "P")

gene_list_B             <- gene_list_B[match(rownames(gene_list_M), 
                                             rownames(gene_list_B)),]

gene_list_out           <- cbind(gene_list_M[,c(1,4,5)], 
                                 gene_list_B[,c(1,4,5)])
colnames(gene_list_out) <- c("M-FC", "M-Pval", "M-Adj.Pval",
                             "B-FC", "B-Pval", "B-Adj.Pval")
gene_list_out           <- cbind(gene_list_out,
                                 annotation[match(rownames(gene_list_out),
                                                  rownames(annotation)),])
gene_list_out           <- gene_list_out[abs(gene_list_out$`B-FC`) > 0.1,]

write.csv(gene_list_out, file=paste0(contrasts[i],
                                     "_gene_list.csv"))


beta_out           <- betaAdj.dendro[rownames(gene_list_out),]
colnames(beta_out) <- paste0(pData(lumi.norm.in)$Sample.ID,
                             "_",
                             pData(lumi.norm.in)$SampleType,
                             "_",
                             pData(lumi.norm.in)$Source,
                             "_",
                             pData(lumi.norm.in)$Dendro_Group) 
write.csv(beta_out, file=paste0(contrasts[i], "_Adjusted_Betas.csv"))

message(paste0(contrasts[i],
               ": ",
               nrow(gene_list_out),
               " Significant Probes"))
}
##'-----------------------------------------------------------------------------------------#

colnames(beta_out) <- paste0(pData(lumi.norm.in)$Dendro_Group, "_", colnames(beta_out))
beta_out           <- beta_out[,grep("NB", colnames(beta_out), invert = T)]

d          <- dist(t(beta_out), 
                   method   = 'euclidian')
h          <- hclust(d, 
                     method = 'complete')
png(paste0("Dendrogram_A.C_B_Probes.png"),
    width  = 12.53, 
    height = 6.98, 
    units  = "in", 
    res    = 600)
ggdendrogram(h,
             leaf_labels = F,
             labels      = T,
             rotate      = T)
dev.off()

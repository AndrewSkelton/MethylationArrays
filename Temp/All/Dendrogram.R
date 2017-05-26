install.packages("ggdendro")
library(ggdendro)

pheno_in <- lumi.norm.filtered %>% pData %>% as.data.frame %>% add_rownames("label") %>% 
            filter(grepl('Preserved_Control_Hip|Preserved_OA_Hip',SampleType)) 
mAdj_in  <- mAdj[match(gene_list_out$Name,rownames(betaAdj)),match(pheno_in$label,colnames(betaAdj))]
# mAdj_in  <- mAdj[,match(pheno_in$label,colnames(betaAdj))]
hc       <- mAdj_in %>% t %>% dist %>% hclust
dd       <- hc %>% dendro_data
labs     <- label(dd) %>% left_join(pheno_in)

gg <- ggdendrogram(hc, rotate = T, size = 2) +
      scale_x_discrete(labels = c()) + 
      scale_y_continuous(expand = c(0.1, 0)) +
      geom_text(data=labs, aes(label=SampleType, x=x, y=y, colour=as.factor(SampleType)), hjust=1.0, size=2.0) +
      theme(legend.text=element_text(size=4), legend.title=element_text(size=6))
png("../26Oct_Query/Hip_samples_DMList.png",width  = 12.53, height = 8, units  = "in", res = 300)
print(gg)
dev.off()

write.csv(dd$labels %>% left_join(pheno_in), file = "../26Oct_Query/Hip_samples_DMList.csv")







library(ComplexHeatmap)
library(circlize)

##'Complex Heatmap - Annotation
##'-----------------------------------------------------------------------------------------#
pheno_in <- pData(lumi.norm.filtered)[match(phenotmp$CpG, rownames(pData(lumi.norm.filtered))),]
pheno_in$Subtype <- phenotmp$SampleType_1
df  <- data.frame(Sample     = pheno_in$Subtype)
ha  <- HeatmapAnnotation(df  = df, 
                         col = list(Sample = c("D"  =  "black", 
                                               "E"       =  "grey")))
bAdj_in  <- betaAdj[match(gene_list_out$Name,rownames(betaAdj)),
                    match(rownames(pheno_in),colnames(betaAdj))] %>% t %>% scale %>% t
##'-----------------------------------------------------------------------------------------#



##'Complex Heatmap - Main Function Call
##'-----------------------------------------------------------------------------------------#
ht2 = Heatmap(bAdj_in, 
              name                     = "Row Z-Score", 
              column_title             = "", 
              top_annotation           = ha,
              col                      = colorRamp2(c(-2, 0, 2), 
                                                    c("blue", "white", "red")),
              show_row_dend            = F,
              cluster_rows             = T,
              show_column_dend         = T,
              show_column_names        = F,
              column_dend_height       = unit(1, "cm"),
              show_row_names = F
)

png("../26Oct_Query/D_and_E_Heatmap.png",
    width  = 12.53, height = 6.98, 
    units  = "in", res    = 600)
draw(ht2)
dev.off()
##'-----------------------------------------------------------------------------------------#


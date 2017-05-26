##'Differential Methylation - Model Fit
##'-----------------------------------------------------------------------------------------#
Treat                  <- pData(lumi.norm.filtered)$Age %>% as.numeric
design                 <- model.matrix(~Treat)
fitB                   <- lmFit(minfi::getBeta(lumi.norm.filtered), design) %>% eBayes
fitM                   <- lmFit(minfi::getM(lumi.norm.filtered),    design) %>% eBayes
##'-----------------------------------------------------------------------------------------#


##'Differential Methylation - Extract Results
##'-----------------------------------------------------------------------------------------#
lvcm.filter    <- topTable(fit           = fitM,
                           coef          = 2,
                           p.value       = 0.05,
                           number        = Inf,
                           adjust.method = "BH",
                           sort.by       = "P")

betaAdj_subset <- getBeta(lumi.norm.filtered)

lvcm.cpgIn  <- rownames(lvcm.filter)[1:nrow(lvcm.filter)]
lvcm.cpgIn  <- rownames(lvcm.filter)[1:3]
lvcm.df     <- melt(betaAdj_subset[lvcm.cpgIn,])

if(length(lvcm.cpgIn) > 1) {
  colnames(lvcm.df)         <- c("CpG", "Sample", "Beta")
} else {
  lvcm.df$CpG               <- lvcm.cpgIn
  lvcm.df$Sample            <- rownames(lvcm.df)
  rownames(lvcm.df)         <- c(1:nrow(lvcm.df))
  lvcm.df                   <- lvcm.df[,c(2,3,1)]
  colnames(lvcm.df)[3]      <- "Beta"
}

lvcm.df     <- lvcm.df %>% left_join(pData(lumi.norm.filtered) %>% as.data.frame, by=c("Sample"="Sentrix_ID"))
lvcm.df$Age <- as.numeric(lvcm.df$Age)

lvcm.gg      <- ggplot(lvcm.df, 
                       aes(x                = Age, 
                           y                = Beta, 
                           colour           = SampleType)) + 
  geom_point() +
  geom_smooth(method        = lm, 
              se            = T, 
              aes(group     = 1)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,100)) +
  #                   scale_x_continuous(limits = c(min(lvcm.df$Age),
  #                                                 max(lvcm.df$Age))) +
  theme_bw()

if(length(lvcm.cpgIn) > 1) {
  lvcm.gg <- lvcm.gg + 
    # facet_grid(SampleType ~ CpG)
    facet_grid(. ~ CpG)
}

png("Top3_cAge.png",
    width  = 8.53,
    height = 8.53,
    units  = "in",
    res    = 600)
print(lvcm.gg)
dev.off()
##'-----------------------------------------------------------------------------------------#


##'cAge Regression - Output
##'-----------------------------------------------------------------------------------------#
betaOut <- betaAdj_subset %>% as.data.frame %>% add_rownames("CpG") %>%  dplyr::slice(match(rownames(lvcm.filter)[1:100], CpG))
wb                     <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Pheno_Data")
openxlsx::writeData(wb, "Pheno_Data", pheno_in)
openxlsx::addWorksheet(wb, "cAgeRegression")
openxlsx::writeData(wb, "cAgeRegression", lvcm.filter %>% add_rownames("CpG"))
openxlsx::addWorksheet(wb, "Top100Beta")
openxlsx::writeData(wb, "Top100Beta", betaOut)
openxlsx::saveWorkbook(wb, "Scripps_cAgeRegression.xlsx", overwrite = T)
##'-----------------------------------------------------------------------------------------#



contrasts              <- fitM$coefficients %>% colnames
names(contrasts)       <- fitM$coefficients %>% colnames
pVal                   <- 0.01
deltaB                 <- 0.1
i                      <- "Treat"

wb                     <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Pheno_Data")
openxlsx::writeData(wb, "Pheno_Data", pheno_in)

  gene_list_B            <- topTable(fitB, coef = contrasts[i], number = Inf, sort.by = "P")
  gene_list_M            <- topTable(fitM, coef = contrasts[i], number = Inf, sort.by = "P",
                                     p.value    = pVal)
  if(nrow(gene_list_M) > 0){
    gene_list_out          <- gene_list_M %>% 
      add_rownames("CpG") %>% 
      left_join(gene_list_B %>% add_rownames("CpG"), by = c("CpG" = "CpG"))
    colnames(gene_list_out)<- colnames(gene_list_out) %>% gsub(".x",".M",.) %>% gsub(".y",".B",.)
    gene_list_out          <- gene_list_out %>% 
      filter(abs(logFC.B) > deltaB) %>% 
      dplyr::select(CpG,logFC.B,P.Value.M,adj.P.Val.M) %>% 
      left_join(annotation, by = c("CpG" = "Name")) %>% 
      as.data.frame
    message(paste0("Differential Methylation in ", names(contrasts)[i], ": ", nrow(gene_list_out)))
    
    # Write xlsx sheet
    openxlsx::addWorksheet(wb, names(contrasts)[i])
    openxlsx::writeData(wb, names(contrasts)[i], gene_list_out)
    betas_in  <- lumi.norm.filtered %>%
      minfi::getBeta(.) %>%
      as.data.frame %>%
      add_rownames("CpG") %>%
      dplyr::slice(match(gene_list_out$CpG, CpG))
    openxlsx::addWorksheet(wb, paste0("Betas ",names(contrasts)[i]))
    openxlsx::writeData(wb, paste0("Betas ",names(contrasts)[i]), betas_in)
    
    # Write Bed Track
    bed_nme <- names(contrasts)[i] %>% gsub(" ", "",.)
    bed_out <- data.frame(chr = gene_list_out$chr) %>%
      mutate(start   = gene_list_out$pos - 1,
             end     = gene_list_out$pos %>% as.numeric,
             name    = paste0(bed_nme, "_DeltaBeta_",
                              gene_list_out$logFC.B %>% round(3)),
             score   = gene_list_out$logFC.B %>% round(3),
             strand  = gene_list_out$strand %>% as.vector,
             box_s   = start,
             box_e   = end,
             colour  = ifelse(score < 0, "0,43,255", "255,0,0"))
    write.table(bed_out,
                sep       = "\t",
                row.names = F,
                col.names = F,
                quote     = F,
                file      = paste0(bed_nme,"_Scripps.bed"))
    exe     <- paste0("sed -i '1s/^/track name=\"",
                      bed_nme, "\" description=\"",
                      bed_nme, "\" visibility=3 ",
                      "itemRgb=\"On\"",
                      "\\n/' ",
                      bed_nme, "_Scripps.bed")
    system(exe)
    message(paste0("Differential Methylation in ", names(contrasts)[i], ": 0"))
openxlsx::saveWorkbook(wb, "DM_450K_Scripps.xlsx", overwrite = T)

##'-----------------------------------------------------------------------------------------#




data_in   <- read_tsv("Query/Louise_Query.txt") %>% as.data.frame
data_out  <- betaAdj_subset %>% as.data.frame %>% add_rownames("CpG") %>% dplyr::slice(match(data_in$CpG, CpG))
wb        <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Pheno_Data")
openxlsx::writeData(wb, "Pheno_Data", pheno_in)
openxlsx::addWorksheet(wb, "Query_Beta")
openxlsx::writeData(wb, "Query_Beta", data_out)
openxlsx::saveWorkbook(wb, "Query/Louise_Query_CpG_17Aug.xlsx", overwrite = T)









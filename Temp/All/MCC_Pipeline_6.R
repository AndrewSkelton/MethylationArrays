#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : DMR Detection                                                              |
#-------------------------------------------------------------------------------------------#



##'Data Setup for DMR Detection - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
# bump.samples.1        <- "OA_Hip"
# bump.samples.2        <- "Hip_Control"
##
# bump.samples.1        <- "OA_Knee"
# bump.samples.2        <- "Hip_Control"
##
bump.samples.1        <- "OA_Knee"
bump.samples.2        <- "OA_Hip"
##

bump.in               <- lumi.norm.filtered[,grep(paste(bump.samples.1,
                                                        bump.samples.2,
                                                        sep="|"), 
                                                  pData(lumi.norm.filtered)$SampleType)]
treatments            <- unique(pData(bump.in)$SampleType)
treatment_arrays      <- factor(pData(bump.in)$SampleType)
batch_arrays          <- factor(pData(bump.in)$Batch)
bump.design           <- model.matrix(~treatment_arrays + batch_arrays)

##Run Once
# gtf_path      <- "~/genome.gtf"
# transcriptdb  <- makeTxDbFromGFF(file     = gtf_path, 
#                                 format   = 'gtf', 
#                                 organism = 'Homo sapiens')
# ensembl    <- useMart("ENSEMBL_MART_ENSEMBL",
#                       dataset="hsapiens_gene_ensembl",
#                       host="grch37.ensembl.org",
#                       path="/biomart/martservice")
##
##'-----------------------------------------------------------------------------------------#



##'Run Bumphunter in bootstap mode for extra model covariates 
##'-----------------------------------------------------------------------------------------#
registerDoParallel(cores = 20)
dmrs <- minfi::bumphunter(object         = bump.in, 
                          coef           = 2,
                          design         = bump.design,
                          cutoff         = 0.1, 
                          B              = 50, 
                          type           = "Beta",
                          smooth         = T,
                          smoothFunction = loessByCluster,
                          nullMethod     = "bootstrap",
                          pickCutoffQ    = 0.99)

bump.table            <- dmrs$table
bump.table            <- bump.table[abs(bump.table$start - bump.table$end) != 0,]
bump.table            <- bump.table[with(bump.table, order(p.value)), ]
bump.table            <- bump.table[bump.table$p.value < 0.05,]

# bump.fit <- lmFit(getBeta(bump.in), bump.design)
write.csv(bump.table,
          file=paste0(bump.samples.1,
                      "_",
                      bump.samples.2,
                      "_DMR.csv"))
##'-----------------------------------------------------------------------------------------#



##'DMR Visualisations - Any range can be used, but this takes the results of a bumphunter
##'                     output table: bump.row 
##'-----------------------------------------------------------------------------------------#
if(nrow(bump.table) > 9) { plot_range <- 10 } else { plot_range <- nrow(bump.table) }
for(i in 1:plot_range) {

bump.offset           <- 10000
bump.row              <- i
bump.pos.offset       <- 200
bump.chr              <- bump.table[bump.row,]$chr
bump.start            <- bump.table[bump.row,]$start - bump.offset
bump.end              <- bump.table[bump.row,]$end   + bump.offset
bump.annotation       <- as.data.frame(getAnnotation(bump.in))
bump.annotation.in    <- bump.annotation[bump.annotation$chr == bump.chr   &
                                           bump.annotation$pos >= bump.start &
                                           bump.annotation$pos <= bump.end,]
bump.subset           <- m2beta(mAdj[match(rownames(bump.annotation.in),
                                           rownames(mAdj)),
                                     grep(paste0(bump.samples.1, 
                                                 "|", 
                                                 bump.samples.2), 
                                          pData(lumi.norm.filtered)$SampleType)])
df                    <- melt(bump.subset)
df                    <- cbind(df,
                               bump.annotation.in[match(df$Var1, 
                                                        rownames(bump.annotation.in)),])
df                    <- cbind(df,
                               pData(bump.in)[match(df$Var2, 
                                                    rownames(pData(bump.in))),])
colnames(df)[1:3]     <- c("CpG", "Sample", "Beta")
df                    <- as.data.frame(df)
bump.pos.limit        <- c(min(df$pos  - bump.pos.offset),
                           max(df$pos))


##Get Transcript Annotation
geneID_tmp             <- crunch(transcriptdb, 
                                which=GRanges(seqnames = as.numeric(gsub("chr", 
                                                                         "", 
                                                                         bump.chr)), 
                                              IRanges(bump.start, 
                                                      bump.end)))

geneID_in              <- subsetByOverlaps(geneID_tmp,
                                           GRanges(seqnames = as.numeric(gsub("chr", 
                                                                              "", 
                                                                              bump.chr)), 
                                                   IRanges(bump.start + bump.offset, 
                                                           bump.end   - bump.offset)))


if(length(geneID_in) > 0) {
  annotation.in         <- getBM(attributes = c("ensembl_gene_id",
                                                "external_gene_name"),
                                 filters    = "ensembl_gene_id",
                                 values     = as.vector(unique(geneID_in@elementMetadata@listData$gene_id)),
                                 mart       = ensembl)
  annotation.in         <- annotation.in$external_gene_name[1]
} else {
  annotation.in <- "No Data"
}
##



##Plot3 - Gene Track
trans         <- crunch(transcriptdb, 
                        which=GRanges(seqnames = as.numeric(gsub("chr", 
                                                                 "", 
                                                                 bump.chr)), 
                                      IRanges(min(df$pos), 
                                              max(df$pos))))
gr1           <- split(trans, 
                       trans$tx_name)
p1            <- ggbio::autoplot(gr1)                  + 
                          theme_bw()                                           +
                          geom_vline(xintercept     = bump.start + bump.offset, 
                                     colour         = "red")                   +
                          geom_vline(xintercept     = bump.end   - bump.offset, 
                                     colour         = "red")                   +
                          scale_x_continuous(limits = bump.pos.limit)   
##


##Plot1 - Betas
gg1 <- ggplot(df, 
              aes(x         = pos, 
                  y         = Beta, 
                  colour    = SampleType))                +
  geom_point()                                            + 
  scale_colour_hue(l        = 50)                         +
  geom_vline(xintercept     = c(bump.start + bump.offset,
                                bump.end   - bump.offset),
             colour         ='red')                       + 
  theme_bw()                                              +
  stat_summary(fun.y        = mean, 
               geom         = "line",
               aes(colour   = SampleType, 
                   group    = SampleType))                +
  scale_x_continuous(limits = bump.pos.limit)             +
  scale_y_continuous(limits = c(0, 1))                    +
  theme(legend.title        = element_blank(),
        legend.position = "top",
        legend.background = element_rect(fill = alpha('white', 0)),
        legend.margin   = unit(-0.6,"cm"))                +
  xlab("")                                                +
  ggtitle(paste0("Detected Bump ",
                 bump.chr,
                 ":",
                 bump.start + bump.offset,
                 "-",
                 bump.end   - bump.offset,
                 " in ",
                 bump.samples.1,
                 " Vs ",
                 bump.samples.2,
                 "\nGene: ",
                 annotation.in))
##


##Plot2 - Strip Plot of Site Type
gg2 <- ggplot(df, aes(x      = pos, 
                      y      = Relation_to_Island))        +
  theme_bw()                                               +
  geom_point(size            = 3)                          +
  ylab("")                                                 +
  theme(legend.title         = element_blank(),
        legend.position      = "bottom",
        legend.margin        = unit(-0.6,"cm"))            +                
  xlab("")                                                 +
  scale_x_continuous(limits = bump.pos.limit)              +
  geom_vline(xintercept      = c(bump.start + bump.offset,
                                 bump.end   - bump.offset), 
             colour          = "red")
##


ggOut <- plot_grid(gg1, 
                   gg2,
                   p1@ggplot,
                   align       = 'v', 
                   labels      = c('A', 'B', 'C'),
                   nrow        = 3,
                   rel_heights = c(4, 1.5, 3))

save_plot(paste0(bump.samples.1,
                 "_",
                 bump.samples.2,
                 "_Idx_",
                 i,
                 "_DMR_.png"),
          ggOut,
          base_height = 10,
          base_width  = 12.53,
          dpi         = 300)
}
##'-----------------------------------------------------------------------------------------#

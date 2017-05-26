#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : R Statistical Programming Language                                         |
#  Study       : Methylation Cartilage Consortium                                           |
#  Data Owner  : Newcastle University - Prof. John Loughlin, Dr. Louise Reynard             |
#  Description : DMR Detection - DMRcate                                                    |
#-------------------------------------------------------------------------------------------#



##'Data Setup for DMR Detection - Subsetting and Model Design
##'-----------------------------------------------------------------------------------------#
# bump.samples.1        <- "Hip_Control"
# bump.samples.2        <- "OA_Hip"
# bump.in               <- 
# treatments            <- unique(pData(bump.in)$SampleType)
# treatment_arrays      <- factor(pData(bump.in)$SampleType)


##Run Once
gtf_path      <- "~/genome.gtf"
transcriptdb  <- makeTxDbFromGFF(file     = gtf_path,
                                 format   = 'gtf',
                                 organism = 'Homo sapiens')
ensembl    <- useMart("ENSEMBL_MART_ENSEMBL",
                      dataset="hsapiens_gene_ensembl",
                      host="grch37.ensembl.org",
                      path="/biomart/martservice")

##'-----------------------------------------------------------------------------------------#

current_dir <- getwd()
dir.create(paste0(current_dir,
                  "/DMR"))
setwd(paste0(current_dir,
             "/DMR"))

##'Run DMRcate
##'-----------------------------------------------------------------------------------------#
for(i in 1:length(contrasts)) {
  dir.create(paste0(current_dir,
                    "/DMR/",
                    gsub(" ", 
                         "",
                         contrasts[i])))
  setwd(paste0(current_dir,
               "/DMR/",
               gsub(" ", 
                    "",
                    contrasts[i])))
  
  bump.design           <- design
  bump.samples          <- strsplit(contrasts[i], " - ")[[1]]
  bump.samples.1        <- bump.samples[1]
  bump.samples.2        <- bump.samples[2]
  
  myannotation       <- cpg.annotate(minfi::getM(lumi.norm.filtered), 
                                     analysis.type = "differential", 
                                     design        = bump.design, 
                                     coef          = contrasts[i],
                                     contrasts     = T, 
                                     cont.matrix   = cont_mat)
  dmrcoutput         <- dmrcate(myannotation)
  dmrresults         <- dmrcoutput$results
  dmrresults$betaAfc <- abs(dmrresults$maxbetafc)
  dmrresults         <- dmrresults[with(dmrresults, order(meanpval)), ]
  dmrresults         <- dmrresults[dmrresults$betaAfc > 0.1,]
  dmrresults         <- dmrresults[with(dmrresults, order(-no.probes)), ]
  dmrresults         <- dmrresults[with(dmrresults, order(-no.probes)), ]
  
  write.csv(dmrresults,
            file=paste0(bump.samples.1,
                        "_",
                        bump.samples.2,
                        "_DMRcate.csv"))

# myannotationv      <- cpg.annotate(getM(bump.in), 
#                                    analysis.type = "variability", 
#                                    design        = bump.design, 
#                                    coef          = 2)
# dmrcoutputv        <- dmrcate(myannotation)
##'-----------------------------------------------------------------------------------------#





##'Run DMRcate
##'-----------------------------------------------------------------------------------------#
if(nrow(dmrresults) > 10) {print_cut <- 10} else {print_cut <- nrow(dmrresults)}
for(j in 1:print_cut) { 

bump.offset           <- 1000
bump.row              <- j
bump.pos.offset       <- 80
bump.coords           <- dmrresults$hg19coord[bump.row]
bump.chr              <- strsplit(bump.coords, ':')[[1]][1]
bump.start            <- as.numeric(strsplit(bump.coords, '-|:')[[1]][2])
bump.end              <- as.numeric(strsplit(bump.coords, '-')[[1]][2])
bump.annotation       <- annotation 
bump.annotation.in    <- bump.annotation[bump.annotation$chr == bump.chr   &
                                         bump.annotation$pos >= (bump.start - bump.offset) &
                                         bump.annotation$pos <= (bump.end   + bump.offset),]
bump.subset           <- betaAdj[match(rownames(bump.annotation.in),
                                       rownames(betaAdj)),
                                 grep(paste(strsplit(contrasts[i], 
                                                     " - ")[[1]], 
                                            collapse="|"),
                                      pData(lumi.norm.filtered)$SampleType)]
df                    <- melt(bump.subset)
df                    <- cbind(df,
                               bump.annotation.in[match(df$Var1, 
                                                        rownames(bump.annotation.in)),])
df                    <- cbind(df,
                               pData(lumi.norm.filtered)[match(df$Var2, 
                                                               rownames(pData(lumi.norm.filtered))),])
colnames(df)[1:3]     <- c("CpG", "Sample", "Beta")
df                    <- as.data.frame(df)
bump.pos.limit        <- c(min(df$pos)  - bump.pos.offset,
                           max(df$pos))

bump.pos.limit        <- c(bump.start  - bump.offset,
                           bump.end    + bump.offset)


##Get Transcript Annotation
# geneID_tmp             <- crunch(transcriptdb, 
#                                  which=GRanges(seqnames = as.numeric(gsub("chr", 
#                                                                           "", 
#                                                                           bump.chr)), 
#                                                IRanges(bump.start - bump.offset, 
#                                                        bump.end   + bump.offset)))

geneID_tmp             <- transcriptsByOverlaps(transcriptdb, 
                                                GRanges(seqnames = as.numeric(gsub("chr", 
                                                                                   "", 
                                                                                   bump.chr)), 
                                                        IRanges(bump.start - bump.offset, 
                                                                bump.end   + bump.offset)))

geneID_in              <- subsetByOverlaps(geneID_tmp,
                                           GRanges(seqnames = as.numeric(gsub("chr", 
                                                                              "", 
                                                                              bump.chr)), 
                                                   IRanges(bump.start, 
                                                           bump.end)))

if(length(geneID_in) > 0) {
  annotation.in         <- getBM(attributes = c("ensembl_gene_id",
                                                "external_gene_name"),
                                 filters    = "ensembl_transcript_id",
                                 values     = as.vector(geneID_in@elementMetadata@listData$tx_name[1]),
                                 mart       = ensembl)
  annotation.in         <- annotation.in$external_gene_name[1]
} else if(length(geneID_tmp) > 0) {
  annotation.in         <- getBM(attributes = c("ensembl_gene_id",
                                                "external_gene_name"),
                                 filters    = "ensembl_transcript_id",
                                 values     = as.vector(geneID_tmp@elementMetadata@listData$tx_name[1]),
                                 mart       = ensembl)
  annotation.in         <- annotation.in$external_gene_name[1]
  annotation.in         <- paste0("Close to ", annotation.in)
} else {
  annotation.in <- "No Data"
}
##
##

if(length(geneID_tmp) == 0) {
  p1_empty <- data.frame(x=bump.pos.limit, y=c(4,5,6,7))
  p1 <- ggplot(p1_empty, aes(x=x, y=y)) + 
    geom_blank() +
    theme_bw() + 
    scale_x_continuous(limits = bump.pos.limit) +
    geom_vline(xintercept     = bump.start, 
               colour         = "red")  +
    geom_vline(xintercept     = bump.end, 
               colour         = "red") +
    labs(x="",
         y="") +
    theme(axis.text.y=element_blank())
} else {
  trans         <- crunch(transcriptdb, which=geneID_tmp)
  gr1           <- split(trans, 
                         trans$tx_name)
  p1            <- ggbio::autoplot(gr1)                  + 
    theme_bw()                                           +
    geom_vline(xintercept     = bump.start, 
               colour         = "red")                   +
    geom_vline(xintercept     = bump.end, 
               colour         = "red")
  p1 <- p1@ggplot
}



##Plot1 - Betas
gg1 <- ggplot(df, 
              aes(x         = pos, 
                  y         = Beta, 
                  colour    = SampleType))                +
  geom_point()                                            + 
  scale_colour_hue(l        = 50)                         +
  geom_vline(xintercept     = c(bump.start,
                                bump.end),
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
  ggtitle(paste0("DMR ",
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
  geom_vline(xintercept      = c(bump.start,
                                 bump.end), 
             colour          = "red")
##


ggOut <- plot_grid(gg1, 
                   gg2,
                   p1,
                   align       = 'v', 
                   labels      = c('A', 'B', 'C'),
                   nrow        = 3,
                   rel_heights = c(4, 1.5, 3))
ggOut

save_plot(paste0(bump.samples.1,
                 "_",
                 bump.samples.2,
                 "_Idx_",
                 j,
                 "_DMR_.png"),
          ggOut,
          base_height = 10,
          base_width  = 12.53,
          dpi         = 300)
}
setwd("../")
}
##'-----------------------------------------------------------------------------------------#
setwd(current_dir)



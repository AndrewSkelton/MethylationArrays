library(dplyr)
library(readr)
# data_in <- c("cg10163220","cg16416184","cg02908942","cg11348345","cg03311459","cg03808158",
#              "cg25130672","cg12411093","cg01583034","cg03313212","cg05119480","cg04578997",
#              "cg21548032","cg11425280","cg07175149","cg12970542","cg16162324","cg11739675",
#              "cg01940855","cg23855505","cg22260952","cg16618104","cg06647068","cg12883014",
#              "cg09038971","cg05900122","cg06893569","cg20278154","cg13591408","cg02575406",
#              "cg03318667","cg24733614","cg22827210","cg16861964","cg10441497","cg19047077",
#              "cg07295775","cg07696842","cg05957190","cg07911905","cg01203340","cg11117177",
#              "cg07161429","cg06930722","cg09237846","cg20094699","cg15404019","cg26537335",
#              "cg01179256","cg20883723","cg02376674","cg22673901","cg14537856","cg11304240",
#              "cg17844339","cg12529671","cg10301869","cg22154485","cg04792380","cg25635544",
#              "cg00608661")
sampletypes    <- c("Preserved_Control_Hip","Preserved_OA_Hip","Preserved_OA_Knee")
# sample_in    <- "Preserved_OA_Hip"
# sample_in    <- "Preserved_OA_Knee"
setwd("Correlation_Tests/")
regions_in   <- read_tsv("../Scripts/Correlation_Regions.txt") %>% 
                as.data.frame %>% 
                mutate(Region = gsub(",","",Region),
                       Region = gsub(" ","",Region),
                       CHR    = gsub(":.*$","",Region),
                       START  = gsub("-.*$", "", Region),
                       START  = sub(".*:", "", START),
                       STOP   = sub(".*-", "", Region)) %>% 
                mutate_each(funs(as.numeric), START, STOP) %>% 
                mutate(LEN = STOP-START)

for(i in c(1:nrow(regions_in))){
  anno_tmp <- annotation %>% 
              filter(chr == regions_in$CHR[i] & pos > regions_in$START[i] & pos < regions_in$STOP[i])
  data_in  <- anno_tmp$Name
  for(j in sampletypes){
    mat_in       <- betaAdj[match(data_in, rownames(betaAdj)),
                            grep(j, pData(lumi.norm.filtered)$SampleType)]
    permutations <- combn(data_in, 2) %>% t %>% as.data.frame %>% mutate(Combo = paste0(V1, "_", V2))
    df_out       <- data.frame()
    
    for(k in 1:nrow(permutations)) {
      cpg_a   <- as.numeric(mat_in[as.vector(permutations[k,1]),])
      cpg_b   <- as.numeric(mat_in[as.vector(permutations[k,2]),])
      cor_out <- cor.test(cpg_a, cpg_b)
      df_out  <- bind_rows(df_out , c(permutations[k,3], 
                                      cor_out$estimate, 
                                      cor_out$p.value) %>% t %>% as.data.frame)
    }
    colnames(df_out) <- c("Pairs", "r", "pval")
    df_out$Sample    <- j
    df_out           <- df_out %>% arrange(pval)
    write.csv(df_out, file = paste0(j,"_",regions_in$Name[i],"_Correlation_Tests.csv"))
    print(j)
  }
  print(regions_in$Name[i])
}

df_out <- c()
anno_lite <- annotation %>% dplyr::select(chr,pos,strand,Name,UCSC_RefGene_Name)
files_in  <- list.files(pattern = "*.csv")
for(i in files_in){
  file_in      <- read_csv(i) %>% 
                  as.data.frame
  file_in$AdjP <- p.adjust(file_in$pval, method = "BH")
  file_in      <- file_in %>% 
                  mutate(CpG_2 = sub(".*_", "", Pairs),
                         CpG_1 = gsub("_.*$","",Pairs)) %>% 
                  left_join(anno_lite, by = c("CpG_1"="Name")) %>% 
                  left_join(anno_lite, by = c("CpG_2"="Name")) 
  colnames(file_in) <- colnames(file_in) %>% gsub(".x",".CpG_1",.)
  colnames(file_in) <- colnames(file_in) %>% gsub(".y",".CpG_2",.)
  file_in           <- file_in %>% mutate(LEN = `pos.CpG_2` - `pos.CpG_1`) %>% arrange(desc(LEN))
  foo    <- paste0(i,": ",file_in$LEN[1])
  df_out <- rbind(df_out, data.frame(File = i, LEN = file_in$LEN[1]))
  write_csv(file_in, path = paste0("Filtered/",i))
}

df_out %>% arrange(desc(LEN))

for(i in c(1:nrow(regions_in))){
  anno_tmp <- annotation %>% 
              filter(chr == regions_in$CHR[i] & pos > regions_in$START[i] & pos < regions_in$STOP[i])
  data_in  <- anno_tmp$Name
  for(j in sampletypes){
    mat_in       <- betaAdj[match(data_in, rownames(betaAdj)),
                            grep(j, pData(lumi.norm.filtered)$SampleType)] %>% 
                    as.data.frame %>% add_rownames("CpG")
    write_csv(mat_in, path = paste0("Betas/",j,"_",regions_in$Name[i],".csv"))
  }
}

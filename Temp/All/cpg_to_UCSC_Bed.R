setwd("~/Loughlin/Collaboration/All/Analysis_7/DM/")

# dm <- read.csv("Lesioned_OA_Knee-Preserved_OA_Knee/Lesioned_OA_Knee-Preserved_OA_Knee_DM.csv")
# dm <- read.csv("Preserved_Control_Hip-Preserved_OA_Hip/Preserved_Control_Hip-Preserved_OA_Hip_DM.csv")
# dm <- read.csv("Preserved_OA_Hip-Lesioned_OA_Hip/Preserved_OA_Hip-Lesioned_OA_Hip_DM.csv")
# dm <- read.csv("Preserved_OA_Hip-Preserved_OA_Knee/Preserved_OA_Hip-Preserved_OA_Knee_DM.csv")
dm         <- list(Preserved_OA_Hip_Lesioned_OA_Hip,
                   Preserved_OA_Hip_Preserved_Control_Hip,
                   Preserved_OA_Hip_Preserved_OA_Knee, 
                   Preserved_OA_Knee_Lesioned_OA_Knee)
dm_name    <- c("Preserved_OA_Hip-Lesioned_OA_Hip",
                "Preserved_OA_Hip-Preserved_Control_Hip",
                "Preserved_OA_Hip-Preserved_OA_Knee",
                "Preserved_OA_Knee-Lesioned_OA_Knee")
bed_name   <- c("Pre Vs Les OA Hip",
                "OA Hip Vs NOF",
                "OA Hip Vs OA Knee",
                "Pre OA Knee Vs Les OA Knee")
for(i in 1:length(dm)) {
  anno_in <- annotation[match(dm[[i]]$Name,
                              rownames(annotation)),]
  
  out        <- c()
  out$chr    <- dm[[i]]$chr
  out        <- as.data.frame(out)
  out$start  <- (dm[[i]]$pos - 1)
  out$end    <- (dm[[i]]$pos)
  out$name   <- paste0(dm[[i]]$Name, "_DeltaBeta_", dm[[i]]$`B-FC` %>% round(3))
  out$score  <- as.numeric(dm[[i]]$`B-FC`)
  out$strand <- dm[[i]]$strand
  out$box_s  <- out$start
  out$box_e  <- out$end
  out$colour <- as.vector("255,0,0")
  out[out$score < 0,]$colour <- "0,43,255" 
  
  write.table(out, 
              sep       = "\t", 
              row.names = F, 
              col.names = F, 
              quote     = F,   
              file      = paste0(dm_name[i],".bed"))
  exe <- paste0("sed -i '1s/^/track name=\"", 
                bed_name[i], "\" description=\"",
                bed_name[i], "\" visibility=3 ",
                "itemRgb=\"On\"",
                "\\n/' ", 
                dm_name[i], ".bed")
  system(exe)
}


library(openxlsx)
wb          <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Annotation")
openxlsx::writeData(wb, "Annotation", annotation %>% select(Name,chr,pos,strand))
openxlsx::saveWorkbook(wb, "Annotation_In.xlsx", overwrite = T)








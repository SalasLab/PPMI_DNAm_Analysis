library(readxl)
library(openxlsx)
library(rasterpdf)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/betas_CHR.RDA")
CHRs <- names(cpg_per_chr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# PD vs HC
parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC"

for (timept in c("BL","V04","V06","V08")){
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("CrossSectional")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
  results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
  results_Many <- results_Many[order(results_Many$adj.P.Val),]
  sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
  colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
  sig_df_summ$Model <- c("Unadj", "NLR", "Many")
  rownames(sig_df_summ) <- sig_df_summ$Mode
  sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
  sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
  sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
  sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
  sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
  sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
  setwd(parent_dir)
  results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
  results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
  results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
  save_data <- list("Summary"=sig_df_summ,
                    "Unadj"=data.frame(results_Unadj),
                    "NLR"=data.frame(results_NLR),
                    "Many"=data.frame(results_Many))
  write.xlsx(x = save_data,
             file = paste0("PDvsHC_",timept,"_full.xlsx"))
}


first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
results_Many <- results_Many[order(results_Many$adj.P.Val),]
sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
sig_df_summ$Model <- c("Unadj", "NLR", "Many")
rownames(sig_df_summ) <- sig_df_summ$Mode
sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
setwd(parent_dir)
results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
save_data <- list("Summary"=sig_df_summ,
                  "Unadj"=data.frame(results_Unadj),
                  "NLR"=data.frame(results_NLR),
                  "Many"=data.frame(results_Many))
write.xlsx(x = save_data,
           file = "PDvsHC_LONGIT_full.xlsx")

# first <- TRUE
# for (CHR in CHRs){
#   print(CHR)
#   setwd(parent_dir)
#   setwd(CHR)
#   setwd("Longit")
#   tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Unadj") 
#   tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "NLR") 
#   tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Many") 
#   if (first == TRUE){
#     results_Unadj <- tmp_results_Unadj
#     results_NLR<- tmp_results_NLR
#     results_Many<- tmp_results_Many
#     first <- FALSE
#   } else {
#     results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
#     results_NLR <- rbind(results_NLR, tmp_results_NLR)
#     results_Many <- rbind(results_Many, tmp_results_Many)
#   }
# }
# results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
# results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
# results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
# results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
# results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
# results_Many <- results_Many[order(results_Many$adj.P.Val),]
# sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
# colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
# sig_df_summ$Model <- c("Unadj", "NLR", "Many")
# rownames(sig_df_summ) <- sig_df_summ$Mode
# sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
# sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
# sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
# sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
# sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
# sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
# setwd(parent_dir)
# results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
# results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
# results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
# save_data <- list("Summary"=sig_df_summ,
#                   "Unadj"=data.frame(results_Unadj),
#                   "NLR"=data.frame(results_NLR),
#                   "Many"=data.frame(results_Many))
# write.xlsx(x = save_data,
#            file = "ProdvsHC_LONGIT_subset.xlsx")
# 


# Prod vs HC
parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC"

for (timept in c("BL","V04","V06","V08")){
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("CrossSectional")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_", timept,"_fulldata_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
  results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
  results_Many <- results_Many[order(results_Many$adj.P.Val),]
  sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
  colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
  sig_df_summ$Model <- c("Unadj", "NLR", "Many")
  rownames(sig_df_summ) <- sig_df_summ$Mode
  sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
  sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
  sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
  sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
  sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
  sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
  setwd(parent_dir)
  results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
  results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
  results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
  save_data <- list("Summary"=sig_df_summ,
                    "Unadj"=data.frame(results_Unadj),
                    "NLR"=data.frame(results_NLR),
                    "Many"=data.frame(results_Many))
  write.xlsx(x = save_data,
             file = paste0("ProdvsHC_",timept,"_full.xlsx"))
}

for (timept in c("BL","V04","V06","V08")){
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("CrossSectional")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_", timept,"_subsetdata_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_", timept,"_subsetdata_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_", timept,"_subsetdata_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
  results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
  results_Many <- results_Many[order(results_Many$adj.P.Val),]
  sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
  colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
  sig_df_summ$Model <- c("Unadj", "NLR", "Many")
  rownames(sig_df_summ) <- sig_df_summ$Mode
  sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
  sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
  sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
  sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
  sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
  sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
  setwd(parent_dir)
  results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
  results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
  results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
  save_data <- list("Summary"=sig_df_summ,
                    "Unadj"=data.frame(results_Unadj),
                    "NLR"=data.frame(results_NLR),
                    "Many"=data.frame(results_Many))
  write.xlsx(x = save_data,
             file = paste0("ProdvsHC_",timept,"_subset.xlsx"))
}

first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
results_Many <- results_Many[order(results_Many$adj.P.Val),]
sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
sig_df_summ$Model <- c("Unadj", "NLR", "Many")
rownames(sig_df_summ) <- sig_df_summ$Mode
sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
setwd(parent_dir)
results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
save_data <- list("Summary"=sig_df_summ,
                  "Unadj"=data.frame(results_Unadj),
                  "NLR"=data.frame(results_NLR),
                  "Many"=data.frame(results_Many))
write.xlsx(x = save_data,
           file = "ProdvsHC_LONGIT_full.xlsx")

first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
results_Many <- results_Many[order(results_Many$adj.P.Val),]
sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
sig_df_summ$Model <- c("Unadj", "NLR", "Many")
rownames(sig_df_summ) <- sig_df_summ$Mode
sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
setwd(parent_dir)
results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
save_data <- list("Summary"=sig_df_summ,
                  "Unadj"=data.frame(results_Unadj),
                  "NLR"=data.frame(results_NLR),
                  "Many"=data.frame(results_Many))
write.xlsx(x = save_data,
           file = "ProdvsHC_LONGIT_subset.xlsx")


# Prod Noncon vs Prod Con
parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdConvsProdNoncon"
first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
results_Many <- results_Many[order(results_Many$adj.P.Val),]
sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
sig_df_summ$Model <- c("Unadj", "NLR", "Many")
rownames(sig_df_summ) <- sig_df_summ$Mode
sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
setwd(parent_dir)
results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
save_data <- list("Summary"=sig_df_summ,
                  "Unadj"=data.frame(results_Unadj),
                  "NLR"=data.frame(results_NLR),
                  "Many"=data.frame(results_Many))
write.xlsx(x = save_data,
           file = "ProdConvsProdNoncon_LONGIT_full.xlsx")

first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_subsetdata_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
results_Many <- results_Many[order(results_Many$adj.P.Val),]
sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
sig_df_summ$Model <- c("Unadj", "NLR", "Many")
rownames(sig_df_summ) <- sig_df_summ$Mode
sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
setwd(parent_dir)
results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
save_data <- list("Summary"=sig_df_summ,
                  "Unadj"=data.frame(results_Unadj),
                  "NLR"=data.frame(results_NLR),
                  "Many"=data.frame(results_Many))
write.xlsx(x = save_data,
           file = "ProdConvsProdNoncon_LONGIT_subset.xlsx")














############## Volcano Plots
library(readxl)
library(openxlsx)
library(rasterpdf)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/betas_CHR.RDA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/CpG_Gene_Names.RDA")
CHRs <- names(cpg_per_chr)
colnames(cpg_gene_names)[1] <- "CpG"

###### Prod vs HC
parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdConvsProdNoncon"
first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj <- left_join(results_Unadj,cpg_gene_names)
results_NLR <- left_join(results_NLR,cpg_gene_names)
results_Many <- left_join(results_Many,cpg_gene_names)
results_Unadj$UCSC_RefGene_Name[results_Unadj$UCSC_RefGene_Name==""] <- results_Unadj$CpG[results_Unadj$UCSC_RefGene_Name==""]
results_NLR$UCSC_RefGene_Name[results_NLR$UCSC_RefGene_Name==""] <- results_NLR$CpG[results_NLR$UCSC_RefGene_Name==""]
results_Many$UCSC_RefGene_Name[results_Many$UCSC_RefGene_Name==""] <- results_Many$CpG[results_Many$UCSC_RefGene_Name==""]

results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj$log10pVal <- -log10(results_Unadj$adj.P.Val)
results_NLR$log10pVal <- -log10(results_NLR$adj.P.Val)
results_Many$log10pVal <- -log10(results_Many$adj.P.Val)

results_Unadj$tmpcat <- ifelse(results_Unadj$adj.P.Val<0.05, "Sig","No Difference")
results_Unadj$tmpcat <- ifelse(results_Unadj$tmpcat == "No Difference", "No Difference",
                               ifelse(results_Unadj$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_Unadj$delta_beta  > 0.1, "Hypermethylated","No Difference")))
results_NLR$tmpcat <- ifelse(results_NLR$adj.P.Val<0.05, "Sig","No Difference")
results_NLR$tmpcat <- ifelse(results_NLR$tmpcat == "No Difference", "No Difference",
                             ifelse(results_NLR$delta_beta  < -0.1, "Hypomethylated",
                                    ifelse(results_NLR$delta_beta  > 0.1, "Hypermethylated","No Difference")))
results_Many$tmpcat <- ifelse(results_Many$adj.P.Val<0.05, "Sig","No Difference")
results_Many$tmpcat <- ifelse(results_Many$tmpcat == "No Difference", "No Difference",
                              ifelse(results_Many$delta_beta  < -0.1, "Hypomethylated",
                                     ifelse(results_Many$delta_beta  > 0.1, "Hypermethylated","No Difference")))


hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")

top_cpg_data <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
p1 <- ggplot(results_Unadj, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,2.5)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "red", size = 0.8, linetype="dashed") +
 geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name),  size = 2)

top_cpg_data <- results_NLR[results_NLR$adj.P.Val<0.05,]
p2 <- ggplot(results_NLR, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,2.5)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "red", size = 0.8, linetype="dashed") +
 geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name),  size = 2)

top_cpg_data <- results_Many[results_Many$adj.P.Val<0.05,]
tmp <- c("cg13892688","SLC6A15","KATNAL1")
top_cpg_data$UCSC_RefGene_Name <- tmp
p3 <- ggplot(results_Many, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,2.5)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "red", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "red", size = 0.8, linetype="dashed") +
 geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)


raster_pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdConvsProdNoncon/LongitVolcanos_raster3.PDF", height = 3.5, width = 4)
p1
p2
p3
dev.off()

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdConvsProdNoncon/LongitVolcanos3.PDF", height = 3.5, width = 4)
p1
p2
p3
dev.off()

###### Prod vs HC
parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC"
first <- TRUE
for (CHR in CHRs){
  print(CHR)
  setwd(parent_dir)
  setwd(CHR)
  setwd("Longit")
  tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Unadj") 
  tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "NLR") 
  tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Many") 
  if (first == TRUE){
    results_Unadj <- tmp_results_Unadj
    results_NLR<- tmp_results_NLR
    results_Many<- tmp_results_Many
    first <- FALSE
  } else {
    results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
    results_NLR <- rbind(results_NLR, tmp_results_NLR)
    results_Many <- rbind(results_Many, tmp_results_Many)
  }
}
results_Unadj <- left_join(results_Unadj,cpg_gene_names)
results_NLR <- left_join(results_NLR,cpg_gene_names)
results_Many <- left_join(results_Many,cpg_gene_names)
results_Unadj$UCSC_RefGene_Name[results_Unadj$UCSC_RefGene_Name==""] <- results_Unadj$CpG[results_Unadj$UCSC_RefGene_Name==""]
results_NLR$UCSC_RefGene_Name[results_NLR$UCSC_RefGene_Name==""] <- results_NLR$CpG[results_NLR$UCSC_RefGene_Name==""]
results_Many$UCSC_RefGene_Name[results_Many$UCSC_RefGene_Name==""] <- results_Many$CpG[results_Many$UCSC_RefGene_Name==""]

results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
results_Unadj$log10pVal <- -log10(results_Unadj$adj.P.Val)
results_NLR$log10pVal <- -log10(results_NLR$adj.P.Val)
results_Many$log10pVal <- -log10(results_Many$adj.P.Val)

results_Unadj$tmpcat <- ifelse(results_Unadj$adj.P.Val<0.05, "Sig","No Difference")
results_Unadj$tmpcat <- ifelse(results_Unadj$tmpcat == "No Difference", "No Difference",
                               ifelse(results_Unadj$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_Unadj$delta_beta  > 0.1, "Hypermethylated","No Difference")))
results_NLR$tmpcat <- ifelse(results_NLR$adj.P.Val<0.05, "Sig","No Difference")
results_NLR$tmpcat <- ifelse(results_NLR$tmpcat == "No Difference", "No Difference",
                               ifelse(results_NLR$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_NLR$delta_beta  > 0.1, "Hypermethylated","No Difference")))
results_Many$tmpcat <- ifelse(results_Many$adj.P.Val<0.05, "Sig","No Difference")
results_Many$tmpcat <- ifelse(results_Many$tmpcat == "No Difference", "No Difference",
                               ifelse(results_Many$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_Many$delta_beta  > 0.1, "Hypermethylated","No Difference")))


hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")

top_cpg_data <- results_Unadj[results_Unadj$tmpcat != "No Difference",]
tmp <- c("cg11173636","PYROXD2","KIRREL3","SHANK2","ACOT1/HEATR4","RAD51B","L3HYPDH","PCNX1","cg06115838","STK4","LSG1","cg18845950","ERICH1","DENND4C")
top_cpg_data$UCSC_RefGene_Name <- tmp
top_cpg_data <- top_cpg_data[order(top_cpg_data$adj.P.Val),]
top_cpg_data <- top_cpg_data[1:10,]
p1 <- ggplot(results_Unadj, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,23)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name),  size = 2)

top_cpg_data <- results_NLR[results_NLR$tmpcat != "No Difference",]
  p2 <- ggplot(results_NLR, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,23)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "#3971B5", size = 0.8, linetype="dashed")# +
 # geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = CpG),  size = 1)
  
  top_cpg_data <- results_Many[results_Many$tmpcat != "No Difference",]
  tmp <- c("cg11173636","PYROXD2","KIRREL3","SHANK2","ACOT1/HEATR4","RAD51B","L3HYPDH","PCNX1","cg06115838","STK4","LSG1","cg18845950","ERICH1","DENND4C")
  top_cpg_data$UCSC_RefGene_Name <- tmp
  top_cpg_data <- top_cpg_data[order(top_cpg_data$adj.P.Val),]
  top_cpg_data <- top_cpg_data[1:10,]
  p3 <- ggplot(results_Many, aes(delta_beta, log10pVal)) +  
  geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
  xlab(expression(paste(Delta,Beta))) + 
  ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
  ylim(c(0,23)) + xlim(c(-0.28, 0.28)) +
  scale_fill_manual(values = hyper_hypo_colors) +
  geom_hline(yintercept = -log10(0.05), color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = -0.1, color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_vline(xintercept = 0.1, color = "#3971B5", size = 0.8, linetype="dashed") +
  geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  
  ProdvsHC_BL_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ProdvsHC_BL_full.xlsx",  sheet = "Summary")
  ProdvsHC_V04_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ProdvsHC_V04_full.xlsx",  sheet = "Summary")
  ProdvsHC_V06_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ProdvsHC_V06_full.xlsx",  sheet = "Summary")
  ProdvsHC_V08_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ProdvsHC_V08_full.xlsx",  sheet = "Summary")
  ProdvsHC_LONGIT_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ProdvsHC_LONGIT_full.xlsx",  sheet = "Summary")
  
  all_data <- rbind(ProdvsHC_BL_full[ProdvsHC_BL_full$Model == "Many",],
                    ProdvsHC_V04_full[ProdvsHC_V04_full$Model == "Many",],
                    ProdvsHC_V06_full[ProdvsHC_V06_full$Model == "Many",],
                    ProdvsHC_V08_full[ProdvsHC_V08_full$Model == "Many",],
                    ProdvsHC_LONGIT_full[ProdvsHC_LONGIT_full$Model == "Many",])
  all_data$`Time Point` <- c("BL","Y1","Y2","Y3","LONGIT")
  all_data$`Time Point` <- factor( all_data$`Time Point`, levels = c("BL","Y1","Y2","Y3","LONGIT"))
  all_data <- melt(all_data, id.vars = c("Model","Time Point"))
  
  hyper_hypo_colors <- c(Hypo = "yellow", Hyper = "blue")
  
  p4 <- ggplot(all_data, aes(fill = variable , y = value, x = `Time Point`)) + geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = hyper_hypo_colors, name = "") + scale_y_log10() +
    ylab("Number of CpGs in which Padj < 0.05") + theme_minimal() 
  
  
  raster_pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/LongitVolcanos_raster4.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/LongitVolcanos4.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/ManyCompare.PDF", height = 3.5, width = 4)
  p4
  dev.off()
  
  
  
  
  ###### PD vs HC
  parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC"
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("Longit")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_fulldata_newbatch_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  results_Unadj <- left_join(results_Unadj,cpg_gene_names)
  results_NLR <- left_join(results_NLR,cpg_gene_names)
  results_Many <- left_join(results_Many,cpg_gene_names)
  results_Unadj$UCSC_RefGene_Name[results_Unadj$UCSC_RefGene_Name==""] <- results_Unadj$CpG[results_Unadj$UCSC_RefGene_Name==""]
  results_NLR$UCSC_RefGene_Name[results_NLR$UCSC_RefGene_Name==""] <- results_NLR$CpG[results_NLR$UCSC_RefGene_Name==""]
  results_Many$UCSC_RefGene_Name[results_Many$UCSC_RefGene_Name==""] <- results_Many$CpG[results_Many$UCSC_RefGene_Name==""]
  
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj$log10pVal <- -log10(results_Unadj$adj.P.Val)
  results_NLR$log10pVal <- -log10(results_NLR$adj.P.Val)
  results_Many$log10pVal <- -log10(results_Many$adj.P.Val)
  
  results_Unadj$tmpcat <- ifelse(results_Unadj$adj.P.Val<0.05, "Sig","No Difference")
  results_Unadj$tmpcat <- ifelse(results_Unadj$tmpcat == "No Difference", "No Difference",
                                 ifelse(results_Unadj$delta_beta  < -0.1, "Hypomethylated",
                                        ifelse(results_Unadj$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  results_NLR$tmpcat <- ifelse(results_NLR$adj.P.Val<0.05, "Sig","No Difference")
  results_NLR$tmpcat <- ifelse(results_NLR$tmpcat == "No Difference", "No Difference",
                               ifelse(results_NLR$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_NLR$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  results_Many$tmpcat <- ifelse(results_Many$adj.P.Val<0.05, "Sig","No Difference")
  results_Many$tmpcat <- ifelse(results_Many$tmpcat == "No Difference", "No Difference",
                                ifelse(results_Many$delta_beta  < -0.1, "Hypomethylated",
                                       ifelse(results_Many$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  
  
  hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")
  
  top_cpg_data <- results_Unadj[results_Unadj$tmpcat != "No Difference",]
  p1 <- ggplot(results_Unadj, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,7)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
  geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  top_cpg_data <- results_NLR[results_NLR$tmpcat != "No Difference",]
  p2 <- ggplot(results_NLR, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,7)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
  geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  top_cpg_data <- results_Many[results_Many$tmpcat != "No Difference",]
  p3 <- ggplot(results_Many, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,7)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
  geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  
  PDvsHC_BL_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/PDvsHC_BL_full.xlsx",  sheet = "Summary")
  PDvsHC_V04_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/PDvsHC_V04_full.xlsx",  sheet = "Summary")
  PDvsHC_V06_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/PDvsHC_V06_full.xlsx",  sheet = "Summary")
  PDvsHC_V08_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/PDvsHC_V08_full.xlsx",  sheet = "Summary")
  PDvsHC_LONGIT_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/PDvsHC_LONGIT_full.xlsx",  sheet = "Summary")
  
  all_data <- rbind(PDvsHC_BL_full[PDvsHC_BL_full$Model == "Many",],
                    PDvsHC_V04_full[PDvsHC_V04_full$Model == "Many",],
                    PDvsHC_V06_full[PDvsHC_V06_full$Model == "Many",],
                    PDvsHC_V08_full[PDvsHC_V08_full$Model == "Many",],
                    PDvsHC_LONGIT_full[PDvsHC_LONGIT_full$Model == "Many",])
  all_data$`Time Point` <- c("BL","Y1","Y2","Y3","LONGIT")
  all_data$`Time Point` <- factor( all_data$`Time Point`, levels = c("BL","Y1","Y2","Y3","LONGIT"))
  all_data <- melt(all_data, id.vars = c("Model","Time Point"))
  
  hyper_hypo_colors <- c(Hypo = "yellow", Hyper = "blue")
  
  p4 <- ggplot(all_data, aes(fill = variable , y = value, x = `Time Point`)) + geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = hyper_hypo_colors, name = "") +
    ylab("Number of CpGs in which Padj < 0.05") + theme_minimal() + scale_y_log10()
  
  
  raster_pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/LongitVolcanos_raster.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/LongitVolcanos.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsHC/ManyCompare.PDF", height = 3.5, width = 4)
  p4
  dev.off()
  
  
  
  ###### PD vs Prod
  parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd"
  
  for (timept in c("BL","V04","V06","V08")){
    first <- TRUE
    for (CHR in CHRs){
      print(CHR)
      setwd(parent_dir)
      setwd(CHR)
      setwd("CrossSectional")
      tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_", timept,"_PDvsProd_EWASResults.xlsx"), sheet = "Unadj") 
      tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_", timept,"_PDvsProd_EWASResults.xlsx"), sheet = "NLR") 
      tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_", timept,"_PDvsProd_EWASResults.xlsx"), sheet = "Many") 
      if (first == TRUE){
        results_Unadj <- tmp_results_Unadj
        results_NLR<- tmp_results_NLR
        results_Many<- tmp_results_Many
        first <- FALSE
      } else {
        results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
        results_NLR <- rbind(results_NLR, tmp_results_NLR)
        results_Many <- rbind(results_Many, tmp_results_Many)
      }
    }
    results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
    results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
    results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
    results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
    results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
    results_Many <- results_Many[order(results_Many$adj.P.Val),]
    sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
    colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
    sig_df_summ$Model <- c("Unadj", "NLR", "Many")
    rownames(sig_df_summ) <- sig_df_summ$Mode
    sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
    sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
    sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
    sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
    sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
    sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
    setwd(parent_dir)
    results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
    results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
    results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
    save_data <- list("Summary"=sig_df_summ,
                      "Unadj"=data.frame(results_Unadj),
                      "NLR"=data.frame(results_NLR),
                      "Many"=data.frame(results_Many))
    write.xlsx(x = save_data,
               file = paste0("PDvsProd_",timept,".xlsx"))
  }
  
      
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("Longit")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj <- results_Unadj[order(results_Unadj$adj.P.Val),]
  results_NLR <- results_NLR[order(results_NLR$adj.P.Val),]
  results_Many <- results_Many[order(results_Many$adj.P.Val),]
  sig_df_summ <- data.frame(matrix(0,nrow=3,ncol=3))
  colnames(sig_df_summ) <- c("Model","Hyper","Hypo")
  sig_df_summ$Model <- c("Unadj", "NLR", "Many")
  rownames(sig_df_summ) <- sig_df_summ$Mode
  sig_df_summ["Unadj", "Hypo"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta < 0)
  sig_df_summ["Unadj", "Hyper"] <- sum(results_Unadj$adj.P.Val < 0.05 & results_Unadj$delta_beta > 0)
  sig_df_summ["NLR", "Hypo"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta < 0)
  sig_df_summ["NLR", "Hyper"] <- sum(results_NLR$adj.P.Val < 0.05 & results_NLR$delta_beta > 0)
  sig_df_summ["Many", "Hypo"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta < 0)
  sig_df_summ["Many", "Hyper"] <- sum(results_Many$adj.P.Val < 0.05 & results_Many$delta_beta > 0)
  setwd(parent_dir)
  results_Unadj <- results_Unadj[results_Unadj$adj.P.Val<0.05,]
  results_NLR <- results_NLR[results_NLR$adj.P.Val<0.05,]
  results_Many <- results_Many[results_Many$adj.P.Val<0.05,]
  save_data <- list("Summary"=sig_df_summ,
                    "Unadj"=data.frame(results_Unadj),
                    "NLR"=data.frame(results_NLR),
                    "Many"=data.frame(results_Many))
  write.xlsx(x = save_data,
             file = "PDvsProd_LONGIT.xlsx")
  
  
  
  
  parent_dir <- "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd"
  first <- TRUE
  for (CHR in CHRs){
    print(CHR)
    setwd(parent_dir)
    setwd(CHR)
    setwd("Longit")
    tmp_results_Unadj <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "Unadj") 
    tmp_results_NLR <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "NLR") 
    tmp_results_Many <- read_excel(paste0(CHR, "_EWAS_LONGIT_PDvsProd_EWASResults.xlsx"), sheet = "Many") 
    if (first == TRUE){
      results_Unadj <- tmp_results_Unadj
      results_NLR<- tmp_results_NLR
      results_Many<- tmp_results_Many
      first <- FALSE
    } else {
      results_Unadj <- rbind(results_Unadj, tmp_results_Unadj)
      results_NLR <- rbind(results_NLR, tmp_results_NLR)
      results_Many <- rbind(results_Many, tmp_results_Many)
    }
  }
  library(dplyr)
  results_Unadj <- left_join(results_Unadj,cpg_gene_names)
  results_NLR <- left_join(results_NLR,cpg_gene_names)
  results_Many <- left_join(results_Many,cpg_gene_names)
  results_Unadj$UCSC_RefGene_Name[results_Unadj$UCSC_RefGene_Name==""] <- results_Unadj$CpG[results_Unadj$UCSC_RefGene_Name==""]
  results_NLR$UCSC_RefGene_Name[results_NLR$UCSC_RefGene_Name==""] <- results_NLR$CpG[results_NLR$UCSC_RefGene_Name==""]
  results_Many$UCSC_RefGene_Name[results_Many$UCSC_RefGene_Name==""] <- results_Many$CpG[results_Many$UCSC_RefGene_Name==""]
  
  results_Unadj$adj.P.Val <- p.adjust(results_Unadj$P.Value, method="fdr")
  results_NLR$adj.P.Val <- p.adjust(results_NLR$P.Value, method="fdr")
  results_Many$adj.P.Val <- p.adjust(results_Many$P.Value, method="fdr")
  results_Unadj$log10pVal <- -log10(results_Unadj$adj.P.Val)
  results_NLR$log10pVal <- -log10(results_NLR$adj.P.Val)
  results_Many$log10pVal <- -log10(results_Many$adj.P.Val)
  
  results_Unadj$tmpcat <- ifelse(results_Unadj$adj.P.Val<0.05, "Sig","No Difference")
  results_Unadj$tmpcat <- ifelse(results_Unadj$tmpcat == "No Difference", "No Difference",
                                 ifelse(results_Unadj$delta_beta  < -0.1, "Hypomethylated",
                                        ifelse(results_Unadj$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  results_NLR$tmpcat <- ifelse(results_NLR$adj.P.Val<0.05, "Sig","No Difference")
  results_NLR$tmpcat <- ifelse(results_NLR$tmpcat == "No Difference", "No Difference",
                               ifelse(results_NLR$delta_beta  < -0.1, "Hypomethylated",
                                      ifelse(results_NLR$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  results_Many$tmpcat <- ifelse(results_Many$adj.P.Val<0.05, "Sig","No Difference")
  results_Many$tmpcat <- ifelse(results_Many$tmpcat == "No Difference", "No Difference",
                                ifelse(results_Many$delta_beta  < -0.1, "Hypomethylated",
                                       ifelse(results_Many$delta_beta  > 0.1, "Hypermethylated","No Difference")))
  
  
  hyper_hypo_colors <- c(`No Difference` = "gray50", Hypomethylated = "yellow", Hypermethylated = "blue")
  
  top_cpg_data <- results_Unadj[results_Unadj$tmpcat != "No Difference",]
  p1 <- ggplot(results_Unadj, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,15.5)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  top_cpg_data <- results_NLR[results_NLR$tmpcat != "No Difference",]
  p2 <- ggplot(results_NLR, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,15.5)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  top_cpg_data <- results_Many[results_Many$tmpcat != "No Difference",]
  p3 <- ggplot(results_Many, aes(delta_beta, log10pVal)) +  
    geom_point(aes(fill = tmpcat),size = 9/5,shape=21, show.legend = FALSE) +
    xlab(expression(paste(Delta,Beta))) + 
    ylab(expression("-log"[10]*"(Padj)")) + theme_minimal() + 
    ylim(c(0,15.5)) + xlim(c(-0.28, 0.28)) +
    scale_fill_manual(values = hyper_hypo_colors) +
    geom_hline(yintercept = -log10(0.05), color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = -0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_vline(xintercept = 0.1, color = "#C56124", size = 0.8, linetype="dashed") +
    geom_label_repel(data = top_cpg_data, mapping = aes(delta_beta, log10pVal, label = UCSC_RefGene_Name), size = 2)
  
  
  PDvsHC_BL_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/PDvsProd_BL.xlsx",  sheet = "Summary")
  PDvsHC_V04_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/PDvsProd_V04.xlsx",  sheet = "Summary")
  PDvsHC_V06_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/PDvsProd_V06.xlsx",  sheet = "Summary")
  PDvsHC_V08_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/PDvsProd_V08.xlsx",  sheet = "Summary")
  PDvsHC_LONGIT_full <- read_excel("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/PDvsProd_LONGIT.xlsx",  sheet = "Summary")
  
  all_data <- rbind(PDvsHC_BL_full[PDvsHC_BL_full$Model == "Many",],
                    PDvsHC_V04_full[PDvsHC_V04_full$Model == "Many",],
                    PDvsHC_V06_full[PDvsHC_V06_full$Model == "Many",],
                    PDvsHC_V08_full[PDvsHC_V08_full$Model == "Many",],
                    PDvsHC_LONGIT_full[PDvsHC_LONGIT_full$Model == "Many",])
  all_data$`Time Point` <- c("BL","Y1","Y2","Y3","LONGIT")
  all_data$`Time Point` <- factor( all_data$`Time Point`, levels = c("BL","Y1","Y2","Y3","LONGIT"))
  all_data <- melt(all_data, id.vars = c("Model","Time Point"))
  
  hyper_hypo_colors <- c(Hypo = "yellow", Hyper = "blue")
  
  p4 <- ggplot(all_data, aes(fill = variable , y = value, x = `Time Point`)) + geom_bar(position="dodge", stat="identity") +
    scale_fill_manual(values = hyper_hypo_colors, name = "") +
    ylab("Number of CpGs in which Padj < 0.05") + theme_minimal() 
  
  
  raster_pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/LongitVolcanos_raster.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/LongitVolcanos.PDF", height = 3.5, width = 4)
  p1
  p2
  p3
  dev.off()
  
  pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/PDvsProd/ManyCompare.PDF", height = 3.5, width = 4)
  p4
  dev.off()
  
  
  
  
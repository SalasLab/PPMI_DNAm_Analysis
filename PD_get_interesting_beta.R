# get CpGs
cell.types <- c("Bas","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Eos","Mono","Neu","NK","Treg")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
pheno <- pheno[pheno$rmv_sample == 0,]
library(readr)
library(dplyr)
PT_status <- read_csv("/dartfs/rc/lab/S/SalasLab/PD/Raw_data/Participant_Status.csv", 
                      col_types = cols(PATNO = col_character()))
PT_status <- PT_status[, c(1,13:23)]
pheno <- left_join(pheno, PT_status, by = "PATNO")
pheno$Prod_Convert <- ifelse(pheno$PHENOCNV == 1, 1, 0)
pheno$Prod_Convert_SampleTime <- ifelse(is.na(pheno$DIAG1VIS), NA ,
                                        ifelse(pheno$DIAG1VIS == "V02", 6, 
                                               ifelse(pheno$DIAG1VIS == "V04", 12, 
                                                      ifelse(pheno$DIAG1VIS == "V06", 24, 
                                                             ifelse(pheno$DIAG1VIS == "V08", 36, 
                                                                    ifelse(pheno$DIAG1VIS == "V10", 48, 
                                                                           ifelse(pheno$DIAG1VIS == "V12", 60, NA)))))))

pheno$Prod_Convert_Relavtive <- ifelse(pheno$PHENOCNV == 1, ifelse(pheno$Sample_Time < pheno$Prod_Convert_SampleTime, "Pre-conversion","Converted" ) , NA)
pheno$Prod_conv_rmv <- ifelse(is.na(pheno$Prod_Convert_Relavtive), 0, 
                              ifelse(pheno$Prod_Convert_Relavtive == "Converted", 1, 0))

pheno <- pheno[pheno$Prod_conv_rmv == 0,]


chr_cpgs <- c("cg13102742", # PD Hit
              "cg02833127", # the rest are prod hit hypo
              "cg02746014",
              "cg21223075",
              "cg00685135",
              "cg03885684",
              "cg04833938",
              "cg18845950",
              "cg13799572",
              "cg14755254",
              "cg18683228",
              "cg02078724",
              "cg09157251",
              "cg06115838",
              "cg00088299",
              "cg02907150",
              "cg13892688",  #COnv hit
              "cg18159740",  #COnv hit
              "cg00154888" ,#COnv hit
              "cg11787544", # rest are prod hit hyper
              "cg16786756",
              "cg10576280",
              "cg16628641",
              "cg11173636",
              "cg00924943",
              "cg11523661",
              "cg22968327",
              "cg26690318",
              "cg01543583",
              "cg26251192",
              "cg11394338",
              "cg06612594"
) 
              
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/CpG_Gene_Names.RDA")
cpg_gene_names2 <- cpg_gene_names[chr_cpgs,]

first <- TRUE
for (i in unique(pheno$processing_batch)){
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i)
  setwd(set_name)
  print(set_name)
  betas <- loadRData(paste0(set_name, "_betas_filtered.RDA"))
  betas <- betas[,which(colnames(betas) %in% pheno$ID)]
  betas <- betas[chr_cpgs,]
  if (first == TRUE){
    agg_betas <- betas
    first <- FALSE
  } else {
    tmp_check <- sum(rownames(agg_betas) != rownames(betas))
    if (tmp_check != 0){
      message("There is some missing data")
      break
    }
    agg_betas <- cbind(agg_betas, betas)
  }
}




df <- pheno
df$DX_new_char_factor <- factor(df$DX_new_char, levels = c("HC","Prod","PD"))
p <- ggplot(df, aes(x = DX_new_char_factor, y = betas_GLIPR2, color = DX_new_char_factor))
p <- p + geom_boxplot(outlier.shape = NA, fill = light_cols, color = dark_cols) + geom_jitter(aes(color = Prod_Convert_Relavtive), size=0.5) + scale_color_manual(values = c("NA" = "black", "Converted" = "red", "Pre-conversion" = "green")) + theme_minimal()
p <- p + geom_hline(yintercept = 0.1, col = "red") + ggtitle(label ="cg00685135 - GLIPR2" )
p
# Updated QC for Parkinsons Disease


# 1. Sex mismatches
# 2. Median beta across entire study
# 3. ewasttools thresholds
# 4. looked at samples that are outliers w cell type

library(dplyr)
library(readr)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Loading in the data
for (i in 1:23){
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i)
  setwd(set_name)
  print(i)
  if (i == 1){
    pheno <- loadRData(paste0(set_name, "_pheno.RDA"))
    ctrl_metric <- loadRData(paste0(set_name, "_EWAStools_metrics.RDA"))
  } else {
    tmppheno <- loadRData(paste0(set_name, "_pheno.RDA"))
    pheno <- rbind(pheno, tmppheno)
    tmpctrl_metric <- loadRData(paste0(set_name, "_EWAStools_metrics.RDA"))
    ctrl_metric <- rbind(ctrl_metric,tmpctrl_metric)
  }
}
rm(tmppheno,tmpctrl_metric)
pheno$rmv_sample <- 0


# 1. Sex mismatches
pheno$rmv_sample_sex <- 0
pheno$is.Female <- ifelse(pheno$Gender_reported == "M", 0, 1)
table(pheno$ewastools_sex, useNA = "ifany")
tmpvector <- ifelse(is.na(pheno$ewastools_sex == "m"), 999, ifelse(pheno$ewastools_sex == "m", 0, 1))
tmpsamps <- which(pheno$is.Female != tmpvector & tmpvector != 999)
pheno[tmpsamps, "rmv_sample"] <- 1
pheno[tmpsamps, "rmv_sample_sex"] <- 1


# 2. Median beta
for (i in 1:23){
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i)
  setwd(set_name)
  print(set_name)
  
  if (i == 1){
    tmpprobestats <- loadRData(paste0(set_name, "_probeStatsfiltered.RDA"))
    median_beta <- tmpprobestats$Median_Beta
    mean_beta <- tmpprobestats$Mean_Beta
    n_cpg_missing <- tmpprobestats$NA_counts
  } else {
    tmpprobestats <- loadRData(paste0(set_name, "_probeStatsfiltered.RDA"))
    tmpmedian_beta <- tmpprobestats$Median_Beta
    tmpmean_beta <- tmpprobestats$Mean_Beta
    tmpn_cpg_missing <- tmpprobestats$NA_counts
    
    median_beta <- c(median_beta, tmpmedian_beta)
    mean_beta <- c(mean_beta, tmpmean_beta)
    n_cpg_missing <- c(n_cpg_missing, tmpn_cpg_missing)
  }
  
}
pheno$Median_beta <- median_beta
pheno$Mean_beta <- mean_beta
pheno$n_cpg_missing <- n_cpg_missing
IQR_median <- IQR(median_beta)
median_beta_val <- median(median_beta)

upper_limit <- median_beta_val + 3.5*IQR_median
lower_limit <- median_beta_val - 3.5*IQR_median

pheno$rmv_sample_median_beta <- ifelse(pheno$Median_beta > upper_limit, 1, ifelse(pheno$Median_beta < lower_limit, 1, 0))
tmpsamps <- which(pheno$rmv_sample_median_beta == 1)
pheno[tmpsamps, "rmv_sample"] <- 1

pheno$rmv_sample_missing_cpg <- ifelse(pheno$n_cpg_missing > 74383, 1,0)
tmpsamps <- which(pheno$rmv_sample_missing_cpg == 1)
pheno[tmpsamps, "rmv_sample"] <- 1


# 3. ewastools
pass_the_test <- ctrl_metric
any_failure <- matrix(data = 0, nrow = 17, ncol = 2)
rownames(any_failure) <- colnames(ctrl_metric)[1:17]
for (ctrl in colnames(pass_the_test)[1:17]){
  for (i in 1:nrow(pass_the_test)){
    if (ctrl %in% c("Staining.Green", "Staining.Red",
                    "Extension.Green", "Extension.Red",
                    "Non.polymorphic.Green", "Non.polymorphic.Red")){
      if (ctrl_metric[i, ctrl] < 5 | is.na(ctrl_metric[i, ctrl])){
        pass_the_test[i, ctrl] <- 1
      } else {
        pass_the_test[i, ctrl] <- 0
      }
      
    } else if (ctrl == "Restoration"){
      if (ctrl_metric[i, ctrl] < 0 | is.na(ctrl_metric[i, ctrl])){
        pass_the_test[i, ctrl] <- 1
      } else {
        pass_the_test[i, ctrl] <- 0
      }
      
    } else {
      if (ctrl_metric[i, ctrl] < 1 | is.na(ctrl_metric[i, ctrl])){ 
        pass_the_test[i, ctrl] <- 1
      } else {
        pass_the_test[i, ctrl] <- 0
      }
      
    }
    
  }
  
  any_failure[rownames(any_failure) == ctrl, 1] <-  sum(pass_the_test[,ctrl])
  if (any_failure[rownames(any_failure) == ctrl, 1] > 0){
    any_failure[rownames(any_failure) == ctrl, 2] <- 1
  }
}
colnames(any_failure) <- c("total.fails", "any.fails")
any_failure <- data.frame(any_failure)
failed_ctrls <- rownames(any_failure[any_failure$total.fails != 0,]) # staining.green, staining.red, specificity.I.Red, nonpolymorphic

pheno$rmv_sample_ctrls_stain_grn <- 0
pheno$rmv_sample_ctrls_stain_red <- 0
pheno$rmv_sample_ctrls_spec_I_red <- 0
pheno$rmv_sample_ctrls_nonpoly_grn <- 0
pheno <- cbind(pheno, ctrl_metric)

tmpsamps <- which(pass_the_test$Staining.Green == 1)
pheno[tmpsamps, "rmv_sample"] <- 1
pheno[tmpsamps, "rmv_sample_ctrls_stain_grn"] <- 1

tmpsamps <- which(pass_the_test$Staining.Red == 1)
pheno[tmpsamps, "rmv_sample"] <- 1
pheno[tmpsamps, "rmv_sample_ctrls_stain_red"] <- 1

tmpsamps <- which(ctrl_metric$Specificity.I.Red < 0.9) # cutoff should be 1 but a few were close enough
pheno[tmpsamps, "rmv_sample"] <- 1
pheno[tmpsamps, "rmv_sample_ctrls_spec_I_red"] <- 1

tmpsamps <- which(pass_the_test$Non.polymorphic.Green == 1)
pheno[tmpsamps, "rmv_sample"] <- 1
pheno[tmpsamps, "rmv_sample_ctrls_nonpoly_grn"] <- 1



#4. Outliers with the deconvolution
cell.types <- c("Bas","Bmem","Bnv","CD4mem","CD4nv","CD8mem","CD8nv","Eos","Mono","Neu","NK","Treg")
cell_outliers <- data.frame(matrix( 0, nrow = nrow(pheno), ncol = length(cell.types)))
colnames(cell_outliers) <- paste0(cell.types, "_outlier")
for (cell in cell.types){
  IQR_median <- IQR(pheno[, cell])
  median_beta_val <- median(pheno[, cell])
  
  upper_limit <- min(median_beta_val + 3.5*IQR_median, 100)
  lower_limit <- max(median_beta_val - 3.5*IQR_median, 0)
  message("For ", cell, " the upper limit is ", upper_limit, " and the lower limit is ", lower_limit)
  for (i in 1:nrow(cell_outliers)){
    if (pheno[i,cell] < lower_limit) {
      cell_outliers[i,cell] <- 1
    } else if (pheno[i,cell] > upper_limit) {
      cell_outliers[i,cell] <- 1
    }
  }
  
}
pheno <- cbind(pheno,cell_outliers )


save(pheno, file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated.RDA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated.RDA")

pheno <- pheno[, c(1:72)]

SCRN <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Screening_Demographics.csv", 
                      col_types = cols(PATNO = col_character()))
SCRN <- SCRN[, c("PATNO","CURRENT_APPRDX")]
colnames(SCRN)[2] <- "subcohort_num"
SCRN$subcohort_char <- ifelse(SCRN$subcohort_num == 1, "PD_OG", 
                              ifelse(SCRN$subcohort_num == 2, "HC_OG", 
                                     ifelse(SCRN$subcohort_num == 3, "SWEDD",
                                            ifelse(SCRN$subcohort_num == 4, "Prod_OG",
                                                   ifelse(SCRN$subcohort_num == 5, "PD_GC",
                                                          ifelse(SCRN$subcohort_num == 6, "HC_GC", "???"
                                                                           ))))))
SCRN$DX_old_char <- ifelse(SCRN$subcohort_num == 1, "PD", 
                              ifelse(SCRN$subcohort_num == 2, "HC", 
                                     ifelse(SCRN$subcohort_num == 3, "SWEDD",
                                            ifelse(SCRN$subcohort_num == 4, "Prod",
                                                   ifelse(SCRN$subcohort_num == 5, "PD",
                                                          ifelse(SCRN$subcohort_num == 6, "HC", "???"
                                                          ))))))

SCRN$DX_old_num <- ifelse(SCRN$DX_old_char == "PD", 2, 
                           ifelse(SCRN$DX_old_char == "Prod", 1, 
                                  ifelse(SCRN$DX_old_char == "HC", 0, "???"
                                                       )))

pheno <- left_join(pheno, SCRN)

pheno$DX_new_char <- ifelse(pheno$CONCOHORT == 1, "PD", 
                            ifelse(pheno$CONCOHORT == 2, "HC", 
                                   ifelse(pheno$CONCOHORT == 4, "Prod", "???"
                                   )))
pheno$DX_new_num <- ifelse(pheno$DX_new_char == "PD", 2, 
                          ifelse(pheno$DX_new_char == "Prod", 1, 
                                 ifelse(pheno$DX_new_char == "HC", 0, "???"
                                 )))

pheno$DX_steve_char <- ifelse(pheno$subcohort_char == "PD_OG", "PD", 
                                  ifelse(pheno$subcohort_char == "HC_OG", "HC", 
                                         ifelse(pheno$subcohort_char == "Prod_OG", "Prod",
                                                ifelse(pheno$subcohort_char == "PD_GC", "PD", "???"
                                                              ))))
pheno$DX_steve_num <- ifelse(pheno$DX_steve_char == "PD", 2, 
                           ifelse(pheno$DX_steve_char == "Prod", 1, 
                                  ifelse(pheno$DX_steve_char == "HC", 0, "???"
                                  )))

pheno$DX_new_num_factor <- as.factor(pheno$DX_new_num)
pheno$DX_old_num_factor <- as.factor(pheno$DX_old_num)
pheno$DX_steve_num_factor <- factor(pheno$DX_steve_num, levels = c("0","1","2","???"))

pheno$mAccel_Hovath <- pheno$mAge_Hovath - pheno$Current_age
pheno$mAccel_Hannum <- pheno$mAge_Hannum - pheno$Current_age
pheno$mAccel_PhenoAge <- pheno$PhenoAge - pheno$Current_age

rownames(pheno) <- paste0(pheno$Slide,"_", pheno$Array)
pheno$ID <- rownames(pheno)

tmp <- pheno[pheno$rmv_sample == 1,]
pheno$run_date2 <- as.character(pheno$run_date)

pheno$run_date_group <- ifelse(startsWith(pheno$run_date2, "2019-06"), "June_2019", 
                               ifelse(startsWith(pheno$run_date2, "2019-07"), "July_2019", 
                                      ifelse(startsWith(pheno$run_date2, "2020-03"), "March_2020", 
                                             ifelse(startsWith(pheno$run_date2, "2020-04"), "March_2020", 
                                                    ifelse(startsWith(pheno$run_date2, "2020-06"), "June_2020",
                                                           ifelse(startsWith(pheno$run_date2, "2021-01"), "Jan_2021", 
                                                                  ifelse(startsWith(pheno$run_date2, "2021-02"), "Feb_2021",
                                                                         ifelse(startsWith(pheno$run_date2, "2021-03"), "March_2021", NA ))))))))

pheno$run_date_group <- factor(pheno$run_date_group, levels = c("June_2019","July_2019","March_2020" ,"June_2020" ,"Jan_2021","Feb_2021","March_2021"))

save(pheno, file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
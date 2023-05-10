load("/dartfs/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/DeconvoEWAS.RDATA")
cell_types <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem","Treg")
cell_types2 <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem")
dark_cols <- c(HC='#1DA072', Prod="#3971B5", PD="#C56124")
light_cols <- c(HC='#A9EFD7', Prod="#76B1E8", PD="#DBA010")
load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/betas_CHR.RDA")
load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/cpg_chr_names.RDA")

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
pheno <- pheno[pheno$rmv_sample == 0,]
setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/")
try(dir.create("EWAS"))
setwd("EWAS")
LoD<-c(Bas=0.478, Bmem=0.328, Bnv=0.483, CD4mem=1.051,
       CD4nv=0.878, CD8mem=0.801, CD8nv=0.758, Eos=0.422,
       Mono=0.394, Neu=0.338, NK=0.459, Treg=0.78)

for (i in 1:ncol(pheno)){
  if (colnames(pheno)[i] %in% cell_types){
    cell <- colnames(pheno)[i]
    pheno[pheno[,cell]<LoD[cell],cell]<- LoD[cell]/sqrt(2)
  }
}
tmp <- pheno[, which(colnames(pheno) %in% cell_types)]
pheno <- pheno[,-which(colnames(pheno) %in% cell_types)]
tmp <- round(tmp, 2)
tmp <- tmp*100/rowSums(tmp)
pheno <- cbind(pheno,tmp)
pheno$NLR <- pheno$Neu/(pheno$Treg + pheno$NK + pheno$Bnv + pheno$Bmem + pheno$CD8nv + pheno$CD8mem + pheno$CD4nv + pheno$CD4mem)
pheno$weird_batch <- ifelse(pheno$processing_batch %in% c(19,20), 1, 0)
pheno <- pheno[pheno$weird_batch == 0,]

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
##### Prod vs HC
setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS")
try(dir.create("ProdvsHC"))
setwd("ProdvsHC")
tmp_cols <- c(HC='#1DA072', PD="#C56124")

my_combos_to_test <- list(Unadj = NULL,
                          NLR = "NLR",
                          Many = c("Neu","Mono","Bnv","CD4mem","Eos"))

for (CHR in names(cpg_chr_names)){
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/")
  try(dir.create(CHR))
  setwd(CHR)
  chr_cpgs <- cpg_chr_names[[CHR]]
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
  
  ##### Cross Sectional
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/")
  setwd(CHR)
  try(dir.create("CrossSectional"))
  setwd("CrossSectional")
  for (timept in c("BL","V04","V06","V08")){
    pheno2 <- pheno[pheno$EVENT_ID == timept,]
    rownames(pheno2) <- pheno2$ID
    tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% pheno2$ID)]
    DeconvoEWAS(df=pheno2, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
      betas=tmp_betas, # beta matrix for dataset
      base_model="DX_new_char + BL_age + is.Female", # the right hand side of the regression model equation
      user_p_adjust_method = "fdr", # adjustment method for multiple testing
      cell_types=cell_types, # list of cell types 
      cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
      model_type = "lm",
      n_cpgs_to_test = cpg_per_chr[CHR], # how many CpGs to do EWAS with
      cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
      cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
      group_var="DX_new_char", # what is the grouping variable
      ref_group="HC", # what is the reference group
      case_group="Prod", # what is the group to compare 
      save_output_csv = TRUE, # save the results in csv table
      save_output_pdf = TRUE, # save a pdf version of the plot
      file_output_prefix = paste0(CHR,"_EWAS_",timept,"_subsetdata"), # string to attach to saved file names
      output_csv_dir = getwd(), # path to folder to save files in
      make_sum_equal_one = FALSE,
      scale_cell_types = TRUE, # scale the cell types
      confint.level=0.95)
  }
  
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/EWAS/ProdvsHC/")
  setwd(CHR)
  timept <- "LONGIT"
  try(dir.create("Longit"))
  setwd("Longit")
  tmp_betas <- agg_betas[, which(colnames(agg_betas) %in% pheno$ID)]
  pheno2 <- pheno
  rownames(pheno2) <- pheno2$ID
  DeconvoEWAS(df=pheno2, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
              betas=tmp_betas, # beta matrix for dataset
              base_model="DX_new_char + BL_age + is.Female + Sample_Time", # the right hand side of the regression model equation
              user_p_adjust_method = "fdr", # adjustment method for multiple testing
              cell_types=cell_types, # list of cell types 
              cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
              model_type = "lmer",
              lmer_randomeffect_var = "PATNO",
              n_cpgs_to_test = cpg_per_chr[CHR], # how many CpGs to do EWAS with
              cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
              cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
              group_var="DX_new_char", # what is the grouping variable
              ref_group="HC", # what is the reference group
              case_group="Prod", # what is the group to compare 
              save_output_csv = TRUE, # save the results in csv table
              save_output_pdf = TRUE, # save a pdf version of the plot
              file_output_prefix = paste0(CHR,"_EWAS_",timept,"_subsetdata"), # string to attach to saved file names
              output_csv_dir = getwd(), # path to folder to save files in
              make_sum_equal_one = FALSE,
              scale_cell_types = TRUE, # scale the cell types
              confint.level=0.95)
  timept <- "LONGIT_inter"
  DeconvoEWAS(df=pheno2, # data frame with all data, colnames should match variable names within reg.model, cell.types, rownames should match colnames of betas 
              betas=tmp_betas, # beta matrix for dataset
              base_model="DX_new_char*Sample_Time + BL_age + is.Female", # the right hand side of the regression model equation
              user_p_adjust_method = "fdr", # adjustment method for multiple testing
              cell_types=cell_types, # list of cell types 
              cell_combos_to_test = my_combos_to_test, # named list of vectors containing names of cell types to try, each entry of list will be one model, leave one entry as NULL to do unadj
              model_type = "lmer",
              lmer_randomeffect_var = "PATNO",
              n_cpgs_to_test = cpg_per_chr[CHR], # how many CpGs to do EWAS with
              cpg_filter_method = "deltaBeta", # method to filter top cpgs, either by total variance or by delta beta
              cpg_user_names=NULL, # if cpg filter method is cpg list, should be names of cpgs user wants to test
              group_var="DX_new_char", # what is the grouping variable
              ref_group="HC", # what is the reference group
              case_group="Prod", # what is the group to compare 
              save_output_csv = TRUE, # save the results in csv table
              save_output_pdf = TRUE, # save a pdf version of the plot
              file_output_prefix = paste0(CHR,"_EWAS_",timept,"_subsetdata"), # string to attach to saved file names
              output_csv_dir = getwd(), # path to folder to save files in
              make_sum_equal_one = FALSE,
              scale_cell_types = TRUE, # scale the cell types
              confint.level=0.95)
}

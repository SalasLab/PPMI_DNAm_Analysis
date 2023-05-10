load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/ggradar2_steve.RDATA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/DeconvoRadar.RDATA")
cell_types <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem","Treg")
cell_types2 <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem")
dark_cols <- c(HC='#1DA072', Prod="#3971B5", PD="#C56124")
light_cols <- c(HC='#A9EFD7', Prod="#76B1E8", PD="#DBA010")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
pheno <- pheno[pheno$rmv_sample == 0,]
pheno <- pheno[pheno$DX_new_char == "PD",]
pheno$rundate_binary <- ifelse(pheno$run_date_group %in% c("Jan_2021", "Feb_2021", "March_2021"),1,0)
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

library(readr)
library(readxl)
DaT <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/DaTScan_Analysis (1).csv", 
                col_types = cols(PATNO = col_character()))
DaT$EVENT_ID <- ifelse(DaT$EVENT_ID == "SC", "BL", DaT$EVENT_ID)
DaT <- DaT[, c(2,3,6,7,8,9)]
DaT$visit_ID <- paste0(DaT$PATNO,"_",DaT$EVENT_ID)
DaT$DATSCAN_PUTAMEN_AI <- 100*abs((DaT$DATSCAN_PUTAMEN_R - DaT$DATSCAN_PUTAMEN_L)/(DaT$DATSCAN_PUTAMEN_R + DaT$DATSCAN_PUTAMEN_L))
DaT$DATSCAN_CAUDATE_AI <- 100*abs((DaT$DATSCAN_CAUDATE_R - DaT$DATSCAN_CAUDATE_L)/(DaT$DATSCAN_CAUDATE_R + DaT$DATSCAN_CAUDATE_L))

DaT_YN <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/DaTScan_Visual_Interpretation_Results.csv", 
                   col_types = cols(PATNO = col_character()))
DaT_YN <- DaT_YN[, c(2,4)]
library(dplyr)
pheno <- pheno %>%
  left_join(DaT) %>%
  left_join(DaT_YN)

PT_status <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Participant_Status.csv", 
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

clin_vars <- c("DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L",
               "DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L",
               "DATSCAN_PUTAMEN_AI","DATSCAN_CAUDATE_AI")
for (clin in clin_vars){
  print(clin)
  print(table(pheno$Sample_Time, is.na(pheno[,clin])))
}


setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
try(dir.create("DaTscan"))
setwd("DaTscan")

##### Cross Sectional
try(dir.create("CrossSectional"))
for (clin in clin_vars){
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
  setwd("DaTscan")
  setwd("CrossSectional")
  for (timept in c("BL","V04","V06")){
    pheno2 <- pheno[pheno$EVENT_ID == timept,]
    pheno2 <- pheno2[!is.na(pheno2[,clin]),]
    pheno3 <- pheno2[pheno2$weird_batch == 0,]
    ## Model1: Univar
    clin_reg_model <- clin
    DeconvoRadarCont(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
      reg_model=clin_reg_model, # the right hand side of the regression model equation
      model_type="lm", # the type of model you want (options: lm, lmer)
      cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
      cont_var = clin,
      group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
      save_output_csv = TRUE, # save the results in csv table
      save_output_pdf = TRUE, # save a pdf version of the plot
      file_output_prefix = paste0("Model1_",clin, " ",timept), # string to attach to saved file names
      output_csv_dir = getwd(), # path to folder to save files in
      grid_label_size_user = 3.5, # size of the ring value labels
      axis_label_size_user = 4,
      plot_title = paste0(clin, " ",timept),
      legend_text_size_user = 10,
      make_sum_equal_one = FALSE, # adjust proportions so they add to 1
      confint.level=0.95)
    
    ## Model3: Univar (subset out bad batch)
    clin_reg_model <- paste0(clin, " rundate_binary")
    DeconvoRadarCont(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                     reg_model=clin_reg_model, # the right hand side of the regression model equation
                     model_type="lm", # the type of model you want (options: lm, lmer)
                     cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                     cont_var = clin,
                     group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                     save_output_csv = TRUE, # save the results in csv table
                     save_output_pdf = TRUE, # save a pdf version of the plot
                     file_output_prefix = paste0("Model3_",clin, " ",timept), # string to attach to saved file names
                     output_csv_dir = getwd(), # path to folder to save files in
                     grid_label_size_user = 3.5, # size of the ring value labels
                     axis_label_size_user = 4,
                     plot_title = paste0(clin, " ",timept),
                     legend_text_size_user = 10,
                     make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                     confint.level=0.95)
    
    ## Model5: Multivar + batch (subset out bad batch)
    clin_reg_model <- paste0(clin, " + rundate_binary + BL_age + is.Female")
    DeconvoRadarCont(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                     reg_model=clin_reg_model, # the right hand side of the regression model equation
                     model_type="lm", # the type of model you want (options: lm, lmer)
                     cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                     cont_var = clin,
                     group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                     save_output_csv = TRUE, # save the results in csv table
                     save_output_pdf = TRUE, # save a pdf version of the plot
                     file_output_prefix = paste0("Model5_",clin, " ",timept), # string to attach to saved file names
                     output_csv_dir = getwd(), # path to folder to save files in
                     grid_label_size_user = 3.5, # size of the ring value labels
                     axis_label_size_user = 4,
                     plot_title = paste0(clin, " ",timept),
                     legend_text_size_user = 10,
                     make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                     confint.level=0.95)
  
  }

  
  ##### Longitudinal
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/DaTscan")
  try(dir.create("Longit"))
  setwd("Longit")
  timept <- "LONGIT"
  
  pheno2 <- pheno[!is.na(pheno[,clin]),]
  pheno3 <- pheno2[pheno2$weird_batch == 0,]
  ## Model1: Univar
  clin_reg_model <- paste0(clin, " + Sample_Time + (1|PATNO)")
  DeconvoRadarCont(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                   reg_model=clin_reg_model, # the right hand side of the regression model equation
                   model_type="lmer", # the type of model you want (options: lm, lmer)
                   cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                   cont_var = clin,
                   group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                   save_output_csv = TRUE, # save the results in csv table
                   save_output_pdf = TRUE, # save a pdf version of the plot
                   file_output_prefix = paste0("Model1_",clin, " ",timept), # string to attach to saved file names
                   output_csv_dir = getwd(), # path to folder to save files in
                   grid_label_size_user = 3.5, # size of the ring value labels
                   axis_label_size_user = 4,
                   plot_title = paste0(clin, " ",timept),
                   legend_text_size_user = 10,
                   make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                   confint.level=0.95)
  
  ## Model2: Univar + batch
  clin_reg_model <- paste0(clin, " + rundate_binary + weird_batch + Sample_Time + (1|PATNO)")
  DeconvoRadarCont(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                   reg_model=clin_reg_model, # the right hand side of the regression model equation
                   model_type="lmer", # the type of model you want (options: lm, lmer)
                   cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                   cont_var = clin,
                   group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                   save_output_csv = TRUE, # save the results in csv table
                   save_output_pdf = TRUE, # save a pdf version of the plot
                   file_output_prefix = paste0("Model2_",clin, " ",timept), # string to attach to saved file names
                   output_csv_dir = getwd(), # path to folder to save files in
                   grid_label_size_user = 3.5, # size of the ring value labels
                   axis_label_size_user = 4,
                   plot_title = paste0(clin, " ",timept),
                   legend_text_size_user = 10,
                   make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                   confint.level=0.95)
  
  ## Model3: Univar (subset out bad batch)
  clin_reg_model <- paste0(clin, " + rundate_binary + Sample_Time + (1|PATNO)")
  DeconvoRadarCont(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                   reg_model=clin_reg_model, # the right hand side of the regression model equation
                   model_type="lmer", # the type of model you want (options: lm, lmer)
                   cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                   cont_var = clin,
                   group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                   save_output_csv = TRUE, # save the results in csv table
                   save_output_pdf = TRUE, # save a pdf version of the plot
                   file_output_prefix = paste0("Model3_",clin, " ",timept), # string to attach to saved file names
                   output_csv_dir = getwd(), # path to folder to save files in
                   grid_label_size_user = 3.5, # size of the ring value labels
                   axis_label_size_user = 4,
                   plot_title = paste0(clin, " ",timept),
                   legend_text_size_user = 10,
                   make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                   confint.level=0.95)
  
  ## Model4: Multivar + batch
  clin_reg_model <- paste0(clin, " + rundate_binary + weird_batch + BL_age + is.Female + Sample_Time + (1|PATNO)")
  DeconvoRadarCont(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                   reg_model=clin_reg_model, # the right hand side of the regression model equation
                   model_type="lmer", # the type of model you want (options: lm, lmer)
                   cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                   cont_var = clin,
                   group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                   save_output_csv = TRUE, # save the results in csv table
                   save_output_pdf = TRUE, # save a pdf version of the plot
                   file_output_prefix = paste0("Model4_",clin, " ",timept), # string to attach to saved file names
                   output_csv_dir = getwd(), # path to folder to save files in
                   grid_label_size_user = 3.5, # size of the ring value labels
                   axis_label_size_user = 4,
                   plot_title = paste0(clin, " ",timept),
                   legend_text_size_user = 10,
                   make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                   confint.level=0.95)
  
  ## Model5: Multivar + batch (subset out bad batch)
  clin_reg_model <- paste0(clin, " + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)")
  DeconvoRadarCont(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                   reg_model=clin_reg_model, # the right hand side of the regression model equation
                   model_type="lmer", # the type of model you want (options: lm, lmer)
                   cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                   cont_var = clin,
                   group_colours_user = c("black","red"), # colors for plot, name the vector the names of the groups
                   save_output_csv = TRUE, # save the results in csv table
                   save_output_pdf = TRUE, # save a pdf version of the plot
                   file_output_prefix = paste0("Model5_",clin, " ",timept), # string to attach to saved file names
                   output_csv_dir = getwd(), # path to folder to save files in
                   grid_label_size_user = 3.5, # size of the ring value labels
                   axis_label_size_user = 4,
                   plot_title = paste0(clin, " ",timept),
                   legend_text_size_user = 10,
                   make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                   confint.level=0.95)
  
  
}

library(readr)
library(dplyr)
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/ggradar2_steve.RDATA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Scripts/Test_Deconvo_Func_Output/DeconvoRadar.RDATA")
cell_types <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem","Treg")
cell_types2 <- c("Bas","Eos","Neu","Mono","Bmem","Bnv","NK","CD4nv","CD8nv","CD8mem","CD4mem")
dark_cols <- c(HC='#1DA072', Prod="#3971B5", PD="#C56124")
light_cols <- c(HC='#A9EFD7', Prod="#76B1E8", PD="#DBA010")

load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
pheno <- pheno[pheno$rmv_sample == 0,]
pheno$rundate_binary <- ifelse(pheno$run_date_group %in% c("Jan_2021", "Feb_2021", "March_2021"),1,0)
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
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
pheno$Prod_conv_rmv <- ifelse(is.na(pheno$Prod_Convert_Relavtive), 0, 
                              ifelse(pheno$Prod_Convert_Relavtive == "Converted", 1, 0))

##### PD vs HC
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
try(dir.create("PDvsHC"))
setwd("PDvsHC")
tmp_cols <- c(HC='#1DA072', PD="#C56124")

##### Cross Sectional
try(dir.create("CrossSectional"))
setwd("CrossSectional")

for (timept in c("BL","V04","V06","V08")){
  pheno2 <- pheno[pheno$EVENT_ID == timept,]
  pheno3 <- pheno2[pheno2$weird_batch == 0,]
  ## Model1: Univar
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
    reg_model="DX_new_char", # the right hand side of the regression model equation
    model_type="lm", # the type of model you want (options: lm, lmer)
    cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
    group_var="DX_new_char", # what is the grouping variable
    ref_group="HC", # what is the reference group
    case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
    group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
    save_output_csv = TRUE, # save the results in csv table
    save_output_pdf = TRUE, # save a pdf version of the plot
    file_output_prefix = paste0("Model1_",timept), # string to attach to saved file names
    output_csv_dir = getwd(), # path to folder to save files in
    grid_label_size_user = 6.2, # size of the ring value labels
    axis_label_size_user = 5.5,
    plot_title = paste0("PDvsHC: ",timept),
    legend_text_size_user = 2,
    make_sum_equal_one = FALSE, # adjust proportions so they add to 1
    confint.level=0.95)
  
  ## Model2: Univar + batch
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + weird_batch + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model2_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("PDvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model3: Univar (subset out bad batch)
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model3_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("PDvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model4: Multivar + batch
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + weird_batch + rundate_binary + BL_age + is.Female", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model4_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("PDvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model5: Multivar + batch (subset out bad batch)
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + BL_age + is.Female + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model5_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("PDvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)

}

##### Longitudinal
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
setwd("PDvsHC")
try(dir.create("Longit"))
setwd("Longit")
timept <- "LONGIT"

pheno3 <- pheno[pheno$weird_batch == 0,]
## Model1: Univar
DeconvoRadar(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model1_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model2: Univar + batch
DeconvoRadar(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + weird_batch + rundate_binary + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model2_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model3: Univar (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + rundate_binary + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model3_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model4: Multivar + batch
DeconvoRadar(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + weird_batch + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model4_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model5: Multivar + batch (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model5_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model6: Multivar + batch + interaction
DeconvoRadar(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char*Sample_Time + rundate_binary + weird_batch + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model6_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model7: Multivar + interaction (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char*Sample_Time + rundate_binary + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model7_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)


load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CMV_status.RDA")
pheno <- left_join(pheno, CMV_status3)

## Model8: try adding CMV status
DeconvoRadar(df=pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="CMV_Status + DX_new_char + weird_batch + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="PD", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model8_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("PDvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)











##### Prod vs HC
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
try(dir.create("ProdvsHC"))
setwd("ProdvsHC")
tmp_cols <- c(HC='#1DA072', Prod="#3971B5")

##### Cross Sectional
try(dir.create("CrossSectional"))
setwd("CrossSectional")

for (timept in c("BL","V04","V06","V08")){
  pheno2 <- pheno[pheno$EVENT_ID == timept,]
  pheno2 <- pheno2[pheno2$Prod_conv_rmv == 0,]
  pheno3 <- pheno2[pheno2$weird_batch == 0,]
  ## Model1: Univar
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model1_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("ProdvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model2: Univar + batch
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + weird_batch + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model2_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("ProdvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model3: Univar (subset out bad batch)
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model3_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("ProdvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model4: Multivar + batch
  DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + weird_batch + BL_age + is.Female + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model4_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("ProdvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  ## Model5: Multivar + batch (subset out bad batch)
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model="DX_new_char + BL_age + is.Female + + rundate_binary", # the right hand side of the regression model equation
               model_type="lm", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="DX_new_char", # what is the grouping variable
               ref_group="HC", # what is the reference group
               case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("Model5_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("ProdvsHC: ",timept),
               legend_text_size_user = 2,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
}

##### Longitudinal
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
setwd("ProdvsHC")
try(dir.create("Longit"))
setwd("Longit")
timept <- "LONGIT"

pheno2 <- pheno[pheno$Prod_conv_rmv == 0,]
pheno3 <- pheno2[pheno2$weird_batch == 0,]
## Model1: Univar
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model1_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model2: Univar + batch
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + weird_batch + rundate_binary + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model2_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model3: Univar (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + Sample_Time + rundate_binary + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model3_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model4: Multivar + batch
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + weird_batch + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model4_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model5: Multivar + batch (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model5_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model6: Multivar + batch + interaction
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char*Sample_Time + rundate_binary + weird_batch + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model6_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model7: Multivar + interaction (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="DX_new_char*Sample_Time + rundate_binary + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="DX_new_char", # what is the grouping variable
             ref_group="HC", # what is the reference group
             case_groups="Prod", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model7_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdvsHC: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)











##### Prod convereted vs Prod nonconverted
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
try(dir.create("ProdConvsProdNoncon"))
setwd("ProdConvsProdNoncon")
tmp_cols <- c(Converted='red', `Pre-conversion`="black")

##### Longitudinal
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp")
setwd("ProdConvsProdNoncon")
try(dir.create("Longit"))
setwd("Longit")
timept <- "LONGIT"

pheno2 <- pheno[!is.na(pheno$Prod_Convert_Relavtive),]
pheno3 <- pheno2[pheno2$weird_batch == 0,]
## Model1: Univar
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model1_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model2: Univar + batch
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive + rundate_binary + weird_batch + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model2_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model3: Univar (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive + rundate_binary + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model3_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model4: Multivar + batch
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive + rundate_binary + weird_batch + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model4_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model5: Multivar + batch (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model5_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model6: Multivar + batch + interaction
DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive*Sample_Time + rundate_binary + weird_batch + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model6_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)

## Model7: Multivar + interaction (subset out bad batch)
DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="Prod_Convert_Relavtive*Sample_Time + rundate_binary + BL_age + is.Female + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="Prod_Convert_Relavtive", # what is the grouping variable
             ref_group="Pre-conversion", # what is the reference group
             case_groups="Converted", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = paste0("Model7_",timept), # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = paste0("ProdConvsProdNoncon: ",timept),
             legend_text_size_user = 2,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)
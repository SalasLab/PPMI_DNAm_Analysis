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
MDS1 <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/MDS-UPDRS_Part_I_Patient_Questionnaire.csv", col_types = cols(PATNO = col_character()))
MDS1$visit_ID <- paste0(MDS1$PATNO,"_",MDS1$EVENT_ID)
MDS1 <- MDS1[, c("visit_ID", "NP1PTOT")]
MDS1 <- MDS1[!is.na(MDS1$NP1PTOT),]

MDS2 <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/MDS_UPDRS_Part_II__Patient_Questionnaire.csv",col_types = cols(PATNO = col_character()))
MDS2$visit_ID <- paste0(MDS2$PATNO,"_",MDS2$EVENT_ID)
MDS2 <- MDS2[, c("visit_ID", "NP2PTOT")]
MDS2 <- MDS2[!is.na(MDS2$NP2PTOT),]

ESS <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Epworth_Sleepiness_Scale.csv", col_types = cols(PATNO = col_character()))
ESS$visit_ID <- paste0(ESS$PATNO,"_",ESS$EVENT_ID)
ESS$ESS_tot <- apply(ESS[, 7:14], 1, sum)
ESS <- ESS[, c("visit_ID", "ESS_tot")]
ESS <- ESS[!is.na(ESS$ESS_tot),]

MCA <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Montreal_Cognitive_Assessment__MoCA_.csv", col_types = cols(PATNO = col_character()))
MCA <- MCA[MCA$EVENT_ID != "BL",]
MCA$EVENT_ID <- ifelse(MCA$EVENT_ID == "SC", "BL", MCA$EVENT_ID)
MCA$visit_ID <- paste0(MCA$PATNO,"_",MCA$EVENT_ID)
MCA <- MCA[, c("visit_ID", "MCATOT")]
MCA <- MCA[!is.na(MCA$MCATOT),]

GDS <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Geriatric_Depression_Scale__Short_Version_.csv", col_types = cols(PATNO = col_character()))
# tmp <- GDS[,c(6:20)]
# switch_ans <- c(1,5,7,11)
# for (i in 1:nrow(tmp)){
#   for (j in switch_ans){
#     if (is.na(tmp[i,j])){next} else{
#     if (tmp[i,j] == 0){tmp[i,j] <- 1} else{
#     if (tmp[i,j] == 1){tmp[i,j] <- 0} else{
#       tmp[i,j] <- NA
#     }}}
#   }
# }
# GDS_calculated <- tmp
# GDS_calculated$GDS_tot <- rowSums(tmp)
# save(GDS_calculated, file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Calculated_GDS.RDA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Calculated_GDS.RDA")
GDS$GDS_tot <- GDS_calculated$GDS_tot
GDS$visit_ID <- paste0(GDS$PATNO,"_",GDS$EVENT_ID)
GDS <- GDS[, c("visit_ID", "GDS_tot")]
GDS <- GDS[!is.na(GDS$GDS_tot),]


STAI <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/State-Trait_Anxiety_Inventory.csv",  col_types = cols(PATNO = col_character()))
# tmp <- STAI[,c(6:45)]
# switch_ans <- c(1,2,5,8,10,11,15,16,19,20,21,23,26,27,30,33,34,36,39)
# for (i in 1:nrow(tmp)){
#   for (j in switch_ans){
#     if (is.na(tmp[i,j])){next} else{
#     if (tmp[i,j] == 4){tmp[i,j] <- 1} else{
#     if (tmp[i,j] == 3){tmp[i,j] <- 2} else{
#     if (tmp[i,j] == 2){tmp[i,j] <- 3} else{
#     if (tmp[i,j] == 1){tmp[i,j] <- 4} else{
#       tmp[i,j] <- NA
#     }}}}}
#   }
# }
# 
# STAI_calculated <- tmp
# STAI_calculated$STAI_tot <- rowSums(tmp)
#save(STAI_calculated, file = "//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Calculated_STAI.RDA")
load("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/Calculated_STAI.RDA")
STAI$STAI_tot <- STAI_calculated$STAI_tot
STAI$visit_ID <- paste0(STAI$PATNO,"_",STAI$EVENT_ID)
STAI <- STAI[, c("visit_ID", "STAI_tot")]
STAI <- STAI[!is.na(STAI$STAI_tot),]


SCOPA <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/SCOPA-AUT.csv",  col_types = cols(PATNO = col_character()))
SCOPA$SCOPA_tot <- rowSums(SCOPA[,c(7:27)])
SCOPA$visit_ID <- paste0(SCOPA$PATNO,"_",SCOPA$EVENT_ID)
SCOPA <- SCOPA[, c("visit_ID", "SCOPA_tot")]
SCOPA <- SCOPA[!is.na(SCOPA$SCOPA_tot),]

REM <- read_csv("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/REM_Sleep_Behavior_Disorder_Questionnaire.csv",  col_types = cols(PATNO = col_character()))
REM$REM_tot <- rowSums(REM[,c(7:27)])
REM$visit_ID <- paste0(REM$PATNO,"_",REM$EVENT_ID)
REM <- REM[, c("visit_ID", "REM_tot")]
REM <- REM[!is.na(REM$REM_tot),]

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
  left_join(MDS1) %>%
  left_join(MDS2) %>%
  left_join(ESS) %>%
  left_join(MCA) %>%
  left_join(STAI) %>%
  left_join(SCOPA) %>%
  left_join(REM) %>%
  left_join(GDS) %>%
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

clin_vars <- c("NP1PTOT","NP2PTOT","ESS_tot","MCATOT","STAI_tot","REM_tot","GDS_tot",
               "DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L",
               "DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L",
               "DATSCAN_PUTAMEN_AI","DATSCAN_CAUDATE_AI")
for (clin in clin_vars){
  print(clin)
  print(table(pheno$Sample_Time, is.na(pheno[,clin])))
}



################# Clustering
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(factoextra)
library(cluster)
BL_pheno <- pheno[pheno$Sample_Time == 0,]
clin_pheno <- BL_pheno
clin_pheno <- clin_pheno[, clin_vars]
rownames(clin_pheno) <- BL_pheno$ID
clin_pheno <- clin_pheno[, -c(8:11)]
clin_pheno <- na.omit(clin_pheno)
clin_pheno2 <- clin_pheno
clin_pheno <- scale(clin_pheno)
fviz_nbclust(clin_pheno, kmeans, method = "wss")

gap_stat <- clusGap(clin_pheno,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)

km <- kmeans(clin_pheno, centers = 2, nstart = 25)

fviz_cluster(km,data=clin_pheno)

table(km$cluster)
aggregate(clin_pheno2, by=list(cluster=km$cluster), mean)
tmp <- aggregate(clin_pheno, by=list(cluster=km$cluster), mean)
dev.off()
pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/ClinVar/BL_clustering.PDF", height = 4, width = 6)
print(pheatmap(t(tmp[,-1]), cluster_cols = FALSE, color = cividis(20)))
dev.off()
cluster_data <- km$cluster
BL_pheno2 <- BL_pheno[which(BL_pheno$ID %in% names(cluster_data)),]
BL_pheno2 <- BL_pheno2[order(BL_pheno2$ID),]
cluster_data <- cluster_data[order(names(cluster_data))]
BL_pheno2 <- cbind(BL_pheno2, cluster =cluster_data)
BL_pheno2$cluster <- ifelse(BL_pheno2$cluster == 1, "C1","C2")
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
try(dir.create("ClinVar"))
setwd("ClinVar")
tmp_cols <- c(C1 = "black", C2="red")

DeconvoRadar(df=BL_pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="cluster + BL_age + is.Female + rundate_binary + weird_batch", # the right hand side of the regression model equation
             model_type="lm", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="cluster", # what is the grouping variable
             ref_group="C1", # what is the reference group
             case_groups="C2", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = "KmeansClinvar_BL", # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = "KmeansClinvar_BL",
             legend_text_size_user = 1,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)
dev.off()




# taking Year 1 estunate and going longitdunally, test interaction
library(lmtest)
library(lmerTest)
tmp_pheno <- pheno[pheno$PATNO %in% BL_pheno2$PATNO,]
tmp_data <- BL_pheno2[, c("PATNO","cluster")]
tmp_pheno <- left_join(tmp_pheno, tmp_data, by = "PATNO")


for (cell in cell_types2){
  tmp_cell_data <- tmp_pheno[, cell]
  tmp_model1 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + cluster + Sample_Time + (1|PATNO),
                     data = tmp_pheno)
  tmp_model2 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + cluster*Sample_Time + (1|PATNO),
                     data = tmp_pheno)
  
  tmp_results <- lrtest(tmp_model2, tmp_model1)
  if(tmp_results[2,5] < 0.05){
    message("Sig Interaction: ", cell, " pval: ", tmp_results[2,5])
  }
}
tmp_cols <- c(C1 = "black", C2="red")

DeconvoRadar(df=tmp_pheno, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="cluster + BL_age + is.Female + rundate_binary + weird_batch + Sample_Time + (1|PATNO)", # the right hand side of the regression model equation
             model_type="lmer", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="cluster", # what is the grouping variable
             ref_group="C1", # what is the reference group
             case_groups="C2", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = "KmeansClinvar_LONGIT", # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = "KmeansClinvar_LONGIT",
             legend_text_size_user = 1,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)








######## Year 2 timepoint
BL_pheno <- pheno[pheno$Sample_Time == 24,]
clin_pheno <- BL_pheno
clin_pheno <- clin_pheno[, clin_vars]
rownames(clin_pheno) <- BL_pheno$ID
clin_pheno <- clin_pheno[, -c(8:11)]
clin_pheno <- na.omit(clin_pheno)
clin_pheno2 <- clin_pheno
clin_pheno <- scale(clin_pheno)
fviz_nbclust(clin_pheno, kmeans, method = "wss")

gap_stat <- clusGap(clin_pheno,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 50)
fviz_gap_stat(gap_stat)

km <- kmeans(clin_pheno, centers = 2, nstart = 25)

fviz_cluster(km,data=clin_pheno)
table(km$cluster)

aggregate(clin_pheno2, by=list(cluster=km$cluster), mean)
tmp <- aggregate(clin_pheno, by=list(cluster=km$cluster), mean)
dev.off()
pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/ClinVar/V06_clustering.PDF", height = 4, width = 6)
print(pheatmap(t(tmp[,-1]), cluster_cols = FALSE, color = cividis(20)))
dev.off()
cluster_data <- km$cluster
BL_pheno2 <- BL_pheno[which(BL_pheno$ID %in% names(cluster_data)),]
BL_pheno2 <- BL_pheno2[order(BL_pheno2$ID),]
cluster_data <- cluster_data[order(names(cluster_data))]
BL_pheno2 <- cbind(BL_pheno2, cluster =cluster_data)
BL_pheno2$cluster <- ifelse(BL_pheno2$cluster == 1, "C1","C2")
setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
try(dir.create("ClinVar"))
setwd("ClinVar")
tmp_cols <- c(C1 = "black", C2="red")

DeconvoRadar(df=BL_pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
             reg_model="cluster + BL_age + is.Female + rundate_binary + weird_batch", # the right hand side of the regression model equation
             model_type="lm", # the type of model you want (options: lm, lmer)
             cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
             group_var="cluster", # what is the grouping variable
             ref_group="C1", # what is the reference group
             case_groups="C2", # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
             group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
             save_output_csv = TRUE, # save the results in csv table
             save_output_pdf = TRUE, # save a pdf version of the plot
             file_output_prefix = "KmeansClinvar_V06", # string to attach to saved file names
             output_csv_dir = getwd(), # path to folder to save files in
             grid_label_size_user = 6.2, # size of the ring value labels
             axis_label_size_user = 5.5,
             plot_title = "KmeansClinvar_V06",
             legend_text_size_user = 1,
             make_sum_equal_one = FALSE, # adjust proportions so they add to 1
             confint.level=0.95)
dev.off()

setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
try(dir.create("ClinVar"))
setwd("ClinVar")

##### Cross Sectional
try(dir.create("CrossSectional"))
for (clin in clin_vars){
  for (timept in c("BL","V04","V06","V08")){
    if (clin %in% c("DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L",
                    "DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L",
                    "DATSCAN_PUTAMEN_AI","DATSCAN_CAUDATE_AI") & timept == "V08"){next}
    setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
    setwd("ClinVar")
    setwd("CrossSectional")
    pheno2 <- pheno[pheno$EVENT_ID == timept,]
    pheno2 <- pheno2[!is.na(pheno2[, clin]),]
    pheno2$clin_quant <- factor(as.character(ntile(pheno2[,clin], n=4)))
    pheno3 <- pheno2[pheno2$weird_batch == 0,]
    
    if (clin %in% c("MCATOT","DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L","DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L")){
        ref_group_clin <- "4"
        case_group_clin <- "1"
        tmp_cols <- c(`4` = "black", `1`="red")
    } else {
        ref_group_clin <- "1"
        case_group_clin <- "4"
        tmp_cols <- c(`1` = "black", `4`="red")
    }
    
    ## Model1: Univar
    clin_reg_model <- "clin_quant"
    DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                     reg_model=clin_reg_model, # the right hand side of the regression model equation
                     model_type="lm", # the type of model you want (options: lm, lmer)
                     cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                     group_var="clin_quant", # what is the grouping variable
                     ref_group=ref_group_clin, # what is the reference group
                     case_groups=case_group_clin, # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
                     group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
                     save_output_csv = TRUE, # save the results in csv table
                     save_output_pdf = TRUE, # save a pdf version of the plot
                     file_output_prefix = paste0("QuantileComp_Model1_", clin,"_",timept), # string to attach to saved file names
                     output_csv_dir = getwd(), # path to folder to save files in
                     grid_label_size_user = 6.2, # size of the ring value labels
                     axis_label_size_user = 5.5,
                     plot_title = paste0("Model1_",clin, " ",timept),
                     legend_text_size_user = 1,
                     make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                     confint.level=0.95)
    
   
    
    ## Model4: Multivar + batch
    clin_reg_model <- "clin_quant + rundate_binary + BL_age + is.Female"
    DeconvoRadar(df=pheno2, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
                     reg_model=clin_reg_model, # the right hand side of the regression model equation
                     model_type="lm", # the type of model you want (options: lm, lmer)
                     cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
                     group_var="clin_quant", # what is the grouping variable
                     ref_group=ref_group_clin, # what is the reference group
                     case_groups=case_group_clin, # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
                     group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
                     save_output_csv = TRUE, # save the results in csv table
                     save_output_pdf = TRUE, # save a pdf version of the plot
                     file_output_prefix = paste0("QuantileComp_Model4_", clin,"_",timept), # string to attach to saved file names
                     output_csv_dir = getwd(), # path to folder to save files in
                     grid_label_size_user = 6.2, # size of the ring value labels
                     axis_label_size_user = 5.5,
                     plot_title = paste0("Model4_",clin, " ",timept),
                     legend_text_size_user = 1,
                     make_sum_equal_one = FALSE, # adjust proportions so they add to 1
                     confint.level=0.95)
    
  }
}

  
  ##### Longitudinal
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/ClinVar")
  try(dir.create("Longit"))
  setwd("Longit")
  timept <- "LONGIT"

for (clin in clin_vars){
  setwd("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023")
  setwd("ClinVar")
  setwd("Longit")
  pheno2 <- pheno[!is.na(pheno[, clin]),]
  PATNO_BL <- pheno2[pheno2$Sample_Time == 0,]
  if (clin %in% c("DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L",
                  "DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L",
                  "DATSCAN_PUTAMEN_AI","DATSCAN_CAUDATE_AI")){
  PATNO_END <- pheno2[pheno2$Sample_Time == 24,]
  } else {
  PATNO_END <- pheno2[pheno2$Sample_Time == 36,]
  }
  PATNOs <- intersect(PATNO_BL$PATNO, PATNO_END$PATNO)
  PATNO_BL <- PATNO_BL[PATNO_BL$PATNO %in% PATNOs,]
  PATNO_END <- PATNO_END[PATNO_END$PATNO %in% PATNOs,]
  PATNO_BL <- PATNO_BL[order(PATNO_BL$PATNO),]
  PATNO_END <- PATNO_END[order(PATNO_END$PATNO),]
  print(identical(PATNO_BL$PATNO,PATNO_END$PATNO))
  PATNO_BL$time_diff <- PATNO_END[, clin] - PATNO_BL[, clin]
  PATNO_BL$clin_quant <- factor(as.character(ntile(PATNO_BL[,clin], n=4)))
  PATNO_BL <- PATNO_BL[, c("PATNO","clin_quant")]
  pheno3 <- pheno2[pheno2$PATNO %in% PATNOs,]
  pheno3 <- left_join(pheno3, PATNO_BL, by = "PATNO")

  if (clin %in% c("MCATOT","DATSCAN_PUTAMEN_R","DATSCAN_PUTAMEN_L","DATSCAN_CAUDATE_R","DATSCAN_CAUDATE_L")){
    ref_group_clin <- "4"
    case_group_clin <- "1"
    tmp_cols <- c(`4` = "black", `1`="red")
  } else {
    ref_group_clin <- "1"
    case_group_clin <- "4"
    tmp_cols <- c(`1` = "black", `4`="red")
  }
  
  ## Model1: Univar
  clin_reg_model <- "clin_quant + Sample_Time + (1|PATNO)"
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model=clin_reg_model, # the right hand side of the regression model equation
               model_type="lmer", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="clin_quant", # what is the grouping variable
               ref_group=ref_group_clin, # what is the reference group
               case_groups=case_group_clin, # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("QuantileComp_Model1_", clin,"_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("Model1_",clin, " ",timept),
               legend_text_size_user = 1,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
  
  
  ## Model4: Multivar + batch
  clin_reg_model <- "clin_quant + rundate_binary + BL_age + is.Female + Sample_Time + (1|PATNO)"
  DeconvoRadar(df=pheno3, # data frame with all covariate data, colnames should match variable names within reg.model and cell.types 
               reg_model=clin_reg_model, # the right hand side of the regression model equation
               model_type="lmer", # the type of model you want (options: lm, lmer)
               cell_types=cell_types2, # which cell types do you want to plot, list in the order you want them
               group_var="clin_quant", # what is the grouping variable
               ref_group=ref_group_clin, # what is the reference group
               case_groups=case_group_clin, # which groups do you want to plot against the reference (reccomend no more than 2 or 3)
               group_colours_user = tmp_cols, # colors for plot, name the vector the names of the groups
               save_output_csv = TRUE, # save the results in csv table
               save_output_pdf = TRUE, # save a pdf version of the plot
               file_output_prefix = paste0("QuantileComp_Model4_", clin,"_",timept), # string to attach to saved file names
               output_csv_dir = getwd(), # path to folder to save files in
               grid_label_size_user = 6.2, # size of the ring value labels
               axis_label_size_user = 5.5,
               plot_title = paste0("Model4_",clin, " ",timept),
               legend_text_size_user = 1,
               make_sum_equal_one = FALSE, # adjust proportions so they add to 1
               confint.level=0.95)
  
}

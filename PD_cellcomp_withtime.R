library(readr)
library(dplyr)
library(ggforestplot)
library(ggplot2)
library(ggpubr)
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
pheno <- pheno[pheno$Prod_conv_rmv == 0,]

pheno$Sample_Time <- pheno$Sample_Time/12
library(lmerTest)
first <- TRUE
for (case_state in c("HC","PD","Prod")){
  for (cell in cell_types2){
    message(paste0(case_state," ",cell))
    tmp_pheno <- pheno[pheno$DX_new_char == case_state,]
    tmp_cell_data <- tmp_pheno[, cell]
    tmp_model <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + Sample_Time + (1|PATNO),
                              data = tmp_pheno)
    tmp_CIs <- confint(tmp_model)
    tmp_model <- summary(tmp_model)$coefficients
    tmp_CIs <- tmp_CIs[rownames(tmp_model),]
    tmp_model <- data.frame(cbind(tmp_model,tmp_CIs))
    tmp_model$cell <- cell
    tmp_model$case <- case_state
    tmp_model$covariate <- rownames(tmp_model)
    rownames(tmp_model) <- NULL
    tmp_model$analysis_rundate <- Sys.time()
    if (first == TRUE){
      agg_data <- tmp_model
      first <- FALSE
    } else {
      agg_data <- rbind(agg_data,tmp_model)
    }
  }
}

colnames(agg_data) <- c("Estimate","SE","df","t.value","p.value","CI_lo","CI_hi","cell","case","covariate","analysis_rundate")
agg_data2 <- agg_data[agg_data$covariate == "Sample_Time",]
agg_data2$case <- factor(agg_data2$case, levels = c("HC","Prod","PD"))
tmp_cols <- c(HC='#1DA072', Prod ="#3971B5",PD="#C56124")

pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp/CelloverTime.pdf", height = 5, width = 6,onefile = TRUE)

cell_agg1 <- agg_data2[which(agg_data2$cell %in% c("Bas","Bmem","Bnv","CD4mem","CD4nv")),]

p <- forestplot(cell_agg1,
                estimate = Estimate,
                logodds = FALSE,
                colour = case,
                name = cell,
                se = SE,
                ci = 0.95,
                xlab = "Change in Cell Composition per Year",
                xlim = c(-0.3, 0.3)) 
p <- ggpar(p, palette = tmp_cols)
print(p)

cell_agg1 <- agg_data2[which(agg_data2$cell %in% c( "CD8mem","CD8nv","Eos","Mono","NK")),]

p <- forestplot(cell_agg1,
                estimate = Estimate,
                logodds = FALSE,
                colour = case,
                name = cell,
                se = SE,
                ci = 0.95,
                xlab = "Change in Cell Composition per Year",
                xlim = c(-0.3, 0.3)) 
p <- ggpar(p, palette = tmp_cols)
print(p)

cell_agg1 <- agg_data2[which(agg_data2$cell %in% c( "Neu")),]

p <- forestplot(cell_agg1,
                estimate = Estimate,
                logodds = FALSE,
                colour = case,
                name = cell,
                se = SE,
                ci = 0.95,
                xlab = "Change in Cell Composition per Year",) 
p <- ggpar(p, palette = tmp_cols)
print(p)

dev.off()


pdf("//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp/CelloverTime_justNeu.pdf", height = 2, width = 6,onefile = TRUE)
cell_agg1 <- agg_data2[which(agg_data2$cell %in% c( "Neu")),]

p <- forestplot(cell_agg1,
                estimate = Estimate,
                logodds = FALSE,
                colour = case,
                name = cell,
                se = SE,
                ci = 0.95,
                xlab = "Change in Cell Composition per Year",) 
p <- ggpar(p, palette = tmp_cols)
print(p)
dev.off()

write.csv(agg_data, file ="//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/CellComp/CelloverTime.csv" )
















######## Testing interaction 
library(lmerTest)
library(lmtest)
tmp_pheno <- pheno[pheno$DX_new_char != "Prod",]
first <- TRUE
  for (cell in cell_types2){
    tmp_cell_data <- tmp_pheno[, cell]
    tmp_model1 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + DX_new_char + Sample_Time + (1|PATNO),
                      data = tmp_pheno)
    tmp_model2 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + DX_new_char*Sample_Time + (1|PATNO),
                       data = tmp_pheno)
    
    tmp_results <- lrtest(tmp_model2, tmp_model1)
    if(tmp_results[2,5] < 0.05){
      message("Sig Interaction: ", cell, " pval: ", tmp_results[2,5])
    }
  }


tmp_pheno <- pheno[pheno$DX_new_char != "PD",]
first <- TRUE
for (cell in cell_types2){
  tmp_cell_data <- tmp_pheno[, cell]
  tmp_model1 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + DX_new_char + Sample_Time + (1|PATNO),
                     data = tmp_pheno)
  tmp_model2 <- lmer(tmp_cell_data ~ BL_age + is.Female + weird_batch + rundate_binary + DX_new_char*Sample_Time + (1|PATNO),
                     data = tmp_pheno)
  
  tmp_results <- lrtest(tmp_model2, tmp_model1)
  if(tmp_results[2,5] < 0.05){
    message("Sig Interaction: ", cell, " pval: ", tmp_results[2,5])
  }
}
    

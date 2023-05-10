
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

library(readr)
library(wateRmelon)
library(minfi)

CMVest <- read_csv("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/CMVest.csv")
CMVest3<-CMVest$Coefficient
names(CMVest3)<-CMVest$Probe

#Functions
load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
pheno <- pheno[pheno$rmv_sample == 0,]

first <- TRUE
for (i in unique(pheno$processing_batch)){
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i)
  setwd(set_name)
  betas <- loadRData(paste0(set_name, "_betas.RDA"))
  betas <- logit2(betas)
  betas <- betas[, which(colnames(betas) %in% pheno$ID)]
  y_pred<-agep(betas = betas, coeff = CMVest3)
  tmp_class_pred <- ifelse(y_pred$custom_age > 0.50, "Yes", "No")
  tmp_class_pred <- factor(tmp_class_pred, c("No", "Yes"))
  names(tmp_class_pred) <- colnames(betas)
  if (first == TRUE){
    class_pred <- tmp_class_pred
    first <- FALSE
  } else {
    class_pred <- c(class_pred,tmp_class_pred)
  }
}



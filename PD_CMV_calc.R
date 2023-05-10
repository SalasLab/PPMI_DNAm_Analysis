library(readr)
library(wateRmelon)
library(minfi)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

CMVest <- read_csv("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/CMVest.csv")
CMVest3<-CMVest$Coefficient
names(CMVest3)<-CMVest$Probe

load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/pheno_annotated2.RDA")
setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/")
first <- TRUE
for(i in unique(pheno$processing_batch)){
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i)
  setwd(set_name)
  print(set_name)
  betas <- loadRData(paste0(set_name, "_betas.RDA"))
  betas <- logit2(betas)
  
  y_pred<-agep(betas = betas, coeff = CMVest3)
  class_pred <- ifelse(y_pred$custom_age > 0.50, "Yes", "No")
  class_pred <- factor(class_pred, c("No", "Yes"))
  names(class_pred) <- colnames(betas)
  
  if (first == TRUE){
    CMV_status <- class_pred
    first <- FALSE
  } else {
    CMV_status <- c(CMV_status,class_pred)
  }
}

CMV_status2 <- CMV_status
CMV_status2 <- CMV_status2[order(names(CMV_status2))]
pheno <- pheno[order(pheno$ID),]

identical(pheno$ID, names(CMV_status2))

CMV_status3 <- data.frame(ID = names(CMV_status2),
                          CMV_Status= CMV_status2)

setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Feb2023/")
save(CMV_status3, file = "CMV_status.RDA")
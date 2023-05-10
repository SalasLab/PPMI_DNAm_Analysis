### re-ran 19 and 20 because of batch effect issue I found

# HPC
library(minfi)
library(dplyr)
library(readr)
library(ENmix)
library(FlowSorted.BloodExtended.EPIC)
library(ewastools)

#Functions
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Raw Data/dataETOC2.Rd") # this loads the CpG information
epiTOC2 <- function(data.m,ages.v=NULL){
  cpgETOC.v <- dataETOC2.l[[2]];
  estETOC2.m <- dataETOC2.l[[1]];
  soloCpG.v <- dataETOC2.l[[3]];
  ### do epiTOC
  common.v <- intersect(rownames(data.m),cpgETOC.v);
  print(paste("Number of represented epiTOC CpGs (max=385)=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  pcgtAge.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
  ### do epiTOC2
  map.idx <- match(rownames(estETOC2.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Number of represented epiTOC2 CpGs (max=163)=",length(rep.idx),sep=""))
  tmp.m <- data.m[map.idx[rep.idx],];
  TNSC.v <- 2*colMeans(diag(1/(estETOC2.m[rep.idx,1]*(1-estETOC2.m[rep.idx,2]))) %*% (tmp.m - estETOC2.m[rep.idx,2]),na.rm=TRUE);
  TNSC2.v <- 2*colMeans(diag(1/estETOC2.m[rep.idx,1]) %*% tmp.m,na.rm=TRUE);
  ### do HypoClock
  common.v <- intersect(rownames(data.m),soloCpG.v);
  print(paste("Number of represented solo-WCGWs (max=678)=",length(common.v),sep=""));
  map.idx <- match(common.v,rownames(data.m));
  hypoSC.v <- colMeans(data.m[map.idx,],na.rm=TRUE);
  
  estIR.v <- NULL; estIR2.v <- NULL;
  estIR <- NULL;  estIR2 <- NULL;
  if(!is.null(ages.v)){
    estIR.v <- TNSC.v/ages.v;
    estIR <- median(estIR.v,na.rm=TRUE);
    estIR2.v <- TNSC2.v/ages.v;
    estIR2 <- median(estIR2.v,na.rm=TRUE);
  }
  
  return(list(tnsc=TNSC.v,tnsc2=TNSC2.v,irS=estIR.v,irS2=estIR2.v,irT=estIR,irT2=estIR2,pcgtAge=pcgtAge.v,hypoSC=hypoSC.v));
}
getProbStats <- function(betas){
  NA_counts <- c()
  for (i in 1:ncol(betas)){
    NA_counts[i] <- sum(is.na(betas[,i]))
    names(NA_counts)[i] <- colnames(betas)[i]
  }
  total_dim <- dim(betas)
  probe_var <- rowVars(betas, na.rm=TRUE)
  names(probe_var) <- rownames(betas)
  probe_mean <- rowMeans(betas, na.rm = TRUE)
  names(probe_var) <- rownames(betas)
  probe_names <- rownames(betas)
  complete_probe_names <- rownames(betas[complete.cases(betas),])
  median_beta <- colMedians(betas, na.rm=TRUE)
  names(median_beta) <- colnames(betas)
  mean_beta <- colMeans(betas, na.rm=TRUE)
  names(median_beta) <- colnames(betas)
  probe_stats <- list(NA_counts, total_dim, probe_var, probe_names, complete_probe_names, median_beta,mean_beta)
  names(probe_stats) <- c("NA_counts", "Dimension","Variance","Probes","Complete_Probes","Median_Beta","Mean_Beta")
  return(probe_stats)
}
load("/dartfs/rc/lab/S/SalasLab/Glioma_UCSF/Raw Data/PredictFetalSignature.RData")
PredictFetalSignature2<-function(Y) {
  projectWBCnew = function (Y, coefCellType, contrastCellType = NULL, nonnegative = TRUE,
                            lessThanOne = T) 
  {
    if (is.null(contrastCellType)) 
      Xmat <- coefCellType
    else Xmat <- tcrossprod(coefCellType, contrastCellType)
    nCol <- dim(Xmat)[2]
    if (nCol == 2) {
      Dmat <- crossprod(Xmat)
      mixCoef <- t(apply(Y, 2, function(x) {
        solve(Dmat, crossprod(Xmat, x))
      }))
      colnames(mixCoef) <- colnames(Xmat)
      return(mixCoef)
    }
    else {
      nSubj <- dim(Y)[2]
      mixCoef <- matrix(0, nSubj, nCol)
      rownames(mixCoef) <- colnames(Y)
      colnames(mixCoef) <- colnames(Xmat)
      if (nonnegative) {
        if (lessThanOne) {
          Amat <- cbind(rep(-1, nCol), diag(nCol))
          b0vec <- c(-1, rep(0, nCol))
        }
        else {
          Amat <- diag(nCol)
          b0vec <- rep(0, nCol)
        }
        for (i in 1:nSubj) {
          obs <- which(!is.na(Y[, i]))
          Dmat <- crossprod(Xmat[obs, ])
          mixCoef[i, ] <- solve.QP(Dmat, crossprod(Xmat[obs, 
          ], Y[obs, i]), Amat, b0vec)$sol
        }
      }
      else {
        for (i in 1:nSubj) {
          obs <- which(!is.na(Y[, i]))
          Dmat <- crossprod(Xmat[obs, ])
          mixCoef[i, ] <- solve(Dmat, t(Xmat[obs, ]) %*% 
                                  Y[obs, i])
        }
      }
      return(mixCoef)
    }
  }
  
  invarCpGs = rownames(meanMethSignature)
  int = intersect(rownames(Y), invarCpGs)
  overlap = length(int)
  print(paste("Of the 27 invariant CpGs, ", overlap, " out of 27 were contained in the supplied data set", 
              sep = ""))
  targetData = Y[int,]
  projData = meanMethSignature[int,]
  pred = projectWBCnew(targetData,  projData)
  pred*100
}

baseDir <- "/dartfs/rc/lab/S/SalasLab/PD/Raw_data/New_Raw_IDAT_Files"
targets <- read.metharray.sheet(base = baseDir, pattern = "samplesheet.csv")


targets$visit_ID <- paste0(targets$PATNO, "_", targets$EVENT_ID)
tmp <- table(targets$visit_ID)
tmp <- tmp[tmp>1]
tmp <- names(tmp)
targets <- targets[-which(targets$visit_ID %in% tmp),]

Demo <- read_csv("/dartfs/rc/lab/S/SalasLab/PD/Raw_data/Demographics_edited.csv", 
                 col_types = cols(PATNO = col_character()))
Demo <- Demo[, c("PATNO","is.Female","BL_age","Gender_reported")]

targets <- left_join(targets, Demo, by = "PATNO")
targets <- targets[!is.na(targets$is.Female),]
targets <- targets[!is.na(targets$BL_age),]

targets$Sample_Time <- NA
for (i in 1:nrow(targets)){
  if (targets$EVENT_ID[i] == "BL"){
    targets$Sample_Time[i] <- 0
  } else if (targets$EVENT_ID[i] == "V03"){
    targets$Sample_Time[i] <- 9
  } else if (targets$EVENT_ID[i] == "V04"){
    targets$Sample_Time[i] <- 12
  } else if (targets$EVENT_ID[i] == "V05"){
    targets$Sample_Time[i] <- 18
  } else if (targets$EVENT_ID[i] == "V06"){
    targets$Sample_Time[i] <- 24
  } else if (targets$EVENT_ID[i] == "V07"){
    targets$Sample_Time[i] <- 30
  } else if (targets$EVENT_ID[i] == "V08"){
    targets$Sample_Time[i] <- 36
  } else if (targets$EVENT_ID[i] == "V09"){
    targets$Sample_Time[i] <- 42
  }
}
targets$Current_age <- targets$BL_age + (targets$Sample_Time/12)



# odd_samps <- targets[startsWith(targets$Basename, "c"),]
# write.csv(odd_samps, file ="//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/New_Raw_IDAT_Files/see_what_isnt_found.csv")


# repeat but load in see_what_isnt_found2
#odd_samps2 <- targets[startsWith(targets$Basename, "c"),]
#(odd_samps2, file ="//dartfs.dartmouth.edu/rc/lab/S/SalasLab/PD/Raw_data/New_Raw_IDAT_Files/see_what_isnt_found2.csv")



PD_status <- read_csv("/dartfs/rc/lab/S/SalasLab/PD/Raw_data/Participant_Status.csv", 
                      col_types = cols(PATNO = col_character()))
PD_status <- PD_status[,c("PATNO","CONCOHORT","CONCOHORT_DEFINITION")]

targets <- left_join(targets, PD_status, by = "PATNO")
targets <- targets[-which(targets$CONCOHORT == 0),]



save(targets, file = "/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/targets.rda")


samps_per_cycle <- 100
stop_at <- nrow(targets)
cycles <- ceiling(stop_at/samps_per_cycle)

load("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/SWEDD_Paper1/FlowSorted.BloodExtended.EPIC.compTable.rda")
FlowSorted.BloodExtended.EPIC.compTable <- FlowSorted.BloodExtended.EPIC.compTable[order(rownames(FlowSorted.BloodExtended.EPIC.compTable)),]

load("/dartfs/rc/lab/S/SalasLab/Annotation files/annotationEPICb5updated2.rda")
rm(annotDF)
annot$probeType <- substr(annot$Name, 1, 2)

bad_rows <- rownames(annot[annot$CHR=="Y" | 
                             annot$CHR=="X" | 
                             annot$probeType=="rs" |
                             annot$mask==TRUE,])
rm(annot)

for (i in 19:22){
  start_row <- (i-1)*samps_per_cycle + 1
  end_row <- min(i*samps_per_cycle, stop_at)
  message("Processing samples ", start_row, " to ", end_row)
  
  setwd("/dartfs/rc/lab/S/SalasLab/PD/Processed_data/Aug2022/")
  set_name <- paste0("all_batch", i,"_retry")
  message("starting ", set_name)
  dir.create(set_name)
  setwd(set_name)
  
  # Read in data with ewastools package to look at QC and sex prediction
  filenames <- targets$Basename[start_row:end_row]
  read_idats_output <- read_idats(filenames, quiet = TRUE)
  ctrl_metric <- data.frame(control_metrics(read_idats_output))
  rownames(ctrl_metric) <- read_idats_output$meta$sample_id
  ctrl_metric$run_date <- read_idats_output$meta$date
  assign(x = paste0(set_name, "_EWAStools_metrics"), value = ctrl_metric)
  save(list = paste0(set_name, "_EWAStools_metrics"), file = paste0(set_name, "_EWAStools_metrics.RDA"))
  rm(list = paste0(set_name, "_EWAStools_metrics"))
  print("Extracted EWAStools QC Data")
  
  tmp <- check_sex(read_idats_output)
  inferred_sex <- predict_sex(tmp$X, tmp$Y)
  rm(tmp)
  rm(read_idats_output)
  print("Predicted Sex EWAStools")
  
  
  # Read in data with minfi
  RGSet_raw <- read.metharray.exp(targets=targets[start_row:end_row, ], extended = TRUE, recursive = TRUE, force = TRUE)
  saveRDS(RGSet_raw, file = paste0(set_name, "_RGset_raw.RDS"))
  print("Extracted Raw RGSet")
  
  ctrls <- getProbeInfo(RGSet_raw, type = "Control")
  ctrls <- ctrls[ctrls$Address %in% featureNames(RGSet_raw), 
  ]
  ctrl_r <- assays(RGSet_raw)$Red[ctrls$Address[!(ctrls$Type %in% 
                                                    "NEGATIVE")], ]
  ctrl_g <- assays(RGSet_raw)$Green[ctrls$Address[!(ctrls$Type %in% 
                                                      "NEGATIVE")], ]
  ctrl_nneg <- rbind(ctrl_r, ctrl_g)
  assign(x = paste0(set_name, "_neg_ctrl_probes"), value = ctrl_nneg)
  save(list = paste0(set_name, "_neg_ctrl_probes"), file = paste0(set_name, "_neg_ctrl_probes.RDA"))
  rm(list = paste0(set_name, "_neg_ctrl_probes"))
  rm(ctrls,ctrl_r,ctrl_g,ctrl_nneg)
  
  pheno_data <- pData(RGSet_raw)
  pheno_data$processing_batch <- i
  pheno_data$within_batch_order <- seq(1,nrow(pheno_data),1)
  
  assign(x = paste0(set_name, "_basepheno"), value = pheno_data)
  save(list = paste0(set_name, "_basepheno"), file = paste0(set_name, "_basepheno.RDA"))
  rm(list = paste0(set_name, "_basepheno"))
  print("Extracted pheno data")
  
  pheno_data$ewastools_sex <- inferred_sex
  rm(inferred_sex)
  
  # ENmix QC
  qc <- QCinfo(RGSet_raw, detPtype="oob", detPthre=0.05)
  assign(x = paste0(set_name, "_qc_data"), value = qc)
  save(list = paste0(set_name, "_qc_data"), file = paste0(set_name, "_qc_data.RDA"))
  rm(list = paste0(set_name, "_qc_data"))
  print("Calculated qc data")
  
  Mset_Noob <- preprocessNoob(RGSet_raw)

  betas <- getBeta(Mset_Noob)
  rm(Mset_Noob)
  assign(x = paste0(set_name, "_betas"), value = betas)
  save(list = paste0(set_name, "_betas"), file = paste0(set_name, "_betas.RDA"))
  rm(list = paste0(set_name, "_betas"))
  betas_filtered <- qcfilter(betas, qcscore = qc, rmcr = FALSE, rthre = 0.999, cthre = 0.15)
  betas_filtered <- betas_filtered[-which(rownames(betas_filtered) %in% bad_rows),]
  assign(x = paste0(set_name, "_betas_filtered"), value = betas_filtered)
  save(list = paste0(set_name, "_betas_filtered"), file = paste0(set_name, "_betas_filtered.RDA"))
  rm(list = paste0(set_name, "_betas_filtered"))
  print("Calculated betas")
  rm(qc)
  
  # Predict blood cell
  cpginter <- intersect(rownames(FlowSorted.BloodExtended.EPIC.compTable), rownames(betas))
  tmp_betas <- betas[which(rownames(betas) %in% cpginter),]
  tmp_flowsrt <- FlowSorted.BloodExtended.EPIC.compTable[which(rownames(FlowSorted.BloodExtended.EPIC.compTable) %in% cpginter),]
  tmp_betas <- tmp_betas[order(rownames(tmp_betas)),]
  tmp_flowsrt <- tmp_flowsrt[order(rownames(tmp_flowsrt)),]
  message(nrow(FlowSorted.BloodExtended.EPIC.compTable) - nrow(tmp_flowsrt)," CpGs were missing in data")
  countsEPIC_ext <- round(projectCellType_CP(tmp_betas, tmp_flowsrt,lessThanOne =T)*100,3)
  pheno_data <- cbind(pheno_data, countsEPIC_ext)
  print("Calculated cell composition")
  rm(countsEPIC_ext,tmp_betas,tmp_flowsrt)
  
  # missingness information
  probeStats <- getProbStats(betas)
  probeStats_filtered <- getProbStats(betas_filtered)
  assign(x = paste0(set_name, "_probeStatsfiltered"), value = probeStats_filtered)
  save(list = paste0(set_name, "_probeStatsfiltered"), file = paste0(set_name, "_probeStatsfiltered.RDA"))
  rm(list = paste0(set_name, "_probeStatsfiltered"))
  assign(x = paste0(set_name, "_probeStats"), value = probeStats)
  save(list = paste0(set_name, "_probeStats"), file = paste0(set_name, "_probeStats.RDA"))
  rm(list = paste0(set_name, "_probeStats"))
  rm(probeStats,probeStats_filtered)
  print("Calculated ProbeStats")
  
  # mitotic clock
  epiTOC2_results <- epiTOC2(betas, pheno_data$Current_age)
  epiTOC2_results <- as.data.frame(do.call(cbind,epiTOC2_results))
  pheno_data <- cbind(pheno_data, epiTOC2_results)
  print("Calculated EpiTOC")
  
  # FCO
  FCO_results <- PredictFetalSignature2(betas)
  colnames(FCO_results) <- paste0("FCO_",colnames(FCO_results))
  pheno_data <- cbind(pheno_data, FCO_results)
  print("Calculated FCO")
  
  # Methylation Age
  mAge <- methyAge(betas, fastImputation = FALSE, normalize = FALSE)
  mAge <- mAge[,-1]
  pheno_data <- cbind(pheno_data, mAge)
  print("Calculated mAge")
  
  
  # phenotype data
  pheno_data <- as.data.frame(pheno_data)
  assign(x = paste0(set_name, "_pheno"), value = pheno_data)
  save(list = paste0(set_name, "_pheno"), file = paste0(set_name, "_pheno.RDA"))
  rm(pheno_data)
  rm(list = paste0(set_name, "_pheno"))
  print("Finished this batch!")
  
}



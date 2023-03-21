#---- Objects and functions ----
#Libraries
library(geepack)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(export)

#paths
path <- ".../"
path_pheno <- ".../Phenotypes/"
path_out <- ".../Results/"

#Obtain subjects IDs from .fam files
omniID <- data.table::fread(paste0(path, "methodOmni.fam"))[,1:2] %>% setNames(., c("IID", "FID"))
gsaID <- data.table::fread(paste0(path, "methodGSA.fam"))[,1:2] %>% setNames(., c("IID", "FID"))

#Output data frame initialization
out <- data.frame()
outColnames <- c("Phenotype", "Methode", "R2", "R2Liab", "AUC", "AUC_lower", "AUC_upper", "OR", "Pvalue", "Upper", "Lower")
out2 <- data.frame()
outColnames2 <- c("Phenotype", "Methode", "Quantile", "OR", "Pvalue", "Upper", "Lower")
out3 <- data.frame()
outColnames3 <- c("Phenotype", "Methode", "Quantile", "OR", "Pvalue", "Upper", "Lower")

#This function generates our phenotypes.
addPheno <- function(omni, GSA, pheno = c("SKZBIP", "BIPavecsans", "SKZaff"), methods, omniOnly = FALSE){
  Loop <- 1:2 #1 us for OmniExpress subject, 2 is for GSA subjects
  liste_naars_elargie <- fread(paste0(path_pheno, "liste_NAARs_elargie.txt"))
  phenos_don <- fread(paste0(path_pheno, "phenos_liabs.txt"), header = TRUE)
  
  if(pheno == "SKZBIP"){
    for(i in Loop){
      if(i == 1){
        prs <- omni
      }else{
        prs <- GSA
      }
      BIP <- fread(paste0(path_pheno, "BPBroad.pre"), header = FALSE)
      SKZ <- fread(paste0(path_pheno, "SZBroad.pre"), header = FALSE)
      
      #Add ng a variable for NAARs
      #SKZ, omni : 239 sujets (59 cas, 180 temoins)
      #   , gsa:    327 sujets (65 cas, 262 temoins)
      prs$NonAtteint <- ifelse(prs$IID %in% liste_naars_elargie$V1, yes = 1, no = 0)
      prs$SKZ <- ifelse(prs$IID %in% SKZ$V2[SKZ$V6 == 2], yes = 2, no = 0)
      prs$SKZ[prs$NonAtteint == 1] <- 1
      prs$SKZ <- prs$SKZ-1
      prs_SKZ <- prs[prs$SKZ>=0,] %>% dplyr::select(., -NonAtteint)
      prs_SKZ <- prs_SKZ[!is.na(prs_SKZ$SKZ),] #si pheno manquant
      
      #BIP, omni : 274 sujets (94  cas, 180 temoins) 
      #   ,   gsa:  373 sujets (111 cas, 262 temoins)
      prs$BIP <- ifelse(prs$IID %in% BIP$V2[BIP$V6 == 2], yes = 2, no = 0)
      prs$BIP[prs$NonAtteint == 1] <- 1
      prs$BIP <- prs$BIP-1
      prs_BIP <- prs[prs$BIP>=0,] %>% dplyr::select(., -SKZ, -NonAtteint)
      prs_BIP <- prs_BIP[!is.na(prs_BIP$BIP),] #si pheno manquant
      
      cond_SKZ_NAAR <- prs_SKZ$SKZ==0
      cond_BIP_NAAR <- prs_BIP$BIP==0
      for(met in methods){
        prs_SKZ[,met] <- (prs_SKZ[,met]-mean(prs_SKZ[cond_SKZ_NAAR, met], na.rm = TRUE))/sd(prs_SKZ[cond_SKZ_NAAR, met], na.rm = TRUE)
        prs_BIP[,met] <- (prs_BIP[,met]-mean(prs_BIP[cond_BIP_NAAR, met], na.rm = TRUE))/sd(prs_BIP[cond_BIP_NAAR, met], na.rm = TRUE)
      }
      
      if(i == 1){
        meriem_prs <- prs
        meriem_prs_SKZ <- prs_SKZ
        meriem_prs_BIP <- prs_BIP
        prs_SKZ <- data.frame(meriem_prs_SKZ, puce = "Omni")
        prs_BIP <- data.frame(meriem_prs_BIP, puce = "Omni")
        prs_SKZ$N_Cas <- sum(prs_SKZ$SKZ == 1)
        prs_SKZ$N_Temoin <- sum(prs_SKZ$SKZ == 0)
        prs_BIP$N_Cas <- sum(prs_BIP$BIP == 1)
        prs_BIP$N_Temoin <- sum(prs_BIP$BIP == 0)
      }else{
        gsa_prs <- prs
        gsa_prs_SKZ <- prs_SKZ
        gsa_prs_BIP <- prs_BIP
        prs_SKZ <- rbind(data.frame(meriem_prs_SKZ, puce = "Omni"),
                         data.frame(gsa_prs_SKZ, puce = "GSA"))
        prs_BIP <- rbind(data.frame(meriem_prs_BIP, puce = "Omni"),
                         data.frame(gsa_prs_BIP, puce = "GSA"))
        prs_SKZ$N_Cas <- sum(prs_SKZ$SKZ == 1)
        prs_SKZ$N_Temoin <- sum(prs_SKZ$SKZ == 0)
        prs_BIP$N_Cas <- sum(prs_BIP$BIP == 1)
        prs_BIP$N_Temoin <- sum(prs_BIP$BIP == 0)
      }
    }
    out <- list("SKZ" = prs_SKZ, "BIP" = prs_BIP)
    
  }else if(pheno == "BIPavecsans"){
    for(i in Loop){
      if(i == 1){
        prs <- omni
      }else{
        prs <- GSA
      }
      BIP <- readRDS("Ind-BIP-avec-sans-psychose.RDS")
      BIPavec <- BIP$V2[BIP$BIP_avec_psychose == 1]
      BIPsans <- BIP$V2[BIP$BIP_sans_psychose == 1]
      
      #BIP with, omni: 226 sujets (46 cas, 180 temoins)
      #        , gsa:    309 sujets (47 cas, 262 temoins)
      prs$NonAtteint <- ifelse(prs$IID %in% liste_naars_elargie$V1, yes = 1, no = 0)
      prs$BIPavec <- ifelse(prs$IID %in% BIPavec, yes = 2, no = 0)
      prs$BIPavec[prs$NonAtteint == 1] <- 1
      prs$BIPavec <- prs$BIPavec-1
      prs_BIPavec <- prs[prs$BIPavec>=0,] %>% dplyr::select(., -NonAtteint)
      prs_BIPavec <- prs_BIPavec[!is.na(prs_BIPavec$BIPavec),] #si pheno manquant
      
      #BIP without, omni : 228 sujets (48 cas, 180 temoins) 
      #        ,      gsa:  326 sujets (64 cas, 262 temoins)
      prs$BIPsans <- ifelse(prs$IID %in% BIPsans, yes = 2, no = 0)
      prs$BIPsans[prs$NonAtteint == 1] <- 1
      prs$BIPsans <- prs$BIPsans-1
      prs_BIPsans <- prs[prs$BIPsans>=0,] %>% dplyr::select(., -NonAtteint, -BIPavec)
      prs_BIPsans <- prs_BIPsans[!is.na(prs_BIPsans$BIPsans),] #si pheno manquant
      
      cond_BIPavec_NAAR <- prs_BIPavec$BIPavec==0
      cond_BIPsans_NAAR <- prs_BIPsans$BIP==0
      for(met in methods){
        prs_BIPavec[,met] <- (prs_BIPavec[,met]-mean(prs_BIPavec[cond_BIPavec_NAAR, met], na.rm = TRUE))/sd(prs_BIPavec[cond_BIPavec_NAAR, met], na.rm = TRUE)
        prs_BIPsans[,met] <- (prs_BIPsans[,met]-mean(prs_BIPsans[cond_BIPsans_NAAR, met], na.rm = TRUE))/sd(prs_BIPsans[cond_BIPsans_NAAR, met], na.rm = TRUE)
      }
      if(i == 1){
        meriem_prs <- prs
        meriem_prs_BIPavec <- prs_BIPavec
        meriem_prs_BIPsans <- prs_BIPsans
        prs_BIPavec <- data.frame(meriem_prs_BIPavec, puce = "Omni")
        prs_BIPsans <- data.frame(meriem_prs_BIPsans, puce = "Omni")
        prs_BIPavec$N_Temoin <- sum(prs_BIPavec$BIPavec == 0)
        prs_BIPavec$N_Cas <- sum(prs_BIPavec$BIPavec == 1)
        prs_BIPsans$N_Temoin <- sum(prs_BIPsans$BIPsans == 0)
        prs_BIPsans$N_Cas <- sum(prs_BIPsans$BIPsans == 1)
      }else{
        gsa_prs <- prs
        gsa_prs_BIPavec <- prs_BIPavec
        gsa_prs_BIPsans <- prs_BIPsans
        prs_BIPavec <- rbind(data.frame(meriem_prs_BIPavec, puce = "Omni"),
                             data.frame(gsa_prs_BIPavec, puce = "GSA"))
        prs_BIPsans <- rbind(data.frame(meriem_prs_BIPsans, puce = "Omni"),
                             data.frame(gsa_prs_BIPsans, puce = "GSA"))
        prs_BIPavec$N_Temoin <- sum(prs_BIPavec$BIPavec == 0)
        prs_BIPavec$N_Cas <- sum(prs_BIPavec$BIPavec == 1)
        prs_BIPsans$N_Temoin <- sum(prs_BIPsans$BIPsans == 0)
        prs_BIPsans$N_Cas <- sum(prs_BIPsans$BIPsans == 1)
      }
    }
    out <- list("BIPavec" = prs_BIPavec, "BIPsans" = prs_BIPsans)
    
  }else if(pheno == "SKZaff"){
    for(i in Loop){
      if(i == 1){
        prs <- omni
      }else{
        prs <- GSA
      }
      GC <- data.table::fread(paste0(path_pheno, "GCbroad.pre"))
      GC <- GC$V2[GC$V6 == 2]
      rm1 <- data.table::fread(paste0(path_pheno, "SZbroad.pre"))
      rm1 <- rm1$V2[rm1$V6 == 2]
      rm2 <- data.table::fread(paste0(path_pheno, "BPBroad.pre"))
      rm2 <- rm2$V2[rm2$V6 == 2]
      rm <- union(rm1,rm2)
      SKZ <- setdiff(GC, rm)
      
      #SAD,      omni: 198 sujets (18 cas, 180 temoins)
      #        , gsa :    288 sujets (26 cas, 262 temoins)
      prs$NonAtteint <- ifelse(prs$IID %in% liste_naars_elargie$V1, yes = 1, no = 0)
      prs$SKZaffectif <- ifelse(prs$IID %in% SKZ, yes = 2, no = 0)
      prs$SKZaffectif[prs$NonAtteint == 1] <- 1
      prs$SKZaffectif <- prs$SKZaffectif-1
      prs_SKZaffectif <- prs[prs$SKZaffectif>=0,] %>% dplyr::select(., -NonAtteint)
      prs_SKZaffectif <- prs_SKZaffectif[!is.na(prs_SKZaffectif$SKZaffectif),] #si pheno manquant
      
      cond_SKZaffectif_NAAR <- prs_SKZaffectif$SKZaffectif==0
      for(met in methods){
        prs_SKZaffectif[,met] <- (prs_SKZaffectif[,met]-mean(prs_SKZaffectif[cond_SKZaffectif_NAAR, met], na.rm = TRUE))/sd(prs_SKZaffectif[cond_SKZaffectif_NAAR, met], na.rm = TRUE)
      }
      if(i == 1){
        meriem_prs <- prs
        meriem_prs_SKZaffectif <- prs_SKZaffectif
        meriem_prs_SKZaffectif <- prs_SKZaffectif
        prs_SKZaffectif <- data.frame(meriem_prs_SKZaffectif, puce = "Omni")
        prs_SKZaffectif <- data.frame(meriem_prs_SKZaffectif, puce = "Omni")
        prs_SKZaffectif$N_Temoin <- sum(prs_SKZaffectif$SKZaffectif == 0)
        prs_SKZaffectif$N_Cas <- sum(prs_SKZaffectif$SKZaffectif == 1)
      }else{
        gsa_prs <- prs
        gsa_prs_SKZaffectif <- prs_SKZaffectif
        gsa_prs_SKZaffectif <- prs_SKZaffectif
        prs_SKZaffectif <- rbind(data.frame(meriem_prs_SKZaffectif, puce = "Omni"),
                                 data.frame(gsa_prs_SKZaffectif, puce = "GSA"))
        prs_SKZaffectif <- rbind(data.frame(meriem_prs_SKZaffectif, puce = "Omni"),
                                 data.frame(gsa_prs_SKZaffectif, puce = "GSA"))
        prs_SKZaffectif$N_Temoin <- sum(prs_SKZaffectif$SKZaffectif == 0)
        prs_SKZaffectif$N_Cas <- sum(prs_SKZaffectif$SKZaffectif == 1)
      }
    }
    out <-  prs_SKZaffectif
  }
  return(out)
}

#This function computes correlation on a liability scale as proposed by lee et al. 2012
corLiability <- function(lmv, pheno){
  if (pheno == "SKZ")        {K <- 0.01; ncase <- 124}
  if (pheno == "BIP")        {K <- 0.02; ncase <- 205}
  if (pheno == "BIPavec")    {K <- 0.009; ncase <- 93}
  if (pheno == "BIPsans")    {K <- 0.011; ncase <- 112}
  if (pheno == "SKZaffectif"){K <- 0.003; ncase <- 35}
  ncont <- 442
  nt <- ncase+ncont
  P <- ncase/nt
  thd  <- qnorm(K,0,1)
  zv  <- dnorm(thd)
  mv  <- zv/K
  R2O  <- var(lmv$fitted.values)/(ncase/nt*ncont/nt)
  theta  <- mv*(P-K)/(1-K)*(mv*(P-K)/(1-K)-thd)
  cv <-K*(1-K)/zv^2*K*(1-K)/(P*(1-P))
  R2 <-R2O*cv/(1+R2O*theta*cv)
  return(R2)
}

#This function builds a geeglm to analyse the association between PRS and given phenotype.
addModel <- function(data, pheno, methods){
  out <- data.frame()
  for(met in methods){
    form <- paste0(pheno, " ~ ", met)
    #Model
    model  <- geeglm(as.formula(form),
                     data = data, id = factor(data$FID),
                     family = binomial(link = "logit"),
                     corstr = "independence",
                     scale.fix = TRUE,
                     scale.value = rep(1, length(data$FID)))
    modelgaussian  <- geeglm(as.formula(form),
                     data = data, id = factor(data$FID),
                     family = gaussian,
                     corstr = "independence",
                     scale.fix = TRUE,
                     scale.value = rep(1, length(data$FID)))
    summaryModel <- summary(model)
    roc <- pROC::roc(predictor = data[, met], response = as.numeric(data[, pheno]), direction = "<", levels = c(0, 1))
    auc <- pROC::ci.auc(roc)
    add <- c(pheno, met,
             round(cor(model$fitted.values, model$y)^2, 5),
             round(corLiability(modelgaussian, pheno), 5),
             round(as.numeric(auc[2]), 5),
             round(as.numeric(auc[1]), 5),
             round(as.numeric(auc[3]), 5),
             round(exp((summaryModel$coef)$Est[2]),2),
             (summaryModel$coef)$Pr[2],
             round(exp((summaryModel$coef)$Est[2]+1.96*(summaryModel$coef)$Std[2]),2),
             round(exp((summaryModel$coef)$Est[2]-1.96*(summaryModel$coef)$Std[2]),2))
    out <- rbind(out, add)
  }
  return(out)
}

#This function builds a geeglm to analyse the association between PRS quantiles and a given phenotype.
quantileOR <- function(data, pheno, methods, prob = seq(0, 1, 0.25), probRef = 1){
  out <- data.frame()
  for(met in methods){
    #Forming quantiles groups
    prs <- data[,met]
    quantValue <- as.numeric(quantile(prs[data[,pheno] == 0], probs = prob))
    quantValue[1] <- min(prs); quantValue[length(quantValue)] <- max(prs)
    data$quantiles <- cut(prs, breaks = quantValue, include.lowest = TRUE)
    quantilesLevels <- levels(data$quantiles)
    
    #Model
    for(quant in setdiff(1:length(quantilesLevels),probRef)){
      quantData <- data[data$quantiles %in% quantilesLevels[c(probRef,quant)],]
      quantData$quantiles <- droplevels(quantData$quantiles)
      quantData$quantiles <- relevel(quantData$quantiles, ref = quantilesLevels[probRef])
      form <- paste0(pheno, " ~ quantiles")
      model  <- geeglm(as.formula(form),
                       data = quantData, id = factor(quantData$FID),
                       family = binomial(link = "logit"),
                       corstr = "independence",
                       scale.fix = TRUE,
                       scale.value = rep(1, length(quantData$FID)))
      summaryModel <- summary(model)
      add <- c(pheno, met, quant,
               round(exp((summaryModel$coef)$Est[2]),2),
               round((summaryModel$coef)$Pr[2],8),
               round(exp((summaryModel$coef)$Est[2]+1.96*(summaryModel$coef)$Std[2]),2),
               round(exp((summaryModel$coef)$Est[2]-1.96*(summaryModel$coef)$Std[2]),2))
      out <- rbind(out, add)
    }
  }
  return(out)
}

#---- PRS ----
#list.files(path)
#Import PRSs - OmniExpress
#In our
prs_ml_SKZ <- readRDS(paste0(path, "CTOmniSKZ.RDS"))
prs_ml_BIP <- readRDS(paste0(path, "CTOmniBIP.RDS"))
prs_ml_SKZ2 <- readRDS(paste0(path, "ThresholdingOmniSKZ.RDS"))
prs_ml_BIP2 <- readRDS(paste0(path, "ThresholdingOmniBIP.RDS"))
prs_ml_SKZ3 <- readRDS(paste0(path, "multiOmniSKZ.RDS"))
prs_ml_BIP3 <- readRDS(paste0(path, "multiOmniBIP.RDS"))
prs_ml_BIP4 <- readRDS(paste0(path, "OGOmniBIP.RDS"))
prs_ml_SKZ4 <- readRDS(paste0(path, "OGOmniSKZ.RDS"))
prs_ml_SKZ5 <- readRDS(paste0(path, "multiGenCovOmniSKZ.RDS"))
prs_ml_BIP5 <- readRDS(paste0(path, "multiGenCovOmniBIP.RDS"))
prs_ml_SKZ8 <- readRDS(paste0(path, "LDpred2OmniSKZ.Rdata")) %>% as.matrix()
prs_ml_BIP8 <- readRDS(paste0(path, "LDpred2OmniBIP.Rdata")) %>% as.matrix()
prs_ml_SKZ9 <- readRDS(paste0(path, "PANPRSOmniSKZ.RDS")) %>% as.matrix()
prs_ml_BIP9 <- readRDS(paste0(path, "PANPRSOmniBIP.RDS")) %>% as.matrix()

#Import PRSs - GSA
gsa_prs_ml_SKZ <- readRDS(paste0(path, "CTGSASKZ.RDS"))
gsa_prs_ml_BIP <- readRDS(paste0(path, "CTGSABIP.RDS"))
gsa_prs_ml_SKZ2 <- readRDS(paste0(path, "ThresholdingGSASKZ.RDS"))
gsa_prs_ml_BIP2 <- readRDS(paste0(path, "ThresholdingGSABIP.RDS"))
gsa_prs_ml_SKZ3 <- readRDS(paste0(path, "multiGSASKZ.RDS"))
gsa_prs_ml_BIP3 <- readRDS(paste0(path, "multiGSABIP.RDS"))
gsa_prs_ml_SKZ4 <- readRDS(paste0(path, "OGGSASKZ.RDS"))
gsa_prs_ml_BIP4 <- readRDS(paste0(path, "OGGSABIP.RDS"))
gsa_prs_ml_SKZ5 <- readRDS(paste0(path, "multiGenCovGSASKZ.RDS"))
gsa_prs_ml_BIP5 <- readRDS(paste0(path, "multiGenCovGSABIP.RDS"))
gsa_prs_ml_SKZ8 <- readRDS(paste0(path, "LDpred2GSASKZ.Rdata")) %>% as.matrix()
gsa_prs_ml_BIP8 <- readRDS(paste0(path, "LDpred2GSABIP.Rdata")) %>% as.matrix()
gsa_prs_ml_SKZ9 <- readRDS(paste0(path, "PANPRSGSASKZ.RDS")) %>% as.matrix()
gsa_prs_ml_BIP9 <- readRDS(paste0(path, "PANPRSGSABIP.RDS")) %>% as.matrix()

#Methods labels
methods <- c("SCORE_CT_SKZ", "SCORE_CT_BIP", "SCORE_THRESHOLD_SKZ", "SCORE_THRESHOLD_BIP",
             "SCORE_ML_SKZ", "SCORE_ML_BIP", "SCORE_OG_SKZ", "SCORE_OG_BIP",
             "SCORE_GC_SKZ", "SCORE_GC_BIP", "SCORE_LDPRED2_SKZ", "SCORE_LDPRED2_BIP", "SCORE_PP_SKZ", "SCORE_PP_BIP"
)
methodsSKZ <- stringr::str_subset(methods, pattern = "_SKZ")
methodsBIP <- stringr::str_subset(methods, pattern = "_BIP")
nMethods <- length(methodsSKZ)

#Methods (better) labels
methodsToNames <- function(x){
  out <- case_when(x %in% c("SCORE_NAIF_SKZ", "SCORE_NAIF_BIP") ~ "Naif",
                   x %in% c("SCORE_CT_SKZ", "SCORE_CT_BIP") ~ "C+T",
                   x %in% c("SCORE_THRESHOLD_SKZ", "SCORE_THRESHOLD_BIP") ~ "Thresholding",
                   x %in% c("SCORE_ML_SKZ", "SCORE_ML_BIP") ~ "mvL",
                   x %in% c("SCORE_GC_SKZ", "SCORE_GC_BIP") ~ "mvL(BLD-X)",
                   x %in% c("SCORE_LDPRED2_SKZ", "SCORE_LDPRED2_BIP") ~ "LDpred2",
                   x == "SCORE_OG_BIP" ~ "OG",
                   x == "SCORE_OG_SKZ" ~ "OG",
                   x %in% c("SCORE_PP_SKZ", "SCORE_PP_BIP") ~ "PANPRS") %>%
    as.factor(.)
  return(out)
}

#---- BIP ET SKZ ----
#En faire un seul jeu de données de PRS
prs <- data.frame(omniID, prs_ml_SKZ[, 1], prs_ml_BIP[, 1], prs_ml_SKZ2[, 1], prs_ml_BIP2[, 1],
                  prs_ml_SKZ3[, 1], prs_ml_BIP3[, 1], prs_ml_SKZ4[, 1], prs_ml_BIP4[, 1],
                  prs_ml_SKZ5[, 1], prs_ml_BIP5[, 1], prs_ml_SKZ8[, 1], prs_ml_BIP8[, 1], prs_ml_SKZ9[, 1], prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methods))
gsa_prs <- data.frame(gsaID, gsa_prs_ml_SKZ[, 1], gsa_prs_ml_BIP[, 1], gsa_prs_ml_SKZ2[, 1], gsa_prs_ml_BIP2[, 1],
                      gsa_prs_ml_SKZ3[, 1], gsa_prs_ml_BIP3[, 1], gsa_prs_ml_SKZ4[, 1], gsa_prs_ml_BIP4[, 1],
                      gsa_prs_ml_SKZ5[, 1], gsa_prs_ml_BIP5[, 1], gsa_prs_ml_SKZ8[, 1], gsa_prs_ml_BIP8[, 1], gsa_prs_ml_SKZ9[, 1], gsa_prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methods))

SKZBIP <- addPheno(omni = prs, GSA = gsa_prs, pheno = "SKZBIP", methods = methods)
out <- rbind(out, 
             addModel(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsSKZ) %>% setNames(., outColnames),
             addModel(data = SKZBIP$BIP, pheno = "BIP", methods = methodsBIP) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsSKZ) %>% setNames(., outColnames2),
              quantileOR(data = SKZBIP$BIP, pheno = "BIP", methods = methodsBIP) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsSKZ, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2),
              quantileOR(data = SKZBIP$BIP, pheno = "BIP", methods = methodsBIP, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- BIP ET SKZ (with inverse PRSs) ----
out <- rbind(out, 
             addModel(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsBIP) %>% setNames(., outColnames),
             addModel(data = SKZBIP$BIP, pheno = "BIP", methods = methodsSKZ) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsBIP) %>% setNames(., outColnames2),
              quantileOR(data = SKZBIP$BIP, pheno = "BIP", methods = methodsSKZ) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = SKZBIP$SKZ, pheno = "SKZ", methods = methodsBIP, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2),
              quantileOR(data = SKZBIP$BIP, pheno = "BIP", methods = methodsSKZ, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- BIP with-without psychosis----
#En faire un seul jeu de données de PRS
prs <- data.frame(omniID, prs_ml_BIP[, 1], prs_ml_BIP2[, 1],
                  prs_ml_BIP3[, 1], prs_ml_BIP4[, 1],
                  prs_ml_BIP5[, 1], prs_ml_BIP8[, 1], prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsBIP))
gsa_prs <- data.frame(gsaID, gsa_prs_ml_BIP[, 1], gsa_prs_ml_BIP2[, 1],
                      gsa_prs_ml_BIP3[, 1], gsa_prs_ml_BIP4[, 1],
                      gsa_prs_ml_BIP5[, 1], gsa_prs_ml_BIP8[, 1], gsa_prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsBIP))

BIPavecsans <- addPheno(omni = prs, GSA = gsa_prs, pheno = "BIPavecsans", methods = methodsBIP)
out <- rbind(out, 
             addModel(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsBIP) %>% setNames(., outColnames),
             addModel(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsBIP) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsBIP) %>% setNames(., outColnames2),
              quantileOR(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsBIP) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsBIP, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2),
              quantileOR(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsBIP, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- SAD ----
#En faire un seul jeu de données de PRS
prs <- data.frame(omniID, prs_ml_SKZ[, 1], prs_ml_SKZ2[, 1],
                  prs_ml_SKZ3[, 1], prs_ml_SKZ4[, 1],
                  prs_ml_SKZ5[, 1], prs_ml_SKZ8[, 1], prs_ml_SKZ9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsSKZ))
gsa_prs <- data.frame(gsaID, gsa_prs_ml_SKZ[, 1], gsa_prs_ml_SKZ2[, 1],
                      gsa_prs_ml_SKZ3[, 1], gsa_prs_ml_SKZ4[, 1],
                      gsa_prs_ml_SKZ5[, 1], gsa_prs_ml_SKZ8[, 1], gsa_prs_ml_SKZ9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsSKZ))

SKZaff <- addPheno(omni = prs, GSA = gsa_prs, pheno = "SKZaff", methods = methodsSKZ)
out <- rbind(out, 
             addModel(data = SKZaff, pheno = "SKZaffectif", methods = methodsSKZ) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = SKZaff, pheno = "SKZaffectif", methods = methodsSKZ) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = SKZaff, pheno = "SKZaffectif", methods = methodsSKZ, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- BIP with-without psychosis (with SZ PRS) ----
#En faire un seul jeu de données de PRS
prs <- data.frame(omniID, prs_ml_SKZ[, 1], prs_ml_SKZ2[, 1],
                  prs_ml_SKZ3[, 1], prs_ml_SKZ4[, 1],
                  prs_ml_SKZ5[, 1], prs_ml_SKZ8[, 1], prs_ml_SKZ9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsSKZ))
gsa_prs <- data.frame(gsaID, gsa_prs_ml_SKZ[, 1], gsa_prs_ml_SKZ2[, 1],
                      gsa_prs_ml_SKZ3[, 1], gsa_prs_ml_SKZ4[, 1],
                      gsa_prs_ml_SKZ5[, 1], gsa_prs_ml_SKZ8[, 1], gsa_prs_ml_SKZ9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsSKZ))

BIPavecsans <- addPheno(omni = prs, GSA = gsa_prs, pheno = "BIPavecsans", methods = methodsSKZ)
out <- rbind(out, 
             addModel(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsSKZ) %>% setNames(., outColnames),
             addModel(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsSKZ) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsSKZ) %>% setNames(., outColnames2),
              quantileOR(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsSKZ) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = BIPavecsans$BIPavec, pheno = "BIPavec", methods = methodsSKZ, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2),
              quantileOR(data = BIPavecsans$BIPsans, pheno = "BIPsans", methods = methodsSKZ, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- SAD (with BD PRS) ----
#En faire un seul jeu de données de PRS
prs <- data.frame(omniID, prs_ml_BIP[, 1], prs_ml_BIP2[, 1],
                  prs_ml_BIP3[, 1], prs_ml_BIP4[, 1],
                  prs_ml_BIP5[, 1], prs_ml_BIP8[, 1], prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsBIP))
gsa_prs <- data.frame(gsaID, gsa_prs_ml_BIP[, 1], gsa_prs_ml_BIP2[, 1],
                      gsa_prs_ml_BIP3[, 1], gsa_prs_ml_BIP4[, 1],
                      gsa_prs_ml_BIP5[, 1], gsa_prs_ml_BIP8[, 1], gsa_prs_ml_BIP9[, 1]
) %>%
  setNames(., c("FID", "IID", methodsBIP))

SKZaff <- addPheno(omni = prs, GSA = gsa_prs, pheno = "SKZaff", methods = methodsBIP)
out <- rbind(out, 
             addModel(data = SKZaff, pheno = "SKZaffectif", methods = methodsBIP) %>% setNames(., outColnames)
)
out2 <- rbind(out2, 
              quantileOR(data = SKZaff, pheno = "SKZaffectif", methods = methodsBIP) %>% setNames(., outColnames2)
)
out3 <- rbind(out3, 
              quantileOR(data = SKZaff, pheno = "SKZaffectif", methods = methodsBIP, prob = c(0,0.1,0.9,1), probRef = 2) %>% setNames(., outColnames2)
)

#---- Results output format ----
out <- out[,c("Phenotype", "Methode", "R2", "R2Liab", "AUC", "AUC_lower", "AUC_upper", "OR", "Lower", "Upper", "Pvalue")]
out$N_Cas <- c(rep(unique(SKZBIP$SKZ$N_Cas), nMethods), rep(unique(SKZBIP$BIP$N_Cas), nMethods),rep(unique(SKZBIP$SKZ$N_Cas), nMethods), rep(unique(SKZBIP$BIP$N_Cas), nMethods), rep(unique(BIPavecsans$BIPavec$N_Cas), nMethods), rep(unique(BIPavecsans$BIPsans$N_Cas), nMethods), rep(unique(SKZaff$N_Cas), nMethods), rep(unique(BIPavecsans$BIPavec$N_Cas), nMethods), rep(unique(BIPavecsans$BIPsans$N_Cas), nMethods), rep(unique(SKZaff$N_Cas), nMethods))
out$N_Temoin <- c(rep(unique(SKZBIP$SKZ$N_Temoin), nMethods), rep(unique(SKZBIP$BIP$N_Temoin), nMethods), rep(unique(SKZBIP$SKZ$N_Temoin), nMethods), rep(unique(SKZBIP$BIP$N_Temoin), nMethods), rep(unique(BIPavecsans$BIPavec$N_Temoin), nMethods), rep(unique(BIPavecsans$BIPsans$N_Temoin), nMethods), rep(unique(SKZaff$N_Temoin), nMethods), rep(unique(BIPavecsans$BIPavec$N_Temoin), nMethods), rep(unique(BIPavecsans$BIPsans$N_Temoin), nMethods), rep(unique(SKZaff$N_Temoin), nMethods))
out$PRS_type <- strsplit(x = out$Methode, split = "_", fixed = TRUE) %>% sapply(., tail, n = 1)
out$Methode <- methodsToNames(out$Methode)
out <- out[order(out$PRS_type, out$Phenotype, out$Methode),]

out2 <- out2[,c("Phenotype", "Methode", "Quantile", "OR", "Lower", "Upper", "Pvalue")]
out2$PRS_type <- strsplit(x = out2$Methode, split = "_", fixed = TRUE) %>% sapply(., tail, n = 1)
out2$Methode <- methodsToNames(out2$Methode)
out2 <- out2[order(out2$PRS_type, out2$Phenotype, out2$Methode),]


out3 <- out3[,c("Phenotype", "Methode", "Quantile", "OR", "Lower", "Upper", "Pvalue")]
out3$PRS_type <- strsplit(x = out3$Methode, split = "_", fixed = TRUE) %>% sapply(., tail, n = 1)
out3$Methode <- methodsToNames(out3$Methode)
out3 <- out3[order(out3$PRS_type, out3$Phenotype, out3$Methode),]

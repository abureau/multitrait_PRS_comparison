library(multivariateLassosum)
library(mvtnorm)
library(fdrtool)
library(matrixcalc)
library(Matrix)
library(Rcpp)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)
library(bigsnpr)
library(bigstatsr)
library(xgboost)
library(SummaryLasso) 
library(reticulate)
library(dplyr)

Data <- ".../Data_Cartagene_imputed"
parsed <- parseselect(Data,extract = NULL, exclude = NULL,keep = NULL, remove = NULL,chr = NULL)
nbr_SNP <- parsed$P
set.seed(42)
rows_randomized <- sample(10139)
keep.1 <- c(rep(FALSE, parsed$N))
keep.1[rows_randomized[1:8139]] <- TRUE
keep.2 <- c(rep(FALSE, parsed$N))
keep.2[rows_randomized[8140:9139]] <- TRUE
keep.3 <- c(rep(FALSE, parsed$N))
keep.3[rows_randomized[9140:parsed$N]] <- TRUE
parsed.1 <- parseselect(Data, keep = keep.1)
parsed.2 <- parseselect(Data, keep = keep.2)
parsed.3 <- parseselect(Data, keep = keep.3)
rows_1 <- sort(rows_randomized[1:8139])
rows_2 <- sort(rows_randomized[8140:9139])
rows_3 <- sort(rows_randomized[9140:parsed$N])
h_obs_SKZ <- 0.47
h_obs_BIP <- 0.45
Var_genetique_SKZ <- h_obs_SKZ
Var_genetique_BIP <- h_obs_BIP
Var_environm_SKZ <- 1 - Var_genetique_SKZ
Var_environm_BIP <- 1 - Var_genetique_BIP
Correlation <- 0.59
Covariance <- Correlation * sqrt(Var_genetique_SKZ) * sqrt(Var_genetique_BIP)
Var_SKZ_beta <- Var_genetique_SKZ / nbr_SNP
Var_BIP_beta <- Var_genetique_BIP / nbr_SNP
Covariance_beta <- Covariance / nbr_SNP
Sigma_b <- matrix(data = c(Var_SKZ_beta, Covariance_beta, Covariance_beta, Var_BIP_beta), nrow = 2, ncol = 2)
Sigma_s <- matrix(data = c(Var_environm_SKZ, 0, 0, Var_environm_BIP), nrow = 2, ncol = 2)
inv_Sb <- solve(Sigma_b)
inv_Ss <- solve(Sigma_s)
inv_Sb <- matrix(as.numeric(format(round(inv_Sb, 3), nsmall = 3)),nrow = 2,ncol = 2)
sigma2 <- 2e-6
maximum <- min(Var_SKZ_beta/sigma2, Var_BIP_beta/sigma2)
minimum <- Covariance_beta/sigma2
pi1 <- 0.35
ro <- Covariance_beta / (pi1 * sigma2)
pi2 <- Var_SKZ_beta / sigma2 - pi1
pi3 <- Var_BIP_beta / sigma2 - pi1
pi4 <- 1 - pi1 - pi2 - pi3
prob <- c(pi1, pi2, pi3, pi4)
sd.1 <- lassosum:::sd.bfile(bfile = Data, keep=keep.1)
sd.2 <- lassosum:::sd.bfile(bfile = Data, keep=keep.2)
sd.3 <- lassosum:::sd.bfile(bfile = Data, keep=keep.3)
weight <- 1/sd.3
weight[!is.finite(weight)] <- 0
weight.2 <- 1/sd.2
weight.2[!is.finite(weight.2)] <- 0
Var_Y_SKZ <- 1
sd_Y_SKZ <- 1
Var_Y_BIP <- 1
sd_Y_BIP <- 1

#Heritabilities and covariance
#To avoid numerical problems, we fixe 0 to 1e-9
heritabilite <- readRDS(".../S-LDXR/h_rho_bysnp.RDS")
minadj <- 1e-9
heritabilite$h_SKZ[heritabilite$h_SKZ < minadj] <- minadj
heritabilite$h_BIP[heritabilite$h_BIP < minadj] <- minadj
heritabilite_estime_SKZ <- sum(heritabilite$h_SKZ)
heritabilite_estime_BIP <- sum(heritabilite$h_BIP)
covariance_estime <- sum(heritabilite$rho)
delta_SKZ <- h_obs_SKZ/(heritabilite_estime_SKZ)
delta_BIP <- h_obs_BIP/(heritabilite_estime_BIP)
delta <- sqrt(delta_SKZ*delta_BIP)

#Var-Covar matrices
Sigma_b_snp <- array(data = NA, dim = c(2,2, nbr_SNP))
Verif_sigma <- array(data = NA, dim = c(2,2, nbr_SNP))
nProb <- 0
for (j in 1:nbr_SNP) {
  h2_BIP <- heritabilite[j, "h_BIP"]
  h2_SKZ <- heritabilite[j, "h_SKZ"]
  rho    <- heritabilite[j, "rho"]
  h2_SKZ <- delta_SKZ*h2_SKZ
  h2_BIP <- delta_BIP*h2_BIP
  rho    <- delta*rho
  detCond<- ((h2_BIP*h2_SKZ)-(rho^2)) <= 0
  if(detCond){rho <- sqrt(max(c((h2_BIP*h2_SKZ)-0.001,0))); nProb <- nProb+1}
  mat_j  <- matrix(data = c(h2_SKZ, rho, rho, h2_BIP), ncol = 2)
  Verif_sigma[,,j] <- mat_j
  Sigma_b_snp[,,j]  <- round(solve(mat_j),0)
}

#Constant var-covar matrices
inv_Sb_snp <- array(data = NA, dim = c(2,2, nbr_SNP))
for (j in 1:nbr_SNP) {
  inv_Sb_snp[,,j] <- inv_Sb
}

#Function
f_lambda <- function(beta_SKZ_lambda, beta_BIP_lambda, r_hat,sd,keep_sujets,beta){
  bXy <- r_hat %*% beta 
  if(!is.null(sd)){
    weight <- 1/sd
    weight[!is.finite(weight)] <- 0
    scaled.beta_SKZ <- as.matrix(Diagonal(x=weight) %*% beta_SKZ_lambda)
    scaled.beta_BIP <- as.matrix(Diagonal(x=weight) %*% beta_BIP_lambda)
    pgs_SKZ <- pgs(Data, keep=keep_sujets, weights=scaled.beta_SKZ)
    pgs_BIP <- pgs(Data, keep=keep_sujets, weights=scaled.beta_BIP)
  } else{
    pgs_SKZ <- pgs(Data, keep=keep_sujets, weights=beta_SKZ_lambda)
    pgs_BIP <- pgs(Data, keep=keep_sujets, weights=beta_BIP_lambda)
  }
  if(ncol(pgs_SKZ)>1){
    pred <- matrix(data = NA,nrow = 2*nrow(pgs_SKZ),ncol = ncol(pgs_SKZ))
    for(i in 1:ncol(pgs_SKZ)){
      pred[,i] <- c(rbind(pgs_SKZ[,i],pgs_BIP[,i]))
    }
  }else{
    pgs_SKZ<- as.vector(pgs_SKZ)
    pgs_BIP<- as.vector(pgs_BIP)
    pred <- c(rbind(pgs_SKZ,pgs_BIP))
  }
  pred2 <- scale(pred, scale=F)
  bXXb <- colSums(pred2^2) / nrow(pgs_SKZ)
  result <- as.vector(bXy / sqrt(bXXb))
  return(result)
}

#Compute LD matrix using PLINK for PANPRS in the method reference bfiles.
system(paste0(".../plink --bfile ", Data, " --keep .../keep.2.txt --r2 --ld-window-kb 250 --ld-window-r2 0 --out .../PANPRSLD"))
plinkLDgenome <- data.table::fread(".../PANPRSLD.ld")
colnames(plinkLDgenome) <- c("CHR_A", "BP_A",  "SNP_A", "CHR_B", "BP_B",  "SNP_B", "R")

for (k in 1:20) {
  #---- Simulated data ----
  print(paste0("Simulation ", k))
  setwd(paste0(".../Simulation_", k, "/"))
  set.seed(k)

  Beta <- readRDS("Beta_simules.Rdata")
  Beta_simule_SKZ <- Beta[1, ]
  Beta_simule_BIP <- Beta[2, ]
  Beta0 <- readRDS("Beta0_simules.Rdata")
  Beta0_simule_SKZ <- Beta0[1, ]
  Beta0_simule_BIP <- Beta0[2, ]

  LDblocks <- "EUR.hg19"
  ref.bim <- read.table2(paste0(Data, ".bim"))
  LDblocks <- read.table2(system.file(paste0("data/Berisa.", LDblocks, ".bed"), package = "lassosum"), header = T)
  LDblocks[, 1] <- as.character(sub("chr", "", LDblocks[, 1], ignore.case = T))
  LDblocks <- splitgenome(CHR = ref.bim$V1, POS = ref.bim$V4, ref.CHR = LDblocks[, 1], ref.breaks = LDblocks[, 3])
  Blocks <- parseblocks(LDblocks)
  Blocks$startvec <- Blocks$startvec + 1
  Blocks$endvec <- Blocks$endvec + 1
  pos.1 <- which(parsed.1$keep) - 1
  keepbytes.1 <- floor(pos.1 / 4)
  keepoffset.1 <- pos.1 %% 4 * 2
  r_SKZ <- readRDS(".../r_SKZ_simules.Rdata")
  r_BIP <- readRDS(".../r_BIP_simules.Rdata")
  r_SKZ_destand <- r_SKZ*(1/sd.1)
  r_BIP_destand <- r_BIP*(1/sd.1)
  r <- rbind(r_SKZ,r_BIP)
  r_SKZ.2 <- as.vector(readRDS(".../r_SKZ_simules_jeu2.Rdata"))
  r_BIP.2 <- as.vector(readRDS(".../r_BIP_simules_jeu2.Rdata"))
  r.2 <- c(rbind(r_SKZ,r_BIP))
  #Number of cases and controls used in the summary statitics GWAS.
  #Values obtained using the references of the summary statistics (see Simulation_Standard.R).
  size_SKZ <- sum(33426, 32541)
  size_BIP <- sum(21524, 20129)
  sample_size <- c(size_SKZ,size_BIP)
  cores <- detectCores()
  pvalue_SKZ <- 2*pnorm(abs(r_SKZ_destand/sqrt(1/size_SKZ)), mean = 0, sd = 1, lower.tail = FALSE)
  pvalue_BIP <- 2*pnorm(abs(r_BIP_destand/sqrt(1/size_BIP)), mean = 0, sd = 1, lower.tail = FALSE)

  #---- Thresholding ----
  NCORES <- nb_cores()
  if(!file.exists(".../Data_Cartagene_imputed.rds")){
    snp_readBed2(".../Data_Cartagene_imputed.bed",
                 backingfile = ".../Data_Cartagene_imputed.rds",
                 ncores = NCORES)
  }
  threshold_param <- c(1, 0.75, 0.5, 0.25, 0.1, 0.05, 0.001, 1e-4)
  obj.bigSNP <- snp_attach(".../Data_Cartagene_imputed.rds")

  #Validation
  ##SKZ
  for(Lam in threshold_param){
    add <- r_SKZ_destand
    add[pvalue_SKZ > Lam] <- 0
    if(which(threshold_param == Lam) == 1) {BETA_SKZ <- data.frame(add) } else {BETA_SKZ <- cbind(BETA_SKZ, add)}
  }
  cl <- makeCluster(cores[1]-2) 
  #We use the pseudovalidate function from the lassosum package for validation as we're using correlation from the same set as the one
  #used to modelize lassosum. (date set #2)
  pv_SKZ <- lassosum:::pseudovalidation(Data, beta = BETA_SKZ, cor = r_SKZ.2, keep = parsed.2$keep, sd = sd.2, cluster = cl)
  stopCluster(cl)
  x_SKZ <- as.vector(pv_SKZ)
  names(x_SKZ) <- as.character(threshold_param)
  saveRDS(x_SKZ, file = "Valeurs_f_lambda_thresholding_SKZ.Rdata")
  
  ##BIP
  for(Lam in threshold_param){
    add <- r_BIP_destand
    add[pvalue_BIP > Lam] <- 0
    if(which(threshold_param == Lam) == 1) {BETA_BIP <- data.frame(add) } else {BETA_BIP <- cbind(BETA_BIP, add)}
  }
  cl <- makeCluster(cores[1]-2) 
  #We use the pseudovalidate function from the lassosum package for validation as we're using correlation from the same set as the one
  #used to modelize lassosum. (date set #2)
  pv_BIP <- lassosum:::pseudovalidation(Data, beta = BETA_BIP, cor = r_BIP.2, keep = parsed.2$keep, sd = sd.2, cluster = cl)
  stopCluster(cl)
  x_BIP <- as.vector(pv_BIP)
  names(x_BIP) <- as.character(threshold_param)
  saveRDS(x_BIP, file = "Valeurs_f_lambda_thresholding_BIP.Rdata")
  
  #PRS
  for(Lam in threshold_param){
    # SKZ
    PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=BETA_SKZ[,which(threshold_param == Lam)])
    saveRDS(PGS_estime_SKZ, file = paste0("PGS_thresholding_SKZ_", Lam, ".Rdata"))
    # BIP 
    PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=BETA_BIP[,which(threshold_param == Lam)])
    saveRDS(PGS_estime_BIP, file = paste0("PGS_thresholding_BIP_", Lam, ".Rdata"))
  }

  #---- MultiLassosum STANDARD ----
  AllLambdas <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.5)
  cl <- makeCluster(cores[1]-2)  
  out_lassosumExt <- multivariateLassosum::lassosum(cor = r, inv_Sb = inv_Sb_snp, inv_Ss = inv_Ss, bfile = Data, keep = keep.2,
                                                    lambda = AllLambdas, shrink = 0.5, maxiter = 100000, trace = 1, blocks = LDblocks,
					            sample_size = sample_size, cluster = cl)
  stopCluster(cl)
  saveRDS(out_lassosumExt,file = "results_lassosumExt_s_0.5.Rdata")
  
  #Validation
  NbrLambdas <- length(AllLambdas)
  mat_Beta_SKZ <- out_lassosumExt$beta[,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
  mat_Beta_BIP <- out_lassosumExt$beta[,seq(from = 2, to = NbrLambdas*2, by = 2)]
  BETA <- data.frame()
  for(idx in 1:NbrLambdas){
    if(idx == 1){
      BETA <- data.frame(c(rbind(out_lassosumExt$beta[,1],rbind(out_lassosumExt$beta[,2]))))
    }else{
      add_idx <- c(rbind(out_lassosumExt$beta[,(2*idx)-1],rbind(out_lassosumExt$beta[,2*idx])))
      BETA <- cbind(BETA, add_idx)
    }
  }
  BETA <- as.matrix(BETA)
  x <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = r.2, sd = sd.2, keep_sujets = parsed.2$keep, beta = BETA)
  x_lassosumExt <- x
  names(x) <- paste0("lamdba_", AllLambdas)
  saveRDS(x, file = "Valeurs_f_lambda_LassosumExtension.Rdata")
  
  #PRS
  for(Lam in AllLambdas){
    # SKZ
    scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weight) %*% mat_Beta_SKZ[,which(AllLambdas == Lam)])
    PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
    saveRDS(PGS_estime_SKZ, file = paste0("PGS_estime_SKZ_", Lam, ".Rdata"))
    # BIP 
    scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weight) %*% mat_Beta_BIP[,which(AllLambdas == Lam)])
    PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
    saveRDS(PGS_estime_BIP, file = paste0("PGS_estime_BIP_", Lam, ".Rdata"))
  }

  #---- MultiLassosum ADAPTIVE (BETA STANDARD) ----
  max_x_lassosumExt <- which.max(x_lassosumExt)
  AW <- t(MultiLassosum$beta[,c((2*max_x_lassosumExt)-1, 2*max_x_lassosumExt)])
  AW[1,] <- 1/abs(AW[1,]*weight.2)
  AW[2,] <- 1/abs(AW[2,]*weight.2)
  saveRDS(AW, file =  "BetaMulti/Adaptive_weights.Rdata")

  AllLambdas <- c(0.0000001, 0.0000005,0.000001, 0.000005, 0.00001, 0.00005, 0.0001)
  cl <- makeCluster(cores[1]-2) 
  out_lassosumExtAdapBM <- multivariateLassosum::lassosum(cor = r, inv_Sb = inv_Sb_snp, inv_Ss = inv_Ss, bfile = Data, weights = AW,
                                                          keep = keep.2, lambda = AllLambdas, shrink = 0.5, maxiter = 100000, trace = 1,
                                                          blocks = LDblocks, sample_size = sample_size, cluster = cl )
  stopCluster(cl)
  saveRDS(out_lassosumExtAdapBM, file = "/BetaMulti/results_lassosumExtAdap_s_0.5.Rdata")

  #Validation
  NbrLambdas <- length(AllLambdas)
  mat_Beta_SKZ <- out_lassosumExtAdapBM$beta[,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
  mat_Beta_BIP <- out_lassosumExtAdapBM$beta[,seq(from = 2, to = NbrLambdas*2, by = 2)]
  BETA <- data.frame()
  for(idx in 1:NbrLambdas){
    if(idx == 1){
      BETA <- data.frame(c(rbind(out_lassosumExtAdapBM$beta[,1],rbind(out_lassosumExtAdapBM$beta[,2]))))
    }else{
      add_idx <- c(rbind(out_lassosumExtAdapBM$beta[,(2*idx)-1],rbind(out_lassosumExtAdapBM$beta[,2*idx])))
      BETA <- cbind(BETA, add_idx)
    }
  }
  BETA <- as.matrix(BETA)
  x <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = r.2, sd = sd.2, keep_sujets = parsed.2$keep, beta = BETA)
  x_lassosumExt <- x
  names(x) <- paste0("lamdba_", AllLambdas)
  saveRDS(x, file = "/BetaMulti/Valeurs_f_lambda_LassosumExtension.Rdata")
  
  #PRS
  for(Lam in AllLambdas){
    # SKZ
    scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weight) %*% mat_Beta_SKZ[,which(AllLambdas == Lam)])
    PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
    saveRDS(PGS_estime_SKZ, file = paste0("BetaMulti/PGS_estime_SKZ_", Lam, ".Rdata"))
    # BIP 
    scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weight) %*% mat_Beta_BIP[,which(AllLambdas == Lam)])
    PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
    saveRDS(PGS_estime_BIP, file = paste0("BetaMulti/PGS_estime_BIP_", Lam, ".Rdata"))
  }

  #---- Genetic_cov ----
  cl <- makeCluster(cores[1]-2)
  #Lambdas ? tester.
  AllLambdas <- c(0.000001, 0.00005, 0.00001, 0.0005, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 5, 10)

  out_lassosumGenCov <- multivariateLassosum::lassosum(cor = r, inv_Sb = Sigma_b_snp, inv_Ss = inv_Ss, bfile = Data,
                                          keep = keep.2, lambda = AllLambdas, shrink = 0.5, maxiter = 100000, trace = 1,
                                          blocks = LDblocks, sample_size = sample_size, cluster = cl)
  stopCluster(cl)
  saveRDS(out_lassosumGenCov,file = "GenCov/results_lassosumExt_s_0.5.Rdata")
  
  #Validation
  NbrLambdas <- length(AllLambdas)
  mat_Beta_SKZ <- out_lassosumGenCov$beta[,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
  mat_Beta_BIP <- out_lassosumGenCov$beta[,seq(from = 2, to = NbrLambdas*2, by = 2)]
  BETA <- data.frame()
  for(idx in 1:NbrLambdas){
    if(idx == 1){
      BETA <- data.frame(c(rbind(out_lassosumGenCov$beta[,1],rbind(out_lassosumGenCov$beta[,2]))))
    }else{
      add_idx <- c(rbind(out_lassosumGenCov$beta[,(2*idx)-1],rbind(out_lassosumGenCov$beta[,2*idx])))
      BETA <- cbind(BETA, add_idx)
    }
  }
  BETA <- as.matrix(BETA)
  x <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = r.2, sd = sd.2, keep_sujets = parsed.2$keep, beta = BETA)
  names(x) <- paste0("lamdba_", AllLambdas)
  x_lassosumGenCov <- x
  saveRDS(x, file = "GenCov/Valeurs_f_lambda_LassosumExtension.Rdata")
  
  #PRS
  for(Lam in AllLambdas){
  ## SKZ :
  scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weight) %*% mat_Beta_SKZ[,which(AllLambdas == Lam)])
  PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
  saveRDS(PGS_estime_SKZ, file = paste0("GenCov/PGS_estime_SKZ_", Lam, ".Rdata"))
  ## BIP : 
  scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weight) %*% mat_Beta_BIP[,which(AllLambdas == Lam)])
  PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
  saveRDS(PGS_estime_BIP, file = paste0("GenCov/PGS_estime_BIP_", Lam, ".Rdata"))
  }

  #---- MultiLassosum ADAPTIVE (BETA GENETIC COV) ----
  max_x_lassosumGenCov <- which.max(x_lassosumGenCov)
  AW <- t(out_lassosumGenCov$beta[,c((2*max_x_lassosumGenCov)-1, 2*max_x_lassosumGenCov)])
  AW[1,] <- 1/abs(AW[1,]*weight.2)
  AW[2,] <- 1/abs(AW[2,]*weight.2)
  saveRDS(AW, file =  "BetaGenCov/Adaptive_weights.Rdata")

  AllLambdas <- c(0.0000001, 0.0000005,0.000001, 0.000005, 0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005)
  cl <- makeCluster(cores[1]-2) 
  out_lassosumExtAdapGenCov <- multivariateLassosum::lassosum(cor = r, inv_Sb = Sigma_b_snp, inv_Ss = inv_Ss, bfile = Data, weights = AW,
                                                          keep = keep.2, lambda = AllLambdas, shrink = 0.5, maxiter = 100000, trace = 1,
                                                          blocks = LDblocks, sample_size = sample_size, cluster = cl )
  stopCluster(cl)
  saveRDS(out_lassosumExtAdapGenCov, file = "BetaGenCov/results_lassosumExtAdap_s_0.5.Rdata")

  #Validation
  NbrLambdas <- length(AllLambdas)
  mat_Beta_SKZ <- out_lassosumExtAdapGenCov$beta[,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
  mat_Beta_BIP <- out_lassosumExtAdapGenCov$beta[,seq(from = 2, to = NbrLambdas*2, by = 2)]
  BETA <- data.frame()
  for(idx in 1:NbrLambdas){
    if(idx == 1){
      BETA <- data.frame(c(rbind(out_lassosumExtAdapGenCov$beta[,1],rbind(out_lassosumExtAdapGenCov$beta[,2]))))
    }else{
      add_idx <- c(rbind(out_lassosumExtAdapGenCov$beta[,(2*idx)-1],rbind(out_lassosumExtAdapGenCov$beta[,2*idx])))
      BETA <- cbind(BETA, add_idx)
    }
  }
  BETA <- as.matrix(BETA)
  x <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = r.2, sd = sd.2, keep_sujets = parsed.2$keep, beta = BETA)
  names(x) <- paste0("lamdba_", AllLambdas)
  saveRDS(x, file = paste0("BetaGenCov/Valeurs_f_lambda_LassosumExtension.Rdata"))
  
  #PRS
  for(Lam in AllLambdas){
    # SKZ
    scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weight) %*% mat_Beta_SKZ[,which(AllLambdas == Lam)])
    PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
    saveRDS(PGS_estime_SKZ, file = paste0("BetaGenCov/PGS_estime_SKZ_", Lam, ".Rdata"))
    # BIP 
    scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weight) %*% mat_Beta_BIP[,which(AllLambdas == Lam)])
    PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
    saveRDS(PGS_estime_BIP, file = paste0("BetaGenCov/PGS_estime_BIP_", Lam, ".Rdata"))
  }

  #---- LASSOSUM ----
  AllLambdas <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05)

  #SKZ
  cl <- makeCluster(cores[1]-2) 
  out_SKZ_lassosum <- lassosum::lassosum(cor = r_SKZ, bfile = Data, keep = keep.2 ,lambda = AllLambdas,shrink = 0.5,
                                         maxiter = 100000, trace = 1, blocks = LDblocks, cluster = cl)
  stopCluster(cl)
  saveRDS(out_SKZ_lassosum, file = "OG/out_SKZ_lassosum_s_0.5.Rdata")

  #BIP
  cl <- makeCluster(cores[1]-2) 
  out_BIP_lassosum<- lassosum::lassosum(cor = r_BIP,bfile = Data, keep = keep.2 ,lambda = AllLambdas, shrink = 0.5,
                                        maxiter = 100000, trace = 1, blocks = LDblocks,cluster = cl)
  stopCluster(cl)
  saveRDS(out_BIP_lassosum, file = "OG/out_BIP_lassosum_s_0.5.Rdata")

  #Validation
  ##SKZ
  BETA_SKZ <- out_SKZ_lassosum$beta
  cl <- makeCluster(cores[1]-2) 
  pv_SKZ <- lassosum:::pseudovalidation(Data, beta = BETA_SKZ, cor = r_SKZ.2, keep = parsed.2$keep, sd = sd.2, cluster = cl)
  stopCluster(cl)
  x_SKZ <- as.vector(pv_SKZ)
  names(x_SKZ) <- as.character(AllLambdas)
  saveRDS(x_SKZ, file = "OG/Valeurs_f_lambda_Lassosum_SKZ.Rdata")
  
  ##BIP
  BETA_BIP <- out_BIP_lassosum$beta
  cl <- makeCluster(cores[1]-2) 
  pv_BIP <- lassosum:::pseudovalidation(Data, beta = BETA_BIP, cor = r_BIP.2, keep = parsed.2$keep, sd = sd.2, cluster = cl)
  stopCluster(cl)
  x_BIP <- as.vector(pv_BIP)
  names(x_BIP) <- as.character(AllLambdas)
  saveRDS(x_BIP, file = "OG/Valeurs_f_lambda_Lassosum_BIP.Rdata")
  
  #PRS 
  for(Lam in AllLambdas){
  ## SKZ :
  scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weight) %*% out_SKZ_lassosum$beta[,which(AllLambdas == Lam)])
  PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
  saveRDS(PGS_estime_SKZ, file = "OG/PGS_estime_SKZ_", Lam, ".Rdata")
  
  ## BIP : 
  scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weight) %*% out_BIP_lassosum$beta[,which(AllLambdas == Lam)])
  PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
  saveRDS(PGS_estime_BIP, file = "OG/PGS_estime_BIP_", Lam, ".Rdata")
  }

  #---- PANPRS ----
  Zmatrix <- matrix(c(sign(r_SKZ_destand)*abs(qnorm(p=(pvalue_SKZ)/2)), sign(r_BIP_destand)*abs(qnorm(p=(pvalue_BIP)/2))), ncol = 2)
  #Some Z are infinite, we replace them by the maximal value by phenotype
  Zmatrix[,1][!is.finite(Zmatrix[,1])] <- max(Zmatrix[,1][is.finite(Zmatrix[,1])])
  Zmatrix[,2][!is.finite(Zmatrix[,2])] <- max(Zmatrix[,2][is.finite(Zmatrix[,2])])
  row.names(Zmatrix) <- ref.bim$V2
  plinkLDgenome <- as.data.frame(plinkLDgenome)

  #generate the set of tuning parameters on a modified version of gsPEN which doesn't fit the model, but only outputs the tuning matrix.
  initialTuning <- SummaryLasso::gsPEN(summaryZ = Zmatrix, Nvec = c(size_SKZ, size_BIP), plinkLD = plinkLDgenome, numChrs = 22, fit = FALSE)

  for(chr in 1:22){
    system(paste0(".../plink --bfile ", Data, " --chr ", chr ," --keep .../keep.2.txt --r2 --ld-window-kb 250 --ld-window-r2 0 --out .../Simulation_", k, "/PANPRSLD"))
    plinkLD <- data.table::fread(paste0(".../Simulation_", k, "/PANPRSLD.ld"))
    colnames(plinkLD) <- c("CHR_A", "BP_A",  "SNP_A", "CHR_B", "BP_B",  "SNP_B", "R")
    plinkLD <- as.data.frame(plinkLD)
    ZmatrixChr <- Zmatrix[ref.bim$V1 == chr,]
    PANPRSchr <- SummaryLasso::gsPEN(summaryZ = ZmatrixChr, Nvec = c(size_SKZ, size_BIP), plinkLD = plinkLD, NumIter = 1000, breaking = 1, 
        numChrs = 1, ChrIndexBeta = 0, Init_summaryBetas = 0, Zscale = 1, 
        RupperVal = NULL, tuningMatrix = initialTuning, penalty = "mixLOG", outputAll = 0, fit = TRUE)
    tuningChr <- apply(PANPRSchr$tuningMatrix, 1, FUN = function(x){paste0(round(x,4), collapse = "-")})
    
    if(chr==1){
      PANPRS <- PANPRSchr
      PANPRS$Numitervec <- NULL
      tuning <- tuningChr
    }else{
      tuningKeep <- intersect(tuning, tuningChr)
      tuningCritChr <- which(tuningChr %in% tuningKeep)
      tuningCrit <- which(tuning %in% tuningKeep)
      PANPRS$BetaMatrix <- cbind(PANPRS$BetaMatrix[tuningCrit,], PANPRSchr$BetaMatrix[tuningCritChr,])
      PANPRS$tuningMatrix <- PANPRS$tuningMatrix[tuningCrit,]
      tuning <- tuningKeep
    }
  }
  saveRDS(PANPRS, "PANPRS/PANPRS.RDS")

  #Validation
  whereSKZ <- which(stringr::str_detect(dimnames(PANPRS$BetaMatrix)[[2]], "trait1"))
  orderSKZ <- dimnames(PANPRS$BetaMatrix)[[2]][whereSKZ]
  orderSKZ <- substr(orderSKZ, 1, nchar(orderSKZ)-7)
  whereBIP <- which(stringr::str_detect(dimnames(PANPRS$BetaMatrix)[[2]], "trait2"))
  orderBIP <- dimnames(PANPRS$BetaMatrix)[[2]][whereBIP]
  orderBIP <- substr(orderBIP, 1, nchar(orderBIP)-7)
  NbrLambdas <- dim(PANPRS$tuningMatrix)[1]
  mat_Beta_SKZ <- PANPRS$BetaMatrix[,whereSKZ]
  mat_Beta_BIP <- PANPRS$BetaMatrix[,whereBIP]
  BETA <- data.frame()
  for(idx in 1:NbrLambdas){
    if(idx == 1){
      BETA <- data.frame(c(rbind(PANPRS$BetaMatrix[1,whereSKZ],rbind(PANPRS$BetaMatrix[1,whereBIP]))))
    }else{
      add_idx <- c(rbind(PANPRS$BetaMatrix[idx,whereSKZ],rbind(PANPRS$BetaMatrix[idx,whereBIP])))
      BETA <- cbind(BETA, add_idx)
    }
  }
  BETA <- as.matrix(BETA)
  x <- f_lambda(beta_SKZ_lambda = t(mat_Beta_SKZ), beta_BIP_lambda = t(mat_Beta_BIP), r_hat = r.2, sd = NULL, keep_sujets = parsed.2$keep, beta = BETA)
  names(x) <- apply(PANPRS$tuningMatrix, 1, FUN = function(x){paste0(round(x,4), collapse = "-")})
  saveRDS(x, file = "PANPRS/Valeurs_f_lambda_PANPRS.Rdata")
  
  #PRS
  AllLambdas <- names(x)
  AllLambdasMax <- which.max(x)
  AllLambdasLoop <- AllLambdas[order(x, decreasing = TRUE)[1:5]]
  for(Lam in AllLambdasLoop){
  # SKZ :
  scaled.beta_estime_SKZ <- as.matrix(t(mat_Beta_SKZ)[,which(AllLambdas == Lam)])
  PGS_estime_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_SKZ)
  saveRDS(PGS_estime_SKZ, file = paste0("PANPRS/PGS_estime_SKZ_", Lam, ".Rdata"))
  # BIP : 
  scaled.beta_estime_BIP <- as.matrix(t(mat_Beta_BIP)[,which(AllLambdas == Lam)])
  PGS_estime_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled.beta_estime_BIP)
  saveRDS(PGS_estime_BIP, file = paste0("PANPRS/PGS_estime_BIP_", Lam, ".Rdata"))
  }
  
  #---- LDPRED2 ----
  G   <- obj.bigSNP$genotypes
  CHR <- as.integer(obj.bigSNP$map$chromosome)
  POS <- obj.bigSNP$map$physical.pos
  y   <- obj.bigSNP$fam$affection - 1
  
  #Correlation computation
  corr <- snp_cor(G, ind.row = rows_2, ncores = NCORES)
  tmp <- tempfile(tmpdir = "corr")
  corr0 <- as_SFBM(corr, tmp)
  n_eff_SKZ <- 4 / (1/size_SKZ[1] + 1/size_SKZ[2])
  n_eff_BIP <- 4 / (1/size_BIP[1] + 1/size_BIP[2])
  Beta_se_SKZ <- (1/sqrt(size_SKZ))*(1/sd.1)
  Beta_se_BIP <- (1/sqrt(size_BIP))*(1/sd.1)
  df_beta_SKZ <- cbind(r_SKZ_destand, Beta_se_SKZ, n_eff_SKZ)
  names(df_beta_SKZ) <- c("beta","beta_se","n_eff")
  df_beta_BIP <- cbind(r_BIP_destand, Beta_se_BIP, n_eff_BIP)
  names(df_beta_BIP) <- c("beta","beta_se","n_eff")
  
  #LDPRED2 on SKZ
  df_beta_SKZ <- as.data.frame(df_beta_SKZ)
  ldsc_SKZ <- snp_ldsc2(corr, df_beta_SKZ)
  #We use the heritability value that was fixed in the simulation
  h2_est_SKZ <- 0.47
  multi_auto_SKZ <- snp_ldpred2_auto(corr0, df_beta_SKZ, h2_init = h2_est_SKZ, vec_p_init = seq_log(1e-4, 0.5, length.out = 30), ncores = NCORES, sparse = TRUE)
  saveRDS(multi_auto_SKZ, file = ".../LDpred2/multi_auto_SKZ.Rdata")
  beta_auto_SKZ_sparse <- sapply(multi_auto_SKZ, function(auto) auto$beta_est_sparse)
  pred_auto_SKZ <- big_prodMat(G, beta_auto_SKZ_sparse, ind.row = rows_3)
  sc <- apply(pred_auto_SKZ, 2, sd)
  keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
  keep[which(is.na(keep))]<-FALSE
  final_beta_auto_SKZ <- rowMeans(beta_auto_SKZ_sparse[, keep])
  saveRDS(final_beta_auto_SKZ, file = ".../LDpred2/Beta_SKZ_LDpred2.Rdata")
  final_pred_auto_SKZ <- rowMeans(pred_auto_SKZ[, keep])
  saveRDS(final_pred_auto_SKZ, file = ".../LDpred2/PGS_estime_SKZ_LDpred2_prodMat.Rdata")
  
  #LDPRED2 on BIP
  df_beta_BIP <- as.data.frame(df_beta_BIP)
  ldsc_BIP <- snp_ldsc2(corr, df_beta_BIP)
  #We use the heritability value that was fixed in the simulation
  h2_est_BIP <- 0.45
  multi_auto_BIP <- snp_ldpred2_auto(corr0, df_beta_BIP, h2_init = h2_est_BIP, vec_p_init = seq_log(1e-4, 0.5, length.out = 30), ncores = NCORES, sparse = TRUE)
  saveRDS(multi_auto_BIP, file = ".../LDpred2/multi_auto_BIP.Rdata")
  beta_auto_BIP_sparse <- sapply(multi_auto_BIP, function(auto) auto$beta_est_sparse)
  pred_auto_BIP <- big_prodMat(G, beta_auto_BIP_sparse, ind.row = rows_3)
  sc <- apply(pred_auto_BIP, 2, sd)
  keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
  keep[which(is.na(keep))]<-FALSE
  final_beta_auto_BIP <- rowMeans(beta_auto_BIP_sparse[, keep])
  saveRDS(final_beta_auto_BIP, file = ".../LDpred2/Beta_BIP_LDpred2.Rdata")
  Beta_BIP_LDpred2 <- final_beta_auto_BIP
  final_pred_auto_BIP <- rowMeans(pred_auto_BIP[, keep])
  saveRDS(final_pred_auto_BIP, file = ".../LDpred2/PGS_estime_BIP_LDpred2_prodMat.Rdata") 
  
}

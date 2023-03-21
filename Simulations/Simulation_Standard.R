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

#CARTaGENE data aren't publicly available.
#This code shows how we generated our simulated data.
#Objects derived from CARTaGENE data are found here:
#https://data.mendeley.com/datasets/jxz9jwssf6/2
#If you wish to use our data, we suggest to download every R objects, 
#from Mendeley, and putting them all in a directory. 
#The codes will be easier to follow.

#Every code suppose that your data are all in the same directory.
#Complete the path to your data
path <- ".../"

#Please specify which kind of scenario you want to simulate.
#This `for` loop will set the parameters needed for each simulation.
#Please enter the simulation type. It needs to be written the same ways as it is in the original paper of this project.
simuType <- "..."
set.seed(42)
Data <- paste0(path, "Data_Cartagene")

if(simuType=="n = 29,330"){
  n.1 <- 23330; n.2 <- 26330; nbr_ind <- 29330
  h_obs_SKZ <- 0.47; h_obs_BIP <- 0.45
  Correlation <- 0.59
  sigma2 <- 2.263e-6
  pi1 <- 0.35
  path <- paste0(path, "29300ind/")
}else if(simuType=="n = 10,139"){
  n.1 <- 8139; n.2 <- 9139; nbr_ind <- 10139
  h_obs_SKZ <- 0.47; h_obs_BIP <- 0.45
  Correlation <- 0.59
  sigma2 <- 2e-6
  pi1 <- 0.35
  path <- paste0(path, "10139ind/")
}else if(simuType=="n = 29,330; Low Polygenicity"){
  n.1 <- 23330; n.2 <- 26330; nbr_ind <- 29330
  h_obs_SKZ <- 0.47; h_obs_BIP <- 0.45
  Correlation <- 0.59
  sigma2 <- 1.06244e-05
  pi1 <- 0.08
  path <- paste0(path, "29300ind_lowPoly/")
}else if(simuType=="n = 29,330; Low Heritability"){
  n.1 <- 23330; n.2 <- 26330; nbr_ind <- 29330
  h_obs_SKZ <- 0.094; h_obs_BIP <- 0.090
  Correlation <- 0.59
  sigma2 <- 4.53e-7
  pi1 <- 0.35
  path <- paste0(path, "29300ind_lowHeri/")
}else if(simuType=="n = 29,330; Moderate Correlation"){
  n.1 <- 23330; n.2 <- 26330; nbr_ind <- 29330
  h_obs_SKZ <- 0.47; h_obs_BIP <- 0.45
  Correlation <- 0.43
  sigma2 <- 2.263e-6
  pi1 <- 0.35
  path <- paste0(path, "29300ind_moderateCor/")
}else{
  warning("Please provide an actual simulation scenario, written as it is in the original paper of this project.")
}

parsed <- parseselect(Data,extract = NULL, exclude = NULL,keep = NULL, remove = NULL,chr = NULL)
nbr_SNP <- parsed$P
rho.overlap <- 0
rows_randomized <- sample(nbr_ind)
keep.1 <- c(rep(FALSE, nbr_ind))
keep.1[rows_randomized[1:n.1]] <- TRUE
keep.2 <- c(rep(FALSE, nbr_ind))
keep.2[rows_randomized[(n.1+1):n.2]] <- TRUE
keep.3 <- c(rep(FALSE, nbr_ind))
keep.3[rows_randomized[(n.2+1):nbr_ind]] <- TRUE
parsed.1 <- parseselect(Data, keep = keep.1)
parsed.2 <- parseselect(Data, keep = keep.2)
parsed.3 <- parseselect(Data, keep = keep.3)
rows_1 <- sort(rows_randomized[1:n.1])
rows_2 <- sort(rows_randomized[(n.1+1):n.2])
rows_3 <- sort(rows_randomized[(n.2+1):nbr_ind])
Var_environm_SKZ <- 1 - h_obs_SKZ; Var_environm_BIP <- 1 - h_obs_BIP
Var_SKZ_beta <- h_obs_SKZ / nbr_SNP; Var_BIP_beta <- h_obs_BIP / nbr_SNP
Covariance <- Correlation * sqrt(h_obs_SKZ) * sqrt(h_obs_BIP)
Covariance_beta <- Covariance / nbr_SNP
Sigma_b <- matrix(data = c(Var_SKZ_beta, Covariance_beta, Covariance_beta, Var_BIP_beta), nrow = 2, ncol = 2)
Sigma_s <- matrix(data = c(Var_environm_SKZ, 0, 0, Var_environm_BIP), nrow = 2, ncol = 2)
inv_Ss <- solve(Sigma_s)
inv_Sb <- matrix(as.numeric(format(round(solve(Sigma_b), 3), nsmall = 3)),nrow = 2,ncol = 2)
ro <- Covariance_beta / (pi1 * sigma2)
pi2 <- Var_SKZ_beta / sigma2 - pi1
pi3 <- Var_BIP_beta / sigma2 - pi1
pi4 <- 1 - pi1 - pi2 - pi3
prob <- c(pi1, pi2, pi3, pi4)
#sd.1 <- lassosum:::sd.bfile(bfile = Data,keep=keep.1)
#Import the ones we generated  from online 
sd.1 <- readRDS(paste0(path, "sd.1.RDS"))
weight <- 1/sd.1
weight[!is.finite(weight)] <- 0
Var_Y_SKZ <- 1
sd_Y_SKZ <- 1
Var_Y_BIP <- 1
sd_Y_BIP <- 1

#Function to simulate phenotypes.
simuPheno <- function(gen_liab, h2, K){
  coeff1 <- sqrt(h2) / stats::sd(gen_liab)
  gen_liab <- (gen_liab * coeff1) - mean(gen_liab * coeff1)
  stopifnot(all.equal(mean(gen_liab), 0))
  stopifnot(all.equal(stats::var(gen_liab), h2))
  env_liab <- stats::rnorm(length(gen_liab), sd = sqrt(1 - h2))
  var_env <- stats::var(env_liab)
  cov_env <- stats::cov(gen_liab, env_liab)
  coeff2 <- (sqrt(cov_env^2 + (1 - h2) * var_env) - cov_env) / var_env
  full_liab <- gen_liab + ((env_liab * coeff2) - mean(env_liab * coeff2))
  stopifnot(all.equal(mean(full_liab), 0))
  stopifnot(all.equal(stats::var(full_liab), 1))

  # possibly make binary outcome using liability threshold model
  pheno <- if (is.null(K)) {
    full_liab
  } else {
    (full_liab > stats::qnorm(K, lower.tail = FALSE)) + 0L
  }
  return(pheno)
}


for (k in 1:20){
  print(paste0("Simulation ", k))
  setwd(paste0(path, "Simulation_", k, "/"))
  set.seed(k)

  #Betas
  Beta <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    c1 <- rmvnorm(1, mean = rep(0, 2), sigma = matrix(data = c(sigma2, sigma2 * ro, sigma2 * ro, sigma2), ncol = 2))
    c2 <- c(rnorm(n = 1, mean = 0, sd = sqrt(sigma2)), 0)
    c3 <- c(0, rnorm(n = 1, mean = 0, sd = sqrt(sigma2)))
    c4 <- c(0, 0)
    vect <- sample(x = c("c1", "c2", "c3", "c4"), size = 1, prob = prob)
    Beta[c(1, 2), j] <- get(vect)
  }
  Beta_simule_SKZ <- Beta[1, ]
  Beta_simule_BIP <- Beta[2, ]
  saveRDS(Beta, file = "Standard/Beta_simules.Rdata")

  #Objects to compute the correlation.
  LDblocks <- "EUR.hg19"
  ref.bim <- read.table2(paste0(Data, ".bim"))
  LDblocks <- read.table2(system.file(paste0(
    "data/Berisa.",
    LDblocks, ".bed"
  ),
  package = "lassosum"
  ), header = T)
  LDblocks[, 1] <- as.character(sub("chr", "", LDblocks[, 1], ignore.case = T))
  LDblocks <- splitgenome(
    CHR = ref.bim$V1,
    POS = ref.bim$V4,
    ref.CHR = LDblocks[, 1],
    ref.breaks = LDblocks[, 3]
  )
  Blocks <- parseblocks(LDblocks)
  Blocks$startvec <- Blocks$startvec + 1
  Blocks$endvec <- Blocks$endvec + 1
  pos.1 <- which(parsed.1$keep) - 1
  keepbytes.1 <- floor(pos.1 / 4)
  keepoffset.1 <- pos.1 %% 4 * 2
  pos.2 <- which(parsed.2$keep) - 1
  keepbytes.2 <- floor(pos.2 / 4)
  keepoffset.2 <- pos.2 %% 4 * 2
  n_SKZ <- sum(33426, 32541)
  n_BIP <- sum(21524, 20129)
  registerDoParallel(cl)
  
  #Correlations in the first data set used as a reference.
  r = foreach (i=1:length(Blocks$startvec), .combine=cbind) %dopar% {
    region <- c(Blocks$startvec[i]:Blocks$endvec[i])
    extract_region <- rep(FALSE,nbr_SNP)
    extract_region[region] <- TRUE
    parsed.ref_region <- multivariateLassosum::parseselect(Data, keep = keep.1, extract = extract_region)
    extract2 <- multivariateLassosum::selectregion(!parsed.ref_region$extract)
    extract2[[1]] <- extract2[[1]] - 1
    genotypeMatrix_region <- multivariateLassosum::genotypeMatrix(
    fileName = paste0(Data, ".bed"), N = parsed.ref_region$N, P = parsed.ref_region$P,
    col_skip_pos = extract2[[1]], col_skip = extract2[[2]],
    keepbytes = keepbytes.1, keepoffset = keepoffset.1, fillmissing = 1
    )
    if (length(region) > 1) {
      R_region <- cor(x = genotypeMatrix_region)
    }
    else {
      R_region <- matrix(data = 1, nrow = 1, ncol = 1)
    }
    Beta_region <- as.matrix(Beta[, region])
    if (rho.overlap>0){
      sigmaB = matrix(c(1/ n_SKZ,rho.overlap/sqrt(n_SKZ*n_BIP),rho.overlap/sqrt(n_SKZ*n_BIP) ,1/ n_BIP),2,2) %x% R_region
      tmp = mvtnorm::rmvnorm(1, mean = c(R_region %*% Beta_region[1, ],R_region %*% Beta_region[2, ]), sigma = sigmaB)
      r_region = rbind(tmp[1:nrow(R_region)],tmp[(nrow(R_region)+1):length(tmp)])
    } else {
      r_SKZ_region <- mvtnorm::rmvnorm(1, mean = R_region %*% Beta_region[1, ], sigma = R_region / n_SKZ)
      r_BIP_region <- mvtnorm::rmvnorm(1, mean = R_region %*% Beta_region[2, ], sigma = R_region / n_BIP)
      r_region = rbind(r_SKZ_region,r_BIP_region)
    }
    r_region
  }
  r_SKZ <- r[1,]
  r_BIP <- r[2,]  
  #Import the ones we generated from online 
  #r_SKZ <- readRDS(paste0("Standard/r_SKZ_simule.RData"))
  #r_BIP <- readRDS(paste0("Standard/r_BIP_simule.RData"))
  saveRDS(r_SKZ, file = "Standard/r_SKZ_simules.Rdata")
  saveRDS(r_BIP, file = "Standard/r_BIP_simules.Rdata")

  #Correlations in the second data set used for validation.
  r = foreach (i=1:length(Blocks$startvec), .combine=cbind) %dopar% {
    region <- c(Blocks$startvec[i]:Blocks$endvec[i])
    extract_region <- rep(FALSE,nbr_SNP)
    extract_region[region] <- TRUE
    parsed.ref_region <- multivariateLassosum::parseselect(Data, keep = keep.2, extract = extract_region)
    extract2 <- multivariateLassosum::selectregion(!parsed.ref_region$extract)
    extract2[[1]] <- extract2[[1]] - 1
    genotypeMatrix_region <- multivariateLassosum::genotypeMatrix(
      fileName = paste0(Data, ".bed"), N = parsed.ref_region$N, P = parsed.ref_region$P,
      col_skip_pos = extract2[[1]], col_skip = extract2[[2]],
      keepbytes = keepbytes.2, keepoffset = keepoffset.2, fillmissing = 1
    )
    if (length(region) > 1) {
      R_region <- cor(x = genotypeMatrix_region)
      R_region[is.na(R_region)] <- 0
    } else {
      R_region <- matrix(data = 1, nrow = 1, ncol = 1)
    }
    Beta_region <- as.matrix(Beta[, region])
    if (rho.overlap>0) {
      sigmaB = matrix(c(1/ n_SKZ,rho.overlap/sqrt(n_SKZ*n_BIP),rho.overlap/sqrt(n_SKZ*n_BIP) ,1/ n_BIP),2,2) %x% R_region
      tmp = mvtnorm::rmvnorm(1, mean = c(R_region %*% Beta_region[1, ],R_region %*% Beta_region[2, ]), sigma = sigmaB)
      r_region = rbind(tmp[1:nrow(R_region)],tmp[(nrow(R_region)+1):length(tmp)])
    } else {
      r_SKZ_region <- mvtnorm::rmvnorm(1, mean = R_region %*% Beta_region[1, ], sigma = R_region / n_SKZ)
      r_BIP_region <- mvtnorm::rmvnorm(1, mean = R_region %*% Beta_region[2, ], sigma = R_region / n_BIP)
      r_region = rbind(r_SKZ_region,r_BIP_region)
    }
    r_region
  }
  r_SKZ <- r[1,]
  r_BIP <- r[2,]
  #Import the ones we generated  from online 
  #r_SKZ <- readRDS(paste0("Standard/Simulation_",k , "/r_SKZ_simule_jeu2.RData"))
  #r_BIP <- readRDS(paste0("Standard/Simulation_",k , "/r_BIP_simule_jeu2.RData"))
  saveRDS(r_SKZ, file = "Standard/r_SKZ_simules_jeu2.Rdata")
  saveRDS(r_BIP, file = "Standard/r_BIP_simules_jeu2.Rdata")
  stopCluster(cl)

  #Simulated PRSs
  scaled_beta_simule_SKZ <- as.matrix(Diagonal(x=weight) %*% Beta_simule_SKZ)
  scaled_beta_simule_BIP <- as.matrix(Diagonal(x=weight) %*% Beta_simule_BIP)
  PGS_simule_SKZ <- pgs(Data, keep=parsed.3$keep, weights=scaled_beta_simule_SKZ)
  PGS_simule_BIP <- pgs(Data, keep=parsed.3$keep, weights=scaled_beta_simule_BIP)
  #Import the ones we generated  from Mendeley 
  #PGS_simule_SKZ <- readRDS(paste0("Standard/PGS_simule_SKZ.Rdata"))
  #PGS_simule_BIP <- readRDS(paste0("Standard/PGS_simule_BIP.Rdata"))
  saveRDS(PGS_simule_SKZ, file = "Standard/PGS_simule_SKZ.Rdata")
  saveRDS(PGS_simule_BIP, file = "Standard/PGS_simule_BIP.Rdata")

  #Simulated Phenos
  gen_liab_SKZ <- as.vector(scale(PGS_simule_SKZ, center = TRUE, scale = FALSE))
  gen_liab_BIP <- as.vector(scale(PGS_simule_BIP, center = TRUE, scale = FALSE))
  simuPheno_SKZ <- simuPheno(gen_liab_SKZ, h2_SKZ = var(gen_liab_SKZ), K = 0.01)
  simuPheno_BIP <- simuPheno(gen_liab_BIP, h2_BIP = var(gen_liab_BIP), K = 0.02)
  #Import the ones we generated  from Mendeley 
  #simuPheno_SKZ <- readRDS(paste0("Standard/pheno_simule_SKZ.Rdata"))
  #simuPheno_BIP <- readRDS(paste0("Standard/pheno_simule_BIP.Rdata"))
  saveRDS(simuPheno_SKZ, file = "Standard/pheno_simule_SKZ.Rdata")
  saveRDS(simuPheno_BIP, file = "Standard/pheno_simule_BIP.Rdata")

}

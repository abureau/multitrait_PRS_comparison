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
#https://data.mendeley.com/datasets/jxz9jwssf6/1
#Except for simulated PRSs. They are found on this GitHub repository.
#If you wish to use our data, we suggest to download every R objects, 
#from Github and Mendeley, and putting them all in a directory. 
#The codes will be easier to follow.

#Every code suppose that your data are all in the same directory.
#Complete the path to your data
path <- ".../"

Data <- paste0(path, "Data_Cartagene")
parsed <- parseselect(Data,extract = NULL, exclude = NULL,keep = NULL, remove = NULL,chr = NULL)
nbr_SNP <- parsed$P
nbr_ind <- parsed$N
set.seed(42)
rows_randomized <- sample(10139)
keep.1 <- c(rep(FALSE, nbr_ind))
keep.1[rows_randomized[1:8139]] <- TRUE
keep.2 <- c(rep(FALSE, nbr_ind))
keep.2[rows_randomized[8140:9139]] <- TRUE
keep.3 <- c(rep(FALSE, nbr_ind))
keep.3[rows_randomized[9140:nbr_ind]] <- TRUE
parsed.1 <- parseselect(Data, keep = keep.1)
parsed.2 <- parseselect(Data, keep = keep.2)
parsed.3 <- parseselect(Data, keep = keep.3)
rows_1 <- sort(rows_randomized[1:8139])
rows_2 <- sort(rows_randomized[8140:9139])
rows_3 <- sort(rows_randomized[9140:nbr_ind])
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
matrixcalc::is.positive.semi.definite(Sigma_b, tol = 1e-8)
Sigma_s <- matrix(data = c(Var_environm_SKZ, 0, 0, Var_environm_BIP), nrow = 2, ncol = 2)
matrixcalc::is.positive.semi.definite(Sigma_s, tol = 1e-8)
inv_Sb <- solve(Sigma_b)
inv_Ss <- solve(Sigma_s)
inv_Sb <- matrix(as.numeric(format(round(inv_Sb, 3), nsmall = 3)),nrow = 2,ncol = 2)
matrixcalc::is.positive.semi.definite(inv_Sb, tol = 1e-8)
sigma2 <- 2e-6
maximum <- min(Var_SKZ_beta/sigma2, Var_BIP_beta/sigma2)
minimum <- Covariance_beta/sigma2
pi1 <- 0.35
ro <- Covariance_beta / (pi1 * sigma2)
pi2 <- Var_SKZ_beta / sigma2 - pi1
pi3 <- Var_BIP_beta / sigma2 - pi1
pi4 <- 1 - pi1 - pi2 - pi3
prob <- c(pi1, pi2, pi3, pi4)
#sd.1 <- lassosum:::sd.bfile(bfile = Data,keep=keep.1)
#sd.2 <- lassosum:::sd.bfile(bfile = Data,keep=keep.2)
#Import the ones we generated  from online 
sd.1 <- readRDS("sd.1.RDS")
sd.2 <- readRDS("sd.2.RDS")
Var_Y_SKZ <- 1
sd_Y_SKZ <- 1
Var_Y_BIP <- 1
sd_Y_BIP <- 1

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
  saveRDS(Beta, file = "Beta_simules.Rdata")

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
  n_SKZ <- sum(33426, 32541)
  n_BIP <- sum(21524, 20129)
  cores=detectCores()
  print("nombre de cores")
  cores
  cl <- makeCluster(cores[1]-10)
  registerDoParallel(cl)
  
  #Standardized betas in the first data set
  Beta0 <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta0) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    Beta0[1, j] <- Beta[1, j] * (sd.1[j] / sd_Y_SKZ)
    Beta0[2, j] <- Beta[2, j] * (sd.1[j] / sd_Y_BIP)
  }
  Beta0_simule_SKZ <- Beta0[1, ]
  Beta0_simule_BIP <- Beta0[2, ]
  saveRDS(Beta0, file = "Beta0_simules.Rdata")

  #Correlations in the first data set used as a reference.
  r_SKZ <- c()
  r_BIP <- c()
  for (i in 1:length(Blocks$startvec)) {
    region <- c(Blocks$startvec[i]:Blocks$endvec[i])
    extract_region <- rep(FALSE,nbr_SNP)
    extract_region[region] <- TRUE
    parsed.ref_region <- parseselect(Data, keep = keep.1, extract = extract_region)
    extract2 <- selectregion(!parsed.ref_region$extract)
    extract2[[1]] <- extract2[[1]] - 1
    genotypeMatrix_region <- genotypeMatrix(
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
    Beta0_region <- as.matrix(Beta0[, region])
    r_SKZ_region <- rmvnorm(1, mean = R_region %*% Beta0_region[1, ], sigma = R_region / n_SKZ)
    r_SKZ <- append(x = r_SKZ, values = r_SKZ_region)
    r_BIP_region <- rmvnorm(1, mean = R_region %*% Beta0_region[2, ], sigma = R_region / n_BIP)
    r_BIP <- append(x = r_BIP, values = r_BIP_region)
  }
  r <- rbind(r_SKZ,r_BIP)
  r_SKZ <- r[1,]
  r_BIP <- r[2,]  
  #Import the ones we generated from online 
  #r_SKZ <- readRDS(paste0("Standard/Simulation_",k , "/r_SKZ_simule.RData"))
  #r_BIP <- readRDS(paste0("Standard/Simulation_",k , "/r_BIP_simule.RData"))
  saveRDS(r_SKZ, file = "r_SKZ_simules.Rdata")
  saveRDS(r_BIP, file = "r_BIP_simules.Rdata")

  #Standardized betas in the second data set
  Beta0 <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta0) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    Beta0[1, j] <- Beta[1, j] * (sd.2[j] / sd_Y_SKZ)
    Beta0[2, j] <- Beta[2, j] * (sd.2[j] / sd_Y_BIP)
  }
  Beta0_simule_SKZ <- Beta0[1, ]
  Beta0_simule_BIP <- Beta0[2, ]
  saveRDS(Beta0, file = "Beta0_simules_jeu2.Rdata")
  
  #Correlations in the second data set used for validation.
  r_SKZ <- c()
  r_BIP <- c()
  for (i in 1:length(Blocks$startvec)) {
    region <- c(Blocks$startvec[i]:Blocks$endvec[i])
    extract_region <- rep(FALSE,nbr_SNP)
    extract_region[region] <- TRUE
    parsed.ref_region <- parseselect(Data, keep = keep.2, extract = extract_region)
    extract2 <- selectregion(!parsed.ref_region$extract)
    extract2[[1]] <- extract2[[1]] - 1
    genotypeMatrix_region <- genotypeMatrix(
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
    Beta0_region <- as.matrix(Beta0[, region])
    r_SKZ_region <- rmvnorm(1, mean = R_region %*% Beta0_region[1, ], sigma = R_region / n_SKZ)
    r_SKZ <- append(x = r_SKZ, values = r_SKZ_region)
    r_BIP_region <- rmvnorm(1, mean = R_region %*% Beta0_region[2, ], sigma = R_region / n_BIP)
    r_BIP <- append(x = r_BIP, values = r_BIP_region)
  }
  r <- rbind(r_SKZ,r_BIP)
  r_SKZ <- r[1,]
  r_BIP <- r[2,]
  #Import the ones we generated  from online 
  #r_SKZ <- readRDS(paste0("Standard/Simulation_",k , "/r_SKZ_simule_jeu2.RData"))
  #r_BIP <- readRDS(paste0("Standard/Simulation_",k , "/r_BIP_simule_jeu2.RData"))
  saveRDS(r_SKZ, file = "GenCovNewSim/r_SKZ_simules_jeu2.Rdata")
  saveRDS(r_BIP, file = "GenCovNewSim/r_BIP_simules_jeu2.Rdata")
  stopCluster(cl)

  #Simulated PRSs
  PGS_simule_SKZ <- pgs(Data, keep=parsed.3$keep, weights=Beta_simule_SKZ)
  PGS_simule_BIP <- pgs(Data, keep=parsed.3$keep, weights=Beta_simule_BIP)
  #Import the ones we generated  from GitHub 
  #PGS_simule_SKZ <- readRDS(paste0("Standard/PGS_simule",k , "_SKZ.Rdata"))
  #PGS_simule_BIP <- readRDS(paste0("Standard/PGS_simule",k , "_BIP.Rdata"))
  saveRDS(PGS_simule_SKZ, file = "PGS_simule_SKZ.Rdata")
  saveRDS(PGS_simule_BIP, file = "PGS_simule_BIP.Rdata")
}

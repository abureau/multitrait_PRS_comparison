library(data.table)
library(dplyr)
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

#---- Summary statistics SNPs matching ----
#Trait 1 - Schizophrenia
#Data are available here : https://pgc.unc.edu/for-researchers/download-results/
ss_SKZ <- fread('PGC3_SCZ_wave3_public.v2.tsv',fill = TRUE)
#There's a problem at line 5328756. Column were not correctly delimited. This is how we corrected it.
ss_SKZ[5328756,1] <- 21 ;ss_SKZ[5328756,2] <- "rs148878475";ss_SKZ[5328756,3]<-9648204;
ss_SKZ[5328756,4] <- "C" ; ss_SKZ[5328756,5] <- "T" ; ss_SKZ[5328756,6] <- 0.9878 ;
ss_SKZ[5328756,7] <- 0.9845 ;ss_SKZ[5328756,8] <- 0.1375 ;ss_SKZ[5328756,9] <- 7.5732 ;
ss_SKZ[5328756,10] <- 1.0014 ;ss_SKZ[5328756,11] <- 0.0432 ;ss_SKZ[5328756,12] <- 0 ;
ss_SKZ$CHR <- as.numeric(ss_SKZ$CHR)

#Trait 2 - Bipolar disorder
#Data are available here : https://pgc.unc.edu/for-researchers/download-results/
ss_BIP <- fread('pgc-bip2021-all.vcf.tsv')
#There is no problem here.
names(ss_BIP)[names(ss_BIP) == '#CHROM'] <- "CHR"
ss_BIP$CHR <- as.numeric(ss_BIP$CHR)

#Matching the two sets of summary statistics
#It was already done in 1000GenomesDataPrep.R, but for those who doesn't use 1000 Genomes, here's the code again.
matchSS <- lassosum:::matchpos(tomatch = ss_SKZ, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "CHR", pos = "BP", ref.pos = "POS", ref = "A1", ref.ref = "A1", alt = "A2", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
ss_SKZ <- ss_SKZ[matchSS$order,]
ss_BIP <- ss_BIP[matchSS$ref.extract,]

#Matching of the reference data set and the two sets of summary statistics. Alleles matching isn't necessary here.
#Reference data from CARTaGENE are not available publicly.
Data <- paste0(path, "Data_Cartagene")
bim <- data.table::fread(paste0(Data, ".bim"))
matchTest <- lassosum:::matchpos(tomatch = ss_BIP, ref.df = bim, chr = "CHR", ref.chr = "V1", pos = "POS", ref.pos = "V4", auto.detect.tomatch = F, auto.detect.ref = F, rm.duplicates = T)
ss_BIP <- ss_BIP[matchTest$order,]
ss_SKZ <- ss_SKZ[matchTest$order,]

#---- Formatting the data for S-LDXR ----
ss_SKZ$BETA <- log(ss_SKZ$OR)
ss_SKZ <- ss_SKZ[,c("SNP", "CHR", "BP", "A1", "A2", "BETA", "P", "Nco", "Nca")]
ss_SKZ$N <- rowSums(ss_SKZ[,c("Nco", "Nca")])
ss_SKZ$Z <- sign(ss_SKZ$BETA)*abs(qnorm(p=(ss_SKZ$P)/2))
ss_SKZ <- ss_SKZ[,c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")]
ss_SKZ <- ss_SKZ[order(S$CHR, S$BP),]
data.table::fwrite(ss_SKZ, paste0(path, "S-LDXR/ss_SKZ.txt"), col.names = TRUE, sep = "\t", row.names = FALSE)

ss_BIP <- ss_BIP[,c(3,1,2,4,5,6,8,15,14)]
colnames(ss_BIP) <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "PVAL", "Nco", "Nca")
ss_BIP$N <- rowSums(ss_BIP[,c("Nco", "Nca")])
ss_BIP$Z <- sign(ss_BIP$BETA)*abs(qnorm(p=(ss_BIP$PVAL)/2))
ss_BIP <- ss_BIP[,c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")]
ss_BIP <- ss_BIP[order(ss_BIP$CHR, ss_BIP$BP),]
data.table::fwrite(ss_BIP, paste0(path, "S-LDXR/ss_BIP.txt"), col.names = TRUE, sep = "\t", row.names = FALSE)
system(paste0("gzip ", path, "S-LDXR/ss_SKZ.txt"))
system(paste0("gzip ", path, "S-LDXR/ss_BIP.txt"))

#Run S-LDXR following this code:
#Simulation_B-LDX.sh
#once it is done, continue following this code.

#---- Heritabilities and covariance computation ----
#Import S-LDXR output.
estimates <- fread(paste0(path, "S-LDXR/out.txt"))
#We want the continuous annotations first, then the binary ones.
estimates <- estimates[c(54,55,56,57,59,61,62,60, 1:53,58),]
bim <- bim %>% setNames(., c("CHR", "SNP", "CM", "POS", "A1", "A2"))

#For every chromosome, we initialize an output by trait
h <- data.frame()
annot_ajout <- 0
for(i in 1:22){
  print(paste0("Chromosome ", i))
  data_i <- bim[bim$CHR == i,]
  comp_i <- data_i
  #Import the annotations of this chromosome
  annotations <- fread(paste0(path, "S-LDXR/annotations/baseline-LD-X.", i, ".annot.gz")) %>% as.data.frame()
  annotations <- distinct(annotations)
  #Be sure that the annotations are in the same order as they are in the files where thetas are estimated.
  annotations <- annotations[, c("CHR", "BP", "SNP", "CM", estimates$ANNOT)]
  
  #We add the missing SNPs to the annotation file.
  #For the continuous annotation, we impute using the mean value.
  match_avant_i <- lassosum:::matchpos(tomatch = data_i, ref.df = annotations,
                                 chr = "CHR", ref.chr = "CHR", pos = "POS",
                                 ref.pos = "BP",# snp = "SNP", ref.snp = "SNP",
		                 silent = TRUE, auto.detect.tomatch = F, auto.detect.ref = F, rm.duplicates = T)
  add_i <- data_i[setdiff(1:nrow(data_i), match_avant_i$order),]
  annotations_i <- annotations[match_avant_i$ref.extract,] %>% as.data.frame()
  data_i <- data_i[match_avant_i$order,]
  data_i <- rbind(data_i, add_i)
  annot_ajout <- annot_ajout + nrow(add_i)
  #Some SNPs are duplicated because they are triallelic, we keep only one copy.
  add2_i <- setdiff(comp_i$SNP, data_i$SNP)
  if(!length(add2_i) == 0 ){
    add2_i <- setdiff(comp_i$SNP, data_comp_i$SNP)
    add2_i <- data_i[data_i$SNP %in% add2_i,]
    annot_ajout <- annot_ajout + nrow(add2_i)
    add_i <- rbind(add_i, add2_i)
    data_i <- rbind(data_i, comp_i[comp_i$SNP %in% add2_i$SNP])
  }
  add_i <- add_i[, c("CHR", "POS", "SNP")]
  add_i <- data.frame(add_i, rep(0, nrow(add_i)),
	         rep(mean(annotations[,5], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,6], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,7], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,8], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,9], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,10], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,11], na.rm = TRUE), nrow(add_i)),
	         rep(mean(annotations[,12], na.rm = TRUE), nrow(add_i)),
	         rep(1,nrow(add_i)),
	         matrix(rep(0,nrow(add_i)*53), ncol = 53))
  colnames(add_i) <- colnames(annotations)
  annotations_i <- rbind(annotations_i, add_i)
  check <- data.frame(annotations_i[,2], data_i[,4], annotations_i[,2]==data_i[,4])
  print(paste0("The dimensions are the same : ", nrow(data_i) == nrow(comp_i)))
  print(paste0("SNPs are in the same order : ", all(check[,3])))

  h_SKZ_i <- list()
  h_BIP_i <- list()
  h_trans_i <- list()
  #We keep every annotations, continuous and binary ones.
  annotations_i <- annotations_i[, 5:66]
  estimates_i <- estimates[, c("ANNOT", "TAU2", "TAU1", "THETA")]
  #we compute heritabilities.
  h_SKZ_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$TAU2)
  h_BIP_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$TAU1)
  h_trans_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$THETA)
  
  h <- rbind(h, data.frame(data_i[,1], data_i[,4], data_i[,2], h_SKZ_i[,1], h_BIP_i[,1], h_trans_i[,1]))
}

names(h) <- c("CHR", "BP", "SNP", "h_SKZ", "h_BIP", "rho")
h <- h[c("h_SKZ", "h_BIP", "rho", "CHR", "BP", "SNP")]

#Go back to the original SNP order, which is necessary for the simulations.
h <- h[match(data$SNP, h$SNP),]
check <- data.frame(data$SNP, h$SNP, data$SNP == h$SNP)
print("The order is ok ", all(check[,3]))
saveRDS(h, paste0(path, "S-LDXR/h_rho_bysnp.RDS"))

#---- Create folders to save objects ----
for(i in 1:20){
  system(paste0("mkdir ", path, "Simulation_", i, "/GenCov"))
}

#---- Modify the simulations ----
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
#Import the ones we generated  from Mendeley  
sd.1 <- readRDS("sd.1.RDS")
sd.2 <- readRDS("sd.2.RDS")
Var_Y_SKZ <- 1
sd_Y_SKZ <- 1
Var_Y_BIP <- 1
sd_Y_BIP <- 1

# Loop over the simulation replicates
for(k in 1:20){
  print(paste0("Simulation : ", k))
  setwd(paste0(path, "Simulation_", k, "/"))
  set.seed(k)

  #---- Simulate the new betas ----
  #Import the estimated heritabilities and covariances. SNPs are in the good order.
  heritabilite <- readRDS(paste0(path, "S-LDXR/h_rho_bysnp.RDS"))

  #Sample which SNP "causes" which trait.
  pot_BOTH <- which(heritabilite$h_BIP>0 & heritabilite$h_SKZ>0)
  BOTH <- sample(pot_BOTH, pi1*nrow(heritabilite))
  pot_SKZ <- c(setdiff(pot_BOTH, BOTH), which(heritabilite$h_BIP<0 & heritabilite$h_SKZ>0))
  SKZ <- sample(pot_SKZ, pi2*nrow(heritabilite))
  pot_BIP <- c(setdiff(setdiff(pot_BOTH, BOTH),SKZ), which(heritabilite$h_BIP>0 & heritabilite$h_SKZ<0))
  BIP <- sample(pot_BIP, pi3*nrow(heritabilite))
  NONE <- setdiff(1:nrow(heritabilite), c(BOTH, SKZ, BIP))
  heritabilite$idx <- 1:nrow(heritabilite)
  heritabilite$Effect <- ""
  heritabilite$Effect <- ifelse(heritabilite$idx %in% BOTH, yes = "BOTH", no = heritabilite$Effect)
  heritabilite$Effect <- ifelse(heritabilite$idx %in% SKZ, yes = "SKZ", no = heritabilite$Effect)
  heritabilite$Effect <- ifelse(heritabilite$idx %in% BIP, yes = "BIP", no = heritabilite$Effect)
  heritabilite$Effect <- ifelse(heritabilite$idx %in% NONE, yes = "NONE", no = heritabilite$Effect)
  #SNPs causing trait 1, then h2 is fixed to 0 and rho to 0
  heritabilite$h_BIP[heritabilite$Effect == "SKZ"] <- 0
  heritabilite$rho[heritabilite$Effect == "SKZ"] <- 0
  #SNPs causing trait 2, then h1 is fixed to 0 and rho to 0
  heritabilite$h_SKZ[heritabilite$Effect == "BIP"] <- 0
  heritabilite$rho[heritabilite$Effect == "BIP"] <- 0
  #SNPs causing nothing, then h1 and h2 is fixed to 0 and rho to 0
  heritabilite$h_BIP[heritabilite$Effect == "NONE"] <- 0
  heritabilite$h_SKZ[heritabilite$Effect == "NONE"] <- 0
  heritabilite$rho[heritabilite$Effect == "NONE"] <- 0

  #Total heritability and correction factors
  heritabilite_estime_SKZ <- sum(heritabilite$h_SKZ)
  heritabilite_estime_BIP <- sum(heritabilite$h_BIP)
  covariance_estime <- sum(heritabilite$rho)
  alpha_SKZ <- h_obs_SKZ/heritabilite_estime_SKZ
  cor_SKZ <- 0.47/0.35
  alpha_BIP <- h_obs_BIP/heritabilite_estime_BIP
  cor_BIP <- 0.45/0.35
  alpha <- sqrt(alpha_SKZ*alpha_BIP)
  cor <- sqrt(cor_SKZ*cor_BIP)

  #Simulate the new betas.
  Beta <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    if(heritabilite$Effect[j] == "BOTH"){
      h2_BIP <- heritabilite$h_BIP[j]
      h2_SKZ <- heritabilite$h_SKZ[j]
      rho <- heritabilite$rho[j]
      h2_SKZ <- cor_SKZ*alpha_SKZ*h2_SKZ
      h2_BIP <- cor_BIP*alpha_BIP*h2_BIP
      rho <- cor*alpha*rho
      #As explained in the paper, when the matrix determinant is not invertible, a correction needs to be done on the rho estimate.
      detCond<- ((h2_BIP*h2_SKZ)-(rho^2)) <= 0
      if(detCond){rho <- sqrt(max(c((h2_BIP*h2_SKZ)-0.001,0)))}
      vect <- rmvnorm(1, mean = rep(0, 2), sigma = matrix(data = c(h2_SKZ, rho, rho, h2_BIP), ncol = 2))
    }else if(heritabilite$Effect[j] == "SKZ"){
      h2_SKZ <- heritabilite$h_SKZ[j]
      vect <- c(rnorm(n = 1, mean = 0, sd = cor_SKZ*alpha_SKZ*h2_SKZ), 0)
    }else if(heritabilite$Effect[j] == "BIP"){
      h2_BIP <- heritabilite$h_BIP[j]
      vect <- c(0, rnorm(n = 1, mean = 0, sd = cor_BIP*alpha_BIP*h2_BIP))
    }else if(heritabilite$Effect[j] == "NONE"){
      vect <- c(0, 0)
    }
    Beta[c(1, 2), j] <- vect
  }
  Beta_simule_SKZ <- Beta[1, ]
  Beta_simule_BIP <- Beta[2, ]
  saveRDS(Beta, file = "GenCov/Beta_simules.Rdata")

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
  cores <- detectCores()
  cl <- makeCluster(cores[1]-10)
  registerDoParallel(cl)
  
  #Standardized betas the first data set
  Beta0 <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta0) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    Beta0[1, j] <- Beta[1, j] * (sd.1[j] / sd_Y_SKZ)
    Beta0[2, j] <- Beta[2, j] * (sd.1[j] / sd_Y_BIP)
  }
  Beta0_simule_SKZ <- Beta0[1, ]
  Beta0_simule_BIP <- Beta0[2, ]
  saveRDS(Beta0, file = "GenCov/Beta0_simules.Rdata")

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
    }else {
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
  #Import the ones we generated  from Mendeley  
  #r_SKZ <- readRDS(paste0("B-LDX/Simulation_",k , "/r_SKZ_simule.RData"))
  #r_BIP <- readRDS(paste0("B-LDX/Simulation_",k , "/r_BIP_simule.RData"))
  saveRDS(r_SKZ, file = "GenCov/r_SKZ_simules.Rdata")
  saveRDS(r_BIP, file = "GenCov/r_BIP_simules.Rdata")

  #Standardized betas the second data set
  Beta0 <- matrix(data = NA, nrow = 2, ncol = nbr_SNP)
  rownames(Beta0) <- c("SKZ", "BIP")
  for (j in 1:nbr_SNP) {
    Beta0[1, j] <- Beta[1, j] * (sd.2[j] / sd_Y_SKZ)
    Beta0[2, j] <- Beta[2, j] * (sd.2[j] / sd_Y_BIP)
  }
  Beta0_simule_SKZ <- Beta0[1, ]
  Beta0_simule_BIP <- Beta0[2, ]
  saveRDS(Beta0, file = "GenCov/Beta0_simules_jeu2.Rdata")
  
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
  #Import the ones we generated  from Mendeley  
  #r_SKZ <- readRDS(paste0("B-LDX/Simulation_",k , "/r_SKZ_simule_jeu2.RData"))
  #r_BIP <- readRDS(paste0("B-LDX/Simulation_",k , "/r_BIP_simule_jeu2.RData"))
  saveRDS(r_SKZ, file = "GenCovNewSim/r_SKZ_simules_jeu2.Rdata")
  saveRDS(r_BIP, file = "GenCovNewSim/r_BIP_simules_jeu2.Rdata")
  stopCluster(cl)

  #Simulated PRSs
  PGS_simule_SKZ <- pgs(Data, keep=parsed.3$keep, weights=Beta_simule_SKZ)
  PGS_simule_BIP <- pgs(Data, keep=parsed.3$keep, weights=Beta_simule_BIP)
  #Import the ones we generated  from GitHub 
  #PGS_simule_SKZ <- readRDS(paste0("B-LDX/PGS_simule",k , "_SKZ.Rdata"))
  #PGS_simule_BIP <- readRDS(paste0("B-LDX/PGS_simule",k , "_BIP.Rdata"))
  saveRDS(PGS_simule_SKZ, file = "GenCov/PGS_simule_SKZ.Rdata")
  saveRDS(PGS_simule_BIP, file = "GenCov/PGS_simule_BIP.Rdata")
}

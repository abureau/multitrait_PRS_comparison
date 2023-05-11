library(multivariateLassosum)
library(Rcpp)
library(data.table)
library(matrixcalc)
library(Matrix)
library(matlib)
library(parallel)
library(bigsnpr)
library(dplyr)
library(purrr)
library(doParallel)

#Please enter the path to the directory where your 1000 genomes data are found:
path <- ".../"
setwd(paste0(path, "RealData"))

#Please enter the path to your plink software.
pathPlink <- ".../plink"

#---- MAF control ----
#Control for MAF<0.1. 3 639 921 SNPs remaining.
#We are using the 1000 Genomes data created in 1000GenomesSumstatsDataPrep.R
system(paste0(pathPlink, " --bfile ", path, "allchrs4 --maf 0.1 --make-bed --out ", path, "Data/methodRef"))

#---- Pseudo summary statistics ----
#Reference panel (1000 Genome)
ref.bfile <- "Data/methodRef"
Xr <- bigsnpr::bed(paste0(ref.bfile,".bed"))

#Test data 1 and test data 2 were manipulated in order for the SNPs to be in the same order as in the reference panel and the sum. stats.
#Alleles are also in the same order.
#Please refer to code 1000GenomesSumstatsDataPrep.R which explains how to do so.
#Test data 1 (OmniExpress)
Omni <- "Data/methodOmni"
omnibim <- data.table::fread(paste0(Omni, ".bim"))

#Test data 2 (GSA)
GSA <- "Data/methodGSA"
GSAbim <- data.table::fread(paste0(GSA, ".bim"))

#Summary statistics
#We are using the data created in 1000GenomesSumstatsDataPrep.R
#SKZ
ss_SKZ <- fread('ss_SKZ.txt',fill = TRUE)
size_SKZ <- sum(94015, 67390)
ss_SKZ$CHR <- as.numeric(ss_SKZ$CHR)

#BIP
ss_BIP <- fread('ss_BIP.txt')
names(ss_BIP)[names(ss_BIP) == '#CHROM'] <- "CHR"
size_BIP <- sum(371549,41917)
ss_BIP$CHR <- as.numeric(ss_BIP$CHR)

#Correlations
cor_BIP <- p2cor(p = ss_BIP$PVAL, n = size_BIP, sign=ss_BIP$BETA)
cor_SKZ <- p2cor(p = ss_SKZ$P, n = size_SKZ, sign=log(ss_SKZ$OR))

#Let's generate the pseudo sumstats
#The reference panel needs to be centered and scaled
#xi=(mi-2pi)/sqrt(2pi(1-pi)) where pi is the minor allele frequency
#of the ith genetic marker and mi is the ith column vector of the allele count matrix M
sample_size <- c(size_SKZ,size_BIP)
nAverage <- mean(c(size_BIP, size_SKZ))
nB <- 0.1*nAverage
nr <- Xr$nrow
set.seed(101)
g <- rnorm(n = nr, mean = 0, sd = 1)
XrCounts <- bigsnpr::bed_counts(Xr)
XrMAF <- vector(mode = "numeric", length = nr)
XrMeans <- vector(mode = "numeric", length = nr)
XrSds <- vector(mode = "numeric", length = nr)
for(i in 1:dim(XrCounts)[2]){
	XrCount0 <- XrCounts[1,i]
	XrCount1 <- XrCounts[2,i]
	XrCount2 <- XrCounts[3,i]
	XrMAF[i] <- (XrCount1*1+XrCount2*2)/(2*nr)
	XrMeans[i] <- 2*XrMAF[i]
	XrSds[i] <- sqrt(XrMeans[i]*(1-XrMAF[i]))
}
Xtrg <- bigsnpr::bed_cprodVec(obj.bed = Xr, y.row = g, center = XrMeans, scale = XrSds)

#SKZ pseudo sumstats
n_SKZ <- size_SKZ
nA_SKZ <- n_SKZ - nB
r_SKZ <- cor_SKZ
varr_SKZ <- var(r_SKZ)
rA_SKZ <- r_SKZ + sqrt((nB*varr_SKZ)/(nA_SKZ*nr))*Xtrg
rB_SKZ <- (1/nB)*(n_SKZ*r_SKZ - nA_SKZ*rA_SKZ)
ssA_SKZ <- ss_SKZ
ssB_SKZ <- ss_SKZ

#BIP pseudo sumstats
n_BIP <- size_BIP
nA_BIP <- n_BIP - nB
r_BIP <- cor_BIP
varr_BIP <- var(r_BIP)
rA_BIP <- r_BIP + sqrt((nB*varr_BIP)/(nA_BIP*nr))*Xtrg
rB_BIP <- (1/nB)*(n_BIP*r_BIP - nA_BIP*rA_BIP)
ssA_BIP <- ss_BIP 
ssB_BIP <- ss_BIP 

#From correlations to p-values.
cor2p <- function(cor,  n){
	Q <- (abs(cor)*sqrt(n-2))/(sqrt(1-abs(cor)^2))
	pval <- pt(q = Q, df = n-2, lower.tail = F)
	return(pval*2)
}

#Reformating the sumstats data. It's easier when they are all in the same format.
ssA_SKZ$Cor <- rA_SKZ
ssB_SKZ$Cor <- rB_SKZ
ssA_BIP$Cor <- rA_BIP
ssB_BIP$Cor <- rB_BIP
ssA_SKZ$P <- cor2p(rA_SKZ, n = size_SKZ)
ssB_SKZ$P <- cor2p(rB_SKZ, n = size_SKZ)
ssA_BIP$PVAL <- cor2p(rA_BIP, n = size_BIP) 
ssB_BIP$PVAL <- cor2p(rB_BIP, n = size_BIP)
ssA_SKZ$OR <- NA
ssB_SKZ$OR <- NA
ssA_BIP$BETA <- NA 
ssB_BIP$BETA <- NA
ssA_SKZ$SE <- NA 
ssB_SKZ$SE <- NA
ssA_BIP$SE <- NA 
ssB_BIP$SE <- NA

#To lower the size of the objects, let's remove unnecessary variables.
ss_SKZ <- dplyr::select(ss_SKZ, -Direction, -INFO, -ngt, -HetISqt, -HetDf, -HetPVa) 
ssA_SKZ <- dplyr::select(ssA_SKZ, -Direction, -INFO, -ngt, -HetISqt, -HetDf, -HetPVa) 
ssB_SKZ <- dplyr::select(ssB_SKZ, -Direction, -INFO, -ngt, -HetISqt, -HetDf, -HetPVa) 
ss_BIP <- dplyr::select(ss_BIP, -DIRE, -NGT, -IMPINFO, - NEFFDIV2) 
ssA_BIP <- dplyr::select(ssA_BIP, -DIRE, -NGT, -IMPINFO, - NEFFDIV2) 
ssB_BIP <- dplyr::select(ssB_BIP, -DIRE, -NGT, -IMPINFO, - NEFFDIV2) 

#Exporattion of the new pseudo sumstats.
data.table::fwrite(ssA_SKZ, "Data/ssA_SKZ.txt")
data.table::fwrite(ssA_BIP, "Data/ssA_BIP.txt")
data.table::fwrite(ssB_SKZ, "Data/ssB_SKZ.txt")
data.table::fwrite(ssB_BIP, "Data/ssB_BIP.txt")
system("mv ss_SKZ.txt Data/ss_SKZ.txt")
system("mv ss_BIP.txt Data/ss_BIP.txt")

#---- S-LDXR format ----
#Trait 1 (Schizophrenia)
ssA_SKZ <- ssA_SKZ[,c("SNP", "CHR", "BP", "A1", "A2", "Cor", "P", "Nco", "Nca")]
ssA_SKZ$N <- rowSums(ssA_SKZ[,c("Nco", "Nca")])
ssA_SKZ$Z <- sign(ssA_SKZ$Cor)*abs(qnorm(p=(ssA_SKZ$P)/2))
ssA_SKZ <- ssA_SKZ[,c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")]
ssA_SKZ <- ssA_SKZ[order(ssA_SKZ$CHR, ssA_SKZ$BP),]
#Trait 2 (Bipolar disorder)
ssA_BIP <- ssA_BIP[,c("ID", "CHR", "POS", "A1", "A2", "Cor", "PVAL", "NCON", "NCAS")]
ssA_BIP$N <- rowSums(ssA_BIP[,c("NCON", "NCAS")])
ssA_BIP$Z <- sign(ssA_BIP$Cor)*abs(qnorm(p=(ssA_BIP$PVAL)/2))
ssA_BIP <- ssA_BIP[,c("ID", "CHR", "POS", "A1", "A2", "Z", "N")]
colnames(ssA_BIP) <- c("SNP", "CHR", "BP", "A1", "A2", "Z", "N")
ssA_BIP <- ssA_BIP[order(ssA_BIP$CHR, ssA_BIP$BP),]
#Exportation the data
data.table::fwrite(ssA_SKZ, "S-LDXR/ssA_SKZ_SLDXR.txt", col.names = TRUE, sep = "\t", row.names = FALSE)
data.table::fwrite(ssA_BIP, "S-LDXR/ssA_BIP_SLDXR.txt", col.names = TRUE, sep = "\t", row.names = FALSE)
system("gzip S-LDXR/ssA_SKZ_SLDXR.txt")
system("gzip S-LDXR/ssA_BIP_SLDXR.txt")

#Run S-LDXR following this code:
#RealData.sh
#Once it is done, continue following the present code.

#---- Heritabilities and covariance estimates ----
#Import S-LDXR output.
estimates <- fread("S-LDXR/out-ssA-SKZ-BIP.txt")
#We want the continuous annotations first, then the binary ones.
estimates <- estimates[c(54,55,56,57,59,61,62,60, 1:53,58),]
#Let's generate a dataframe containing the SNP information in the order of our sumstats.
#SKZ as much as BIP sumstats could be used.
SNPinfo <- fread("S-LDXR/ssA_SKZ_SLDXR.txt", select = c("SNP", "BP", "CHR", "A1", "A2"))
data <- setNames(SNPinfo, c("SNP", "POS", "CHR", "A1", "A2"))

#For every chromosome, we initialize an output by trait
h <- data.frame()
annot_ajout <- 0
for(i in 1:22){
  print(paste0("Traitement du chromosome ", i, " en cours"))
  data_i <- data[data$CHR == i,]
  comp_i <- data_i
  #Import the annotations of this chromosome.
  annotations <- fread(paste0("/home/jricard/Projet_Meriem/S-LDXR/annotations/baseline-LD-X.", i, ".annot.gz")) %>% as.data.frame()
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
  check <- data.frame(annotations_i[,2], data_i[,2], annotations_i[,2]==data_i[,2])
  print(paste0("Les dimensions sont les memes : ", nrow(data_i) == nrow(comp_i)))
  print(paste0("Les SNPs sont dans le bon ordre : ", all(check[,3])))

  h_SKZ_i <- list()
  h_BIP_i <- list()
  h_trans_i <- list()
  #We keep every annotations, continuous and binary ones.
  annotations_i <- annotations_i[, 5:66]
  estimates_i <- estimates[, c("ANNOT", "TAU2", "TAU1", "THETA")]
  #Ee compute heritabilities.
  h_SKZ_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$TAU2)
  h_BIP_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$TAU1)
  h_trans_i <- as.matrix(annotations_i) %*% as.matrix(estimates_i$THETA)
  
  h <- rbind(h, data.frame(data_i[,3], data_i[,2], data_i[,1], data_i[,4], data_i[,5], h_SKZ_i[,1], h_BIP_i[,1], h_trans_i[,1]))
}

#Export the output in the initial SNPs order.
names(h) <- c("CHR", "BP", "SNP", "A1", "A2", "h_SKZ", "h_BIP", "rho")
h <- h[c("h_SKZ", "h_BIP", "rho", "CHR", "BP", "SNP", "A1", "A2")]
h <- h[match(data$SNP, h$SNP),]
saveRDS(h, "S-LDXR/h_rho_bysnp.RDS")

#---- Methods ----
#Test data OmniExpress.
omnifam <- data.table::fread(paste0(Omni, ".fam"))

#Sumstats SKZ.
corA_SKZ <- ssA_SKZ$Cor
corB_SKZ <- ssB_SKZ$Cor
size_SKZ <- sum(94015, 67390)
ss_SKZ <- ss_SKZ[ss_SKZ$SNP %in% ssA_SKZ$SNP]

#Sumstats BIP.
corA_BIP <- ssA_BIP$Cor
corB_BIP <- ssB_BIP$Cor
size_BIP <- sum(371549, 41917)
ss_BIP <- ss_BIP[ss_BIP$ID %in% ssA_BIP$ID]

#Weights and LD blocks.
LDblocks <- "EUR.hg19" 
sdOmni <- lassosum:::sd.bfile(bfile = omni)
weightOmni <- 1/sdOmni
weightOmni[!is.finite(weightOmni)] <- 0
sdGSA <- lassosum:::sd.bfile(bfile = GSA)
weightGSA <- 1/sdGSA 
weightGSA[!is.finite(weightGSA)] <- 0

#Heritabilities and covariances.
#As mentionned in the paper, to avoid numeric problems, we fix negative and little estimates to 1e-9.
heritabilite <- readRDS("/S-LDXR/h_rho_bysnp.RDS")
heritabilite$h_SKZ[heritabilite$h_SKZ < 1e-9] <- 1e-9
heritabilite$h_BIP[heritabilite$h_BIP < 1e-9] <- 1e-9

#Let's validate the order of the heritability estimates.
#We use ssA_BIP as a reference, but ssA_SKZ, GSAbim or omnibim could be used as the SNPs order is the same.
matchHeritabilite <- lassosum:::matchpos(tomatch = heritabilite, ref.df = ssA_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  					chr = "CHR", ref.chr = "CHR", pos = "BP", ref.pos = "POS",
					ref = "A1", ref.ref = "A1", alt = "A2", ref.alt = "A2",
  					exclude.ambiguous = F, silent = F, rm.duplicates = F)
heritabilite <- heritabilite[matchHeritabilite$order,]

#Var-Covar matrices
nbr_SNP <- nrow(heritabilite)
phenotypic.genetic.Var.Cov.matrix.gencov <- array(data = NA, dim = c(2,2, nbr_SNP))
for (j in 1:nbr_SNP) {
  h2_BIP <- heritabilite[j, "h_BIP"]
  h2_SKZ <- heritabilite[j, "h_SKZ"]
  rho    <- heritabilite[j, "rho"]
  detCond<- ((h2_BIP*h2_SKZ)-(rho^2)) <= 0
  if(detCond){rho <- sqrt(max(c((h2_BIP*h2_SKZ)-0.001,0)))}
  mat_j  <- matrix(data = c(h2_SKZ, rho, rho, h2_BIP), ncol = 2)
  phenotypic.genetic.Var.Cov.matrix.gencov[,,j] <- mat_j
}
totalMatrix <- apply(phenotypic.genetic.Var.Cov.matrix.gencov, 1:2, sum)

#Constant var-covar matrices
phenotypic.genetic.Var.Cov.matrix <- array(data = NA, dim = c(2,2, nbr_SNP))
bySNP <- 1/nbr_SNP
for (j in 1:nbr_SNP) {
  phenotypic.genetic.Var.Cov.matrix[,,j] <- matrix(data = c(bySNP*totalMatrix[1,1], bySNP*totalMatrix[1,2], bySNP*totalMatrix[2,1], bySNP*totalMatrix[2,2]), nrow = 2, ncol = 2)
}

#---- multivariateLassosum STANDARD ----
cl<-makeCluster(ncores[1]-5)
AllLambdas <- exp(seq(from = log(1e-2), to = log(2500), length.out = 20))
outMulti <- lassosum.pipeline(cor = list(corA_SKZ,corA_BIP),
                         phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix,
                         Var.phenotypic = Var.phenotypic, chr = list(ssA_SKZ$CHR,ssA_BIP$CHR),
                         pos = list(ssA_SKZ$BP,ssA_BIP$POS), cluster = cl,
                         A1 =  list(ssA_SKZ$A1,ssA_BIP$A1), A2 = list(ssA_SKZ$A2,ssA_BIP$A2),
                         ref.bfile = omni, #test.bfile = testGSA.file, destandardize = F,
 			 LDblocks = LDblocks, lambda = AllLambdas, s = 0.5, sample_size = c(nA_SKZ, nA_BIP))
stopCluster(cl)
saveRDS(outMulti , file = "Results/multi.RDS")

#Validation
NbrS <- length(outMulti$s)
NbrLambdas <- length(outMulti$lambda)
NbrTrait <- length(outMulti$sumstats)
mat_Beta_SKZ <- outMulti$beta[[NbrS]][,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
mat_Beta_BIP <- outMulti$beta[[NbrS]][,seq(from = 2, to = NbrLambdas*2, by = 2)]
array_Beta <- array(c(mat_Beta_SKZ, mat_Beta_BIP), dim = c(dim(mat_Beta_SKZ)[1], dim(mat_Beta_SKZ)[2], 2))
# We use the pseudovalidation function for validation as we're using correlation (computed using pseudo-summary statistics) instead of phenotypes
xMulti <- multivariateLassosum::pseudovalidation(omni, r = cbind(corB_SKZ, corB_BIP), sd = sdOmni, beta = array_Beta)
names(xMulti) <- paste0("lamdba_", AllLambdas)
maxMulti <- which.max(xMulti)
saveRDS(xMulti, file = "Results/multiFlambda.RDS")

#Final model on real sumstats
cl<-makeCluster(ncores[1]-5)
outMulti <- lassosum.pipeline(cor = list(cor_SKZ,cor_BIP),
                              phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix,
                              Var.phenotypic = Var.phenotypic, chr = list(ss_SKZ$CHR,ss_BIP$CHR),
                              pos = list(ss_SKZ$BP,ss_BIP$POS), cluster = cl,
                              A1 =  list(ss_SKZ$A1,ss_BIP$A1), A2 = list(ss_SKZ$A2,ss_BIP$A2),
                              ref.bfile = testOmni.bfile, LDblocks = LDblocks,
                              lambda = AllLambdas[maxMulti], s = 0.5, sample_size = c(size_SKZ, size_BIP))
stopCluster(cl)
BETA_SKZ <- outMulti$beta$`0.5`[,1]
BETA_BIP <- outMulti$beta$`0.5`[,2]


#PRS (OmniExpress data set)
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightOmni*BETA_SKZ)
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiOmniSKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightOmni*BETA_BIP)
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiOmniBIP.RDS"))

#PRS (GSA data set)
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightGSA*BETA_SKZ)
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGSASKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightGSA*BETA_BIP)
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGSABIP.RDS"))

#---- multivariateLassosum B-LDX ----
cl<-makeCluster(ncores[1]-5)
AllLambdas <- exp(seq(from = log(1e-2), to = log(2500), length.out = 20))
outMultiGenCov <- lassosum.pipeline(cor = list(corA_SKZ,corA_BIP),
                         phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix.gencov,
                         Var.phenotypic = Var.phenotypic, chr = list(ssA_SKZ$CHR,ssA_BIP$CHR),
                         pos = list(ssA_SKZ$BP,ssA_BIP$POS), cluster = cl,
                         A1 =  list(ssA_SKZ$A1,ssA_BIP$A1), A2 = list(ssA_SKZ$A2,ssA_BIP$A2),
                         ref.bfile = omni,# test.bfile = testGSA.file, destandardize = F,
 			 LDblocks = LDblocks, lambda = AllLambdas, s = 0.5, sample_size = c(nA_SKZ, nA_BIP))
stopCluster(cl)
saveRDS(outMultiGenCov, file = "Results/multiGenCov.RDS")

#Validation
NbrLambdas <- length(outMultiGenCov$lambda)
mat_Beta_SKZ <- outMultiGenCov$beta[[1]][,seq(from = 1, to = (NbrLambdas*2)-1, by = 2)]
mat_Beta_BIP <- outMultiGenCov$beta[[1]][,seq(from = 2, to = NbrLambdas*2, by = 2)]
array_Beta <- array(c(mat_Beta_SKZ, mat_Beta_BIP), dim = c(dim(mat_Beta_SKZ)[1], dim(mat_Beta_SKZ)[2], 2))
# We use the pseudovalidation function for validation as we're using correlation (computed using pseudo-summary statistics) instead of phenotypes
xMultiGenCov <- multivariateLassosum::pseudovalidation(omni, r = cbind(corB_SKZ, corB_BIP), sd = sdOmni, beta = array_Beta)
names(xMultiGenCov) <- paste0("lamdba_", AllLambdas)
maxMultiGenCov <- which.max(xMultiGenCov)
saveRDS(xMultiGenCov, file = paste0("Results/multiGenCovFlambda.RDS"))

#Final model on real sumstats
cl<-makeCluster(ncores[1]-5)
outMultiGenCov <- lassosum.pipeline(cor = list(cor_SKZ,cor_BIP),
                                    phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix.gencov,
                                    Var.phenotypic = Var.phenotypic, chr = list(ss_SKZ$CHR,ss_BIP$CHR),
                                    pos = list(ss_SKZ$BP,ss_BIP$POS), cluster = cl,
                                    A1 =  list(ss_SKZ$A1,ss_BIP$A1), A2 = list(ss_SKZ$A2,ss_BIP$A2),
                                    ref.bfile = testOmni.bfile, LDblocks = LDblocks, 
                                    lambda = AllLambdas[maxMultiGenCov], s = 0.5, sample_size = c(size_SKZ, size_BIP))
stopCluster(cl)
BETA_SKZ <- outMultiGenCov$beta$`0.5`[,1]
BETA_BIP <- outMultiGenCov$beta$`0.5`[,2]

#PRS (OmniExpress data set)
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightOmni*BETA_SKZ)
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGenCovOmniSKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightOmni*BETA_BIP)
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGenCovOmniBIP.RDS"))

#PRS (GSA data set)
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightGSA*BETA_SKZ)
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGenCovGSASKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightGSA*BETA_BIP)
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGenCovGSABIP.RDS"))

#---- Lassosum ----
cl <- makeCluster(ncores[1]-5) 
#We are using the lambdas proposed by the authors.
AllLambdas <- exp(seq(log(0.001), log(0.1), length.out = 20))

#Lassosum for SKZ
OGSKZ <- lassosum::lassosum.pipeline(cor = corA_SKZ, ref.bfile = omni, lambda = AllLambdas,
		              chr = ssA_SKZ$CHR, pos = ssA_SKZ$BP, A1 = ssA_SKZ$A1, A2 = ssA_SKZ$A2, s = 0.5,
                              trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)

#Lassosum for BIP
cl <- makeCluster(ncores[1]-5) 
OGBIP <- lassosum::lassosum.pipeline(cor = corA_BIP, ref.bfile = omni, lambda = AllLambdas,
		              chr = ssA_BIP$CHR, pos = ssA_BIP$POS, A1 = ssA_BIP$A1, A2 = ssA_BIP$A2, s = 0.5,
                              trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)

#SKZ Validation
BETA_SKZ <- OGSKZ$beta
cl <- makeCluster(ncores[1]-5) 
# We use the pseudovalidation function for validation as we're using correlation (computed using pseudo-summary statistics) instead of phenotypes
xOGSKZ <- lassosum:::pseudovalidation(bfile = omni, beta = BETA_SKZ[[1]], cor = corB_SKZ, sd = sdOmni, cluster = cl)
stopCluster(cl)
xOGSKZ <- as.vector(xOGSKZ)
names(xOGSKZ) <- as.character(AllLambdas)
maxOGSKZ <- which.max(xOGSKZ)
saveRDS(xOGSKZ, file = "Results/OGSKZFlambda.RDS")

#Final model
cl<-makeCluster(ncores[1]-5)
OGSKZ <- lassosum::lassosum.pipeline(cor = cor_SKZ, ref.bfile = testOmni.bfile, lambda = AllLambdas[maxOGSKZ],
                                     chr = ss_SKZ$CHR, pos = ss_SKZ$BP, A1 = ss_SKZ$A1, A2 = ss_SKZ$A2, s = 0.5,
                                     trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)
saveRDS(OGBIP, file = "Results/OGSKZ.RDS")
BETA_SKZ <- OGSKZ$beta

#BIP Validation
BETA_BIP <- OGBIP$beta
cl <- makeCluster(ncores[1]-5) 
# We use the pseudovalidation function for validation as we're using correlation (computed using pseudo-summary statistics) instead of phenotypes
xOGBIP <- lassosum:::pseudovalidation(bfile = omni, beta = BETA_BIP[[1]], cor = corB_BIP, sd = sdOmni, cluster = cl)
stopCluster(cl)
xOGBIP <- as.vector(xOGBIP)
names(xOGBIP) <- as.character(AllLambdas)
maxOGBIP <- which.max(xOGBIP)
saveRDS(xOGBIP, file = "Results/OGBIPFlambda.RDS")
  
#Final model
cl <- makeCluster(ncores[1]-5) 
OGBIP <- lassosum::lassosum.pipeline(cor = cor_BIP, ref.bfile = testOmni.bfile, lambda = AllLambdas[maxOGBIP],
                                     chr = ss_BIP$CHR, pos = ss_BIP$BP, A1 = ss_BIP$A1, A2 = ss_BIP$A2, s = 0.5,
                                     trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)
saveRDS(OGBIP, file = "Results/OGBIP.RDS")
BETA_BIP <- OGBIP$beta

#PRS Omni
#On destandardise avec le sd du jeu de test respectif
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightOmni*BETA_SKZ)
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/OGOmniSKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightOmni*BETA_BIP)
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP )
saveRDS(PGS_estime_BIP, file = paste0("Results/OGOmniBIP.RDS"))

#On destandardise avec le sd du jeu de test respectif
#PRS GSA
# SKZ
scaled.beta_estime_SKZ <- as.matrix(weightGSA*BETA_SKZ)
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/OGGSASKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(weightGSA*BETA_BIP)
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/OGGSABIP.RDS"))


#---- thresholding and Clumping + thresholding ----
#To avoid parallelism issues
options(bigstatsr.check.parallel.blas = FALSE)

#Create the bigsnpr .rds objects if they are not created yet.
ncores <- nb_cores()-5 
if(!file.exists(paste0("Data/methodGSABigsnprObj.rds"))){
  snp_readBed2(paste0(GSA, ".bed"),
               backingfile = "Data/methodGSABigsnprObj",
               ncores = ncores)
}
if(!file.exists(paste0("Data/methodOmniBigsnprObj.rds"))){
  snp_readBed2(paste0(omni, ".bed"),
               backingfile = "Data/methodOmniBigsnprObj",
               ncores = ncores)
}
stopCluster(cl)

#Import the bigsnpr objects.
obj.bigSNPOmni <- snp_attach("Data/methodOmniBigsnprObj.rds")
GOmni   <- obj.bigSNPOmni$genotypes
CHROmni <- as.integer(obj.bigSNPOmni$map$chromosome)
POSOmni <- obj.bigSNPOmni$map$physical.pos
obj.bigSNPGSA <- snp_attach("Data/methodGSABigsnprObj.rds")
GGSA   <- obj.bigSNPGSA$genotypes
CHRGSA <- as.integer(obj.bigSNPGSA$map$chromosome)
POSGSA <- obj.bigSNPGSA$map$physical.pos

#Let's use the sumstats format that the bigsnpr authors demands.
#Here we use the sumstats A.
ss_SKZCT <- setNames(ssA_SKZ[,c("CHR", "SNP", "BP", "A1", "A2", "beta", "P")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
ss_BIPCT <- setNames(ssA_BIP[,c("CHR", "ID", "POS", "A1", "A2", "beta", "PVAL")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
map <- setNames(obj.bigSNPOmni$map[,-(2:3)], c("chr", "pos", "a0", "a1"))
info_snpSKZ <- snp_match(ss_SKZCT, map)
lpSSKZ <- -log10(info_snpSKZ$p)
info_snpBIP <- snp_match(ss_BIPCT, map)
lpSBIP <- -log10(info_snpBIP$p)


#Thresholding
p <- data.frame(p = c(1e-4, 0.001, seq(from = 0.05, to = 1, length.out = 20)))
array_BetaSKZ <- array(NA, dim = c(nrow(info_snpSKZ), nrow(p), 1))
array_BetaBIP <- array(NA, dim = c(nrow(info_snpBIP), nrow(p), 1))
for(t in p$p){
  ind.keepSKZ <- which(info_snpSKZ$p < t)
  betaSKZ_t <- rep(0, nrow(info_snpSKZ))
  betaSKZ_t[ind.keepSKZ] <- info_snpSKZ$beta[ind.keepSKZ]
  idx <- which(p$p == t)
  array_BetaSKZ[,idx,1] <- betaSKZ_t
  
  ind.keepBIP <- which(info_snpBIP$p < t)
  betaBIP_t <- rep(0, nrow(info_snpBIP))
  betaBIP_t[ind.keepBIP] <- info_snpBIP$beta[ind.keepBIP]
  array_BetaBIP[,idx,1] <- betaBIP_t
}
p$xSKZ <- multivariateLassosum::pseudovalidation(omni, beta = array_BetaSKZ, cor = cbind(corB_SKZ), destandardize = FALSE)
xt.max_SKZ <- which.max(p$xSKZ)
p$xBIP <- multivariateLassosum::pseudovalidation(omni, beta = array_BetaBIP, cor = cbind(corB_BIP), destandardize = FALSE)
xt.max_BIP <- which.max(p$xBIP)
saveRDS(p, file = "Results/T_flambda.RDS")


#Clumping + thresholding
p2 <- c(1e-4, 0.001, seq(from = 0.05, to = 1, length.out = 20))
r2 <- c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.95)
pr2 <- data.frame(r2 = rep(r2, each = length(p2)), p = p2)
array_BetaSKZ <- array(NA, dim = c(nrow(info_snpSKZ), nrow(pr2), 1))
array_BetaBIP <- array(NA, dim = c(nrow(info_snpBIP), nrow(pr2), 1))
for(c in r2){
  ind.keepSKZ_c <- snp_clumping(GOmni, infos.chr = info_snpSKZ$chr, S = lpSSKZ, thr.r2 = c, size = 250, infos.pos = info_snpSKZ$pos, ncores = ncores)
  ind.keepBIP_c <- snp_clumping(GOmni, infos.chr = info_snpBIP$chr, S = lpSBIP, thr.r2 = c, size = 250, infos.pos = info_snpBIP$pos, ncores = ncores)
  for(t in p2){
    ind.keepSKZ_t <- which(info_snpSKZ$p < t)
    ind.keepSKZ <- intersect(ind.keepSKZ_c, ind.keepSKZ_t)
    betaSKZ_ct <- rep(0, nrow(info_snpSKZ))
    betaSKZ_ct[ind.keepSKZ] <- info_snpSKZ$beta[ind.keepSKZ]
    idx <- which(pr2$r2 == c & pr2$p == t)
    array_BetaSKZ[,idx,1] <- betaSKZ_ct
    
    ind.keepBIP_t <- which(info_snpBIP$p < t)
    ind.keepBIP <- intersect(ind.keepBIP_c, ind.keepBIP_t)
    betaBIP_ct <- rep(0, nrow(info_snpBIP))
    betaBIP_ct[ind.keepBIP] <- info_snpBIP$beta[ind.keepBIP]
    array_BetaBIP[,idx,1] <- betaBIP_ct
  }
}
pr2$xSKZ <- multivariateLassosum::pseudovalidation(omni, beta = array_BetaSKZ, cor = cbind(corB_SKZ), destandardize = FALSE)
x.max_SKZ <- which.max(pr2$xSKZ)
pr2$xBIP <- multivariateLassosum::pseudovalidation(omni, beta = array_BetaBIP, cor = cbind(corB_BIP), destandardize = FALSE)
x.max_BIP <- which.max(pr2$xBIP)
saveRDS(pr2, file = "Results/CT_flambda_plus.RDS")

#Here we use the real sumstats to run the final model
ss_SKZCT <- setNames(ss_SKZ[,c("CHR", "SNP", "BP", "A1", "A2", "OR", "P")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
ss_SKZCT$beta <- log(ss_SKZCT$beta)
ss_BIPCT <- setNames(ss_BIP[,c("CHR", "ID", "POS", "A1", "A2", "BETA", "PVAL")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
map <- setNames(obj.bigSNPOmni$map[,-(2:3)], c("chr", "pos", "a0", "a1"))
info_snpSKZ <- snp_match(ss_SKZCT, map)
lpSSKZ <- -log10(info_snpSKZ$p)
info_snpBIP <- snp_match(ss_BIPCT, map)
lpSBIP <- -log10(info_snpBIP$p)

#Thresholding
# PRS (OmniExpress data set)
#SKZ
prs <- snp_PRS(GOmni, betas.keep = info_snpSKZ$beta, lpS.keep = lpSSKZ, thr.list = -log10(p$p[xt.max_SKZ]))
saveRDS(prs, file = paste0("Results/ThresholdingOmniSKZ.RDS"))
#BIP
prs <- snp_PRS(GOmni, betas.keep = info_snpBIP$beta, lpS.keep = lpSBIP, thr.list = -log10(p$p[xt.max_BIP]))
saveRDS(prs, file = paste0("Results/ThresholdingOmniBIP.RDS"))

# PRS (GSA data set)
#SKZ
prs <- snp_PRS(GGSA, betas.keep = info_snpSKZ$beta, lpS.keep = lpSSKZ, thr.list = -log10(p$p[xt.max_SKZ]))
saveRDS(prs, file = paste0("Results/ThresholdingGSASKZ.RDS"))
#BIP
prs <- snp_PRS(GGSA, betas.keep = info_snpBIP$beta, lpS.keep = lpSBIP, thr.list = -log10(p$p[xt.max_BIP]))
saveRDS(prs, file = paste0("Results/ThresholdingGSABIP.RDS"))

#Clumping + thresholding
ind.keepSKZ <- snp_clumping(GOmni, infos.chr = info_snpSKZ$chr, S = lpSSKZ, thr.r2 = pr2$r2[x.max_SKZ], size = 250, infos.pos = info_snpSKZ$pos, ncores = ncores)
ind.keepBIP <- snp_clumping(GOmni, infos.chr = info_snpBIP$chr, S = lpSBIP, thr.r2 = pr2$r2[x.max_BIP], size = 250, infos.pos = info_snpBIP$pos, ncores = ncores)

# PRS (OmniExpress data set)
#SKZ
prs <- snp_PRS(GOmni, betas.keep = info_snpSKZ$beta[ind.keepSKZ], ind.keep = ind.keepSKZ, lpS.keep = lpSSKZ[ind.keepSKZ], thr.list = -log10(pr2$p[x.max_SKZ]))
saveRDS(prs, file = paste0("Results/CTOmniSKZ.RDS"))
#BIP
prs <- snp_PRS(GOmni, betas.keep = info_snpBIP$beta[ind.keepBIP], ind.keep = ind.keepBIP, lpS.keep = lpSBIP[ind.keepBIP], thr.list = -log10(pr2$p[x.max_BIP]))
saveRDS(prs, file = paste0("Results/CTOmniBIP.RDS"))

# PRS (GSA data set)
#SKZ
prs <- snp_PRS(GGSA, betas.keep = info_snpSKZ$beta[ind.keepSKZ], ind.keep = ind.keepSKZ, lpS.keep = lpSSKZ[ind.keepSKZ], thr.list = -log10(pr2$p[x.max_SKZ]))
saveRDS(prs, file = paste0("Results/CTGSASKZ.RDS"))
#BIP
prs <- snp_PRS(GGSA, betas.keep = info_snpBIP$beta[ind.keepBIP], ind.keep = ind.keepBIP, lpS.keep = lpSBIP[ind.keepBIP], thr.list = -log10(pr2$p[x.max_BIP]))
saveRDS(prs, file = paste0("Results/CTGSABIP.RDS"))

#---- LDPRED2 ----
#Correlations 
tmp <- tempfile(tmpdir = "corr")
for (chr in 1:22) {
  print(chr)
  ind.chr <- which(CHROmni == chr)
  corr <- snp_cor(GOmni, ncores = ncores, ind.col = ind.chr, size = 500/2)
  if (chr == 1) {
    ld <- Matrix::colSums(corr^2)
    corr0 <- as_SFBM(corr, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr^2))
    corr0$add_columns(corr, nrow(corr0))
  }
}

#We chose to test new values of polygenicity, as explain in the original paper.
vec_p_init <- c(c(0.0001, 0.0005, 0.001, 0.005, 0.01), seq(0.05, 0.5, length.out = 25))

#let's use the sumstats A in the format that the bigsnpr authors demands
#to run LDpred2-grid-sp.
ssA_SKZ <- ssA_SKZ[, c("SNP", "CHR", "BP", "A1", "A2", "beta", "SE", "P")]
colnames(ssA_SKZ) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
ssA_SKZ$n_eff <- 4 / (1 / sizeA_SKZ[2] + 1 / sizeA_SKZ[1])
ssA_SKZ$beta_se <- ssA_SKZ$beta/(abs(qnorm(p = ssA_SKZ$p/2))*sign(ssA_SKZ$beta))
dfA_beta_SKZ <- snp_match(ssA_SKZ, mapOmni)

ssA_BIP <- ssA_BIP[, c("ID", "CHR", "POS", "A1", "A2", "beta", "SE", "PVAL")]
colnames(ssA_BIP) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
ssA_BIP$n_eff <- 4 / (1 / sizeA_BIP[2] + 1 / sizeA_BIP[1])
ssA_BIP$beta_se <- ssA_BIP$beta/(abs(qnorm(p = ssA_BIP$p/2))*sign(ssA_BIP$beta))
dfA_beta_BIP <- snp_match(ssA_BIP, mapOmni)

#Method
#SKZ
dfA_beta_SKZ <- as.data.frame(dfA_beta_SKZ)[,c("beta","beta_se","n_eff")]
ldsc_SKZ <- snp_ldsc(ld, length(ld), chi2 = (dfA_beta_SKZ$beta / dfA_beta_SKZ$beta_se)^2, sample_size = dfA_beta_SKZ$n_eff, blocks = NULL, ncores = ncores)
h2_est_SKZ <- ldsc_SKZ[[2]]
h2_seq_SKZ <- round(h2_est_SKZ * c(0.3, 0.7, 1, 1.4), 4)
paramsA_SKZ <- expand.grid(p = vec_p_init, h2 = h2_seq_SKZ, sparse = TRUE)
ldpred2_grid_SKZ <- snp_ldpred2_grid(corr0, dfA_beta_SKZ, paramsA_SKZ, ncores = ncores )
saveRDS(ldpred2_grid_SKZ, file = "Results/LDpred2grid_SKZ.RDS")
#Let's keep models using parameters that converge only.
any_nas_SKZ <- apply(ldpred2_grid_SKZ, 2, function(x) any(is.na(x)))
array_Beta <- array(c(ldpred2_grid_SKZ[,!any_nas_SKZ]), dim = c(nrow(ldpred2_grid_SKZ), sum(!any_nas_SKZ), 1))
x <- multivariateLassosum::pseudovalidation(omni, beta = array_Beta, cor = cbind(corB_SKZ), destandardize = FALSE)
x.max_SKZ <- which.max(x)
paramsA_SKZ <- paramsA_SKZ[!any_nas_SKZ,]
saveRDS(x, file = "Results/LDpred2grid_SKZ_flambda.RDS")

#BIP
dfA_beta_BIP <- as.data.frame(dfA_beta_BIP)[,c("beta","beta_se","n_eff")]
ldsc_BIP <- snp_ldsc(ld, length(ld), chi2 = (dfA_beta_BIP$beta / dfA_beta_BIP$beta_se)^2, sample_size = dfA_beta_BIP$n_eff, blocks = NULL, ncores = ncores)
h2_est_BIP <- ldsc_BIP[[2]]
h2_seq_BIP <- round(h2_est_BIP * c(0.3, 0.7, 1, 1.4), 4)
paramsA_BIP <- expand.grid(p = vec_p_init, h2 = h2_seq_BIP, sparse = TRUE)
ldpred2_grid_BIP <- snp_ldpred2_grid(corr0, dfA_beta_BIP, paramsA_BIP, ncores = ncores )
saveRDS(ldpred2_grid_BIP, file = "Results/LDpred2grid_BIP.RDS")
#Let's keep models using parameters that converge only.
any_nas_BIP <- apply(ldpred2_grid_BIP, 2, function(x) any(is.na(x)))
array_Beta <- array(c(ldpred2_grid_BIP[,!any_nas_BIP]), dim = c(nrow(ldpred2_grid_BIP), sum(!any_nas_BIP), 1))
x <- multivariateLassosum::pseudovalidation(omni, beta = array_Beta, cor = cbind(corB_BIP), destandardize = FALSE)
x.max_BIP <- which.max(x)
paramsA_BIP <- paramsA_BIP[!any_nas_BIP,]
saveRDS(x, file = "Results/LDpred2grid_BIP_flambda.RDS")

#Now, let's use the sumstats format that the bigsnpr authors demands.
#We run the final -grid-sp models and the -auto models.
ss_SKZ$beta <- log(ss_SKZ$OR)
ss_SKZ <- ss_SKZ[, c("SNP", "CHR", "BP", "A1", "A2", "beta", "SE", "P")]
colnames(ss_SKZ) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
size_SKZ <- c(94015, 67390)
ss_SKZ$n_eff <- 4 / (1 / size_SKZ[2] + 1 / size_SKZ[1])
df_beta_SKZ <- snp_match(ss_SKZ, mapOmni)
ss_BIP <- ss_BIP[, c("ID", "CHR", "POS", "A1", "A2", "BETA", "SE", "PVAL")]
colnames(ss_BIP) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
size_BIP <- c(371549, 41917)
ss_BIP$n_eff <- 4 / (1 / size_BIP[2] + 1 / size_BIP[1])
df_beta_BIP <- snp_match(ss_BIP, mapOmni)
 
#Final model grid-sp models
#SKZ
params_SKZ <- expand.grid(p = paramsA_SKZ[x.max_SKZ,"p"], h2 = paramsA_SKZ[x.max_SKZ,"h2"], sparse = TRUE)
ldpred2_grid_SKZ <- snp_ldpred2_grid(corr0, df_beta_SKZ, params_SKZ, ncores = ncores)
saveRDS(ldpred2_grid_SKZ, file = "Results/LDpred2grid_SKZ.RDS")
final_beta_SKZ <- as.matrix(ldpred2_grid_SKZ)
#BIP
params_BIP <- expand.grid(p = paramsA_BIP[x.max_BIP,"p"], h2 = paramsA_BIP[x.max_BIP,"h2"], sparse = TRUE)
ldpred2_grid_BIP <- snp_ldpred2_grid(corr0, df_beta_BIP, params_BIP, ncores = ncores)
saveRDS(ldpred2_grid_BIP, file = "Results/LDpred2grid_BIP.RDS")
final_beta_BIP <- as.matrix(ldpred2_grid_BIP)

# PRS (GSA data set)
#SKZ
final_pred_SKZ <- big_prodMat(GOmni, as.matrix(final_beta_SKZ))
saveRDS(final_pred_SKZ, file = paste0("Results/LDpred2gridOmniSKZ.Rdata"))
#BIP
final_pred_BIP <- big_prodMat(GOmni, as.matrix(final_beta_BIP))
saveRDS(final_pred_BIP, file = paste0("Results/LDpred2gridOmniBIP.Rdata")) 

# PRS (Omni data set)
#SKZ
final_pred_SKZ <- big_prodMat(GGSA, as.matrix(final_beta_SKZ))
saveRDS(final_pred_SKZ, file = paste0("Results/LDpred2gridGSASKZ.Rdata"))
#BIP
final_pred_BIP <- big_prodMat(GGSA, as.matrix(final_beta_BIP))
saveRDS(final_pred_BIP, file = paste0("Results/LDpred2gridGSABIP.Rdata")) 


#LDpred2-auto for SKZ
df_beta_SKZ <- as.data.frame(df_beta_SKZ)[,c("beta","beta_se","n_eff")]
ldsc_SKZ <- snp_ldsc(ld, length(ld), chi2 = (df_beta_SKZ$beta / df_beta_SKZ$beta_se)^2, sample_size = df_beta_SKZ$n_eff, blocks = NULL, ncores = ncores)
h2_est_SKZ <- ldsc_SKZ[[2]]
multi_auto_SKZ <- snp_ldpred2_auto(corr0, df_beta_SKZ, h2_init = h2_est_SKZ, vec_p_init,
				   ncores = NCORES, sparse = TRUE, allow_jump_sign = FALSE)
saveRDS(multi_auto_SKZ, file = "Results/LDpred2SKZ.RDS") 
beta_auto_SKZ <- sapply(multi_auto_SKZ, function(auto) auto$beta_est)
beta_auto_SKZ_sparse <- sapply(multi_auto_SKZ, function(auto) auto$beta_est_sparse)

#PRS (OmniExpress data set)
pred_auto_SKZ <- big_prodMat(GOmni, beta_auto_SKZ_sparse)
sc <- apply(pred_auto_SKZ, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_SKZ_LDpred2 <- rowMeans(as.data.frame(beta_auto_SKZ_sparse)[, keep])
final_pred_auto_SKZ <- big_prodVec(GOmni, Beta_SKZ_LDpred2)
saveRDS(final_pred_auto_SKZ, file = paste0("Results/LDpred2OmniSKZ.Rdata"))
 
#PRS (GSA data set)
pred_auto_SKZ <- big_prodMat(GGSA, beta_auto_SKZ_sparse)
sc <- apply(pred_auto_SKZ, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_SKZ_LDpred2 <- rowMeans(as.data.frame(beta_auto_SKZ_sparse)[, keep])
final_pred_auto_SKZ <- big_prodVec(GGSA, Beta_SKZ_LDpred2)
saveRDS(final_pred_auto_SKZ, file = paste0("Results/LDpred2GSASKZ.Rdata"))

#LDpred2-auto for BIP
df_beta_BIP <- as.data.frame(df_beta_BIP)[,c("beta","beta_se","n_eff")]
ldsc_BIP <- snp_ldsc(ld, length(ld), chi2 = (df_beta_BIP$beta / df_beta_BIP$beta_se)^2, sample_size = df_beta_BIP$n_eff, blocks = NULL, ncores = ncores)
h2_est_BIP <- ldsc_BIP[[2]]
multi_auto_BIP <- snp_ldpred2_auto(corr0, df_beta_BIP, h2_init = h2_est_BIP, vec_p_init,
				   ncores = ncores, sparse = TRUE, allow_jump_sign = FALSE)
saveRDS(multi_auto_BIP, file = "Results/LDpred2BIP.RDS")
beta_auto_BIP <- sapply(multi_auto_BIP, function(auto) auto$beta_est)
beta_auto_BIP_sparse <- sapply(multi_auto_BIP, function(auto) auto$beta_est_sparse)
 
#PRS (OmniExpress data set)
pred_auto_BIP <- big_prodMat(GOmni, beta_auto_BIP_sparse)
sc <- apply(pred_auto_BIP, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_BIP_LDpred2 <- rowMeans(as.data.frame(beta_auto_BIP_sparse)[, keep])
final_pred_auto_BIP <- big_prodVec(GOmni, Beta_BIP_LDpred2)
saveRDS(final_pred_auto_BIP, file = paste0("Results/LDpred2OmniBIP.Rdata")) 

#PRS (GSA data set)
pred_auto_BIP <- big_prodMat(GGSA, beta_auto_BIP_sparse)
sc <- apply(pred_auto_BIP, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_BIP_LDpred2 <- rowMeans(as.data.frame(beta_auto_BIP_sparse)[, keep])
final_pred_auto_BIP <- big_prodVec(GGSA, Beta_BIP_LDpred2)
saveRDS(final_pred_auto_BIP, file = paste0("Results/LDpred2GSABIP.Rdata")) 
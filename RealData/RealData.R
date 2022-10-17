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
setwd(".../RealData")

#---- MAF control ----
#We keep the SNPs where the minor allele frequency is higher than 10% when calculated among European subject 
#in our 1000 Genomes reference data.
sujetsRef <- data.table::fread("all_chrs_with_cms4.fam")
ancestry <- data.table::fread("igsr_samples.tsv")
out <- ancestry[which(ancestry$`Sample name` %in% sujetsRef$V1 & ancestry$`Superpopulation code` == "EUR"), c("Sample name")]
out$V2 <- out[,1]
data.table::fwrite(data.table::data.table(out), "/Data/sujetsEUR.txt", col.names = FALSE, sep = "\t", row.names = FALSE)
#Control for MAF<0.1. 3 639 921 SNPs remaining.
system(".../plink --bfile .../all_chrs_with_cms4 --extract .../Data/communMarkerRef.txt --keep .../Data/sujetsEUR.txt --maf 0.1 --make-bed --out .../Data/methodRef")

#---- Pseudo summary statistics ----
#Reference panel (1000 Genome)
ref.bfile <- "Data/methodRef"
Xr <- bigsnpr::bed(paste0(ref.bfile,".bed"))

#Test data 1 (OmniExpress)
Omni <- "Data/methodOmni"
omnibim <- data.table::fread(paste0(Omni, ".bim"))

#Test data 2 (GSA)
GSA <- "Data/methodGSA"
GSAbim <- data.table::fread(paste0(GSA, ".bim"))

#Summary statistics
#SKZ
ss_SKZ <- fread('PGC3_SCZ_wave3_public.v2.tsv',fill = TRUE)
size_SKZ <- sum(94015, 67390)
ss_SKZ$CHR <- as.numeric(ss_SKZ$CHR)

#BIP
ss_BIP <- fread('pgc-bip2021-all.vcf.tsv')
names(ss_BIP)[names(ss_BIP) == '#CHROM'] <- "CHR"
size_BIP <- sum(371549,41917)
ss_BIP$CHR <- as.numeric(ss_BIP$CHR)

#Matching of the sumstats SKZ and BIP (number of SNPs in common = 7 540 573)
matchSS <- lassosum:::matchpos(tomatch = ss_SKZ, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "CHR", pos = "BP", ref.pos = "POS", ref = "A1", ref.ref = "A1", alt = "A2", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
ss_SKZ <- ss_SKZ[matchSS$order,]
ss_BIP <- ss_BIP[matchSS$ref.extract,]

#Matching of the OmniExpress data set with the sumstats (ss_BIP and ss_SKZ are already in the same order) (number of SNPs in common = 6 398 846)
matchTest <- lassosum:::matchpos(tomatch = omnibim, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
testOrder <- matchTest$order
ss_BIP <- ss_BIP[matchTest$ref.extract,]
ss_SKZ <- ss_SKZ[matchTest$ref.extract,]

#Matching of the GSA data set with the sumstats (ss_BIP, ss_SKZ, OmniExpress data set are already in the same order) (number of SNPs in common = 6 398 846)
matchTestGSA <- lassosum:::matchpos(tomatch = GSAbim, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
testOrderGSA <- matchTestGSA$order
ss_BIP <- ss_BIP[matchTestGSA$ref.extract,]
ss_SKZ <- ss_SKZ[matchTestGSA$ref.extract,]
testOrder <- testOrder[matchTestGSA$ref.extract]

#Matching of the reference data set with the sumstats (ss_BIP, ss_SKZ, OmniExpress and GSA data are already in the same order) (number of SNPs in common = 3 639 921)
matchXr <- lassosum:::matchpos(tomatch = Xr$map, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "CHR", pos = "physical.pos", ref.pos = "POS", ref = "allele1", ref.ref = "A1", alt = "allele2", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
XrOrder <- matchXr$order
ss_BIP <- ss_BIP[matchXr$ref.extract,]
ss_SKZ <- ss_SKZ[matchXr$ref.extract,]
testOrder <- testOrder[matchXr$ref.extract]
testOrderGSA <- testOrderGSA[matchXr$ref.extract]

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
Xtrg <- bigsnpr::bed_cprodVec(obj.bed = Xr, y.row = g, ind.col = XrOrder, center = XrMeans[XrOrder], scale = XrSds[XrOrder])

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

#Exportation of the sumstats where only the SNPs in common are found.
data.table::fwrite(ss_SKZ, "Data/ss_SKZ.txt")
data.table::fwrite(ss_BIP, "Data/ss_BIP.txt")

#Exporattion of the new pseudo sumstats.
data.table::fwrite(ssA_SKZ, "Data/ssA_SKZ.txt")
data.table::fwrite(ssA_BIP, "Data/ssA_BIP.txt")
data.table::fwrite(ssB_SKZ, "Data/ssB_SKZ.txt")
data.table::fwrite(ssB_BIP, "Data/ssB_BIP.txt")
XrMap <- Xr$map[XrOrder,]
data.table::fwrite(data.table::data.table(XrMap$marker.ID), "Data/communMarkerRef.txt", col.names = FALSE)

#Only keep these SNPs in common in the test data bfiles.
GSAbim <- GSAbim[testOrderGSA,]
omnibim <- omnibim[testOrder,]
data.table::fwrite(data.table::data.table(omnibim$V2), "Data/communMarkerOmni.txt", col.names = FALSE, row.names = FALSE)
data.table::fwrite(data.table::data.table(GSAbim$V2), "Data/communMarkerGSA.txt", col.names = FALSE, row.names = FALSE)
system("/mnt-biostats/Plink/plink --bfile .../Data/methodOmni --extract .../Data/communMarkerOmni.txt --make-bed --out .../Data/methodOmni")
system("/mnt-biostats/Plink/plink --bfile .../Data/methodGSA --extract .../Data/communMarkerGSA.txt --make-bed --out .../Data/methodGSA")

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
omni <- "/home/jricard/Projet_Meriem/DonneesReellesMAF/Data/methodOmni"
omnibim <- data.table::fread(paste0(Omni, ".bim"))
omnifam <- data.table::fread(paste0(Omni, ".fam"))

#Test data GSA.
GSA <- "/home/jricard/Projet_Meriem/DonneesReellesMAF/Data/methodGSA"
GSAbim <- data.table::fread(paste0(GSA, ".bim"))

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

#Just to validate that sumstats are in the same order as the OmniExpress and GSA test data sets.
#If not, be sure that every data sets are in the same SNPs order before using the methods.
#Objects that result from lassosum:::matchpos will then be used to reverse effect size if alleles are reversed between sets.
matchGSA_SKZ <- lassosum:::matchpos(tomatch = ss_SKZ, ref.df = GSAbim, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "BP", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchOmni_SKZ <- lassosum:::matchpos(tomatch = ss_SKZ, ref.df = omnibim, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "BP", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSA_BIP <- lassosum:::matchpos(tomatch = ss_BIP, ref.df = GSAbim, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "POS", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchOmni_BIP <- lassosum:::matchpos(tomatch = ss_BIP, ref.df = omnibim, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "POS", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

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
nProb <- 0
phenotypic.genetic.Var.Cov.matrix.gencov <- array(data = NA, dim = c(2,2, nbr_SNP))
for (j in 1:nbr_SNP) {
  h2_BIP <- heritabilite[j, "h_BIP"]
  h2_SKZ <- heritabilite[j, "h_SKZ"]
  rho    <- heritabilite[j, "rho"]
  detCond<- ((h2_BIP*h2_SKZ)-(rho^2)) <= 0
  if(detCond){rho <- sqrt(max(c((h2_BIP*h2_SKZ)-0.001,0))); nProb <- nProb+1}
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

#flambda function for validation
f_lambda <- function(beta_SKZ_lambda, beta_BIP_lambda, r_hat, sd, beta){
  bXy <- r_hat %*% beta 
  weight <- 1/sd
  weight[!is.finite(weight)] <- 0
  scaled.beta_SKZ <- as.matrix(Diagonal(x=weight) %*% beta_SKZ_lambda)
  scaled.beta_BIP <- as.matrix(Diagonal(x=weight) %*% beta_BIP_lambda)
  pgs_SKZ <- pgs(bfile = omni, weights=scaled.beta_SKZ)
  pgs_BIP <- pgs(bfile = omni, weights=scaled.beta_BIP)
  n <- nrow(omnifam)
  if(ncol(pgs_SKZ)>1){
    pred <- matrix(data = NA,nrow = 2*nrow(pgs_SKZ),ncol = ncol(pgs_SKZ))
    for(i in 1:ncol(pgs_SKZ)){pred[,i] <- c(rbind(pgs_SKZ[,i],pgs_BIP[,i]))}
  }else{
    pgs_SKZ<- as.vector(pgs_SKZ)
    pgs_BIP<- as.vector(pgs_BIP)
    pred <- c(rbind(pgs_SKZ,pgs_BIP))
  }
  pred2 <- scale(pred, scale=F)
  bXXb <- colSums(pred2^2) / n
  result <- as.vector(bXy / sqrt(bXXb))
  return(result)
}

#---- Thresholding ----
#Thresholds are fixed using Trubetskoy (2022) and Mullins (2021) recommendations for, respectively, SZ and BP.
betaSKZ <- log(ss_SKZ$OR)*matchGSA_SKZ$rev
betaBIP <- ss_BIP$BETA*matchGSA_BIP$rev
betaSKZ[ss_SKZ$P > 0.05] <- 0
betaBIP[ss_BIP$PVAL > 0.10] <- 0
PRSthresholdGSASKZ <- pgs(GSA, weights=betaSKZ)
PRSthresholdGSABIP <- pgs(GSA, weights=betaBIP)
saveRDS(PRSthresholdGSASKZ, file = "Results/ThresholdingGSASKZ.RDS")
saveRDS(PRSthresholdGSABIP, file = "Results/ThresholdingGSABIP.RDS")

betaSKZ<- log(ss_SKZ$OR)*matchOmni$rev
betaBIP <- ss_BIP$BETA*matchOmni$rev
betaSKZ[ss_SKZ$P > 0.05] <- 0
betaBIP[ss_BIP$PVAL > 0.10] <- 0
PRSthresholdOmniSKZ <- pgs(omni, weights=betaSKZ)
PRSthresholdOmniBIP <- pgs(omni, weights=betaBIP)
saveRDS(PRSthresholdOmniSKZ, file = "Results/ThresholdingOmniSKZ.RDS")
saveRDS(PRSthresholdOmniBIP, file = "Results/ThresholdingOmniBIP.RDS")


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
BETA <- data.frame()
for(idx in 1:NbrLambdas){
  if(idx == 1){
    BETA <- data.frame(c(rbind(outMulti$beta[[1]][,1],rbind(outMulti$beta[[1]][,2]))))
  }else{
    add_idx <- c(rbind(outMulti$beta[[1]][,(2*idx)-1],rbind(outMulti$beta[[1]][,2*idx])))
    BETA <- cbind(BETA, add_idx)
  }
}
BETA <- as.matrix(BETA)
xMulti <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = c(rbind(corB_SKZ, corB_BIP)), sd = sdOmni, beta = BETA)
names(xMulti) <- paste0("lamdba_", AllLambdas)
maxMulti <- which.max(xMulti)
saveRDS(xMulti, file = "Results/multiFlambda.RDS")

#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
omni$V1 <- as.factor(omni$V1)
matchOmni <- lassosum:::matchpos(tomatch = outMulti$sumstats[[1]], ref.df = omni, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
GSA$V1 <- as.factor(GSA$V1)
matchGSA <- lassosum:::matchpos(tomatch = outMulti$sumstats[[1]], ref.df = GSA, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)


#PRS (OmniExpress data set)
#Reverse alleles if necessary
Lam <- AllLambdas[maxMulti]
# SKZ
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightOmni) %*% (mat_Beta_SKZ[,which(AllLambdas == Lam)]*matchOmni$rev))
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiOmniSKZ_", Lam, ".RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightOmni) %*% (mat_Beta_BIP[,which(AllLambdas == Lam)]*matchOmni$rev))
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiOmniBIP_", Lam, ".RDS"))

#PRS (GSA data set)
#Reverse alleles if necessary
# SKZ
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightGSA) %*% (mat_Beta_SKZ[,which(AllLambdas == Lam)]*matchGSA$rev))
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGSASKZ_", Lam, ".RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightGSA) %*% (mat_Beta_BIP[,which(AllLambdas == Lam)]*matchGSA$rev))
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGSABIP_", Lam, ".RDS"))

#---- Genetic_cov ----
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
BETA <- data.frame()
for(idx in 1:NbrLambdas){
  if(idx == 1){
    BETA <- data.frame(c(rbind(outMultiGenCov$beta[[1]][,1],rbind(outMultiGenCov$beta[[1]][,2]))))
  }else{
    add_idx <- c(rbind(outMultiGenCov$beta[[1]][,(2*idx)-1],rbind(outMultiGenCov$beta[[1]][,2*idx])))
    BETA <- cbind(BETA, add_idx)
  }
}
BETA <- as.matrix(BETA)
xMultiGenCov <- f_lambda(beta_SKZ_lambda = mat_Beta_SKZ, beta_BIP_lambda = mat_Beta_BIP, r_hat = c(rbind(corB_SKZ, corB_BIP)), sd = sdOmni, beta = BETA)
names(xMultiGenCov) <- paste0("lamdba_", AllLambdas)
maxMultiGenCov <- which.max(xMultiGenCov)
saveRDS(xMultiGenCov, file = paste0("Results/multiGenCovFlambda.RDS"))

#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
matchOmni <- lassosum:::matchpos(tomatch = outMultiGenCov$sumstats[[1]], ref.df = omni, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSA <- lassosum:::matchpos(tomatch = outMultiGenCov$sumstats[[1]], ref.df = GSA, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

#PRS (OmniExpress data set)
#Reverse alleles if necessary
Lam <- AllLambdas[maxMultiGenCov]
# SKZ
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightOmni) %*% (mat_Beta_SKZ[,which(AllLambdas == Lam)]*matchOmni$rev))
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGenCovOmniSKZ_", Lam, ".RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightOmni) %*% (mat_Beta_BIP[,which(AllLambdas == Lam)]*matchOmni$rev))
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGenCovOmniBIP_", Lam, ".RDS"))

#PRS (GSA data set)
#Reverse alleles if necessary
# SKZ
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightGSA) %*% (mat_Beta_SKZ[,which(AllLambdas == Lam)]*matchGSA$rev))
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/multiGenCovGSASKZ_", Lam, ".RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightGSA) %*% (mat_Beta_BIP[,which(AllLambdas == Lam)]*matchGSA$rev))
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/multiGenCovGSABIP_", Lam, ".RDS"))

#---- LASSOSUM ----
cl <- makeCluster(ncores[1]-5) 
#We are using the lambdas proposed by the authors.
AllLambdas <- exp(seq(log(0.001), log(0.1), length.out = 20))

#Lassosum for SKZ
OGSKZ <- lassosum::lassosum.pipeline(cor = corA_SKZ, ref.bfile = omni, lambda = AllLambdas, #test.bfile = testGSA.file, destandardize = FALSE,
		              chr = ssA_SKZ$CHR, pos = ssA_SKZ$BP, A1 = ssA_SKZ$A1, A2 = ssA_SKZ$A2, s = 0.5,
                              trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)
saveRDS(OGSKZ, file = "Results/OGSKZ.RDS")

#Lassosum for BIP
cl <- makeCluster(ncores[1]-5) 
OGBIP <- lassosum::lassosum.pipeline(cor = corA_BIP, ref.bfile = omni, lambda = AllLambdas, #test.bfile = testGSA.file, destandardize = FALSE,
		              chr = ssA_BIP$CHR, pos = ssA_BIP$POS, A1 = ssA_BIP$A1, A2 = ssA_BIP$A2, s = 0.5,
                              trace = 1, LDblocks = LDblocks, cluster = cl)
stopCluster(cl)
saveRDS(OGBIP, file = "Results/OGBIP.RDS")

#SKZ Validation
BETA_SKZ <- OGSKZ$beta
cl <- makeCluster(ncores[1]-5) 
xOGSKZ <- lassosum:::pseudovalidation(bfile = omni, beta = BETA_SKZ[[1]], cor = corB_SKZ, sd = sdOmni, cluster = cl)
stopCluster(cl)
xOGSKZ <- as.vector(xOGSKZ)
names(xOGSKZ) <- as.character(AllLambdas)
maxOGSKZ <- which.max(xOGSKZ)
saveRDS(xOGSKZ, file = "Results/OGSKZFlambda.RDS")

#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
matchOmniSKZ <- lassosum:::matchpos(tomatch = OGSKZ$sumstats, ref.df = omni, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSASKZ <- lassosum:::matchpos(tomatch = OGSKZ$sumstats, ref.df = GSA, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

#BIP Validation
BETA_BIP <- OGBIP$beta
cl <- makeCluster(ncores[1]-5) 
xOGBIP <- lassosum:::pseudovalidation(bfile = omni, beta = BETA_BIP[[1]], cor = corB_BIP, sd = sdOmni, cluster = cl)
stopCluster(cl)
xOGBIP <- as.vector(xOGBIP)
names(xOGBIP) <- as.character(AllLambdas)
maxOGBIP <- which.max(xOGBIP)
saveRDS(xOGBIP, file = "Results/OGBIPFlambda.RDS")
  
#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
matchOmniBIP <- lassosum:::matchpos(tomatch = OGBIP$sumstats, ref.df = omni, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSABIP <- lassosum:::matchpos(tomatch = OGBIP$sumstats, ref.df = GSA, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chr", ref.chr = "V1", pos = "pos", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

#PRS (OmniExpress data set)
#Reverse alleles if necessary
# SKZ
Lam <- AllLambdas[maxOGSKZ]
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightOmni) %*% (BETA_SKZ[[1]][,which(AllLambdas == Lam)]*matchOmniSKZ$rev))
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/OGOmniSKZ_", Lam, ".RDS"))
# BIP 
Lam <- AllLambdas[maxOGBIP]
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightOmni) %*% (BETA_BIP[[1]][,which(AllLambdas == Lam)]*matchOmniBIP$rev))
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/OGOmniBIP_", Lam, ".RDS"))

#PRS (GSA data set)
#Reverse alleles if necessary
# SKZ
Lam <- AllLambdas[maxOGSKZ]
scaled.beta_estime_SKZ <- as.matrix(Diagonal(x=weightGSA) %*% (BETA_SKZ[[1]][,which(AllLambdas == Lam)]*matchGSASKZ$rev))
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("Results/OGGSASKZ_", Lam, ".RDS"))
# BIP 
Lam <- AllLambdas[maxOGBIP]
scaled.beta_estime_BIP <- as.matrix(Diagonal(x=weightGSA) %*% (BETA_BIP[[1]][,which(AllLambdas == Lam)]*matchGSABIP$rev))
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("Results/OGGSABIP_", Lam, ".RDS"))

#---- Clumping and thresholding ----
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
#Here we use the real sumstats, same as the thresholding method.
ss_SKZCT <- setNames(ss_SKZ[,c("CHR", "SNP", "BP", "A1", "A2", "OR", "P")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
ss_SKZCT$beta <- log(ss_SKZCT$beta)
ss_BIPCT <- setNames(ss_BIP[,c("CHR", "ID", "POS", "A1", "A2", "BETA", "PVAL")], c("chr", "rsid", "pos", "a0", "a1", "beta", "p"))
map <- setNames(obj.bigSNPOmni$map[,-(2:3)], c("chr", "pos", "a0", "a1"))
info_snpSKZ <- snp_match(ss_SKZCT, map)
lpSSKZ <- -log10(info_snpSKZ$p)
info_snpBIP <- snp_match(ss_BIPCT, map)
lpSBIP <- -log10(info_snpBIP$p)

#Clumping
#Thresholds are fixed using Trubetskoy (2022) and Mullins (2021) recommendations for, respectively, SZ and BP.
ind.keepSKZ <- snp_clumping(GOmni, infos.chr = info_snpSKZ$chr, S = info_snpSKZ$p, thr.r2 = 0.5, size = 250, infos.pos = info_snpSKZ$pos, ncores = NCORES)
ind.keepBIP <- snp_clumping(GOmni, infos.chr = info_snpBIP$chr, S = info_snpBIP$p, thr.r2 = 0.5, size = 250, infos.pos = info_snpBIP$pos, ncores = NCORES)

#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
matchGSA <- lassosum:::matchpos(tomatch = obj.bigSNPGSA$map, ref.df = obj.bigSNPOmni$map, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "chromosome", pos = "physical.pos", ref.pos = "physical.pos", ref = "allele1", ref.ref = "allele1", alt = "allele2", ref.alt = "allele2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

# PRS (OmniExpress data set)
#SKZ
prs <- snp_PRS(GOmni, betas.keep = info_snpSKZ$beta[ind.keepSKZ], ind.keep = ind.keepSKZ, lpS.keep = lpSSKZ[ind.keepSKZ], thr.list = -log10(0.05))
saveRDS(prs, file = paste0("Results/CTOmniSKZ.RDS"))
#BIP
prs <- snp_PRS(GOmni, betas.keep = info_snpBIP$beta[ind.keepBIP], ind.keep = ind.keepBIP, lpS.keep = lpSBIP[ind.keepBIP], thr.list = -log10(0.1))
saveRDS(prs, file = paste0("Results/CTOmniBIP.RDS"))

# PRS (GSA data set)
#Reverse alleles if necessary
#SKZ
prs <- snp_PRS(GGSA, betas.keep = info_snpSKZ$beta[ind.keepSKZ], ind.keep = ind.keepSKZ, lpS.keep = lpSSKZ[ind.keepSKZ], thr.list = -log10(0.05), same.keep = revGSA[ind.keepSKZ])
saveRDS(prs, file = paste0("Results/CTGSASKZ.RDS"))
#BIP
prs <- snp_PRS(GGSA, betas.keep = info_snpBIP$beta[ind.keepBIP], ind.keep = ind.keepBIP, lpS.keep = lpSBIP[ind.keepBIP], thr.list = -log10(0.1), same.keep = revGSA[ind.keepBIP])
saveRDS(prs, file = paste0("Results/CTGSABIP.RDS"))

#---- LDPRED2 ----
#Correlations 
corr <- snp_cor(GOmni, ncores = NCORES, size = 500/2)
tmp <- tempfile(tmpdir = "corr")
corr0 <- as_SFBM(corr, tmp)
saveRDS(corr, "Data/LDpred2Corr.RDS")
saveRDS(corr0, "Data/LDpred2Corr0.RDS")

#Again, let's use the sumstats format that the bigsnpr authors demands.
ss_SKZ$beta <- log(ss_SKZ$OR)
ss_SKZ <- ss_SKZ[, c("SNP", "CHR", "BP", "A1", "A2", "beta", "SE", "P")]
colnames(ss_SKZ) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
ss_SKZ$n_eff <- 4 / (1 / size_SKZ[2] + 1 / size_SKZ[1])
df_beta_SKZ <- snp_match(ss_SKZ, mapOmni)
ss_BIP <- ss_BIP[, c("ID", "CHR", "POS", "A1", "A2", "BETA", "SE", "PVAL")]
colnames(ss_BIP) <- c("rsid", "chr", "pos", "a1", "a0", "beta", "beta_se", "p")
ss_BIP$n_eff <- 4 / (1 / size_BIP[2] + 1 / size_BIP[1])
df_beta_BIP <- snp_match(ss_BIP, mapOmni)

#Always check the order of the SNPs and the order of the alleles in the object resulting from the method. Everything should be good.
matchOmniSKZ <- lassosum:::matchpos(tomatch = obj.bigSNPOmni$map, ref.df = df_beta_SKZ, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "chr", pos = "physical.pos", ref.pos = "pos", ref = "allele2", ref.ref = "a0", alt = "allele1", ref.alt = "a1",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchOmniBIP <- lassosum:::matchpos(tomatch = obj.bigSNPOmni$map, ref.df = df_beta_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "chr", pos = "physical.pos", ref.pos = "pos", ref = "allele2", ref.ref = "a0", alt = "allele1", ref.alt = "a1",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSASKZ <- lassosum:::matchpos(tomatch = obj.bigSNPGSA$map, ref.df = df_beta_SKZ, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "chr", pos = "physical.pos", ref.pos = "pos", ref = "allele2", ref.ref = "a0", alt = "allele1", ref.alt = "a1",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSABIP <- lassosum:::matchpos(tomatch = obj.bigSNPGSA$map, ref.df = df_beta_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "chromosome", ref.chr = "chr", pos = "physical.pos", ref.pos = "pos", ref = "allele2", ref.ref = "a0", alt = "allele1", ref.alt = "a1",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

#LDpred2 for SKZ
df_beta_SKZ <- as.data.frame(df_beta_SKZ)[,c("beta","beta_se","n_eff")]
ldsc_SKZ <- snp_ldsc2(corr, df_beta_SKZ)
#We use the heritability value estimated by SLDX-R and we reduce the burn-in as we conclude that they are not necessary for convergence here: 
h2_est_SKZ <- totalMatrix[1,1]
multi_auto_SKZ <- snp_ldpred2_auto(corr0, df_beta_SKZ, h2_init = h2_est_SKZ, vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
				   ncores = NCORES, sparse = TRUE, allow_jump_sign = FALSE, burn_in = 0, num_iter = 100)
saveRDS(multi_auto_SKZ, file = "Results/LDpred2SKZ.RDS") 
beta_auto_SKZ <- sapply(multi_auto_SKZ, function(auto) auto$beta_est)
beta_auto_SKZ_sparse <- sapply(multi_auto_SKZ, function(auto) auto$beta_est_sparse)

#PRS (OmniExpress data set)
#Reverse alleles if necessary
pred_auto_SKZ <- big_prodMat(GOmni, beta_auto_SKZ_sparse*matchOmniSKZ$rev)
sc <- apply(pred_auto_SKZ, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_SKZ_LDpred2 <- rowMeans(as.data.frame(beta_auto_SKZ_sparse*matchOmniSKZ$rev)[, keep])
final_pred_auto_SKZ <- big_prodVec(GOmni, Beta_SKZ_LDpred2)
saveRDS(final_pred_auto_SKZ, file = paste0("Results/LDpred2OmniSKZ.Rdata"))
 
#PRS (GSA data set)
#Reverse alleles if necessary
pred_auto_SKZ <- big_prodMat(GGSA, beta_auto_SKZ_sparse*matchGSASKZ$rev)
sc <- apply(pred_auto_SKZ, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_SKZ_LDpred2 <- rowMeans(as.data.frame(beta_auto_SKZ_sparse*matchGSASKZ$rev)[, keep])
final_pred_auto_SKZ <- big_prodVec(GGSA, Beta_SKZ_LDpred2)
saveRDS(final_pred_auto_SKZ, file = paste0("Results/LDpred2GSASKZ.Rdata"))

#LDpred2 for BIP
df_beta_BIP <- as.data.frame(df_beta_BIP)[,c("beta","beta_se","n_eff")]
ldsc_BIP <- snp_ldsc2(corr, df_beta_BIP)
#We use the heritability value estimated by SLDX-R and we reduce the burn-in as we conclude that they are not necessary for convergence here: 
h2_est_BIP <- totalMatrix[2,2]
multi_auto_BIP <- snp_ldpred2_auto(corr0, df_beta_BIP, h2_init = h2_est_BIP, vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
				   ncores = NCORES, sparse = TRUE, allow_jump_sign = FALSE, burn_in = 0, num_iter = 100)
saveRDS(multi_auto_BIP, file = "Results/LDpred2BIP.RDS")
beta_auto_BIP <- sapply(multi_auto_BIP, function(auto) auto$beta_est)
beta_auto_BIP_sparse <- sapply(multi_auto_BIP, function(auto) auto$beta_est_sparse)
 
#PRS (OmniExpress data set)
#Reverse alleles if necessary
pred_auto_BIP <- big_prodMat(GOmni, beta_auto_BIP_sparse*matchOmniBIP$rev)
sc <- apply(pred_auto_BIP, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_BIP_LDpred2 <- rowMeans(as.data.frame(beta_auto_BIP_sparse*matchOmniBIP$rev)[, keep])
final_pred_auto_BIP <- big_prodVec(GOmni, Beta_BIP_LDpred2)
saveRDS(final_pred_auto_BIP, file = paste0("Results/LDpred2OmniBIP.Rdata")) 

#PRS (GSA data set)
#Reverse alleles if necessary
pred_auto_BIP <- big_prodMat(GGSA, beta_auto_BIP_sparse*matchGSABIP$rev)
sc <- apply(pred_auto_BIP, 2, sd)
keep <- abs(sc - median(sc,na.rm = T)) < 3 * mad(sc,na.rm = T)
keep[which(is.na(keep))]<-FALSE
Beta_BIP_LDpred2 <- rowMeans(as.data.frame(beta_auto_BIP_sparse*matchGSABIP$rev)[, keep])
final_pred_auto_BIP <- big_prodVec(GGSA, Beta_BIP_LDpred2)
saveRDS(final_pred_auto_BIP, file = paste0("Results/LDpred2GSABIP.Rdata")) 
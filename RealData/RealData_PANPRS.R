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

#Follow RealData.R first.

#Control for MAF<0.1 on European subjects is already done. Data are found in
#.../RealData/Data/methodRef bfiles.

#---- Clumping, pvalue minimum ----
#Reference panel (1000 Genome).
ref.bfile <- "Data/methodRef"
Xr <- bigsnpr::bed(paste0(ref.bfile,".bed"))

#Sumstats
#SKZ.
ss_SKZ <- fread('PGC3_SCZ_wave3_public.v2.tsv',fill = TRUE)
size_SKZ <- sum(94015, 67390)
ss_SKZ$CHR <- as.numeric(ss_SKZ$CHR)

#BIP.
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

#Matching of the reference data set with the sumstats (ss_BIP and ss_SKZ are already in the same order) (number of SNPs in common = 3 639 921)
#To be sure that every SNP in the reference panel has a pvalue.
matchXr <- lassosum:::matchpos(tomatch = Xr, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
Xr <- Xr[matchRef$order,]
ss_BIP <- ss_BIP[matchXr$ref.extract,]
ss_SKZ <- ss_SKZ[matchXr$ref.extract,]

#Minimum p-value between BIP and SKZ symstats.
out <- data.frame("SNP" = Xr$V2, "pval" = apply(data.frame(ss_SKZ$P, ss_BIP$PVAL), 1, min))
out <- data.table::data.table(out)
data.table::fwrite(out, "Data/pvalueClumping.txt", row.names = FALSE, sep = "\t", col.names = TRUE)

#Clumping (number of resulting SNPs = 365 527)
system("/mnt-biostats/Plink/plink --bfile .../Data/clumpingRef --clump .../Data/pvalueClumping.txt --clump-p1 1 --clump-r2 0.5 --clump-kb 250 --clump-snp-field SNP --clump-field pval --out .../PANPRS/Data/clumpingOut")

#Only keep these SNPs in the reference data set.
system("/mnt-biostats/Plink/plink --bfile .../Data/clumpingRef --extract .../PANPRS/Data/clumpingOut.clumped --make-bed --out .../PANPRS/Data/methodRef")

#---- Pseudo summary statistics ----
#Ref panel
ref.bfile <- "PANPRS/Data/methodRef"
refbim <- data.table::fread(paste0(ref.bfile,".bim"))
Xr <- bigsnpr::bed(paste0(ref.bfile,".bed"))

#Test data 1 (OmniExpress)
Omni <- "Data/methodOmni"
omnibim <- data.table::fread(paste0(Omni, ".bim"))

#Test data 2 (GSA)
GSA <- "Data/methodGSA"
GSAbim <- data.table::fread(paste0(GSA, ".bim"))

#Matching of the reference data set after clumping with the sumstats (number of SNPs in common = 365 527)
matchRef <- lassosum:::matchpos(tomatch = refbim, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
XrOrder <- matchRef$order
ss_BIP <- ss_BIP[matchRef$ref.extract,]
ss_SKZ <- ss_SKZ[matchRef$ref.extract,]

#Matching of the OmniExpress data set with the sumstats (number of SNPs in common = 365 527)
matchTest <- lassosum:::matchpos(tomatch = omnibim, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
testOrder <- matchTest$order
ss_BIP <- ss_BIP[matchTest$ref.extract,]
ss_SKZ <- ss_SKZ[matchTest$ref.extract,]
XrOrder <- XrOrder[matchTest$ref.extract]

#Matching of the GSA data set with the sumstats (number of SNPs in common = 365 527)
matchTestGSA <- lassosum:::matchpos(tomatch = GSAbim, ref.df = ss_BIP, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "POS", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
testOrderGSA <- matchTestGSA$order
ss_BIP <- ss_BIP[matchTestGSA$ref.extract,]
ss_SKZ <- ss_SKZ[matchTestGSA$ref.extract,]
testOrder <- testOrder[matchTestGSA$ref.extract]
XrOrder <- XrOrder[matchTestGSA$ref.extract]

#From p-values to correlations.
cor_BIP <- lassosum::p2cor(p = ss_BIP$PVAL, n = size_BIP, sign=ss_BIP$BETA)
cor_SKZ <- lassosum::p2cor(p = ss_SKZ$P, n = size_SKZ, sign=log(ss_SKZ$OR))

#Let's generate the pseudo sumstats
#The reference panel needs to be centered and scaled
#xi=(mi-2pi)/sqrt(2pi(1-pi)) where pi is the minor allele frequency
#of the ith genetic marker and mi is the ith column vector of the allele count matrix M.
nAverage <- mean(c(size_BIP, size_SKZ))
nB <- 0.1*nAverage
nr <- Xr$nrow
set.seed(1)
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

#SKZ pseudo sumstats.
n_SKZ <- size_SKZ
nA_SKZ <- n_SKZ - nB
r_SKZ <- cor_SKZ
varr_SKZ <- var(r_SKZ)
rA_SKZ <- r_SKZ + sqrt((nB*varr_SKZ)/(nA_SKZ*nr))*Xtrg
rB_SKZ <- (1/nB)*(n_SKZ*r_SKZ - nA_SKZ*rA_SKZ)
ssA_SKZ <- ss_SKZ
ssB_SKZ <- ss_SKZ

#BIP pseudo sumstats.
n_BIP <- size_BIP
nA_BIP <- n_BIP - nB
r_BIP <- cor_BIP
varr_BIP <- var(r_BIP)
rA_BIP <- r_BIP + sqrt((nB*varr_BIP)/(nA_BIP*nr))*Xtrg
rB_BIP <- (1/nB)*(n_BIP*r_BIP - nA_BIP*rA_BIP)
ssA_BIP <- ss_BIP 
ssB_BIP <- ss_BIP 

#From correlation to p-values.
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

#Exportation of the new pseudo sumstats.
data.table::fwrite(ssA_SKZ, "PANPRS/Data/ssA_SKZ.txt")
data.table::fwrite(ssA_BIP, "PANPRS/Data/ssA_BIP.txt")
data.table::fwrite(ssB_SKZ, "PANPRS/Data/ssB_SKZ.txt")
data.table::fwrite(ssB_BIP, "PANPRS/Data/ssB_BIP.txt")

#Keep these SNPs only in the test data bfiles.
GSAbim <- GSAbim[testOrderGSA,]
omnibim <- omnibim[testOrder,]
data.table::fwrite(data.table::data.table(omnibim$V2), "PANPRS/Data/markerClumpedOmni.txt", col.names = FALSE, row.names = FALSE)
data.table::fwrite(data.table::data.table(GSAbim$V2), "PANPRS/Data/markerClumpedGSA.txt", col.names = FALSE, row.names = FALSE)
system(".../plink --bfile .../PANPRS/Data/methodOmni --extract .../PANPRS/Data/markerClumpedOmni.txt --make-bed --out .../PANPRS/Data/methodOmni")
system(".../plink --bfile .../PANPRS/Data/methodGSA --extract .../PANPRS/Data/markerClumpedGSA.txt --make-bed --out .../PANPRS/Data/methodGSA")

#---- PANPRS ----
#Test data OmniExpress and GSA
omni <- "PANPRS/Data/methodOmni"
omnibim <- data.table::fread(paste0(omni, ".bim"))
GSA <- "PANPRS/Data/methodGSA"

#Sumstats SKZ
corA_SKZ <- ssA_SKZ$Cor
size_SKZ <- sum(94015, 67390)
ss_SKZ <- ss_SKZ[ss_SKZ$SNP %in% ssA_SKZ$SNP]
corB_SKZ <- ssB_SKZ$Cor

#Sumstats BIP
corA_BIP <- ssA_BIP$Cor
size_BIP <- sum(371549, 41917)
ss_BIP <- ss_BIP[ss_BIP$ID %in% ssA_BIP$ID]
corB_BIP <- ssB_BIP$Cor

#Just to validate that sumstats are in the same order as the OmniExpress and GSA test data sets.
#If not, be sure that every data sets are in the same SNPs order before using the methods.
matchOmni <- lassosum:::matchpos(tomatch = omnibim, ref.df = ssA_SKZ, auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "V1", ref.chr = "CHR", pos = "V4", ref.pos = "BP", ref = "V5", ref.ref = "A1", alt = "V6", ref.alt = "A2",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)

#PANPRS
Zmatrix <- matrix(c(sign(ssA_SKZ$Cor)*abs(qnorm(p=(ssA_SKZ$P)/2)), sign(ssA_BIP$Cor)*abs(qnorm(p=(ssA_BIP$PVAL)/2))), ncol = 2)
row.names(Zmatrix) <- omni$V2[matchOmni$order]

#Compute LDscore using PLINK in the method reference bfiles, OmniExpress.
system(paste0(".../plink --bfile ", omni, " --r2 0.2 --ld-window-kb 250 --ld-window-r2 0 --out .../PANPRS/Data/PANPRSLD-RefOmni"))
plinkLDgenome <- fread("PANPRS/Data/PANPRSLD-RefOmni.ld")
colnames(plinkLDgenome) <- c("CHR_A", "BP_A",  "SNP_A", "CHR_B", "BP_B",  "SNP_B", "R")
#Generate the set of tuning parameters on a modified version of gsPEN which doesn't fit the model, but only outputs the tuning matrix.
initialTuning <- SummaryLasso::gsPEN(summaryZ = Zmatrix, Nvec = c(size_SKZ, size_BIP), plinkLD = plinkLDgenome, numChrs = 2, fit = FALSE) 

for(chr in 1:22){
  plinkLD <- plinkLDgenome[plinkLDgenome$CHR_A == chr & plinkLDgenome$CHR_B == chr,]
  plinkLD <- as.data.frame(plinkLD)
  ZmatrixChr <- Zmatrix[ssA_SKZ$CHR == chr,]
  PANPRSchr <- summaryLasso::gsPEN(summaryZ = ZmatrixChr, Nvec = c(size_SKZ, size_BIP), plinkLD = plinkLD, NumIter = 1000, breaking = 1, 
    numChrs = 1, ChrIndexBeta = 0, Init_summaryBetas = 0, Zscale = 1, 
    RupperVal = NULL, tuningMatrix = initialTuning, penalty = c("mixLOG"), outputAll = 0, fit = TRUE)
  saveRDS(PANPRSchr, paste0("PANPRS/Results/PANPRS-chr", chr, ".RDS"))
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
  saveRDS(PANPRS, paste0("PANPRS/Results/PANPRS.RDS"))
}

#Validation
whereSKZ <- which(stringr::str_detect(dimnames(PANPRS$BetaMatrix)[[2]], "trait1"))
orderSKZ <- dimnames(PANPRS$BetaMatrix)[[2]][whereSKZ]
orderSKZ <- substr(orderSKZ, 1, nchar(orderSKZ)-7)
whereBIP <- which(stringr::str_detect(dimnames(PANPRS$BetaMatrix)[[2]], "trait2"))
orderBIP <- dimnames(PANPRS$BetaMatrix)[[2]][whereBIP]
orderBIP <- substr(orderBIP, 1, nchar(orderBIP)-7)
#BIP and SKZ betas are in the same order
NbrLambdas <- dim(PANPRS$tuningMatrix)[1]
mat_Beta_SKZ <- PANPRS$BetaMatrix[,whereSKZ]
mat_Beta_BIP <- PANPRS$BetaMatrix[,whereBIP]

#Be sure that the betas are in the same SNP order as they are in the bim file of the OmniExpress data set.
order <- match(omnibim$V2, orderSKZ)
mat_Beta_SKZ <- mat_Beta_SKZ[,order]
mat_Beta_BIP <- mat_Beta_BIP[,order]
array_Beta <- array(c(mat_Beta_SKZ, mat_Beta_BIP), dim = c(dim(mat_Beta_SKZ)[1], dim(mat_Beta_SKZ)[2], 2))
x <- multivariateLassosum::pseudovalidation(r = cbind(corB_SKZ, corB_BIP), keep_sujets = parsed.2$keep, beta = array_Beta, destandardize = FALSE)
names(x) <- apply(PANPRS$tuningMatrix[DiffZeroCrit,], 1, FUN = function(x){paste0(round(x,4), collapse = "-")})
saveRDS(x, file = "PANPRS/Results/Valeurs_f_lambda_PANPRS.Rdata")

#Finally, watch out for allele order.
matchOmni <- lassosum:::matchpos(tomatch = ssA_SKZ, ref.df = omnibim , auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "BP", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)
matchGSA <- lassosum:::matchpos(tomatch = ssA_SKZ, ref.df = testGSA.bim , auto.detect.tomatch = F, auto.detect.ref = F,
  chr = "CHR", ref.chr = "V1", pos = "BP", ref.pos = "V4", ref = "A1", ref.ref = "V5", alt = "A2", ref.alt = "V6",
  exclude.ambiguous = F, silent = F, rm.duplicates = F)


#PRS Omni
#Reverse alleles if necessary
AllLambdas <- names(x)
AllLambdasMax <- which.max(x)
Lam <- AllLambdas[AllLambdasMax]
# SKZ
scaled.beta_estime_SKZ <- as.matrix(t(mat_Beta_SKZ)[,which(AllLambdas == Lam)]*matchOmni$rev)
PGS_estime_SKZ <- pgs(omni, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("PANPRS/Results/PANPRSOmniSKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(t(mat_Beta_BIP)[,which(AllLambdas == Lam)]*matchOmni$rev)
PGS_estime_BIP <- pgs(omni, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("PANPRS/Results/PANPRSOmniBIP.RDS"))

#PRS GSA
#Reverse alleles if necessary
# SKZ
scaled.beta_estime_SKZ <- as.matrix(t(mat_Beta_SKZ)[,which(AllLambdas == Lam)]*matchGSA$rev)
PGS_estime_SKZ <- pgs(GSA, weights=scaled.beta_estime_SKZ)
saveRDS(PGS_estime_SKZ, file = paste0("PANPRS/Results/PANPRSGSASKZ.RDS"))
# BIP 
scaled.beta_estime_BIP <- as.matrix(t(mat_Beta_BIP)[,which(AllLambdas == Lam)]*matchGSA$rev)
PGS_estime_BIP <- pgs(GSA, weights=scaled.beta_estime_BIP)
saveRDS(PGS_estime_BIP, file = paste0("PANPRS/Results/PANPRSGSABIP.RDS"))


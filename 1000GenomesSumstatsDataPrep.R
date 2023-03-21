#1000Genomes phase 3 data are available here:
#https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/

#Using 1000 Genomes data, we followed these steps:
#1- merge the chr1 to chr22 files;
#2- extract only the SNPs for which summary statistics are available for both traits;
#3- remove SNPs where call rate is lower than 0.01;
#4- remove SNPs where minor allele frequency is lower than 0.001;

#FIRST download the data.
#Our codes suppose that every data are in the same directory.
#Please enter the path to the directory where your 1000 genomes data are found:
path <- ".../"
setwd(path)

#Please enter the path to your plink and vcftools software (to use vcf-concat).
pathPlink <- ".../plink"
pathVcf <- ".../vcftools"

library(lassosum)
library(Rcpp)
library(data.table)
library(matrixcalc)
library(Matrix)

#---- Merging chr1 to chr22 ----
#This step requires vcftools software. Please feel free to change the path to you vcftools software.
system(paste0(pathVcf, " ", path, "ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | gzip -c > ", path, "ALL_allchr.vcg.gz"))

#Convert to PLINK bfiles. Please feel free to change the path to you PLINK software.
system(paste0(pathPlink, " --vcf ", path, "ALL_allchr.vcf.gz --make-bed --out ", path, "allchrs"))

#---- Extract SNPs with sumstats ----
ref.bfile <- "allchrs"

#Import summary statistics
##SKZ
#Data are available here : https://pgc.unc.edu/for-researchers/download-results/
ss_SKZ <- fread('PGC3_SCZ_wave3_public.v2.tsv',fill = TRUE)
#There's a problem at line 5328756. Column were not correctly delimited. This is how we corrected it.
ss_SKZ[5328756,1] <- 21 ;ss_SKZ[5328756,2] <- "rs148878475";ss_SKZ[5328756,3]<-9648204;
ss_SKZ[5328756,4] <- "C" ; ss_SKZ[5328756,5] <- "T" ; ss_SKZ[5328756,6] <- 0.9878 ;
ss_SKZ[5328756,7] <- 0.9845 ;ss_SKZ[5328756,8] <- 0.1375 ;ss_SKZ[5328756,9] <- 7.5732 ;
ss_SKZ[5328756,10] <- 1.0014 ;ss_SKZ[5328756,11] <- 0.0432 ;ss_SKZ[5328756,12] <- 0 ;
ss_SKZ$CHR <- as.numeric(ss_SKZ$CHR)

##BIP
#Data are available here : https://pgc.unc.edu/for-researchers/download-results/
ss_BIP <- fread('pgc-bip2021-all.vcf.tsv')
names(ss_BIP)[names(ss_BIP) == '#CHROM'] <- "CHR"
ss_BIP$CHR <- as.numeric(ss_BIP$CHR)

#SNPs found in BIP and SKZ summary stats:
#To do so, we prepare lists of SNPs information by trait.
sums_SKZ <- list(chr=ss_SKZ$CHR, pos=ss_SKZ$BP, A1=ss_SKZ$A1, A2=ss_SKZ$A2, cor=cor_SKZ, snp=NULL)
sums_SKZ <- as.data.frame(sums_SKZ)

sums_BIP <- list(chr=ss_BIP$CHR, pos=ss_BIP$POS, A1=ss_BIP$A1, A2=ss_BIP$A2, cor=cor_BIP, snp=NULL)
sums_BIP <- as.data.frame(sums_BIP)

#We compare the datasets using matchpos from the lassosum package. BIP is the reference.
m.commun_SKZ_BIP <- matchpos(tomatch = sums_SKZ, ref.df = sums_BIP, auto.detect.ref = F, 
                             ref.chr = "chr", chr = "chr",
                             ref.pos="pos", pos="pos",
                             ref.alt="A1", ref.ref="A2", alt="A1", ref="A2",
                             rm.duplicates = T, exclude.ambiguous = T, silent=T)

#The initial summary statistics, where we only keep the matching SNPs.
ss_SKZ.BIP <- ss_SKZ[m.commun_SKZ_BIP$order,]
ss_BIP.SKZ <- ss_BIP[m.commun_SKZ_BIP$ref.extract,]

#If alleles order doesn't match, modify it and reverse the association with the phenotype.
ss_SKZ.BIP[m.commun_SKZ_BIP$rev == -1, c("A1", "A2")] <- ss_SKZ.BIP[m.commun_SKZ_BIP$rev == -1, c("A2", "A1")]
ss_SKZ.BIP[m.commun_SKZ_BIP$rev == -1, "BETA"] <- ss_SKZ.BIP[m.commun_SKZ_BIP$rev == -1, "BETA"]*(-1)

#Save this change.
fwrite(data.table(ss_SKZ.BIP), "ss_SKZ.txt", col.names = TRUE, row.names = FALSE, sep = "\t")
fwrite(data.table(ss_BIP.SKZ), "ss_BIP.txt", col.names = TRUE, row.names = FALSE, sep = "\t")

#SNPs found in BIP, SKZ summary stats and 1000 genomes reference panel:
sums_BIP <- sums_BIP[m.commun_SKZ_BIP$ref.extract,]
ref.bim <- read.table2(paste0(ref.bfile, ".bim"))
m.ref <- matchpos(tomatch = ref.bim, ref.df = sums_BIP, auto.detect.ref = F, 
                  ref.chr = "chr", chr = "V1",
                  ref.pos="pos", pos="V4", 
                  ref.alt="A1", ref.ref="A2", alt="V5", ref="V6", 
                  rm.duplicates = T, exclude.ambiguous = T, silent=T)

ref.bim <- data.table(ref.bim[m.ref$order, "V2"])
fwrite(ref.bim$V2, "snpToExtract.txt", row.names = FALSE, col.names = FALSE)

toSwitch <- dat.table(ref.bim[m.ref$rev == -1, c("V2", "V6")])
fwrite(toSwitch, "snpToSwitch.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#Only keep these matching SNPs in our 1000 genomes bfiles
system(paste0(pathPlink, " --bfile ", path, "allchrs --extract ", path, "snpToExtract.txt --reference-allele ", path, "snpToSwitch.txt --make-bed --out ", path, "allchrs2"))

#---- MAF and call rate on European subjects QC ----
#We keep European subject.
#Ancestry are found here: https://www.internationalgenome.org/data-portal/sample
sujetsRef <- data.table::fread("allchrs4.fam")
ancestry <- data.table::fread("igsr_samples.tsv")
ancestry <- ancestry[which(ancestry$`Sample name` %in% sujetsRef$V1 & ancestry$`Superpopulation code` == "EUR"), c("Sample name")]
ancestry$V2 <- ancestry[,1]
data.table::fwrite(data.table::data.table(ancestry), "sujetsEUR.txt", col.names = FALSE, sep = "\t", row.names = FALSE)
system(paste0(pathPlink, " --bfile ", path, "allchrs2 --geno 0.01 --maf 0.001 --keep ", path, "sujetsEUR.txt --make-bed --keep-allele-order --out ", path, "allchrs4"))


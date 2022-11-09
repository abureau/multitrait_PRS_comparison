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
setwd(".../")

library(lassosum)
library(Rcpp)
library(data.table)
library(matrixcalc)
library(Matrix)

#---- Merging chr1 to chr22 ----
#This step requires vcftools software. Please feel free to change the path to you vcftools software.
system(".../vcf-concat concat ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | gzip -c > ALL_allchr.vcg.gz")

#Convert to PLINK bfiles. Please feel free to change the path to you PLINK software.
system(".../plink --vcf ALL_allchr.vcf.gz --make-bed --out allchrs")

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
Sums_SKZ <- list(chr=ss_SKZ$CHR, pos=ss_SKZ$BP, A1=ss_SKZ$A1, A2=ss_SKZ$A2, cor=cor_SKZ, snp=NULL)
Sums_SKZ[sapply(Sums_SKZ, is.null)] <- NULL
Sums_SKZ <- as.data.frame(Sums_SKZ)
Sums_SKZ[,"chr"] = as.numeric(Sums_SKZ[,"chr"])
Sums_BIP <- list(chr=ss_BIP$CHR, pos=ss_BIP$POS, A1=ss_BIP$A1, A2=ss_BIP$A2, cor=cor_BIP, snp=NULL)
Sums_BIP[sapply(Sums_BIP, is.null)] <- NULL
Sums_BIP <- as.data.frame(Sums_BIP)

#We compare the datasets using matchpos from the lassosum package.
m.commun_SKZ_BIP <- matchpos(tomatch = Sums_SKZ,ref.df = Sums_BIP, auto.detect.ref = F, 
                             ref.chr = "chr", ref.pos="pos",
                             ref.alt="A1", ref.ref="A2", 
                             rm.duplicates = T, exclude.ambiguous = T, silent=T)

ss_SKZ.BIP<-Sums_SKZ[m.commun_SKZ_BIP$order,]
ss_BIP.SKZ<-Sums_BIP[m.commun_SKZ_BIP$ref.extract,]

#SNPs found in BIP, SKZ summary stats and 1000 genomes reference panel:
ref.bim <- read.table2(paste0(ref.bfile, ".bim"))
m.ref <- matchpos(ss_SKZ.BIP, ref.bim, auto.detect.ref = F, 
                  ref.chr = "V1", ref.snp="V2", 
                  ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                  rm.duplicates = T, exclude.ambiguous = T, 
                  silent=T)

ref.bim <- ref.bim[m.ref$ref.extract,]
fwrite(ref.bim$V2, "snpToExtract.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

#Only keep these matching SNPs in our 1000 genomes bfiles
#Please feel free to change the path to you PLINK software.
system(".../plink --bfile allchrs --extract snpToExtract.txt --make-bed --out allchrs2")

#---- Call rate QC ----
#Please feel free to change the path to you PLINK software.
system(".../plink --bfile allchrs2 --geno 0.01 --make-bed --out allchrs3")

#---- MAF QC ----
#Please feel free to change the path to you PLINK software.
system(".../plink --bfile allchrs3 --maf 0.001 --make-bed --out allchrs4")

#---- Only keep European subjects ----
#We keep European subject.
#Ancestry are found here: https://www.internationalgenome.org/data-portal/sample
sujetsRef <- data.table::fread("allchrs4.fam")
ancestry <- data.table::fread("igsr_samples.tsv")
ancestry <- ancestry[which(ancestry$`Sample name` %in% sujetsRef$V1 & ancestry$`Superpopulation code` == "EUR"), c("Sample name")]
ancestry$V2 <- ancestry[,1]
data.table::fwrite(data.table::data.table(ancestry), "sujetsEUR.txt", col.names = FALSE, sep = "\t", row.names = FALSE)
#Please feel free to change the path to you PLINK software.
system(".../plink --bfileallchrs4 --keep sujetsEUR.txt --make-bed --out allchrs4")


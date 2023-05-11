library(pROC)

#Every code suppose that your data are all in the same directory.
#Complete the path to your data
pathBase <- ".../"

#Please specify which kind of scenario you want to simulate.
#This `for` loop will set the parameters needed for each simulation.
#Please enter the simulation type. It needs to be written the same ways as it is in the original paper of this project:
#"n = 29,330", "n = 10,139", "n = 29,330; Low Polygenicity", "n = 29,330; Moderate Correlation" or "29300ind_moderateCor/".#If an overlap between both trait was simulated, you can now use "n = 29,330; Overlap"
#If an overlap between traits was simulated, you can now use "n = 29,330; Overlap"
simuType <- "..."

if(simuType=="n = 29,330"){
  path <- paste0(pathBase, "29300ind/")
}else if(simuType=="n = 10,139"){
  path <- paste0(pathBase, "10139ind/")
}else if(simuType=="n = 29,330; Low Polygenicity"){
  path <- paste0(pathBase, "29300ind_lowPoly/")
}else if(simuType=="n = 29,330; Low Heritability"){
  path <- paste0(pathBase, "29300ind_lowHeri/")
}else if(simuType=="n = 29,330; Overlap"){
  path <- paste0(pathBase, "29300ind_overlap/")
}else{
  warning("Please provide an actual simulation scenario, written as it is in the original paper of this project.")
}

Cor <- data.frame()
for(k in 1:20){
        setwd(paste0(path, "Simulation_", k, "/B-LDX/"))

	#Simulated PRS
	simule_SKZ <- readRDS("PGS_simule_SKZ.Rdata")
	simule_BIP <- readRDS("PGS_simule_BIP.Rdata")

	#Simulated Phenotypes
	pheno_SKZ <- readRDS("pheno_simule_SKZ.Rdata")
	pheno_BIP <- readRDS("pheno_simule_BIP.Rdata")

	#Thresholding
	valid_SKZ <- readRDS("Valeurs_f_lambda_thresholding_SKZ.Rdata")
	valid_BIP <- readRDS("Valeurs_f_lambda_thresholding_BIP.Rdata")
	order_SKZ <- sort(valid_SKZ, decreasing=T)
	order_BIP <- sort(valid_BIP, decreasing=T)
	for(j in 1:length(valid_SKZ)){
		order_SKZ_j <- which(order_SKZ == valid_SKZ[j])
		order_BIP_j <- which(order_BIP == valid_BIP[j])
		val_SKZ <- names(valid_SKZ[j])
		val_BIP <- names(valid_BIP[j])
		th_SKZ <- readRDS(paste0("PGS_thresholding_SKZ_", val_SKZ,".Rdata"))
		th_BIP <- readRDS(paste0("PGS_thresholding_BIP_", val_BIP,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Thresholding", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = cor(simule_SKZ, th_SKZ), "BIP" = cor(simule_BIP, th_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(th_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(th_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "Thresholding", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}

	#MultivariateLassosum Standard
	valid <- readRDS("Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid, decreasing=T)
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		multi_SKZ <- readRDS(paste0("PGS_estime_SKZ_", val,".Rdata"))
		multi_BIP <- readRDS(paste0("PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "multiLassosum", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, multi_SKZ), "BIP" = cor(simule_BIP, multi_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(multi_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(multi_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "multiLassosum", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}

	#multivariateLassosum B-LDX
	valid <- readRDS("GenCov/Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid, decreasing=T)
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		gencov_SKZ <- readRDS(paste0("GenCov/PGS_estime_SKZ_", val,".Rdata"))
		gencov_BIP <- readRDS(paste0("GenCov/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Genetic_cov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, gencov_SKZ), "BIP" = cor(simule_BIP, gencov_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(gencov_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(gencov_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "Genetic_cov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}


	#multivariateLassosum adaptive Standard
	valid <- readRDS("BetaMulti/Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid, decreasing=T)
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		bm_SKZ <- readRDS(paste0("BetaMulti/PGS_estime_SKZ_", val,".Rdata"))
		bm_BIP <- readRDS(paste0("BetaMulti/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Adap_BetaMulti", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, bm_SKZ), "BIP" = cor(simule_BIP, bm_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(bm_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(bm_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "Adap_BetaMulti", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}

	#multivariateLassosum adaptive B-LDX
	valid <- readRDS("BetaGenCov/Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid, decreasing=T)
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		adap_gencov_SKZ <- readRDS(paste0("BetaGenCov/PGS_estime_SKZ_", val,".Rdata"))
		adap_gencov_BIP <- readRDS(paste0("BetaGenCov/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Adap_GenCov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, adap_gencov_SKZ), "BIP" = cor(simule_BIP, adap_gencov_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(adap_gencov_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(adap_gencov_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "Adap_GenCov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}

	#Lassosum
	valid_SKZ <- readRDS("OG/Valeurs_f_lambda_Lassosum_SKZ.Rdata")
	valid_BIP <- readRDS("OG/Valeurs_f_lambda_Lassosum_BIP.Rdata")
	order_SKZ <- sort(valid_SKZ, decreasing=T)
	order_BIP <- sort(valid_BIP, decreasing=T)
	for(j in 1:length(valid_SKZ)){
		order_SKZ_j <- which(order_SKZ == valid_SKZ[j])
		order_BIP_j <- which(order_BIP == valid_BIP[j])
		val_SKZ <- names(valid_SKZ[j])
		val_BIP <- names(valid_BIP[j])
		og_SKZ <- readRDS(paste0("OG/PGS_estime_SKZ_", val_SKZ,".Rdata"))
		og_BIP <- readRDS(paste0("OG/PGS_estime_BIP_", val_BIP,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Lassosum_og", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = cor(simule_SKZ, og_SKZ), "BIP" = cor(simule_BIP, og_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(og_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(og_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "Lassosum_og", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))

	}

	#PANPRS
	valid <- readRDS("PANPRS/Valeurs_f_lambda_PANPRS.Rdata")
	validmax <- which.max(valid)
	#There is a lot of parameters, we keep the 3 best.
	order <- sort(valid, decreasing=T)[1:3]
	for(j in 1:length(order)){
                val <- names(order[j])
		order_j <- which(names(order) == val)
		PANPRS_SKZ <- readRDS(paste0("PANPRS/PGS_estime_SKZ_", val,".Rdata"))
		PANPRS_BIP <- readRDS(paste0("PANPRS/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "PANPRS", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, PANPRS_SKZ), "BIP" = cor(simule_BIP, PANPRS_BIP)))
		roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(PANPRS_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
		roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(PANPRS_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
		AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "PANPRS", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
	}

	#LDpred2
	LDpred2_SKZ <- readRDS("LDpred2/PGS_estime_SKZ_LDpred2_prodMat.Rdata")
	LDpred2_BIP <- readRDS("LDpred2/PGS_estime_BIP_LDpred2_prodMat.Rdata")
	Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "LDpred2", "Param" = NA, "RankPseudoValSKZ" = 1,"RankPseudoValBIP" = 1, "SKZ" = cor(simule_SKZ, LDpred2_SKZ), "BIP" = cor(simule_BIP, LDpred2_BIP)))
	roc_SKZ <- roc(direction = "<", levels = c(0, 1), pheno_SKZ, as.numeric(LDpred2_SKZ)); auc_SKZ <- as.numeric(auc(roc_SKZ))
	roc_BIP <- roc(direction = "<", levels = c(0, 1), pheno_BIP, as.numeric(LDpred2_BIP)); auc_BIP <- as.numeric(auc(roc_BIP))
	AUC <- rbind(AUC, data.frame("Simulation" = k, "Methode" = "LDpred2", "Param" = NA, "RankPseudoValSKZ" = 1,"RankPseudoValBIP" = 1, "SKZ" = auc_SKZ, "BIP" = auc_BIP))
}

saveRDS(Cor, paste0(path, "Correlations_B-LDX.Rdata"))
saveRDS(AUC, paste0(path, "AUC_B-LDX.Rdata"))

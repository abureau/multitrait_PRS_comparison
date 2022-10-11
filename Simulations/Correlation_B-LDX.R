Cor <- data.frame()
for(k in 1:20){
        setwd(paste0(".../Simulation_", k, "/GenCov/"))

	#Simulated PRS
	simule_SKZ <- readRDS("PGS_simule_SKZ.Rdata")
	simule_BIP <- readRDS("PGS_simule_BIP.Rdata")

	#Thresholding
	valid_SKZ <- readRDS("Valeurs_f_lambda_thresholding_SKZ.Rdata")
	valid_BIP <- readRDS("Valeurs_f_lambda_thresholding_BIP.Rdata")
	order_SKZ <- sort(valid_SKZ)
	order_SKZ <- order_SKZ[length(order_SKZ):1]
	order_BIP <- sort(valid_BIP)
	order_BIP <- order_BIP[length(order_BIP):1]
	for(j in 1:length(valid_SKZ)){
		order_SKZ_j <- which(order_SKZ == valid_SKZ[j])
		order_BIP_j <- which(order_BIP == valid_BIP[j])
		val_SKZ <- names(valid_SKZ[j])
		val_BIP <- names(valid_BIP[j])
		th_SKZ <- readRDS(paste0("PGS_thresholding_SKZ_", val_SKZ,".Rdata"))
		th_BIP <- readRDS(paste0("PGS_thresholding_BIP_", val_BIP,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Thresholding", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = cor(simule_SKZ, th_SKZ), "BIP" = cor(simule_BIP, th_BIP)))
	}

	#MultivariateLassosum
	valid <- readRDS("Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid)
	order <- order[length(order):1]
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		multi_SKZ <- readRDS(paste0("PGS_estime_SKZ_", val,".Rdata"))
		multi_BIP <- readRDS(paste0("PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "multiLassosum", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, multi_SKZ), "BIP" = cor(simule_BIP, multi_BIP)))
	}

	#multivariateLassosum Genetic cov
	valid <- readRDS("GenCov/Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid)
	order <- order[length(order):1]
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		gencov_SKZ <- readRDS(paste0("GenCov/PGS_estime_SKZ_", val,".Rdata"))
		gencov_BIP <- readRDS(paste0("GenCov/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Genetic_cov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, gencov_SKZ), "BIP" = cor(simule_BIP, gencov_BIP)))
	}


	#multivariateLassosum adaptive (Beta from multivariateLassosum Genetic Cov)
	valid <- readRDS("BetaMulti/Valeurs_f_lambda_LassosumExtension.Rdata")
	order <- sort(valid)
	order <- order[length(order):1]
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		bm_SKZ <- readRDS(paste0("BetaMulti/PGS_estime_SKZ_", val,".Rdata"))
		bm_BIP <- readRDS(paste0("BetaMulti/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Adap_BetaMulti", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, bm_SKZ), "BIP" = cor(simule_BIP, bm_BIP)))
	}

	#multivariateLassosum adaptive (Beta from multivariateLassosum standard)
	valid <- readRDS("BetaGenCov/Valeurs_f_lambda_LassosumExtension.Rdata"))
	order <- sort(valid)
	order <- order[length(order):1]
	for(j in 1:length(valid)){
		order_j <- which(order == valid[j])
		val <- strsplit(x = names(valid[j]), split = "_", fixed = TRUE)[[1]][2]
		adap_gencov_SKZ <- readRDS(paste0("BetaGenCov/PGS_estime_SKZ_", val,".Rdata"))
		adap_gencov_BIP <- readRDS(paste0("BetaGenCov/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Adap_GenCov", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, adap_gencov_SKZ), "BIP" = cor(simule_BIP, adap_gencov_BIP)))
	}

	#Lassosum
	valid_SKZ <- readRDS("OG/Valeurs_f_lambda_Lassosum_SKZ.Rdata")
	valid_BIP <- readRDS("OG/Valeurs_f_lambda_Lassosum_BIP.Rdata")
	order_SKZ <- sort(valid_SKZ)
	order_SKZ <- order_SKZ[length(order_SKZ):1]
	order_BIP <- sort(valid_BIP)
	order_BIP <- order_BIP[length(order_BIP):1]
	for(j in 1:length(valid_SKZ)){
		order_SKZ_j <- which(order_SKZ == valid_SKZ[j])
		order_BIP_j <- which(order_BIP == valid_BIP[j])
		val_SKZ <- names(valid_SKZ[j])
		val_BIP <- names(valid_BIP[j])
		og_SKZ <- readRDS(paste0("OG/PGS_estime_SKZ_", val_SKZ,".Rdata"))
		og_BIP <- readRDS(paste0("OG/PGS_estime_BIP_", val_BIP,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "Lassosum_og", "Param" = val_SKZ, "RankPseudoValSKZ" = order_SKZ_j, "RankPseudoValBIP" = order_BIP_j, "SKZ" = cor(simule_SKZ, og_SKZ), "BIP" = cor(simule_BIP, og_BIP)))
	}

	#PANPRS
	valid <- readRDS("PANPRS/Valeurs_f_lambda_PANPRS.Rdata")
	validmax <- which.max(valid)
	valid <- valid[!is.nan(valid)]
	order <- sort(valid)
	order <- order[length(order):1]
	order <- order[1:3]
	for(j in 1:length(order)){
		val <- order[j]
                val <- names(val)
		order_j <- which(names(order) == val)
		PANPRS_SKZ <- readRDS(paste0("PANPRS/PGS_estime_SKZ_", val,".Rdata"))
		PANPRS_BIP <- readRDS(paste0("PANPRS/PGS_estime_BIP_", val,".Rdata"))
		Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "PANPRS", "Param" = val, "RankPseudoValSKZ" = order_j, "RankPseudoValBIP" = order_j, "SKZ" = cor(simule_SKZ, PANPRS_SKZ), "BIP" = cor(simule_BIP, PANPRS_BIP)))
	}

	#LDpred2
	LDpred2_SKZ <- readRDS("LDpred2/PGS_estime_SKZ_LDpred2_prodMat.Rdata")
	LDpred2_BIP <- readRDS("LDpred2/PGS_estime_BIP_LDpred2_prodMat.Rdata")
	Cor <- rbind(Cor, data.frame("Simulation" = i, "Methode" = "LDpred2", "Param" = NA, "RankPseudoValSKZ" = 1,"RankPseudoValBIP" = 1, "SKZ" = cor(simule_SKZ, LDpred2_SKZ), "BIP" = cor(simule_BIP, LDpred2_BIP)))
}

saveRDS(Cor, ".../Correlations_B-LDX.Rdata")


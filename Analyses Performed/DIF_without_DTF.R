# =====================================================================
# Monte Carlo Analysis: Non-TIF Conditions — Variable Importance + ANOVAs
# Author: Catherine M. Bain
# Date: 2024-12-02
# Note: This version excludes ALL plotting code.
# =====================================================================

rm(list = ls(all = TRUE))

#################### Working directory ################################
# NOTE: Hard-coded path. Comment for shared/HPC use and rely on project root.
setwd("~/Documents/research/currentProjects/difml")

################## Packages ###########################################
# Keep only packages used for data wrangling and stats.
# Avoid plyr to prevent masking dplyr verbs.
library(dplyr)
library(purrr)
library(tidyr)
library(effectsize)  # for eta_squared()
library(clipr)       # to copy ANOVA effect sizes to clipboard (optional)

########################### Utility ###################################
# cbind.fill(): column-bind objects with different row counts by padding
# NA rows to shorter ones. Useful when merging uneven summaries.
cbind.fill <- function(...) {
	nm <- list(...)
	nm <- lapply(nm, as.matrix)
	n  <- max(sapply(nm, nrow))
	do.call(cbind, lapply(nm, function(x) rbind(x, matrix(, n - nrow(x), ncol(x)))))
}

########################## Load data ##################################
# Expects an .Rda with an object `fullDataExpanded` created within allConditionsAnalysis.R
load("~/Documents/research/currentProjects/difml/results/fullDataExpanded.Rda")

######################### Condition sets ##############################
# "Visible TIF" conditions from prior work. We will analyze the complement.
conditions_visible_tif <- c(
	56,57,58,59,60,69,71,72,73,74,75,76,77,78,110,111,112,113,114,123,125,126,127,128,129,130,131,132,
	164,165,166,167,168,177,179,180,181,182,183,184,185,186,218,219,220,221,222,231,233,234,235,236,237,238,239,240,
	272,273,274,275,276,285,287,288,289,290,291,292,293,294,326,327,328,329,330,339,341,342,343,344,345,346,347,348,
	380,381,382,383,384,393,395,396,397,398,399,400,401,402,434,435,436,437,438,447,449,450,451,452,453,454,455,456,
	488,489,490,491,492,501,503,504,505,506,507,508,509,510,542,543,544,545,546,555,557,558,559,560,561,562,563,564,
	596,597,598,599,600,609,611,612,613,614,615,616,617,618,650,651,652,653,654,663,665,666,667,668,669,670,671,672,
	704,705,706,707,708,717,719,720,721,722,723,724,725,726,758,759,760,761,762,771,773,774,775,776,777,778,779,780,
	812,813,814,815,816,825,827,828,829,830,831,832,833,834,866,867,868,869,870,879,881,882,883,884,885,886,887,888,
	920,921,922,923,924,933,935,936,937,938,939,940,941,942,974,975,976,977,978,987,989,990,991,992,993,994,995,996,
	1028,1029,1030,1031,1032,1041,1043,1044,1045,1046,1047,1048,1049,1050,1082,1083,1084,1085,1086,1095,1097,1098,1099,
	1100,1101,1102,1103,1104,1136,1137,1138,1139,1140,1149,1151,1152,1153,1154,1155,1156,1157,1158,1190,1191,1192,1193,
	1194,1203,1205,1206,1207,1208,1209,1210,1211,1212,1244,1245,1246,1247,1248,1257,1259,1260,1261,1262,1263,1264,1265,
	1266,1298,1299,1300,1301,1302,1311,1313,1314,1315,1316,1317,1318,1319,1320,1331,1332,1337,1338,1349,1350,1355,1356,
	1367,1368,1373,1374,1385,1386,1391,1392,1403,1404,1409,1410,1421,1422,1427,1428,1439,1440,1445,1446,1457,1458,1463,
	1464,1475,1476,1481,1482,1493,1494,1499,1500,1511,1512,1517,1518,1529,1530,1535,1536,1547,1548,1553,1554,1565,1566,
	1571,1572,1583,1584,1589,1590,1601,1602,1607,1608,1619,1620,1625,1626,1637,1638,1643,1644,1655,1656,1661,1662,1673,
	1674,1679,1680,1691,1692,1697,1698,1709,1710,1715,1716,1727,1728,1733,1734,1745,1746,1751,1752,1770,1788,1806,1824,
	1842,1860,1878,1896,1914,1932,1950,1968,1986,2004,2022,2040,2058,2076,2094,2112,2130,2148,2166,2184
)

# Non-TIF conditions to analyze here
conditions_non_tif <- setdiff(25:2184, conditions_visible_tif)

######################### Variable importance #########################
# For each non-TIF condition:
# 1) Read condition-level results to get DIF item indices (from rf_class row).
# 2) Read variable-importance matrix (columns V1..VJ).
# 3) Rank variables, compute top 5% list, normalize importances,
#    then average normalized importance across DIF vs non-DIF items.

allImportance <- vector()
for (i in conditions_non_tif) {
	message("Processing file: ", i)
	
	# ---- Results: DIF list for this condition (rf_class row only) -----
	file_path_results <- sprintf("results/final/condition%dFinalResults.csv", i)
	results_data <- read.csv(file_path_results, header = TRUE)
	colnames(results_data)[11:12] <- c("difItems", "method")
	results_data <- results_data %>%
		dplyr::select(difItems, method) %>%
		dplyr::filter(method == "rf_class") %>%
		mutate(condition = i) %>%
		dplyr::select(-method)
	
	# ---- Variable importance: columns V1..VJ, one row per model fit ----
	file_path_var <- sprintf("results/finalVarImportance/condition%dFinalVarImp.csv", i)
	var_ranking <- read.csv(file_path_var, header = TRUE)
	
	# Rank variables and compute top-5% set; keep per-row totals for normalization
	var_sorted <- var_ranking %>%
		mutate(
			rankedVarImp = pmap_chr(., ~{
				vars <- c(...)
				ranked_vars <- gsub("V", "", names(vars)[order(-as.numeric(vars))])
				paste(ranked_vars, collapse = ", ")
			}),
			topVars = pmap_chr(., ~{
				vars <- c(...)
				ranked_vars <- gsub("V", "", names(vars)[order(-as.numeric(vars))])
				n_vars <- length(ranked_vars)
				top_n <- ceiling(n_vars * 0.05)
				paste(ranked_vars[seq_len(top_n)], collapse = ", ")
			})
		) %>%
		mutate(
			topVars_list  = lapply(strsplit(as.character(topVars), ",\\s*"), as.numeric),
			difItems_list = lapply(strsplit(as.character(results_data$difItems), ",\\s*"), as.numeric),
			propTopDIF = purrr::map2_dbl(
				topVars_list, difItems_list,
				~ if (length(.x) == 0) NA_real_ else sum(.x %in% .y) / length(.x)
			),
			totalImp = rowSums(dplyr::select(., matches("^V\\d+$")), na.rm = TRUE)
		)
	
	# Normalize row-wise variable importance and average by DIF vs non-DIF
	prop_df <- var_sorted %>% dplyr::select(matches("^V\\d+$"))
	prop_df_normalized <- prop_df / var_sorted$totalImp
	colnames(prop_df_normalized) <- gsub("V", "", colnames(prop_df_normalized))
	
	dif_average    <- numeric(nrow(prop_df_normalized))
	nonDIF_average <- numeric(nrow(prop_df_normalized))
	
	for (r in seq_len(nrow(prop_df_normalized))) {
		cols_to_include <- var_sorted$difItems_list[[r]]
		if (length(cols_to_include) == 0 || all(is.na(cols_to_include))) {
			# No DIF items listed; record NA for DIF mean and overall mean for non-DIF
			dif_average[r]    <- NA_real_
			nonDIF_average[r] <- mean(as.numeric(prop_df_normalized[r, , drop = FALSE]), na.rm = TRUE)
			next
		}
		dif_cols <- intersect(colnames(prop_df_normalized), as.character(cols_to_include))
		non_cols <- setdiff(colnames(prop_df_normalized), dif_cols)
		
		dif_subset     <- prop_df_normalized[r, dif_cols,  drop = FALSE]
		non_dif_subset <- prop_df_normalized[r, non_cols, drop = FALSE]
		
		dif_average[r]    <- mean(as.numeric(dif_subset),     na.rm = TRUE)
		nonDIF_average[r] <- mean(as.numeric(non_dif_subset), na.rm = TRUE)
	}
	
	var_sorted <- cbind(var_sorted, dif_average, nonDIF_average) %>%
		dplyr::select(topVars, propTopDIF, dif_average, nonDIF_average)
	
	# Combine DIF list with summarized importance metrics
	file_data <- cbind.fill(results_data, var_sorted)
	allImportance <- rbind(allImportance, file_data)
}

# Tidy numeric types for downstream use
fullDataNumeric <- data.frame(allImportance) %>%
	mutate(
		propTopDIF     = as.numeric(propTopDIF),
		dif_average    = as.numeric(dif_average),
		nonDIF_average = as.numeric(nonDIF_average),
		condition      = as.numeric(condition)
	)

########################### Accuracy variables ########################
# Keep only non-TIF (non-visible) DIF conditions for accuracy models.
nonVisibleDIF <- subset(fullDataExpanded, Condition %in% conditions_non_tif)

# Convenience fields and method filter for ANOVA models
noRFprob <- nonVisibleDIF %>%
	filter(method != "RF Probabilities") %>%
	mutate(
		totalSample = RefGroupSize + FocalGroupSize,
		FocRefRatio = FocalGroupSize / RefGroupSize
	)

# Long format for marginal means (no plots produced; keep for QA/export)
fullDataLong <- nonVisibleDIF %>%
	dplyr::select(-c(starts_with("TN"), starts_with("TP"), starts_with("FP"), starts_with("FN"))) %>%
	dplyr::select(method, Condition, ends_with("full"), ends_with("ref"), ends_with("focal")) %>%
	tidyr::pivot_longer(
		cols = -c(method, Condition),
		names_to = c("metric", "group"),
		names_pattern = "(.*)_(full|ref|focal)$",
		values_to = "Value"
	)

marginalMeans <- fullDataLong %>%
	group_by(method, group, metric) %>%
	dplyr::summarise(
		mean_out = mean(Value, na.rm = TRUE),
		se_out   = sd(Value,   na.rm = TRUE),
		.groups  = "drop"
	)

############################## ANOVAs: Full Sample ####################
# Important: The original formula referenced `Sample`, which is not a field
# in `noRFprob`. Remove `Sample` terms so models run on wide-format data.

# Accuracy
accANOVA <- lm(
	accuracy_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(accANOVA, partial = TRUE); clipr::write_clip(e)

# Sensitivity
senANOVA_full <- lm(
	sensitivity_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(senANOVA_full, partial = TRUE); clipr::write_clip(e)

# Specificity
specANOVA_full <- lm(
	specificity_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(specANOVA_full, partial = TRUE); clipr::write_clip(e)

# Precision
preANOVA_full <- lm(
	precision_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(preANOVA_full, partial = TRUE); clipr::write_clip(e)

# Negative predictive value
npvANOVA_full <- lm(
	npv_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(npvANOVA_full, partial = TRUE); clipr::write_clip(e)

# F1 score
f1ANOVA_full <- lm(
	f1score_full ~ (method + totalSample + nItems + prevRate)^2 +
		method * (FocRefRatio + difRate + balanced + aStrength + bStrength)^2,
	data = noRFprob
)
e <- eta_squared(f1ANOVA_full, partial = TRUE); clipr::write_clip(e)

# ============================== Outputs ==============================
# Objects created for downstream reporting or export (no plots):
# - fullDataNumeric: DIF vs non-DIF variable-importance summaries per condition
# - marginalMeans:   mean and SE by method × group × metric
# - accANOVA, senANOVA_full, specANOVA_full, preANOVA_full, npvANOVA_full, f1ANOVA_full
#   with partial eta-squared copied to clipboard via clipr::write_clip()
# =====================================================================

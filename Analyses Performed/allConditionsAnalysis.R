# =====================================================================
# Monte Carlo Analysis: Summarize and Model Simulation Results
# Author: Catherine M. Bain
# Date: 2024-12-02
# =====================================================================

rm(list = ls(all = TRUE))

#################### Working directory ################################
# Set to project root containing `results/` and `allConditions.csv`.

################## Packages ###########################################
# Visualization, wrangling, models, parallel, and reporting
library(ggplot2)
library(dplyr)
library(ggpubr)
library(effectsize)
library(tidyr)
library(foreach)
library(doParallel)
library(jmv)          # for ANOVA()
library(interactions) # not used below; keep if you add interaction plots
# NOTE: car and clipr are used with :: so no library() needed.

########################### Helper: interaction line plot ##############
# plot_int()
# Inputs:
#   data   : data.frame
#   iv     : x-axis variable (unquoted tidy-eval)
#   group  : linetype grouping variable (unquoted)
#   outcome: numeric outcome (unquoted)
# Output: ggplot with mean ± SE by iv × group
plot_int <- function(data, iv, group, outcome) {
	plot_data <- data %>%
		group_by({{ group }}, {{ iv }}) %>%
		summarise(
			mean_out = mean({{ outcome }}, na.rm = TRUE),
			se_out   = sd({{ outcome }},   na.rm = TRUE) / sqrt(n()),
			.groups = "drop"
		)
	
	ggplot(
		plot_data,
		aes(x = factor({{ iv }}), y = mean_out, linetype = {{ group }}, group = {{ group }})
	) +
		geom_point(size = 3) +
		geom_line() +
		geom_errorbar(aes(ymin = mean_out - se_out, ymax = mean_out + se_out), width = 0.2) +
		theme_minimal()
}

################### Load all condition result files in parallel ########
# Expects files:
#   results/final/condition<id>FinalResults.csv for id = 1..2184
# Returns a single data.frame bound across conditions.
fullData <- vector(mode = "list", length = 0)

# Configure cluster
n.cores <- 9
cl <- makeCluster(n.cores, type = "SOCK")
registerDoParallel(cl, cores = n.cores)

fullData <- foreach(
	i = 1:2184,
	.packages = c("dplyr"),
	.combine = rbind,
	.inorder = FALSE
) %dopar% {
	file_path <- sprintf("results/final/condition%dFinalResults.csv", i)
	
	if (!file.exists(file_path)) return(NULL)
	
	# Read and standardize column names per file, tag with `condition = i`
	tryCatch(
		{
			data <- read.csv(file_path, header = TRUE)
			# Align to expected 41-column schema produced by the simulation
			colnames(data) <- c(
				"RefGroupSize", "FocalGroupSize", "nItems", "prevRate", "difRate",
				"difType", "balanced", "aStrength", "bStrength", "rep", "difItems",
				"method",
				"TP_full", "TN_full", "FP_full", "FN_full", "accuracy_full",
				"sensitivity_full", "specificity_full", "precision_full", "npv_full",
				"f1score_full",
				"method_focal",
				"TP_focal", "TN_focal", "FP_focal", "FN_focal", "accuracy_focal",
				"sensitivity_focal", "specificity_focal", "precision_focal", "npv_focal",
				"f1score_focal",
				"method_ref",
				"TP_ref", "TN_ref", "FP_ref", "FN_ref", "accuracy_ref",
				"sensitivity_ref", "specificity_ref", "precision_ref", "npv_ref",
				"f1score_ref"
			)
			dplyr::mutate(data, condition = i)
		},
		error = function(e) NULL
	)
}

stopCluster(cl)

################### Column pruning and renaming ########################
# Drop duplicated method tag columns carried for focal/ref blocks.
# Original code used positional indices 11,23,34. Keep but document.
# NOTE: This assumes the 41-col schema above. If upstream changes, adjust.
fullData <- fullData[, -c(11, 23, 34)]

colnames(fullData) <- c(
	"RefGroupSize", "FocalGroupSize", "nItems", "prevRate", "difRate",
	"difType", "balanced", "aStrength", "bStrength", "rep",
	"method",
	"TP_full", "TN_full", "FP_full", "FN_full", "accuracy_full",
	"sensitivity_full", "specificity_full", "precision_full", "npv_full",
	"f1score_full",
	"TP_focal", "TN_focal", "FP_focal", "FN_focal", "accuracy_focal",
	"sensitivity_focal", "specificity_focal", "precision_focal", "npv_focal",
	"f1score_focal",
	"TP_ref", "TN_ref", "FP_ref", "FN_ref", "accuracy_ref",
	"sensitivity_ref", "specificity_ref", "precision_ref", "npv_ref",
	"f1score_ref",
	"Condition"
)

################### Derived metrics and basic QC flags #################
fullDataExpanded <- fullData %>%
	mutate(
		# Reference minus focal; positive means ref performance > focal
		relAccuracy     = accuracy_ref     - accuracy_focal,
		relSensitivity  = sensitivity_ref  - sensitivity_focal,
		relSpecificity  = specificity_ref  - specificity_focal,
		relPrecision    = precision_ref    - precision_focal,
		relNPV          = npv_ref          - npv_focal,
		relF1           = f1score_ref      - f1score_focal
	) %>%
	mutate(
		# QC: any predicted positive in the overall set?
		AtLeastOneInClass1 = (TP_full + FP_full) > 0
	) %>%
	group_by(Condition) %>%
	mutate(
		# Keep conditions where at least half of reps yielded any positives
		Keep = mean(AtLeastOneInClass1) >= 0.5
	) %>%
	ungroup()

# Convert some design columns to factors.
# NOTE: Using positions (6,7,11) is fragile; explicit names are safer.
fullDataExpanded <- fullDataExpanded %>%
	mutate(across(c(difType, balanced, method), as.factor))

################### Tidy labels for method #############################
fullDataExpanded$method <- factor(
	fullDataExpanded$method,
	levels = c("IRT", "rf_probs", "rf_class"),
	labels = c("IRT", "RF Probabilities", "RF Classification")
)

########################### Baseline (no DIF) visuals ##################
# Filter no-DIF subset
fullDataExpandedNoDIF <- subset(fullDataExpanded, difRate == 0)

# Keep only overall metrics for plotting
fullDataExpandedUseful <- fullDataExpandedNoDIF %>%
	dplyr::select(-c(FN_full, FP_full, TN_full, TP_full))

# Method means (optional quick summary)
rfprobMeans <- subset(fullDataExpandedUseful, method == "RF Probabilities") %>%
	summarise(
		accuracy    = mean(accuracy_full,    na.rm = TRUE),
		sensitivity = mean(sensitivity_full, na.rm = TRUE),
		specificity = mean(specificity_full, na.rm = TRUE),
		precision   = mean(precision_full,   na.rm = TRUE),
		npv         = mean(npv_full,         na.rm = TRUE),
		f1score     = mean(f1score_full,     na.rm = TRUE)
	)

# Long format for facetting
fullDataLong <- fullDataExpandedUseful %>%
	select(method, ends_with("full")) %>%
	pivot_longer(cols = -method, names_to = "metric", values_to = "Overall")

# Faceted violin+boxplots for all overall metrics
allMetricsPlot <- ggplot(fullDataLong, aes(x = method, y = Overall, fill = method)) +
	geom_violin(alpha = 0.8, aes(color = method)) +
	geom_boxplot(alpha = 0, width = 0.1) +
	facet_wrap(~ metric, scales = "free") +
	labs(title = "Performance Metrics (Overall, No DIF)", y = "Overall") +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(allMetricsPlot)

############################ ANOVAs: no DIF ###########################
# If the next CSV is not precomputed, uncomment the write step below.
# noRFprobnoDIF <- subset(fullDataExpandedNoDIF, method != "RF Probabilities")
# write.csv(noRFprobnoDIF, file = "results/combinedConditions/noDIFnoRFprobs.csv", row.names = FALSE)

noRFprobnoDIF <- read.csv("results/combinedConditions/noDIFnoRFprobs.csv", header = TRUE)

# Factorize design fields and add Focal/Ref ratio
noRFprobnoDIF <- noRFprobnoDIF %>%
	mutate(across(c(difType, balanced, method), as.factor)) %>%
	mutate(FocRefRatio = FocalGroupSize / RefGroupSize)

# ---------- Focal group: linear models with 3-way interactions --------
accANOVA_focalnoDIF <- lm(accuracy_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(accANOVA_focalnoDIF, partial = TRUE))

senANOVA_focalnoDIF <- lm(sensitivity_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(senANOVA_focalnoDIF, partial = TRUE))

specANOVA_focalnoDIF <- lm(specificity_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(specANOVA_focalnoDIF, partial = TRUE))

precisionANOVA_focalnoDIF <- lm(precision_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(precisionANOVA_focalnoDIF, partial = TRUE))

npvANOVA_focalnoDIF <- lm(npv_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(npvANOVA_focalnoDIF, partial = TRUE))

f1ANOVA_focalnoDIF <- lm(f1score_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(f1ANOVA_focalnoDIF, partial = TRUE))

# jmv::ANOVA examples with partial eta and EMMs
senANOVA_focalnoDIF <- ANOVA(
	formula = sensitivity_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
senANOVA_focalnoDIF

precANOVA_focalnoDIF <- ANOVA(
	formula = precision_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
precANOVA_focalnoDIF

npvANOVA_focalnoDIF <- ANOVA(
	formula = npv_focal ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method:prevRate, emmTables = TRUE, ss = "2"
)
npvANOVA_focalnoDIF

# Focal specificity interaction plot
spec_int <- plot_int(noRFprobnoDIF, prevRate, method, specificity_focal) +
	scale_linetype_manual(values = c("solid", "dashed"), labels = c("IRT", "RF Classification")) +
	labs(
		x = "Prevalence Rate",
		y = "Mean Specificity of Focal Group",
		linetype = "Classification Method",
		title = "Interaction: Prevalence × Method on Specificity (Focal)"
	) +
	ylim(0, 1)

######################## Reference group analyses ######################
accANOVA_refnoDIF <- lm(accuracy_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(accANOVA_refnoDIF, partial = TRUE))

senANOVA_refnoDIF <- lm(sensitivity_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(senANOVA_refnoDIF, partial = TRUE))

specANOVA_refnoDIF <- lm(specificity_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(specANOVA_refnoDIF, partial = TRUE))

precisionANOVA_refnoDIF <- lm(precision_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(precisionANOVA_refnoDIF, partial = TRUE))

npvANOVA_refnoDIF <- lm(npv_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(npvANOVA_refnoDIF, partial = TRUE))

f1ANOVA_refnoDIF <- lm(f1score_ref ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(f1ANOVA_refnoDIF, partial = TRUE))

# jmv::ANOVA examples (reduced formulas)
senANOVA_refnoDIF <- ANOVA(
	formula = sensitivity_ref ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
senANOVA_refnoDIF

precANOVA_refnoDIF <- ANOVA(
	formula = precision_ref ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
precANOVA_refnoDIF

npvANOVA_refnoDIF <- ANOVA(
	formula = npv_ref ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
npvANOVA_refnoDIF

# Reference specificity interaction plot
spec_int <- plot_int(noRFprobnoDIF, prevRate, method, specificity_ref) +
	scale_linetype_manual(values = c("solid", "dashed"), labels = c("IRT", "RF Classification")) +
	labs(
		x = "Prevalence Rate",
		y = "Mean Specificity of Reference Group",
		linetype = "Classification Method",
		title = "Interaction: Prevalence × Method on Specificity (Reference)"
	) +
	ylim(0, 1)

######################## Overall analyses: no DIF ######################
accANOVA_fullnoDIF <- lm(accuracy_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(accANOVA_fullnoDIF, partial = TRUE))

senANOVA_fullnoDIF <- lm(sensitivity_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(senANOVA_fullnoDIF, partial = TRUE))

specANOVA_fullnoDIF <- lm(specificity_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(specANOVA_fullnoDIF, partial = TRUE))

precisionANOVA_fullnoDIF <- lm(precision_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(precisionANOVA_fullnoDIF, partial = TRUE))

npvANOVA_fullnoDIF <- lm(npv_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(npvANOVA_fullnoDIF, partial = TRUE))

f1ANOVA_fullnoDIF <- lm(f1score_full ~ (FocalGroupSize + RefGroupSize + FocRefRatio + nItems + prevRate + method)^3, data = noRFprobnoDIF)
clipr::write_clip(eta_squared(f1ANOVA_fullnoDIF, partial = TRUE))

# jmv::ANOVA examples (reduced formulas)
senANOVA_fullnoDIF <- ANOVA(
	formula = sensitivity_full ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
senANOVA_fullnoDIF

precANOVA_fullnoDIF <- ANOVA(
	formula = precision_full ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
precANOVA_fullnoDIF

npvANOVA_fullnoDIF <- ANOVA(
	formula = npv_full ~ (FocalGroupSize + nItems + prevRate + method)^3,
	data = noRFprobnoDIF,
	effectSize = "partEta", emMeans = ~ method, emmTables = TRUE, ss = "2"
)
npvANOVA_fullnoDIF

# Overall specificity interaction plot
spec_int <- plot_int(noRFprobnoDIF, prevRate, method, specificity_full) +
	scale_linetype_manual(values = c("solid", "dashed"), labels = c("IRT", "RF Classification")) +
	labs(
		x = "Prevalence Rate",
		y = "Mean Specificity (Overall)",
		linetype = "Classification Method",
		title = "Interaction: Prevalence × Method on Specificity (Overall)"
	) +
	ylim(0, 1)

########################### DIF-only violin panels ####################
# FIX: treat difRate as numeric (was compared to character '0')
graphDataDIF <- subset(fullData, difRate != 0)

accuracyDIF <- ggplot(graphDataDIF, aes(x = method, y = accuracy_focal, color = method, fill = method)) +
	geom_violin(alpha = 0.8) +
	geom_boxplot(alpha = 0) +
	ylim(0, 1) +
	theme_minimal()

sensitivityDIF <- ggplot(graphDataDIF, aes(x = method, y = sensitivity_focal, color = method, fill = method)) +
	geom_boxplot(alpha = 0.8) +
	ylim(0, 1) +
	theme_minimal()

specificityDIF <- ggplot(graphDataDIF, aes(x = method, y = specificity_focal, color = method, fill = method)) +
	geom_violin(alpha = 0.8) +
	ylim(0, 1) +
	theme_minimal()

precisionDIF <- ggplot(graphDataDIF, aes(x = method, y = precision_focal, color = method, fill = method)) +
	geom_violin(alpha = 0.8) +
	ylim(0, 1) +
	theme_minimal()

npvDIF <- ggplot(graphDataDIF, aes(x = method, y = npv_focal, color = method, fill = method)) +
	geom_violin(alpha = 0.8) +
	ylim(0, 1) +
	theme_minimal()

f1DIF <- ggplot(graphDataDIF, aes(x = method, y = f1score_focal, color = method, fill = method)) +
	geom_violin(alpha = 0.8) +
	ylim(0, 1) +
	theme_minimal()

ggarrange(accuracyDIF, sensitivityDIF, specificityDIF, precisionDIF, npvDIF, f1DIF,
		  ncol = 3, nrow = 2, legend = "none")

########################### Combined boxplots by group #################
# Wide-to-long using name capture:
#  metric_full -> metric + Group="full", etc.
longData <- fullData %>%
	pivot_longer(
		cols = ends_with("_full") | ends_with("_focal") | ends_with("_ref"),
		names_to = c(".value", "Group"),
		names_pattern = "^(.*)_(full|focal|ref)$"
	)

# Harmonize group labels
longData$Group <- factor(
	longData$Group,
	levels = c("focal", "ref", "full"),
	labels = c("focal", "reference", "overall")
)

# DIF-only slice for method × group comparisons
difOnly <- subset(longData, difRate != 0)

# Median overlays (optional)
median_data <- difOnly %>%
	group_by(method, Group, difRate) %>%
	summarize(median_accuracy = median(accuracy, na.rm = TRUE), .groups = "drop")

# Accuracy
a <- ggplot(difOnly, aes(x = Group, y = accuracy, color = method)) +
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

# Sensitivity
sen <- ggplot(difOnly, aes(x = Group, y = sensitivity, color = method)) +
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

# Specificity
spec <- ggplot(difOnly, aes(x = Group, y = specificity, color = method)) +
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

# Precision
pre <- ggplot(difOnly, aes(x = Group, y = precision, color = method)) +
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

# NPV
npvAll <- ggplot(difOnly, aes(x = Group, y = npv, color = method)) +  # FIX: define as npvAll
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

# F1
f1All <- ggplot(difOnly, aes(x = Group, y = f1score, color = method)) +
	geom_violin(alpha = 0.8, draw_quantiles = 0.5) +
	facet_wrap(~ difRate, scales = "free") +
	ylim(0, 1) +
	theme_minimal() +
	theme(strip.text = element_text(size = 14, face = "bold"))

ggarrange(a, sen, spec, pre, npvAll, f1All,
		  ncol = 3, nrow = 2, legend = "bottom", common.legend = TRUE)

################################## Accuracy ANOVA ######################
# Motivation: compare methods while adjusting for design factors.
# Use Type-II/III SS when order should not matter.
noRFprob <- subset(fullData, method != "rf_probs")

# Overall accuracy with additive predictors
accANOVA <- lm(accuracy_full ~ FocalGroupSize + nItems + prevRate + difRate +
			   	balanced + aStrength + bStrength + method,
			   data = noRFprob)
summary(accANOVA)
eta_squared(accANOVA, partial = TRUE)

# Focal accuracy; Type-II SS via car::Anova
accANOVA_focal <- car::Anova(
	lm(accuracy_focal ~ FocalGroupSize + nItems + prevRate + difRate +
	   	balanced + aStrength + bStrength + method, data = fullData),
	type = 2
)
eta_squared(accANOVA_focal, partial = TRUE)

# Restrict to DIF>0 and drop RF probabilities
noRFprobDIF <- subset(fullData, method != "rf_probs" & difRate != 0)

accANOVA_focalDIF <- lm(accuracy_focal ~ FocalGroupSize + nItems + prevRate +
							difRate + balanced + aStrength + bStrength + method,
						data = noRFprobDIF)
eta_squared(accANOVA_focalDIF, partial = TRUE)

senANOVA_focalDIF <- lm(sensitivity_focal ~ FocalGroupSize + nItems + prevRate +
							difRate + balanced + aStrength + bStrength + method,
						data = noRFprobDIF)
eta_squared(senANOVA_focalDIF, partial = TRUE)

specANOVA_focalDIF <- lm(specificity_focal ~ FocalGroupSize + nItems + prevRate +
						 	difRate + balanced + aStrength + bStrength + method,
						 data = noRFprobDIF)
eta_squared(specANOVA_focalDIF, partial = TRUE)

precisionANOVA_focalDIF <- lm(precision_focal ~ FocalGroupSize + nItems + prevRate +
							  	difRate + balanced + aStrength + bStrength + method,
							  data = noRFprobDIF)
eta_squared(precisionANOVA_focalDIF, partial = TRUE)

npvANOVA_focalDIF <- lm(npv_focal ~ FocalGroupSize + nItems + prevRate +
							difRate + balanced + aStrength + bStrength + method,
						data = noRFprobDIF)
eta_squared(npvANOVA_focalDIF, partial = TRUE)

f1ANOVA_focalDIF <- lm(f1score_focal ~ FocalGroupSize + nItems + prevRate +
					   	difRate + balanced + aStrength + bStrength + method,
					   data = noRFprobDIF)
eta_squared(f1ANOVA_focalDIF, partial = TRUE)


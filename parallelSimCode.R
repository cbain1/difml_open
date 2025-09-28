# ====================================================================
# Monte Carlo Simulation: DIF Effects on Test Information Function
# --------------------------------------------------------------------
# Purpose
#   Simulate polytomous IRT data with DIF, fit IRT and Random Forest
#   classifiers for a binary diagnosis, and export performance metrics.
#
# Assumptions
#   - Response categories: 0..4 (5-category graded response).
#   - diag_classes is binary {0,1}; downstream models use factor.
#   - allConditions.csv exists with columns described below.
#
# Parallel
#   - Outer: foreach dopar over replications.
#   - Inner: optional small clusters in item-prob and response generation.
#
# Outputs
#   - results/condition<id>FinalResults.csv
#   - results/condition<id>FinalVarImp.csv
#
# Dependencies
#   mirt, mclust, pROC, randomForest, dplyr, MASS, caret, foreach,
#   doParallel, parallel
# ====================================================================

rm(list = ls(all = TRUE))

################## Packages ###########################################
suppressPackageStartupMessages({
	library(mirt)
	library(mclust)
	library(pROC)
	library(randomForest)
	library(dplyr)
	library(MASS)
	library(caret)
	library(foreach)
	library(doParallel)
	library(parallel)
})

################## Reproducibility ####################################
# Set a single top-level seed for reproducibility across runs.
# Reviewers can change this one constant to re-generate all.
TOP_LEVEL_SEED <- 20250928
set.seed(TOP_LEVEL_SEED)

################## Helpers and Inner Functions ########################

#' Generate cumulative category probabilities for graded response items
#'
#' @param N Integer. Number of persons.
#' @param I Integer. Number of items.
#' @param theta_ass Numeric vector length N. Assessment ability.
#' @param a Numeric vector length I. Discrimination.
#' @param bMat Numeric I x 4 matrix of ordered thresholds b1<b2<b3<b4.
#' @param cores Integer >=1. Worker count for inner parallel.
#'
#' @return data.frame with columns:
#'   person, item, c0..c4 (cumulative probabilities for scores 0..4).
#' @details
#'   Uses the graded response model. For each person-item, returns
#'   cumulative category probabilities to enable fast sampling.
generate_cprobs_par <- function(N, I, theta_ass, a, bMat, cores = 1) {
	# Inner cluster is optional; keep small to avoid oversubscription.
	cl <- makeCluster(cores)
	on.exit(stopCluster(cl), add = TRUE)
	registerDoParallel(cl)
	
	combos <- expand.grid(j = seq_len(N), i = seq_len(I))
	
	out <- foreach(idx = seq_len(nrow(combos)), .combine = 'rbind') %dopar% {
		j <- combos$j[idx]; i <- combos$i[idx]
		z1 <- a[i] * (theta_ass[j] - bMat[i, 1])
		z2 <- a[i] * (theta_ass[j] - bMat[i, 2])
		z3 <- a[i] * (theta_ass[j] - bMat[i, 3])
		z4 <- a[i] * (theta_ass[j] - bMat[i, 4])
		
		# Cumulative GRM category probabilities
		p1 <- 1 - plogis(z1)
		p2 <- plogis(z1) - plogis(z2)
		p3 <- plogis(z2) - plogis(z3)
		p4 <- plogis(z3) - plogis(z4)
		p5 <- plogis(z4)
		
		c0 <- p1
		c1 <- p1 + p2
		c2 <- p1 + p2 + p3
		c3 <- p1 + p2 + p3 + p4
		c4 <- p1 + p2 + p3 + p4 + p5
		c(j, i, c0, c1, c2, c3, c4)
	}
	
	out <- as.data.frame(out)
	colnames(out) <- c("person", "item", "c0", "c1", "c2", "c3", "c4")
	# Ensure integer indices
	out$person <- as.integer(out$person)
	out$item   <- as.integer(out$item)
	out
}

#' Simulate item responses from cumulative probabilities
#'
#' @param N Integer. Number of persons.
#' @param I Integer. Number of items.
#' @param cprobs data.frame from generate_cprobs_par.
#' @param seed Numeric. RNG seed for U draws.
#' @param cores Integer >=1. Worker count for inner parallel.
#'
#' @return N x I integer matrix of responses 0..4.
generate_item_responses <- function(N, I, cprobs, seed, cores = 1) {
	set.seed(seed)
	
	cl <- makeCluster(cores)
	on.exit(stopCluster(cl), add = TRUE)
	registerDoParallel(cl)
	
	U <- matrix(runif(N * I, 0, 1), nrow = N, ncol = I)
	
	idx_people <- seq_len(N)
	resp <- foreach(j = idx_people, .combine = 'rbind') %dopar% {
		person_probs <- cprobs[cprobs$person == j, 3:7]
		# Ensure item order 1..I
		person_probs <- person_probs[order(cprobs$item[cprobs$person == j]), , drop = FALSE]
		out_j <- integer(I)
		for (it in seq_len(I)) {
			# First cumulative prob > U gives category index 1..5; subtract 1 -> 0..4
			out_j[it] <- min(which(person_probs[it, ] > U[j, it])) - 1L
		}
		out_j
	}
	colnames(resp) <- paste0("x", seq_len(I))
	resp
}

#' Validate per-item response distribution
#'
#' @param dat data.frame or matrix with item columns only.
#' @param min_count Integer. Minimum count per category 0..4.
#' @param I Integer. Number of items.
#' @return TRUE if all items meet min_count in all categories, else FALSE.
validate_response_distribution <- function(dat, min_count = 5, I) {
	invalid <- character(0)
	for (i in seq_len(I)) {
		counts <- table(factor(dat[, i], levels = 0:4))
		if (!all(counts >= min_count)) invalid <- c(invalid, colnames(dat)[i])
	}
	if (length(invalid) > 0) {
		cat("Columns failing count requirement:", paste(invalid, collapse = ", "), "\n")
		return(FALSE)
	} else {
		cat("All columns meet the count requirement.\n")
		return(TRUE)
	}
}

#' Weighted Random Forest wrapper
#'
#' @param data data.frame with predictors and target_var factor {0,1}.
#' @param target_var Character. Name of target column.
#' @param weights Optional numeric vector of per-row weights.
#' @return randomForest model object.
weighted_random_forest <- function(data, target_var, weights = NULL) {
	if (is.null(weights)) {
		cc <- table(data[[target_var]])
		weights <- 1 / cc[as.character(data[[target_var]])]
	}
	form <- as.formula(paste(target_var, "~ ."))
	randomForest(
		form,
		data = data,
		weights = weights,
		importance = TRUE,
		keep.forest = TRUE
	)
}

#' Maximize Youden's J for a scalar score
#'
#' @param scores Numeric vector of scores.
#' @param true_classes Numeric or factor {0,1}.
#' @return list(cutpoint, youdens_index)
maximize_youdens_index <- function(scores, true_classes) {
	scores <- as.numeric(scores)
	y <- as.integer(as.factor(true_classes)) - 1L  # map to {0,1}
	sc <- sort(unique(scores))
	best_cp <- NA_real_; best_J <- -Inf
	for (cp in sc) {
		pred <- ifelse(scores >= cp, 1L, 0L)
		tp <- sum(pred == 1L & y == 1L); fn <- sum(pred == 0L & y == 1L)
		tn <- sum(pred == 0L & y == 0L); fp <- sum(pred == 1L & y == 0L)
		sens <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
		spec <- ifelse(tn + fp > 0, tn / (tn + fp), 0)
		J <- sens + spec - 1
		if (J > best_J) { best_J <- J; best_cp <- cp }
	}
	list(cutpoint = best_cp, youdens_index = best_J)
}

#' Evaluate predictions vs truth
#'
#' @param method Character tag.
#' @param predictions Vector or factor of predicted classes {0,1}.
#' @param test_data data.frame with truth in last column.
#' @return single-row data.frame of metrics.
evaluate_model <- function(method, predictions, test_data) {
	truth <- as.factor(test_data[, ncol(test_data)])
	pred  <- as.factor(predictions)
	cm <- confusionMatrix(pred, truth)
	data.frame(
		method = method,
		TP = cm$table[2, 2],
		TN = cm$table[1, 1],
		FP = cm$table[2, 1],
		FN = cm$table[1, 2],
		accuracy = cm$overall["Accuracy"],
		sensitivity = cm$byClass["Sensitivity"],
		specificity = cm$byClass["Specificity"],
		precision = cm$byClass["Precision"],
		npv = cm$table[1, 1] / (cm$table[1, 1] + cm$table[2, 1]),
		f1_score = cm$byClass["F1"],
		row.names = NULL
	)
}

################## Load design ########################################
# condition       = integer ID for the condition row in allConditions.csv
# 
# refGroupSize    = number of simulees in the reference group (N_ref)
# 
# focalGroupSize  = number of simulees in the focal group (N_foc)
# 
# nItems          = number of items (I)
# 
# prevRate        = prevalence of the diagnosed class (Pr[Y=1])
# 
# difRate         = proportion of items with DIF (|D|/I), where D is the DIF item set
# 
# difType         = which parameters receive DIF: "none", "a", "b", or "both"
# (the script currently applies both when DIF is present)
# 
# balanced        = DIF direction:
# 	TRUE  → half of DIF items increase difficulty/discrimination,
# half decrease.
# FALSE → all DIF items move in the same direction (increase).
# 
# aStrength       = multiplicative change to discrimination on DIF items
# For item i ∈ D:
# 	if increase: a_i^(foc) = a_i^(ref) * (1 + aStrength)
# if decrease: a_i^(foc) = a_i^(ref) * (1 - aStrength)
# 
# bStrength       = additive shift to all thresholds on DIF items
# For item i ∈ D and category thresholds k=1..4:
# 	if increase: b_{ik}^(foc) = b_{ik}^(ref) + bStrength
# if decrease: b_{ik}^(foc) = b_{ik}^(ref) - bStrength

design <- read.csv("allConditions.csv", header = TRUE, stringsAsFactors = FALSE)

# Create results folder if missing
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

################## Main loop over design ##############################
conds <- unique(design$condition)

for (c_id in conds) {
	conditions <- subset(design, condition == c_id)
	
	# -------- Outer parallel sizing -----------------------------------
	total_cores <- detectCores()
	outer_cores <- max(1L, total_cores - 2L)   # keep system responsive
	inner_cores <- 1L                          # safe default for nested work
	
	outer_cl <- makeCluster(outer_cores)
	registerDoParallel(outer_cl)
	
	# -------- Read dynamic condition fields ----------------------------
	rg          <- conditions$refGroupSize
	fg          <- conditions$focalGroupSize
	I           <- conditions$nItems
	prevRate    <- conditions$prevRate
	difRate     <- conditions$difRate
	difType     <- conditions$difType
	balanced    <- conditions$balanced
	a_strength  <- conditions$aStrength
	b_strength  <- conditions$bStrength
	
	# -------- Fixed hyperparameters ------------------------------------
	n_cat     <- 5
	a_mean    <- 1.7
	a_sd      <- 0.3
	b1_mean   <- -1.5
	b1_sd     <- 0.5
	diff_par  <- 1
	diff_sd   <- 0.2
	
	clusterExport(
		outer_cl,
		c("n_cat", "a_mean", "a_sd",
		  "b1_mean", "b1_sd", "diff_par", "diff_sd", "inner_cores",
		  "generate_cprobs_par", "generate_item_responses",
		  "validate_response_distribution", "weighted_random_forest",
		  "evaluate_model", "maximize_youdens_index",
		  "rg", "fg", "I", "prevRate", "difRate", "difType",
		  "balanced", "a_strength", "b_strength", "TOP_LEVEL_SEED"),
		envir = environment()
	)
	
	# -------- Replications ---------------------------------------------
	# .combine='c' collects a flat list of length 2*reps: [results, importances, ...]
	allResults <- foreach(
		rep = 1:30,
		.combine = 'c',
		.errorhandling = "remove",
		.packages = c('mirt', 'pROC', 'randomForest', 'dplyr',
					  'MASS', 'caret', 'mclust', 'foreach', 'doParallel', 'parallel')
	) %dopar% {
		tryCatch({
			# Independent seed per worker and attempt
			rep_seed <- TOP_LEVEL_SEED + c_id * 1e5 + rep * 1009
			set.seed(rep_seed)
			
			local_results <- NULL
			local_importance <- NULL
			
			valid_data <- FALSE
			attempts <- 0L
			max_attempts <- 20L
			
			while (!valid_data && attempts < max_attempts) {
				seed <- rep_seed + attempts
				attempts <- attempts + 1L
				
				# ----- Generate Theta Vector --------------
				# Generate a single latent variable for all persons
				theta_diag <- rnorm(rg + fg, mean = 0, sd = 1)
				
				# Split by group
				thetaRef_diag   <- theta_diag[1:rg]
				thetaFocal_diag <- theta_diag[(rg + 1):(rg + fg)]
				
				
				# ----- Create binary diagnosis by prevalence ------------------
				cut_point_diag <- quantile(theta_diag, probs = 1 - prevRate)
				diag_classes <- ifelse(theta_diag > cut_point_diag, 1L, 0L)
				
				# ----- Reference group item parameters ------------------------
				a <- rnorm(I, a_mean, a_sd)
				bMat <- matrix(0, I, n_cat - 1)
				for (j in 1:I) {
					bMat[j, 1] <- rnorm(1, b1_mean, b1_sd)
					bMat[j, 2] <- bMat[j, 1] + rnorm(1, diff_par, diff_sd)
					bMat[j, 3] <- bMat[j, 2] + rnorm(1, diff_par, diff_sd)
					bMat[j, 4] <- bMat[j, 3] + rnorm(1, diff_par, diff_sd)
				}
				colnames(bMat) <- c("b1", "b2", "b3", "b4")
				
				# ----- Reference group responses ------------------------------
				cprobs_ref <- generate_cprobs_par(rg, I, thetaRef_diag, a, bMat, cores = inner_cores)
				resp_ref   <- generate_item_responses(rg, I, cprobs_ref, seed, cores = inner_cores)
				
				ref_params <- cbind(a = a, bMat)
				
				# ----- Focal DIF application ----------------------------------
				n_dif <- round(I * difRate)
				dif_items_idx <- if (n_dif > 0) sample.int(I, n_dif) else integer(0)
				dif_items_str <- toString(dif_items_idx)
				
				a_foc   <- a
				bMat_f  <- bMat
				
				if (isTRUE(as.logical(balanced))) {
					n_half <- floor(n_dif / 2)
					inc_idx <- dif_items_idx[seq_len(n_half)]
					dec_idx <- setdiff(dif_items_idx, inc_idx)
					
					if (length(inc_idx)) {
						bMat_f[inc_idx, ] <- bMat[inc_idx, ] + b_strength
						a_foc[inc_idx]    <- a[inc_idx] * (1 + a_strength)
					}
					if (length(dec_idx)) {
						bMat_f[dec_idx, ] <- bMat[dec_idx, ] - b_strength
						a_foc[dec_idx]    <- a[dec_idx] * (1 - a_strength)
					}
				} else {
					if (length(dif_items_idx)) {
						bMat_f[dif_items_idx, ] <- bMat[dif_items_idx, ] + b_strength
						a_foc[dif_items_idx]    <- a[dif_items_idx] * (1 + a_strength)
					}
				}
				
				foc_params <- cbind(a = a_foc, bMat_f)
				
				# ----- Focal group responses ---------------------------------
				cprobs_foc <- generate_cprobs_par(fg, I, thetaFocal_diag, a_foc, bMat_f, cores = inner_cores)
				resp_foc   <- generate_item_responses(fg, I, cprobs_foc, seed, cores = inner_cores)
				
				# ----- Bind data and attach classes ---------------------------
				dat <- rbind(resp_ref, resp_foc)
				dat <- cbind(dat, diag_classes = diag_classes)
				rownames(dat) <- seq_len(nrow(dat))
				
				gen_item_params <- dplyr::bind_rows(
					data.frame(ref_params,  group = "ref",   item = seq_len(I)),
					data.frame(foc_params,  group = "focal", item = seq_len(I))
				)
				
				# ----- Distribution validity ----------------------------------
				valid_data <- validate_response_distribution(dat[, 1:I, drop = FALSE], I = I)
				
				if (valid_data) {
					# ----- Train/Test split with validity checks ----------------
					ok_splits <- 0L
					while (ok_splits < 2L) {
						tr_idx <- sample(seq_len(nrow(dat)), size = 0.5 * nrow(dat))
						train  <- as.data.frame(dat[tr_idx, , drop = FALSE])
						test   <- as.data.frame(dat[-tr_idx, , drop = FALSE])
						
						ok_splits <- 0L
						ok_splits <- ok_splits + as.integer(
							validate_response_distribution(train[, 1:I, drop = FALSE], min_count = 1, I = I)
						)
						ok_splits <- ok_splits + as.integer(
							validate_response_distribution(test[,  1:I, drop = FALSE], min_count = 1, I = I)
						)
					}
					
					focal_rows_test <- which(as.integer(rownames(test)) > rg)
					
					# ==================== IRT analysis ==========================
					items_only <- train[, 1:I, drop = FALSE]
					irtModel <- mirt(items_only, model = 1, itemtype = 'graded')
					
					items_test <- test[, 1:I, drop = FALSE]
					theta_test <- fscores(irtModel, response.pattern = items_test, method = "EAP")[, 1]
					theta_tr   <- fscores(irtModel)[, 1]
					
					cp <- maximize_youdens_index(theta_tr, train$diag_classes)$cutpoint
					pred_class <- ifelse(theta_test >= cp, 1L, 0L)
					
					eval_full  <- evaluate_model("IRT", pred_class, test)
					eval_focal <- evaluate_model("IRT", pred_class[focal_rows_test], test[focal_rows_test, ])
					eval_ref   <- evaluate_model("IRT", pred_class[-focal_rows_test], test[-focal_rows_test, ])
					
					irt_row <- cbind(
						RefGroupSize = rg, FocalGroupSize = fg, nItems = I,
						prevRate = prevRate, difRate = difRate, difType = difType,
						balanced = balanced, aStrength = a_strength, bStrength = b_strength,
						rep = rep, dif_items = dif_items_str,
						eval_full, eval_focal, eval_ref
					)
					local_results <- rbind(local_results, irt_row)
					
					# ==================== Random Forest =========================
					n1 <- nrow(dat) * prevRate
					n0 <- nrow(dat) - n1
					w1 <- ifelse(n1 > 0, (1 * n0) / n1, 1)
					
					train$diag_classes <- factor(train$diag_classes, levels = c(0, 1))
					class_w <- ifelse(train$diag_classes == 1, w1, 1)
					
					rfModel <- weighted_random_forest(train, target_var = "diag_classes", weights = class_w)
					
					# Use probability of class "1" for Youden
					tr_prob <- predict(rfModel, train[, 1:I, drop = FALSE], type = "prob")[, "1"]
					cp_rf   <- maximize_youdens_index(tr_prob, train$diag_classes)$cutpoint
					
					# Variable importance
					imp <- importance(rfModel)
					imp_df <- as.data.frame(imp)
					imp_df$Variable <- rownames(imp_df)
					local_importance <- imp_df$MeanDecreaseAccuracy
					
					te_prob <- predict(rfModel, test[, 1:I, drop = FALSE], type = "prob")[, "1"]
					pred_rf_prob <- ifelse(te_prob >= cp_rf, 1L, 0L)
					
					eval_full  <- evaluate_model("rf_probs", pred_rf_prob, test)
					eval_focal <- evaluate_model("rf_probs", pred_rf_prob[focal_rows_test], test[focal_rows_test, ])
					eval_ref   <- evaluate_model("rf_probs", pred_rf_prob[-focal_rows_test], test[-focal_rows_test, ])
					
					rf_prob_row <- cbind(
						RefGroupSize = rg, FocalGroupSize = fg, nItems = I,
						prevRate = prevRate, difRate = difRate, difType = difType,
						balanced = balanced, aStrength = a_strength, bStrength = b_strength,
						rep = rep, dif_items = dif_items_str,
						eval_full, eval_focal, eval_ref
					)
					local_results <- rbind(local_results, rf_prob_row)
					
					# RF direct class predictions
					pred_rf_cls <- predict(rfModel, test[, 1:I, drop = FALSE])
					eval_full  <- evaluate_model("rf_class", pred_rf_cls, test)
					eval_focal <- evaluate_model("rf_class", pred_rf_cls[focal_rows_test], test[focal_rows_test, ])
					eval_ref   <- evaluate_model("rf_class", pred_rf_cls[-focal_rows_test], test[-focal_rows_test, ])
					
					rf_class_row <- cbind(
						RefGroupSize = rg, FocalGroupSize = fg, nItems = I,
						prevRate = prevRate, difRate = difRate, difType = difType,
						balanced = balanced, aStrength = a_strength, bStrength = b_strength,
						rep = rep, dif_items = dif_items_str,
						eval_full, eval_focal, eval_ref
					)
					local_results <- rbind(local_results, rf_class_row)
				}
			} # while attempts
			
			if (attempts >= max_attempts) {
				local_results <- matrix(NA, nrow = 1, ncol = 41)
				local_importance <- NA
			}
			
			# Return pair [results, importance] for this rep
			list(local_results, local_importance)
		}, error = function(e) {
			e_file_name <- paste0("error_log_condition", c_id, ".txt")
			cat(format(Sys.time()), " Error in rep", rep, ":", e$message, "\n",
				file = e_file_name, append = TRUE)
			NULL  # removed due to .errorhandling="remove"
		})
	} # foreach reps
	
	# -------- Collect outputs -----------------------------------------
	imp_elems    <- allResults[seq(2, length(allResults), by = 2)]
	results_elems<- allResults[seq(1, length(allResults), by = 2)]
	
	allImportanceData <- do.call(rbind, imp_elems)
	allResultsData    <- do.call(rbind, results_elems)
	
	colnames(allResultsData) <- c(
		"RefGroupSize", "FocalGroupSize", "nItems", "prevRate", "difRate",
		"difType", "balanced", "aStrength", "bStrength", "rep", "dif_items",
		"method",
		"TP_full", "TN_full", "FP_full", "FN_full", "accuracy_full",
		"sensitivity_full", "specificity_full", "precision_full", "npv_full",
		"f1score_full",
		"TP_focal", "TN_focal", "FP_focal", "FN_focal", "accuracy_focal",
		"sensitivity_focal", "specificity_focal", "precision_focal", "npv_focal",
		"f1score_focal",
		"TP_ref", "TN_ref", "FP_ref", "FN_ref", "accuracy_ref",
		"sensitivity_ref", "specificity_ref", "precision_ref", "npv_ref",
		"f1score_ref"
	)
	
	# -------- Write to disk -------------------------------------------
	write.csv(allResultsData,
			  file = file.path("results", paste0("condition", c_id, "FinalResults.csv")),
			  row.names = FALSE)
	
	write.csv(allImportanceData,
			  file = file.path("results", paste0("condition", c_id, "FinalVarImp.csv")),
			  row.names = FALSE)
	
	stopCluster(outer_cl)
} # end for conditions

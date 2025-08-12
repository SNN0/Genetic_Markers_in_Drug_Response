library("TCGAbiolinks")
library("glmnet")
library("factoextra") 
library("FactoMineR")
library("caret")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library('tidyverse')
library(ggplot2)
library(ggfortify)
library(glmnet)
library(boot)
library(dplyr)
library(stringr)
library(readxl)
library(survcomp)
library(patchwork)
library(gridExtra)
library(grid)
library(data.table)





# Function to process two TSV files (CNA file and gene expression file)
process_tsv_files <- function(cna_file, gene_expression_file) {
        # Read the CNA file
        cna_data <- read.table(cna_file, header = TRUE, sep = "\t")
        
        # Read the gene expression file
        gene_expression_data <- read.table(gene_expression_file, header = TRUE, sep = "\t")
        
        # Extract column names from the CNA file
        cna_columns <- colnames(cna_data)
        
        # Filter the gene expression file columns to match those in the CNA file
        gene_expression_filtered <- gene_expression_data %>%
                select(all_of(cna_columns))
        
        # Write the filtered data back to the same gene expression file
        write.table(gene_expression_filtered, file = gene_expression_file, sep = "\t", row.names = FALSE, quote = FALSE)
}


# Function to process drug sensitivity and resistance data
process_drug_data <- function(expr_sensitive_file, cna_sensitive_file = NULL, mut_sensitive_file = NULL, 
                              expr_resistant_file, cna_resistant_file = NULL, mut_resistant_file = NULL, 
                              conversion_table = NULL) {
        
        # Use gene symbol conversion table if provided as a dataframe
        if (!is.null(conversion_table)) {
                old_to_new <- setNames(conversion_table$Converted.Symbol, conversion_table$Old.Symbol)
        }
        
        # Load sensitive expression data
        expr_sensitive <- read.table(expr_sensitive_file, header = TRUE)
        
        # Convert gene symbols in expression data if conversion table is provided
        if (!is.null(conversion_table)) {
                expr_sensitive$symbol <- ifelse(expr_sensitive$symbol %in% names(old_to_new), 
                                                old_to_new[expr_sensitive$symbol], 
                                                expr_sensitive$symbol)
        }
        
        # Set rownames and clean column names for sensitive expression data
        expr_sensitive <- column_to_rownames(expr_sensitive, var = 'symbol')
        colnames(expr_sensitive) <- gsub("\\.\\d+[A-Z]$", "", colnames(expr_sensitive))
        colnames(expr_sensitive) <- gsub("\\.", "-", colnames(expr_sensitive))
        
        if (!is.null(cna_sensitive_file)) {
                # Load CNA data for sensitive samples
                cna_sensitive <- read.table(cna_sensitive_file, header = TRUE)
                
                # Convert gene symbols in CNA data if conversion table is provided
                if (!is.null(conversion_table)) {
                        cna_sensitive$symbol <- ifelse(cna_sensitive$symbol %in% names(old_to_new), 
                                                       old_to_new[cna_sensitive$symbol], 
                                                       cna_sensitive$symbol)
                }
                
                # Prepare CNA data for sensitive samples
                cna_sensitive <- column_to_rownames(cna_sensitive, var = 'symbol')
                colnames(cna_sensitive) <- gsub("\\.\\d+[A-Z]$", "", colnames(cna_sensitive))
                colnames(cna_sensitive) <- gsub("\\.", "-", colnames(cna_sensitive))
                cna_sensitive_t <- t(cna_sensitive)
                colnames(cna_sensitive_t) <- paste0("cna_", colnames(cna_sensitive_t))
        } else {
                # Create an all-zero CNA matrix for sensitive samples
                cna_sensitive_t <- matrix(0, nrow = ncol(expr_sensitive), ncol = nrow(expr_sensitive))
                colnames(cna_sensitive_t) <- paste0("cna_", rownames(expr_sensitive))
                rownames(cna_sensitive_t) <- colnames(expr_sensitive)
        }
        
        # Prepare transposed and prefixed expression data for sensitive samples
        expr_sensitive_t <- t(expr_sensitive)
        colnames(expr_sensitive_t) <- paste0("expr_", colnames(expr_sensitive_t))
        
        # Check if mutation data file for sensitive samples exists
        if (!is.null(mut_sensitive_file)) {
                # Load mutation data for sensitive samples
                mut_sensitive <- read.table(mut_sensitive_file, header = TRUE)
                
                # Convert gene symbols in mutation data if conversion table is provided
                if (!is.null(conversion_table)) {
                        mut_sensitive$symbol <- ifelse(mut_sensitive$symbol %in% names(old_to_new), 
                                                       old_to_new[mut_sensitive$symbol], 
                                                       mut_sensitive$symbol)
                }
                
                # Prepare mutation data for sensitive samples
                mut_sensitive <- column_to_rownames(mut_sensitive, var = 'symbol')
                colnames(mut_sensitive) <- gsub("\\.\\d+[A-Z]$", "", colnames(mut_sensitive))
                colnames(mut_sensitive) <- gsub("\\.", "-", colnames(mut_sensitive))
                mut_sensitive_t <- t(mut_sensitive)
                colnames(mut_sensitive_t) <- paste0("mut_", colnames(mut_sensitive_t))
        } else {
                # Create an all-zero mutation matrix for sensitive samples
                mut_sensitive_t <- matrix(0, nrow = nrow(expr_sensitive_t), ncol = ncol(expr_sensitive_t))
                colnames(mut_sensitive_t) <- paste0("mut_", rownames(expr_sensitive))
                rownames(mut_sensitive_t) <- rownames(expr_sensitive_t)
        }
        
        # Combine sensitive data
        sensitive_data <- cbind(expr_sensitive_t, mut_sensitive_t, cna_sensitive_t)
        sensitive_data <- as.data.frame(sensitive_data)
        sensitive_data$condition <- 'sensitive'
        
        # Load resistant expression data
        expr_resistant <- read.table(expr_resistant_file, header = TRUE)
        
        # Convert gene symbols in expression data if conversion table is provided
        if (!is.null(conversion_table)) {
                expr_resistant$symbol <- ifelse(expr_resistant$symbol %in% names(old_to_new), 
                                                old_to_new[expr_resistant$symbol], 
                                                expr_resistant$symbol)
        }
        
        # Set rownames and clean column names for resistant expression data
        expr_resistant <- column_to_rownames(expr_resistant, var = 'symbol')
        colnames(expr_resistant) <- gsub("\\.\\d+[A-Z]$", "", colnames(expr_resistant))
        colnames(expr_resistant) <- gsub("\\.", "-", colnames(expr_resistant))
        
        if (!is.null(cna_resistant_file)) {
                # Load CNA data for resistant samples
                cna_resistant <- read.table(cna_resistant_file, header = TRUE)
                
                # Convert gene symbols in CNA data if conversion table is provided
                if (!is.null(conversion_table)) {
                        cna_resistant$symbol <- ifelse(cna_resistant$symbol %in% names(old_to_new), 
                                                       old_to_new[cna_resistant$symbol], 
                                                       cna_resistant$symbol)
                }
                
                # Prepare CNA data for resistant samples
                cna_resistant <- column_to_rownames(cna_resistant, var = 'symbol')
                colnames(cna_resistant) <- gsub("\\.\\d+[A-Z]$", "", colnames(cna_resistant))
                colnames(cna_resistant) <- gsub("\\.", "-", colnames(cna_resistant))
                cna_resistant_t <- t(cna_resistant)
                colnames(cna_resistant_t) <- paste0("cna_", colnames(cna_resistant_t))
        } else {
                # Create an all-zero CNA matrix for resistant samples
                cna_resistant_t <- matrix(0, nrow = ncol(expr_resistant), ncol = nrow(expr_resistant))
                colnames(cna_resistant_t) <- paste0("cna_", rownames(expr_resistant))
                rownames(cna_resistant_t) <- colnames(expr_resistant)
        }
        
        # Prepare transposed and prefixed expression data for resistant samples
        expr_resistant_t <- t(expr_resistant)
        colnames(expr_resistant_t) <- paste0("expr_", colnames(expr_resistant_t))
        
        # Check if mutation data file for resistant samples exists
        if (!is.null(mut_resistant_file)) {
                # Load mutation data for resistant samples
                mut_resistant <- read.table(mut_resistant_file, header = TRUE)
                
                # Convert gene symbols in mutation data if conversion table is provided
                if (!is.null(conversion_table)) {
                        mut_resistant$symbol <- ifelse(mut_resistant$symbol %in% names(old_to_new), 
                                                       old_to_new[mut_resistant$symbol], 
                                                       mut_resistant$symbol)
                }
                
                # Prepare mutation data for resistant samples
                mut_resistant <- column_to_rownames(mut_resistant, var = 'symbol')
                colnames(mut_resistant) <- gsub("\\.\\d+[A-Z]$", "", colnames(mut_resistant))
                colnames(mut_resistant) <- gsub("\\.", "-", colnames(mut_resistant))
                mut_resistant_t <- t(mut_resistant)
                colnames(mut_resistant_t) <- paste0("mut_", colnames(mut_resistant_t))
        } else {
                # Create an all-zero mutation matrix for resistant samples
                mut_resistant_t <- matrix(0, nrow = nrow(expr_resistant_t), ncol = ncol(expr_resistant_t))
                colnames(mut_resistant_t) <- paste0("mut_", rownames(expr_resistant))
                rownames(mut_resistant_t) <- rownames(expr_resistant_t)
        }
        
        # Combine resistant data
        resistant_data <- cbind(expr_resistant_t, mut_resistant_t, cna_resistant_t)
        resistant_data <- as.data.frame(resistant_data)
        resistant_data$condition <- 'resistant'
        
        # Combine sensitive and resistant datasets
        full_data <- rbind(sensitive_data, resistant_data)
        
        # Set rownames for the combined dataset
        rownames(full_data) <- rownames(full_data)
        
        return(full_data)
}


# Function to retrieve and process clinical data for multiple TCGA projects
get_multiple_clinical_data <- function(tcga_projects) {
        # Create an empty list to store clinical data for each project
        clinical_data_list <- list()
        
        # Loop through each project to fetch and process clinical data
        for (i in seq_along(tcga_projects)) {
                tcga_project <- tcga_projects[i]
                
                # Fetch clinical data for the current TCGA project
                clinical_data <- GDCquery_clinic(tcga_project)
                
                # Ensure the required columns are present in the clinical data
                if (!all(c('vital_status', 'days_to_last_follow_up', 'days_to_death') %in% colnames(clinical_data))) {
                        stop(paste("The clinical data for project", tcga_project, "does not contain the required columns."))
                }
                
                # Create a 'deceased' column based on 'vital_status'
                clinical_data$deceased <- ifelse(clinical_data$vital_status == 'Alive', FALSE, TRUE)
                
                # Create an 'overall_survival' column using survival days and vital status
                clinical_data$overall_survival <- ifelse(clinical_data$vital_status == 'Alive',
                                                         clinical_data$days_to_last_follow_up,
                                                         clinical_data$days_to_death)
                
                # Filter and select only the required columns
                clinical_data_filtered <- clinical_data[, c('submitter_id', 'deceased', 'overall_survival')]
                
                # Append the processed clinical data to the list
                clinical_data_list[[i]] <- clinical_data_filtered
        }
        
        # Combine all processed clinical data into a single dataframe
        combined_clinical_data <- do.call(rbind, clinical_data_list)
        
        return(combined_clinical_data)
}




drug_small_sample_models_cond_penalized_with_pvalue_optimized <- function(data, surv_obj) {
        
        #' Penalized Cox Model Function using glmnet package
        #' This function performs penalized Cox regression for genes in the provided dataset
        #' and calculates HR, p-values, CIs, and other metrics for two models:
        #' 1. Model with CNA 
        #' 2. Model with Mutation data.
        
        # Extract gene names
        gene_expr_cols <- grep("^expr_", colnames(data), value = TRUE)
        gene_names <- gsub("expr_", "", gene_expr_cols)
        
        cat("Number of genes to process:", length(gene_names), "\n") # Print the number of genes
        
        # Initialize results for both models with predefined columns
        model_1_results <- data.frame(
                Gene = character(),
                Expression_HR = numeric(),
                Expression_pvalue = numeric(),
                Expression_LowerCI = numeric(),
                Expression_UpperCI = numeric(),
                CNA_HR = numeric(),
                CNA_pvalue = numeric(),
                CNA_LowerCI = numeric(),
                CNA_UpperCI = numeric(),
                Condition_HR = numeric(),
                Condition_pvalue = numeric(),
                Condition_LowerCI = numeric(),
                Condition_UpperCI = numeric(),
                Interaction_cna_condition_HR = numeric(),
                Interaction_cna_condition_pvalue = numeric(),
                Interaction_cna_condition_LowerCI = numeric(),
                Interaction_cna_condition_UpperCI = numeric(),
                Concordance_Index = numeric(),
                AIC = numeric(),
                stringsAsFactors = FALSE
        )
        
        model_2_results <- data.frame(
                Gene = character(),
                Expression_HR = numeric(),
                Expression_pvalue = numeric(),
                Expression_LowerCI = numeric(),
                Expression_UpperCI = numeric(),
                Mut_HR = numeric(),
                Mut_pvalue = numeric(),
                Mut_LowerCI = numeric(),
                Mut_UpperCI = numeric(),
                Condition_HR = numeric(),
                Condition_pvalue = numeric(),
                Condition_LowerCI = numeric(),
                Condition_UpperCI = numeric(),
                Interaction_mut_condition_HR = numeric(),
                Interaction_mut_condition_pvalue = numeric(),
                Interaction_mut_condition_LowerCI = numeric(),
                Interaction_mut_condition_UpperCI = numeric(),
                Concordance_Index = numeric(),
                AIC = numeric(),
                stringsAsFactors = FALSE
        )
        
        # Load parallel library
        library(doParallel)
        library(survcomp)
        
        # Set up parallel backend
        num_cores <- detectCores() - 1
        cl <- makeCluster(num_cores)
        registerDoParallel(cl)
        
        # Loop through each gene and build models
        for (gene in gene_names) {
                cat("Processing gene:", gene, "\n") # Progress message
                
                # Check if mut and cna are present for the gene
                mut_included <- paste0("mut_", gene) %in% colnames(data) && length(unique(data[[paste0("mut_", gene)]])) > 1
                cna_included <- paste0("cna_", gene) %in% colnames(data) && length(unique(data[[paste0("cna_", gene)]])) > 1
                
                cat("mut_included:", mut_included, ", cna_included:", cna_included, "\n") # Print inclusion status
                
                # Skip gene if neither mut nor cna is present
                if (!cna_included && !mut_included) {
                        next
                }
                
                # Model 1: exp + cna + condition + cna:condition
                if (cna_included) {
                        model_matrix <- model.matrix(~ data[[paste0("expr_", gene)]] + factor(data[[paste0("cna_", gene)]]) + 
                                                             factor(data$condition) + factor(data[[paste0("cna_", gene)]]) : factor(data$condition))
                        fit <- tryCatch({
                                cat("Fitting penalized Cox model using glmnet...\n")
                                set.seed(73)
                                cv.glmnet(model_matrix, surv_obj, family = "cox", alpha = 0)
                        }, error = function(e) {
                                cat("Error during model fitting:\n", e$message, "\n")
                                NULL
                        })
                        
                        if (!is.null(fit)) {
                                best_lambda <- fit$lambda.min
                                if (!is.null(best_lambda)) {
                                        cox_fit <- glmnet(model_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                        
                                        # Extract coefficients and calculate HR
                                        coef_info <- coef(cox_fit)
                                        exp_coef <- exp(coef_info)
                                        
                                        # Calculate concordance index using glmnet predicted values
                                        risk_score <- predict(cox_fit, newx = model_matrix, type = "link")
                                        cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                        
                                        # Calculate AIC
                                        log_likelihood <- sum(dpois(surv_obj[, 1], lambda = exp(risk_score), log = TRUE))
                                        aic_value <- -2 * log_likelihood + 2 * length(coef_info)
                                        
                                        cat("Starting p-value estimation using bootstrapping...")
                                        error_count <- 0
                                        valid_bootstrap_coefs <- list()
                                        set.seed(73)
                                        
                                        # Use foreach for parallel bootstrapping
                                        boot_results <- foreach(iteration = 1:1000, .combine = 'c', .packages = c("glmnet")) %dopar% {
                                                sample_idx <- sample(1:nrow(data), replace = TRUE)
                                                boot_matrix <- model_matrix[sample_idx, ]
                                                boot_surv <- surv_obj[sample_idx]
                                                tryCatch({
                                                        boot_fit <- glmnet(boot_matrix, boot_surv, family = "cox", alpha = 0, lambda = best_lambda)
                                                        coef(boot_fit)
                                                }, error = function(e) {
                                                        error_count <<- error_count + 1
                                                        return(NULL)
                                                })
                                        }
                                        
                                        # Filter out NULL results
                                        valid_bootstrap_coefs <- Filter(Negate(is.null), boot_results)
                                        
                                        if (length(valid_bootstrap_coefs) > 0) {
                                                valid_bootstrap_coefs <- do.call(cbind, valid_bootstrap_coefs)
                                                
                                                # Create empty p-values and CI vectors
                                                p_values <- numeric(nrow(valid_bootstrap_coefs))
                                                lower_ci <- numeric(nrow(valid_bootstrap_coefs))
                                                upper_ci <- numeric(nrow(valid_bootstrap_coefs))
                                                
                                                #  Calculate p-value and CI for each coefficient
                                                for (i in 1:nrow(valid_bootstrap_coefs)) {
                                                        bootstrap_values <- valid_bootstrap_coefs[i, ]
                                                        p_values[i] <- mean(abs(bootstrap_values) >= abs(coef_info[i]), na.rm = TRUE)
                                                        lower_ci[i] <- exp(quantile(bootstrap_values, 0.025, na.rm = TRUE))
                                                        upper_ci[i] <- exp(quantile(bootstrap_values, 0.975, na.rm = TRUE))
                                                }
                                        } else {
                                                p_values <- rep(NA, length(coef_info))
                                                lower_ci <- rep(NA, length(coef_info))
                                                upper_ci <- rep(NA, length(coef_info))
                                        }
                                        
                                        cat("Number of errors during bootstrapping:", error_count, "\n")
                                        
                                        # Store results
                                        gene_results <- data.frame(
                                                Gene = gene,
                                                Expression_HR = exp_coef[2],
                                                Expression_pvalue = p_values[2],
                                                Expression_LowerCI = lower_ci[2],
                                                Expression_UpperCI = upper_ci[2],
                                                CNA_HR = exp_coef[3],
                                                CNA_pvalue = p_values[3],
                                                CNA_LowerCI = lower_ci[3],
                                                CNA_UpperCI = upper_ci[3],
                                                Condition_HR = exp_coef[4],
                                                Condition_pvalue = p_values[4],
                                                Condition_LowerCI = lower_ci[4],
                                                Condition_UpperCI = upper_ci[4],
                                                Interaction_cna_condition_HR = exp_coef[5],
                                                Interaction_cna_condition_pvalue = p_values[5],
                                                Interaction_cna_condition_LowerCI = lower_ci[5],
                                                Interaction_cna_condition_UpperCI = upper_ci[5],
                                                Concordance_Index = cindex,
                                                AIC = aic_value
                                        )
                                        missing_cols <- setdiff(names(model_1_results), names(gene_results))
                                        gene_results[missing_cols] <- NA
                                        model_1_results <- rbind(model_1_results, gene_results)
                                }
                        }
                }
                
                # Model 2: exp + mut + condition + mut:condition
                if (mut_included) {
                        model_matrix <- model.matrix(~ data[[paste0("expr_", gene)]] + factor(data[[paste0("mut_", gene)]]) + 
                                                             factor(data$condition) + factor(data[[paste0("mut_", gene)]]) : factor(data$condition))
                        fit <- tryCatch({
                                cat("Fitting penalized Cox model using glmnet...\n")
                                set.seed(73)
                                cv.glmnet(model_matrix, surv_obj, family = "cox", alpha = 0)
                        }, error = function(e) {
                                cat("Error during model fitting:\n", e$message, "\n")
                                NULL
                        })
                        
                        if (!is.null(fit)) {
                                best_lambda <- fit$lambda.min
                                if (!is.null(best_lambda)) {
                                        cox_fit <- glmnet(model_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                        
                                        # Extract coefficients and calculate HR
                                        coef_info <- coef(cox_fit)
                                        exp_coef <- exp(coef_info)
                                        
                                        # Calculate concordance index using glmnet predicted values
                                        risk_score <- predict(cox_fit, newx = model_matrix, type = "link")
                                        cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                        
                                        # Calculate AIC
                                        log_likelihood <- sum(dpois(surv_obj[, 1], lambda = exp(risk_score), log = TRUE))
                                        aic_value <- -2 * log_likelihood + 2 * length(coef_info)
                                        
                                        cat("Starting p-value estimation using bootstrapping...")
                                        error_count <- 0
                                        valid_bootstrap_coefs <- list()
                                        set.seed(73)
                                        
                                        # Use foreach for parallel bootstrapping
                                        boot_results <- foreach(iteration = 1:1000, .combine = 'c', .packages = c("glmnet")) %dopar% {
                                                sample_idx <- sample(1:nrow(data), replace = TRUE)
                                                boot_matrix <- model_matrix[sample_idx, ]
                                                boot_surv <- surv_obj[sample_idx]
                                                tryCatch({
                                                        boot_fit <- glmnet(boot_matrix, boot_surv, family = "cox", alpha = 0, lambda = best_lambda)
                                                        coef(boot_fit)
                                                }, error = function(e) {
                                                        error_count <<- error_count + 1
                                                        return(NULL)
                                                })
                                        }
                                        
                                        # Filter out NULL results
                                        valid_bootstrap_coefs <- Filter(Negate(is.null), boot_results)
                                        
                                        if (length(valid_bootstrap_coefs) > 0) {
                                                valid_bootstrap_coefs <- do.call(cbind, valid_bootstrap_coefs)
                                                
                                                # Create empty p-values and CI vectors
                                                p_values <- numeric(nrow(valid_bootstrap_coefs))
                                                lower_ci <- numeric(nrow(valid_bootstrap_coefs))
                                                upper_ci <- numeric(nrow(valid_bootstrap_coefs))
                                                
                                                #  Calculate p-value and CI for each coefficient
                                                for (i in 1:nrow(valid_bootstrap_coefs)) {
                                                        bootstrap_values <- valid_bootstrap_coefs[i, ]
                                                        p_values[i] <- mean(abs(bootstrap_values) >= abs(coef_info[i]), na.rm = TRUE)
                                                        lower_ci[i] <- exp(quantile(bootstrap_values, 0.025, na.rm = TRUE))
                                                        upper_ci[i] <- exp(quantile(bootstrap_values, 0.975, na.rm = TRUE))
                                                }
                                        } else {
                                                p_values <- rep(NA, length(coef_info))
                                                lower_ci <- rep(NA, length(coef_info))
                                                upper_ci <- rep(NA, length(coef_info))
                                        }
                                        
                                        cat("Number of errors during bootstrapping:", error_count, "\n")
                                        
                                        # Store results
                                        gene_results <- data.frame(
                                                Gene = gene,
                                                Expression_HR = exp_coef[2],
                                                Expression_pvalue = p_values[2],
                                                Expression_LowerCI = lower_ci[2],
                                                Expression_UpperCI = upper_ci[2],
                                                Mut_HR = exp_coef[3],
                                                Mut_pvalue = p_values[3],
                                                Mut_LowerCI = lower_ci[3],
                                                Mut_UpperCI = upper_ci[3],
                                                Condition_HR = exp_coef[4],
                                                Condition_pvalue = p_values[4],
                                                Condition_LowerCI = lower_ci[4],
                                                Condition_UpperCI = upper_ci[4],
                                                Interaction_mut_condition_HR = exp_coef[5],
                                                Interaction_mut_condition_pvalue = p_values[5],
                                                Interaction_mut_condition_LowerCI = lower_ci[5],
                                                Interaction_mut_condition_UpperCI = upper_ci[5],
                                                Concordance_Index = cindex,
                                                AIC = aic_value
                                        )
                                        missing_cols <- setdiff(names(model_2_results), names(gene_results))
                                        gene_results[missing_cols] <- NA
                                        model_2_results <- rbind(model_2_results, gene_results)
                                }
                        }
                }
        }
        
        # Stop parallel backend
        stopCluster(cl)
        
        # Return the results
        return(list(
                Model_1 = model_1_results,
                Model_2 = model_2_results
        ))
}




# Penalized Cox Model Function using glmnet package for multiple genes
drug_multigene_models_cond_penalized <- function(data, surv_obj, pathway_file, single_gene_models, combination = TRUE) {
        # Load gene list from the pathway file
        pathway_data <- read_excel(pathway_file)
        colnames(pathway_data) <- pathway_data[1, ]
        pathway_data <- pathway_data[-1, ]
        
        # Initialize results for both models with predefined columns
        model_1_results <- data.frame(
                Pathway = character(),
                Genes = character(),
                Concordance_Index = numeric(),
                AIC = numeric(),
                stringsAsFactors = FALSE
        )
        
        model_2_results <- data.frame(
                Pathway = character(),
                Genes = character(),
                Concordance_Index = numeric(),
                AIC = numeric(),
                stringsAsFactors = FALSE
        )
        
        best_pathways_model_1 <- data.frame(
                Pathway = character(),
                Best_Gene = character(),
                Best_Gene_C_Index = numeric(),
                Multi_C_Index = numeric(),
                C_Index_Improvement = numeric(),
                Best_Gene_AIC = numeric(),
                Multi_AIC = numeric(),
                AIC_Improvement = numeric(),
                stringsAsFactors = FALSE
        )
        
        best_pathways_model_2 <- data.frame(
                Pathway = character(),
                Best_Gene = character(),
                Best_Gene_C_Index = numeric(),
                Multi_C_Index = numeric(),
                C_Index_Improvement = numeric(),
                Best_Gene_AIC = numeric(),
                Multi_AIC = numeric(),
                AIC_Improvement = numeric(),
                stringsAsFactors = FALSE
        )
        
        # HR results for each pathway
        model_1_hr <- list()
        model_2_hr <- list()
        
        # Iterate over each pathway
        for (i in 1:nrow(pathway_data)) {
                pathway <- pathway_data$Term_Description[i]
                up_genes <- unlist(strsplit(as.character(pathway_data$Up_regulated[i]), ", "))
                down_genes <- unlist(strsplit(as.character(pathway_data$Down_regulated[i]), ", "))
                gene_list <- unique(c(up_genes, down_genes))
                
                # Extract gene names from the provided list
                gene_list <- intersect(gene_list, gsub("expr_", "", grep("^expr_", colnames(data), value = TRUE)))
                
                cat("Number of genes to process for pathway:", pathway, " - ", length(gene_list), "\n") # Print the number of genes in the list
                
                # Generate gene combinations if combination is TRUE, otherwise use the full list
                if (combination) {
                        gene_combinations <- unlist(lapply(1:length(gene_list), function(x) combn(gene_list, x, simplify = FALSE)), recursive = FALSE)
                } else {
                        gene_combinations <- list(gene_list)
                }
                
                # Loop through each combination of genes
                for (genes in gene_combinations) {
                        # Model 1: exp + cna + condition + cna:condition for multiple genes
                        cat("Building Model 1 for pathway:", pathway, " with genes:", paste(genes, collapse = ", "), "\n")
                        model_1_matrix <- NULL
                        selected_genes <- c()
                        for (gene in genes) {
                                cna_included <- paste0("cna_", gene) %in% colnames(data) && length(unique(data[[paste0("cna_", gene)]])) > 1
                                if (cna_included) {
                                        selected_genes <- c(selected_genes, gene)
                                        model_1_matrix <- cbind(model_1_matrix, model.matrix(~ data[[paste0("expr_", gene)]] + factor(data[[paste0("cna_", gene)]]) + 
                                                                                                     factor(data$condition) + factor(data[[paste0("cna_", gene)]]) : factor(data$condition))[, -1])
                                }
                        }
                        
                        if (!is.null(model_1_matrix)) {
                                fit <- tryCatch({
                                        cat("Fitting penalized Cox model using glmnet for Model 1...\n")
                                        set.seed(73)
                                        cv.glmnet(model_1_matrix, surv_obj, family = "cox", alpha = 0)
                                }, error = function(e) {
                                        cat("Error during model fitting for Model 1:\n", e$message, "\n")
                                        NULL
                                })
                                
                                if (!is.null(fit)) {
                                        best_lambda <- fit$lambda.min
                                        if (!is.null(best_lambda)) {
                                                cox_fit <- glmnet(model_1_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                                
                                                # Calculate concordance index using glmnet predicted values
                                                set.seed(73)
                                                risk_score <- predict(cox_fit, newx = model_1_matrix, type = "link")
                                                cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                                
                                                # Calculate AIC
                                                log_likelihood <- sum(dpois(surv_obj[, 1], lambda = exp(risk_score), log = TRUE))
                                                aic_value <- -2 * log_likelihood + 2 * length(coef(cox_fit))
                                                
                                                # Extract HRs
                                                hr_values <- exp(as.numeric(coef(cox_fit)))
                                                hr_names <- colnames(model_1_matrix)
                                                model_1_hr[[pathway]] <- data.frame(Term = hr_names, HR = hr_values, stringsAsFactors = FALSE)
                                                
                                                # Store results
                                                model_1_results <- rbind(model_1_results, data.frame(
                                                        Pathway = pathway,
                                                        Genes = paste(selected_genes, collapse = ", "),
                                                        Concordance_Index = cindex,
                                                        AIC = aic_value
                                                ))
                                                
                                                # Best gene analysis for Model 1 using single_gene_models input
                                                single_gene_model <- single_gene_models$Filtered_Models$Model_1
                                                single_gene_cindices <- single_gene_model[single_gene_model$Gene %in% selected_genes, ]
                                                best_gene_row <- single_gene_cindices[which.max(single_gene_cindices$Concordance_Index), ]
                                                best_gene <- best_gene_row$Gene
                                                best_gene_cindex <- best_gene_row$Concordance_Index
                                                best_gene_aic <- best_gene_row$AIC
                                                
                                                # Store best pathway analysis
                                                c_index_improvement <- cindex - best_gene_cindex
                                                aic_improvement <- best_gene_aic - aic_value
                                                best_pathways_model_1 <- rbind(best_pathways_model_1, data.frame(
                                                        Pathway = pathway,
                                                        Best_Gene = best_gene,
                                                        Best_Gene_C_Index = best_gene_cindex,
                                                        Multi_C_Index = cindex,
                                                        C_Index_Improvement = c_index_improvement,
                                                        Best_Gene_AIC = best_gene_aic,
                                                        Multi_AIC = aic_value,
                                                        AIC_Improvement = aic_improvement,
                                                        stringsAsFactors = FALSE
                                                ))
                                        }
                                }
                        }
                        
                        # Model 2: exp + mut + condition + mut:condition for multiple genes
                        cat("Building Model 2 for pathway:", pathway, " with genes:", paste(genes, collapse = ", "), "\n")
                        model_2_matrix <- NULL
                        selected_genes <- c()
                        for (gene in genes) {
                                mut_included <- paste0("mut_", gene) %in% colnames(data) && length(unique(data[[paste0("mut_", gene)]])) > 1
                                if (mut_included) {
                                        selected_genes <- c(selected_genes, gene)
                                        model_2_matrix <- cbind(model_2_matrix, model.matrix(~ data[[paste0("expr_", gene)]] + factor(data[[paste0("mut_", gene)]]) + 
                                                                                                     factor(data$condition) + factor(data[[paste0("mut_", gene)]]) : factor(data$condition))[, -1])
                                }
                        }
                        
                        if (!is.null(model_2_matrix)) {
                                fit <- tryCatch({
                                        cat("Fitting penalized Cox model using glmnet for Model 2...\n")
                                        set.seed(73)
                                        cv.glmnet(model_2_matrix, surv_obj, family = "cox", alpha = 0)
                                }, error = function(e) {
                                        cat("Error during model fitting for Model 2:\n", e$message, "\n")
                                        NULL
                                })
                                
                                if (!is.null(fit)) {
                                        best_lambda <- fit$lambda.min
                                        if (!is.null(best_lambda)) {
                                                cox_fit <- glmnet(model_2_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                                
                                                # Calculate concordance index using glmnet predicted values
                                                set.seed(73)
                                                risk_score <- predict(cox_fit, newx = model_2_matrix, type = "link")
                                                cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                                
                                                # Calculate AIC
                                                log_likelihood <- sum(dpois(surv_obj[, 1], lambda = exp(risk_score), log = TRUE))
                                                aic_value <- -2 * log_likelihood + 2 * length(coef(cox_fit))
                                                
                                                # Extract HRs
                                                hr_values <- exp(as.numeric(coef(cox_fit)))
                                                hr_names <- colnames(model_2_matrix)
                                                model_2_hr[[pathway]] <- data.frame(Term = hr_names, HR = hr_values, stringsAsFactors = FALSE)
                                                
                                                # Store results
                                                model_2_results <- rbind(model_2_results, data.frame(
                                                        Pathway = pathway,
                                                        Genes = paste(selected_genes, collapse = ", "),
                                                        Concordance_Index = cindex,
                                                        AIC = aic_value
                                                ))
                                                
                                                # Best gene analysis for Model 2 using single_gene_models input
                                                single_gene_model <- single_gene_models$Filtered_Models$Model_2
                                                single_gene_cindices <- single_gene_model[single_gene_model$Gene %in% selected_genes, ]
                                                best_gene_row <- single_gene_cindices[which.max(single_gene_cindices$Concordance_Index), ]
                                                best_gene <- best_gene_row$Gene
                                                best_gene_cindex <- best_gene_row$Concordance_Index
                                                best_gene_aic <- best_gene_row$AIC
                                                
                                                # Store best pathway analysis
                                                c_index_improvement <- cindex - best_gene_cindex
                                                aic_improvement <- best_gene_aic - aic_value
                                                best_pathways_model_2 <- rbind(best_pathways_model_2, data.frame(
                                                        Pathway = pathway,
                                                        Best_Gene = best_gene,
                                                        Best_Gene_C_Index = best_gene_cindex,
                                                        Multi_C_Index = cindex,
                                                        C_Index_Improvement = c_index_improvement,
                                                        Best_Gene_AIC = best_gene_aic,
                                                        Multi_AIC = aic_value,
                                                        AIC_Improvement = aic_improvement,
                                                        stringsAsFactors = FALSE
                                                ))
                                        }
                                }
                        }
                }
        }
        
        # Return the results, including best pathways and HR values
        return(list(
                Model_1 = model_1_results,
                Model_2 = model_2_results,
                Best_Pathways_Model_1 = best_pathways_model_1,
                Best_Pathways_Model_2 = best_pathways_model_2,
                Model_1_HR = model_1_hr,
                Model_2_HR = model_2_hr
        ))
}



#---------------------------- GRAPH ---------------------------------------------

kaplan_meier_analysis <- function(model1_data, model2_data, surv_obj, processed_data, hr_threshold = 1.1, plot = TRUE) {
        analyze_model <- function(data, model_name) {
                if (is.null(data) || nrow(data) == 0) {
                        cat("Skipping", model_name, "because it is empty.\n")
                        return(list(Results = data.frame(), Plots = list()))
                }
                
                results <- data.frame(
                        Gene = character(),
                        HR_Type = character(),
                        P_Value = numeric(),
                        stringsAsFactors = FALSE
                )
                
                plot_list <- list()  # To store Kaplan-Meier plots
                
                # Mapping for shorter HR type names
                hr_type_labels <- list(
                        Expression_HR = "Expression",
                        CNA_HR = "CNA",
                        Condition_HR = "Response",
                        Interaction_cna_condition_HR = "CNA x Response Interaction",
                        Mut_HR = "Mutation",
                        Interaction_mut_condition_HR = "Mutation x Response Interaction"
                )
                
                for (i in 1:nrow(data)) {
                        gene <- data$Gene[i]
                        if (is.na(gene) || gene == "") next  # Skip if gene is NA or empty
                        cat("Analyzing gene:", gene, "in", model_name, "\n")
                        
                        # Determine HR types dynamically based on columns in the dataframe
                        hr_types <- intersect(
                                c("Expression_HR", "CNA_HR", "Condition_HR", "Interaction_cna_condition_HR", 
                                  "Mut_HR", "Interaction_mut_condition_HR"),
                                colnames(data)
                        )
                        
                        # Loop through all HR types
                        for (hr_type in hr_types) {
                                if (!is.na(data[[hr_type]][i]) && data[[hr_type]][i] > hr_threshold) {
                                        cat("HR type exceeding threshold:", hr_type, "\n")
                                        
                                        # Get HR value for the title
                                        hr_value <- round(data[[hr_type]][i], 2)
                                        
                                        # Determine grouping based on HR type
                                        if (hr_type == "Expression_HR") {
                                                median_expr <- median(processed_data[[paste0("expr_", gene)]], na.rm = TRUE)
                                                group <- ifelse(processed_data[[paste0("expr_", gene)]] > median_expr, "High", "Low")
                                                group <- factor(group, levels = c("High", "Low"))  # Ensure group is a factor
                                        } else if (hr_type == "CNA_HR") {
                                                group <- factor(processed_data[[paste0("cna_", gene)]], labels = c("CNA 0", "CNA 1"))
                                        } else if (hr_type == "Condition_HR") {
                                                group <- factor(processed_data$condition, labels = c("Resistant", "Sensitive"))
                                        } else if (hr_type == "Interaction_cna_condition_HR") {
                                                group <- rep(NA, nrow(processed_data))
                                                group[processed_data[[paste0("cna_", gene)]] == 1 & processed_data$condition == "sensitive"] <- "CNA 1: Response 1"
                                                group[processed_data[[paste0("cna_", gene)]] == 1 & processed_data$condition == "resistant"] <- "CNA 1: Response 0"
                                                group[processed_data[[paste0("cna_", gene)]] == 0 & processed_data$condition == "sensitive"] <- "CNA 0: Response 1"
                                                group[processed_data[[paste0("cna_", gene)]] == 0 & processed_data$condition == "resistant"] <- "CNA 0: Response 0"
                                                group <- factor(group, levels = c("CNA 1: Response 1", "CNA 1: Response 0", "CNA 0: Response 1", "CNA 0: Response 0"))
                                        } else if (hr_type == "Mut_HR") {
                                                group <- factor(processed_data[[paste0("mut_", gene)]], labels = c("Mutation 0", "Mutation 1"))
                                        } else if (hr_type == "Interaction_mut_condition_HR") {
                                                group <- rep(NA, nrow(processed_data))
                                                group[processed_data[[paste0("mut_", gene)]] == 1 & processed_data$condition == "sensitive"] <- "Mutation 1: Response 1"
                                                group[processed_data[[paste0("mut_", gene)]] == 1 & processed_data$condition == "resistant"] <- "Mutation 1: Response 0"
                                                group[processed_data[[paste0("mut_", gene)]] == 0 & processed_data$condition == "sensitive"] <- "Mutation 0: Response 1"
                                                group[processed_data[[paste0("mut_", gene)]] == 0 & processed_data$condition == "resistant"] <- "Mutation 0: Response 0"
                                                group <- factor(group, levels = c("Mutation 1: Response 1", "Mutation 1: Response 0", "Mutation 0: Response 1", "Mutation 0: Response 0"))
                                        }
                                        
                                        # Remove unused levels
                                        group <- droplevels(group)
                                        
                                        # Fit Kaplan-Meier survival curves
                                        surv_fit <- survfit(Surv(surv_obj[, 1], surv_obj[, 2]) ~ group)
                                        surv_diff <- survdiff(Surv(surv_obj[, 1], surv_obj[, 2]) ~ group)
                                        p_value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
                                        
                                        # Determine the maximum time for x-axis limit
                                        max_time <- max(surv_obj[, 1], na.rm = TRUE)
                                        
                                        # Add results to dataframe
                                        results <- rbind(results, data.frame(
                                                Gene = gene,
                                                HR_Type = hr_type_labels[[hr_type]],
                                                P_Value = p_value
                                        ))
                                        
                                        # Plot if p-value < 0.05 and plot = TRUE
                                        if (p_value < 0.05) {
                                                cat("Plotting Kaplan-Meier for gene:", gene, "with", hr_type, "\n")
                                                
                                                # Dynamically adjust legend labels
                                                legend_labels <- levels(group)
                                                
                                                km_plot <- ggsurvplot(
                                                        surv_fit, 
                                                        data = data.frame(surv_obj, group), 
                                                        title = paste(gene, "-", hr_type_labels[[hr_type]], "(HR =", hr_value, ")"),
                                                        legend.title = "Group",
                                                        legend.labs = legend_labels,
                                                        xlim = c(0, max_time),  # Adjust x-axis limit
                                                        ggtheme = theme_bw() + 
                                                                theme(
                                                                        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
                                                                        aspect.ratio = 0.612
                                                                ),
                                                        risk.table = F,
                                                        pval = TRUE
                                                )
                                                
                                                # Save plot to list
                                                plot_list[[paste(gene, hr_type, sep = "_")]] <- km_plot
                                                
                                                # Optionally plot if plot = TRUE
                                                if (plot) {
                                                        print(km_plot)
                                                }
                                        }
                                }
                        }
                }
                
                return(list(Results = results, Plots = plot_list))
        }
        
        # Analyze Model 1
        model1_results <- analyze_model(model1_data, "Model 1")
        
        # Analyze Model 2
        model2_results <- analyze_model(model2_data, "Model 2")
        
        return(list(Model1 = model1_results, Model2 = model2_results))
}



# Function to filter Kaplan-Meier results based on survival data and enrichment table
filter_kaplan_meier_results <- function(model1_results, model2_results, enrichment_table_path = NULL) {
        # Load and process enrichment table if provided
        enrichment_table <- NULL
        if (!is.null(enrichment_table_path)) {
                cat("Reading enrichment table...\n")
                enrichment_table <- read_excel(enrichment_table_path)
                colnames(enrichment_table) <- enrichment_table[1, ]  # Set the first row as column names
                enrichment_table <- enrichment_table[-1, ]  # Remove the first row from the dataset
                enrichment_table$Up_regulated <- strsplit(as.character(enrichment_table$Up_regulated), ", ")
                enrichment_table$Down_regulated <- strsplit(as.character(enrichment_table$Down_regulated), ", ")
                cat("Enrichment table successfully loaded and processed.\n")
        }
        
        # Helper function to find the group with the worst survival
        find_worst_group <- function(surv_data) {
                # Calculate the minimum survival probability for each group
                cat("Calculating worst group...\n")
                min_surv <- aggregate(surv ~ group, data = surv_data, FUN = min)
                cat("Minimum survival per group:\n")
                print(min_surv)
                
                # Determine the group(s) with the worst survival
                worst_groups <- min_surv$group[min_surv$surv == min(min_surv$surv)]
                
                if (length(worst_groups) > 1) {
                        cat("Multiple groups with the same minimum survival, resolving by time...\n")
                        relevant_data <- surv_data[surv_data$group %in% worst_groups, ]
                        relevant_data <- relevant_data[relevant_data$surv == min(relevant_data$surv), ]
                        worst_group <- relevant_data$group[which.min(relevant_data$time)]
                } else {
                        worst_group <- worst_groups
                }
                
                cat("Worst group determined:", worst_group, "\n")
                return(worst_group)
        }
        
        # Helper function to process each model's results
        process_model <- function(model_results, hr_type_filter, enrichment_table = NULL) {
                if (is.null(model_results) || length(model_results$Plots) == 0) {
                        cat("No plots found in the model. Skipping...\n")
                        return(list(Results = data.frame(), Plots = list()))
                }
                
                filtered_results <- data.frame()
                filtered_plots <- list()
                
                # Define HR types
                hr_types <- c("Expression_HR", "CNA_HR", "Condition_HR", 
                              "Interaction_cna_condition_HR", "Mut_HR", "Interaction_mut_condition_HR")
                
                for (plot_name in names(model_results$Plots)) {
                        cat("Processing plot:", plot_name, "\n")
                        km_plot <- model_results$Plots[[plot_name]]
                        surv_data <- km_plot$data.survplot
                        
                        # Extract gene and HR type from the plot name
                        parts <- strsplit(plot_name, "_", fixed = TRUE)[[1]]
                        gene <- parts[1]  # Gene name is the first part before "_"
                        hr_type <- paste(parts[-1], collapse = "_")  # HR type is the part after "_"
                        cat("Gene:", gene, "- HR Type:", hr_type, "\n")
                        
                        # Skip if HR type is not valid
                        if (!(hr_type %in% hr_types)) {
                                cat("Invalid HR type. Skipping...\n")
                                next
                        }
                        
                        # Identify the worst group
                        worst_group <- find_worst_group(surv_data)
                        
                        # Apply filtering logic based on HR type and enrichment table
                        if (hr_type == "Expression_HR" && !is.null(enrichment_table)) {
                                upregulated_genes <- unlist(enrichment_table$Up_regulated)
                                downregulated_genes <- unlist(enrichment_table$Down_regulated)
                                
                                if (worst_group == "High" && gene %in% upregulated_genes) {
                                        cat("Gene is upregulated and High group has worst survival. Adding to results.\n")
                                        filtered_results <- rbind(filtered_results, model_results$Results[model_results$Results$Gene == gene, ])
                                        filtered_plots[[plot_name]] <- km_plot
                                } else if (worst_group == "Low" && gene %in% downregulated_genes) {
                                        cat("Gene is downregulated and Low group has worst survival. Adding to results.\n")
                                        filtered_results <- rbind(filtered_results, model_results$Results[model_results$Results$Gene == gene, ])
                                        filtered_plots[[plot_name]] <- km_plot
                                }
                        }
                        
                        # Filter for CNA or Mutation HR types
                        if ((hr_type == "CNA_HR" || hr_type == "Mut_HR") && grepl("1", worst_group)) {
                                cat("Filtering CNA/Mutation HR. Group 1 has worst survival. Adding to results.\n")
                                filtered_results <- rbind(filtered_results, model_results$Results[model_results$Results$Gene == gene, ])
                                filtered_plots[[plot_name]] <- km_plot
                        }
                        
                        # Filter for Condition HR
                        if (hr_type == "Condition_HR" && grepl("Resistant", worst_group)) {
                                cat("Filtering Condition HR. Resistant group has worst survival. Adding to results.\n")
                                filtered_results <- rbind(filtered_results, model_results$Results[model_results$Results$Gene == gene, ])
                                filtered_plots[[plot_name]] <- km_plot
                        }
                        
                        # Filter for Interaction HR types
                        if ((hr_type == "Interaction_cna_condition_HR" && grepl("CNA 1: Response 0", worst_group)) || 
                            (hr_type == "Interaction_mut_condition_HR" && grepl("Mut 1: Response 0", worst_group))) {
                                cat("Filtering Interaction HR. Group has worst survival. Adding to results.\n")
                                filtered_results <- rbind(filtered_results, model_results$Results[model_results$Results$Gene == gene, ])
                                filtered_plots[[plot_name]] <- km_plot
                        }
                }
                
                return(list(Results = filtered_results, Plots = filtered_plots))
        }
        
        # Process results for Model 1
        cat("Processing Model 1...\n")
        model1_filtered <- process_model(
                model1_results,
                hr_type_filter = "CNA_HR",
                enrichment_table = enrichment_table
        )
        
        # Process results for Model 2
        cat("Processing Model 2...\n")
        model2_filtered <- process_model(
                model2_results,
                hr_type_filter = "Mut_HR",
                enrichment_table = enrichment_table
        )
        
        return(list(Model1 = model1_filtered, Model2 = model2_filtered))
}





# Function to create a horizontal bar plot for C-Index values with enhanced visualization
cindex_barplot_horizontal_v5_improved <- function(model_results, model_name, title_main) {
        if (nrow(model_results) == 0) {
                cat("No results to plot for", model_name, "\n")
                return(NULL)
        }
        
        library(RColorBrewer)
        library(ggplot2)
        
        # Extract pathway genes
        model_results$Pathway_Genes <- NA
        model_results$Pathway_Genes[model_results$Type == "Pathway"] <- gsub(".*\\((.*?)\\)", "\\1", model_results$Label[model_results$Type == "Pathway"])
        
        # Remove parentheses from Labels
        model_results$Label <- gsub("\\s*\\(.*?\\)", "", model_results$Label)
        
        # Separate data by type
        single_gene_results <- model_results[model_results$Type == "Single Gene", , drop = FALSE]
        combination_results <- model_results[model_results$Type == "Combination", , drop = FALSE]
        pathway_results <- model_results[model_results$Type == "Pathway", , drop = FALSE]
        
        # Assign colors to single genes if available
        single_gene_colors <- NULL
        if (nrow(single_gene_results) > 0) {
                single_gene_colors <- setNames(
                        colorRampPalette(brewer.pal(8, "Set1"))(nrow(single_gene_results)),
                        single_gene_results$Label
                )
        }
        
        # Assign colors to combinations if available
        combination_colors <- NULL
        if (nrow(combination_results) > 0) {
                combination_colors <- setNames(
                        colorRampPalette(brewer.pal(8, "Dark2"))(nrow(combination_results)),
                        combination_results$Label
                )
        }
        
        # Add Color and ColorGroup columns
        model_results$Color <- NA
        model_results$ColorGroup <- NA
        
        # Single Genes
        if (!is.null(single_gene_colors)) {
                model_results$Color[model_results$Type == "Single Gene"] <- single_gene_colors[model_results$Label[model_results$Type == "Single Gene"]]
                model_results$ColorGroup[model_results$Type == "Single Gene"] <- model_results$Label[model_results$Type == "Single Gene"]
        }
        
        # Combinations
        if (!is.null(combination_colors) && nrow(combination_results) > 0) {
                model_results$Color[model_results$Type == "Combination"] <- combination_colors[model_results$Label[model_results$Type == "Combination"]]
                model_results$ColorGroup[model_results$Type == "Combination"] <- model_results$Label[model_results$Type == "Combination"]
        }
        
        # Pathways
        if (nrow(pathway_results) > 0) {
                for (i in 1:nrow(pathway_results)) {
                        pathway_genes <- unlist(strsplit(pathway_results$Pathway_Genes[i], ", "))
                        
                        matched_combinations <- character(0)
                        
                        # Match combinations if available
                        if (!is.null(combination_colors) && nrow(combination_results) > 0) {
                                matched_combinations <- combination_results$Label[
                                        sapply(combination_results$Label, function(combo) {
                                                all(unlist(strsplit(combo, " \\+ ")) %in% pathway_genes)
                                        })
                                ]
                                
                                matched_combinations <- matched_combinations[order(-nchar(matched_combinations))]
                        }
                        
                        if (length(matched_combinations) > 0) {
                                # Use combination color
                                chosen_comb <- matched_combinations[1]
                                pathway_results$Color[i] <- combination_colors[chosen_comb]
                                pathway_results$ColorGroup[i] <- chosen_comb
                        } else {
                                # Fallback to a single gene color if available
                                if (!is.null(single_gene_colors)) {
                                        matched_genes <- intersect(pathway_genes, names(single_gene_colors))
                                        if (length(matched_genes) > 0) {
                                                chosen_gene <- matched_genes[1]
                                                pathway_results$Color[i] <- single_gene_colors[chosen_gene]
                                                pathway_results$ColorGroup[i] <- chosen_gene
                                        } else {
                                                pathway_results$Color[i] <- "grey70"
                                                pathway_results$ColorGroup[i] <- "Unknown"
                                        }
                                } else {
                                        pathway_results$Color[i] <- "grey70"
                                        pathway_results$ColorGroup[i] <- "Unknown"
                                }
                        }
                }
                
                # Merge pathway colors back
                model_results$Color[model_results$Type == "Pathway"] <- pathway_results$Color
                model_results$ColorGroup[model_results$Type == "Pathway"] <- pathway_results$ColorGroup
        }
        
        # Sort data
        sg <- model_results[model_results$Type == "Single Gene", , drop = FALSE]
        co <- model_results[model_results$Type == "Combination", , drop = FALSE]
        pa <- model_results[model_results$Type == "Pathway", , drop = FALSE]
        
        if (nrow(sg) > 0) {
                sg <- sg[order(sg$C_Index, decreasing = FALSE), ]
        }
        if (nrow(co) > 0) {
                co <- co[order(co$C_Index, decreasing = FALSE), ]
        }
        if (nrow(pa) > 0) {
                pa <- pa[order(pa$ColorGroup, pa$C_Index), ]
        }
        
        model_results <- rbind(sg, co, pa)
        
        model_results$Type <- factor(model_results$Type, levels = c("Single Gene", "Combination", "Pathway"))
        model_results$Label <- factor(model_results$Label, levels = model_results$Label)
        
        # Facet labels
        facet_labels <- c(
                "Single Gene" = "Single gene\nmodels",
                "Combination" = "Combined gene\nmodels",
                "Pathway" = "Pathway models"
        )
        
        unique_groups <- unique(model_results[, c("ColorGroup", "Color")])
        unique_groups <- unique_groups[!is.na(unique_groups$Color), , drop = FALSE]
        
        fill_values <- NULL
        if (nrow(unique_groups) > 0) {
                fill_values <- setNames(unique_groups$Color, unique_groups$ColorGroup)
        }
        
        # Create plot
        modern_plot <- ggplot(model_results, aes(x = C_Index, y = Label, fill = ColorGroup)) +
                geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.2) +
                geom_text(aes(label = sprintf("%.2f", C_Index), x = C_Index - 0.01),
                          hjust = 1, color = "white", size = 3, fontface = "bold") +
                (if (!is.null(fill_values)) scale_fill_manual(values = fill_values) else scale_fill_identity()) +
                scale_x_continuous(expand = c(0, 0)) +
                facet_grid(Type ~ ., scales = "free", space = "free_y", labeller = as_labeller(facet_labels)) +
                labs(
                        title = title_main,
                        x = "C-Index",
                        y = NULL
                ) +
                theme_minimal(base_size = 14) +
                theme(
                        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "#34495e"),
                        strip.text.y = element_text(size = 12, face = "bold", color = "#34495e", angle = 270, lineheight = 1.1),
                        strip.background = element_blank(),
                        axis.text.x = element_text(size = 10, face = "bold", color = "#2c3e50"),
                        axis.text.y = element_text(size = 10, face = "bold", color = "#2c3e50", margin = margin(r = 0)),
                        panel.grid.major.x = element_line(color = "gray85", size = 0.5),
                        panel.grid.major.y = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_blank(),
                        axis.ticks = element_blank(),
                        panel.border = element_blank(),
                        legend.position = "none"
                ) +
                coord_cartesian(clip = "off")
        
        return(modern_plot)
}



# Function to create a vertical bar plot for C-Index values with enhanced visualization
cindex_barplot_vertical_v5_improved <- function(model_results, model_name, title_main) {
        if (nrow(model_results) == 0) {
                cat("No results to plot for", model_name, "\n")
                return(NULL)
        }
        
        library(RColorBrewer)
        library(ggplot2)
        
        # Extract pathway genes
        model_results$Pathway_Genes <- NA
        model_results$Pathway_Genes[model_results$Type == "Pathway"] <- gsub(".*\\((.*?)\\)", "\\1", model_results$Label[model_results$Type == "Pathway"])
        
        # Clean parentheses from labels
        model_results$Label <- gsub("\\s*\\(.*?\\)", "", model_results$Label)
        
        # Separate data by type
        single_gene_results <- model_results[model_results$Type == "Single Gene", , drop = FALSE]
        combination_results <- model_results[model_results$Type == "Combination", , drop = FALSE]
        pathway_results <- model_results[model_results$Type == "Pathway", , drop = FALSE]
        
        # Assign colors to single genes
        single_gene_colors <- NULL
        if (nrow(single_gene_results) > 0) {
                single_gene_colors <- setNames(
                        colorRampPalette(brewer.pal(8, "Set1"))(nrow(single_gene_results)),
                        single_gene_results$Label
                )
        }
        
        # Assign colors to combinations
        combination_colors <- NULL
        if (nrow(combination_results) > 0) {
                combination_colors <- setNames(
                        colorRampPalette(brewer.pal(8, "Dark2"))(nrow(combination_results)),
                        combination_results$Label
                )
        }
        
        # Add Color and ColorGroup columns
        model_results$Color <- NA
        model_results$ColorGroup <- NA
        
        # Single Genes
        if (!is.null(single_gene_colors)) {
                model_results$Color[model_results$Type == "Single Gene"] <- single_gene_colors[model_results$Label[model_results$Type == "Single Gene"]]
                model_results$ColorGroup[model_results$Type == "Single Gene"] <- model_results$Label[model_results$Type == "Single Gene"]
        }
        
        # Combinations
        if (!is.null(combination_colors) && nrow(combination_results) > 0) {
                model_results$Color[model_results$Type == "Combination"] <- combination_colors[model_results$Label[model_results$Type == "Combination"]]
                model_results$ColorGroup[model_results$Type == "Combination"] <- model_results$Label[model_results$Type == "Combination"]
        }
        
        # Pathways
        if (nrow(pathway_results) > 0) {
                for (i in 1:nrow(pathway_results)) {
                        pathway_genes <- unlist(strsplit(pathway_results$Pathway_Genes[i], ", "))
                        
                        matched_combinations <- character(0)
                        if (!is.null(combination_colors) && nrow(combination_results) > 0) {
                                matched_combinations <- combination_results$Label[
                                        sapply(combination_results$Label, function(combo) {
                                                all(unlist(strsplit(combo, " \\+ ")) %in% pathway_genes)
                                        })
                                ]
                                
                                matched_combinations <- matched_combinations[order(-nchar(matched_combinations))]
                        }
                        
                        if (length(matched_combinations) > 0) {
                                chosen_comb <- matched_combinations[1]
                                pathway_results$Color[i] <- combination_colors[chosen_comb]
                                pathway_results$ColorGroup[i] <- chosen_comb
                        } else {
                                matched_genes <- character(0)
                                if (!is.null(single_gene_colors)) {
                                        matched_genes <- intersect(pathway_genes, names(single_gene_colors))
                                }
                                if (length(matched_genes) > 0) {
                                        chosen_gene <- matched_genes[1]
                                        pathway_results$Color[i] <- single_gene_colors[chosen_gene]
                                        pathway_results$ColorGroup[i] <- chosen_gene
                                } else {
                                        pathway_results$Color[i] <- "grey70"
                                        pathway_results$ColorGroup[i] <- "Unknown"
                                }
                        }
                }
                
                # Merge pathway colors back
                model_results$Color[model_results$Type == "Pathway"] <- pathway_results$Color
                model_results$ColorGroup[model_results$Type == "Pathway"] <- pathway_results$ColorGroup
        }
        
        # Sort data
        sg <- model_results[model_results$Type == "Single Gene", , drop = FALSE]
        co <- model_results[model_results$Type == "Combination", , drop = FALSE]
        pa <- model_results[model_results$Type == "Pathway", , drop = FALSE]
        
        if (nrow(sg) > 0) {
                sg <- sg[order(sg$C_Index, decreasing = TRUE), ]
        }
        if (nrow(co) > 0) {
                co <- co[order(co$C_Index, decreasing = TRUE), ]
        }
        if (nrow(pa) > 0) {
                pa <- pa[order(pa$ColorGroup, -pa$C_Index), ]
        }
        
        model_results <- rbind(sg, co, pa)
        
        model_results$Type <- factor(model_results$Type, levels = c("Single Gene", "Combination", "Pathway"))
        model_results$Label <- factor(model_results$Label, levels = model_results$Label)
        
        # Facet labels
        facet_labels <- c(
                "Single Gene" = "Single gene\nmodels",
                "Combination" = "Combined gene\nmodels",
                "Pathway" = "Pathway Models"
        )
        
        unique_groups <- unique(model_results[, c("ColorGroup", "Color")])
        unique_groups <- unique_groups[!is.na(unique_groups$Color), , drop = FALSE]
        
        fill_values <- NULL
        if (nrow(unique_groups) > 0) {
                fill_values <- setNames(unique_groups$Color, unique_groups$ColorGroup)
        }
        
        # Plot
        modern_plot <- ggplot(model_results, aes(x = Label, y = C_Index, fill = ColorGroup)) +
                geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.2) +
                geom_text(aes(label = sprintf("%.2f", C_Index), y = C_Index + 0.03),
                          vjust = 0, color = "black", size = 2.5, fontface = "bold") +
                (if (!is.null(fill_values)) scale_fill_manual(values = fill_values) else scale_fill_identity()) +
                scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                facet_grid(. ~ Type, scales = "free", space = "free_x", labeller = as_labeller(facet_labels)) +
                labs(
                        title = title_main,
                        x = NULL,
                        y = "C-Index"
                ) +
                theme_minimal(base_size = 14) +
                theme(
                        plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#34495e"),
                        strip.text = element_text(size = 9, face = "bold", color = "#34495e"),
                        strip.background = element_blank(),
                        axis.text.x = element_text(size = 8, face = "bold", color = "#2c3e50", angle = 30, hjust = 1),
                        axis.text.y = element_text(size = 6, face = "bold", color = "#2c3e50"),
                        panel.grid.major.y = element_line(color = "gray85", size = 0.5),
                        panel.grid.major.x = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(color = "gray50"),
                        axis.ticks = element_line(color = "gray50"),
                        panel.border = element_blank(),
                        legend.position = "none",
                        panel.spacing = unit(0.5, "lines") # Adjust facet spacing
                )
        
        return(modern_plot)
}


# Function to analyze genes and pathways for multiple survival models
analyze_genes_and_pathways_for_models <- function(
                gene_list,
                excluded_pathways,
                processed_data,
                surv_obj,
                single_gene_model,
                multi_cox_model,
                enrichment_table_path
) {
        # Read enrichment table
        enrichment_table <- read_excel(enrichment_table_path)
        colnames(enrichment_table) <- enrichment_table[1, ]
        enrichment_table <- enrichment_table[-1, ]
        
        # Helper function to process each model
        process_model <- function(model_name, gene_list, single_gene_model, multi_cox_model, surv_obj, cna_or_mut_col) {
                cat("Processing", model_name, "\n")
                
                # Single gene C-index values
                single_gene_cindex <- single_gene_model[[model_name]]
                
                # Filter gene_list to include only genes present in single_gene_model
                gene_list <- gene_list[gene_list %in% single_gene_cindex$Gene]
                if (length(gene_list) == 0) {
                        cat("No valid genes for", model_name, ". Skipping...\n")
                        return(data.frame(Label = character(), C_Index = numeric(), Type = character()))
                }
                cat("Filtered gene list for", model_name, ":", paste(gene_list, collapse = ", "), "\n")
                
                single_gene_cindex <- single_gene_cindex[single_gene_cindex$Gene %in% gene_list, c("Gene", "Concordance_Index")]
                
                # Pathways containing the genes
                multi_cox_pathways <- multi_cox_model[[model_name]]
                relevant_pathways <- multi_cox_pathways[
                        sapply(multi_cox_pathways$Genes, function(genes) any(gene_list %in% unlist(strsplit(genes, ", ")))) &
                                !(multi_cox_pathways$Pathway %in% excluded_pathways), 
                ]
                
                # Extract pathway C-index values and relevant genes
                relevant_pathways <- relevant_pathways[, c("Pathway", "Genes", "Concordance_Index")]
                relevant_pathways$Relevant_Genes <- sapply(relevant_pathways$Genes, function(genes) {
                        intersect(gene_list, unlist(strsplit(genes, ", ")))
                })
                
                # Handle single-gene cases
                if (length(gene_list) < 2) {
                        cat("Not enough genes for combinations in", model_name, ". Only single genes will be used.\n")
                        single_gene_cindex$Type <- "Single Gene"
                        colnames(single_gene_cindex) <- c("Label", "C_Index", "Type")
                        
                        relevant_pathways$Type <- "Pathway"
                        relevant_pathways$Label <- paste0(relevant_pathways$Pathway, " (", sapply(relevant_pathways$Relevant_Genes, paste, collapse = ", "), ")")
                        relevant_pathways <- relevant_pathways[, c("Label", "Concordance_Index", "Type")]
                        colnames(relevant_pathways) <- c("Label", "C_Index", "Type")
                        
                        return(rbind(single_gene_cindex, relevant_pathways))
                }
                
                # Generate combinations of genes for 2 to length(gene_list)
                all_combinations <- list()
                for (k in 2:length(gene_list)) {
                        all_combinations <- c(all_combinations, combn(gene_list, k, simplify = FALSE))
                }
                
                # Penalized Cox model for combinations
                comb_cindex <- data.frame(Combination = character(), C_Index = numeric(), stringsAsFactors = FALSE)
                
                for (comb in all_combinations) {
                        # Check if CNA/Mutation columns are included for all genes in comb
                        # Also build model matrix step by step
                        model_matrix <- NULL
                        valid_comb <- TRUE
                        
                        # Create formula components dynamically
                        expr_terms <- character()
                        cna_mut_terms <- character()
                        interaction_terms <- character()
                        
                        for (gene in comb) {
                                expr_gene <- paste0("expr_", gene)
                                cna_mut_col_name <- paste0(cna_or_mut_col, "_", gene)
                                if (!(expr_gene %in% colnames(processed_data)) || !(cna_mut_col_name %in% colnames(processed_data))) {
                                        valid_comb <- FALSE
                                        break
                                }
                                
                                # Add terms
                                expr_terms <- c(expr_terms, paste0("processed_data[['", expr_gene, "']]"))
                                cna_mut_terms <- c(cna_mut_terms, paste0("factor(processed_data[['", cna_mut_col_name, "']])"))
                                interaction_terms <- c(interaction_terms, 
                                                       paste0("factor(processed_data[['", cna_mut_col_name, "']]) : factor(processed_data$condition)"))
                        }
                        
                        if (!valid_comb) next
                        
                        # Build formula string
                        # Base: expr terms + cna_mut terms + condition + interactions
                        full_formula_str <- paste("~", 
                                                  paste(c(expr_terms, cna_mut_terms, "factor(processed_data$condition)", interaction_terms), collapse = " + "))
                        
                        # Construct model matrix
                        model_matrix <- model.matrix(as.formula(full_formula_str))
                        
                        # Fit penalized Cox model
                        fit <- tryCatch({
                                set.seed(73)
                                cv.glmnet(model_matrix, surv_obj, family = "cox", alpha = 0)
                        }, error = function(e) {
                                NULL
                        })
                        
                        if (!is.null(fit)) {
                                best_lambda <- fit$lambda.min
                                if (!is.null(best_lambda)) {
                                        cox_fit <- glmnet(model_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                        risk_score <- predict(cox_fit, newx = model_matrix, type = "link")
                                        cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                        
                                        comb_cindex <- rbind(comb_cindex, data.frame(Combination = paste(comb, collapse = " + "), C_Index = cindex))
                                }
                        }
                }
                
                # Combine results for plotting
                single_gene_cindex$Type <- "Single Gene"
                colnames(single_gene_cindex) <- c("Label", "C_Index", "Type")
                
                relevant_pathways$Type <- "Pathway"
                relevant_pathways$Label <- paste0(relevant_pathways$Pathway, " (", sapply(relevant_pathways$Relevant_Genes, paste, collapse = ", "), ")")
                relevant_pathways <- relevant_pathways[, c("Label", "Concordance_Index", "Type")]
                colnames(relevant_pathways) <- c("Label", "C_Index", "Type")
                
                if (nrow(comb_cindex) > 0) {
                        comb_cindex$Type <- "Combination"
                        colnames(comb_cindex) <- c("Label", "C_Index", "Type")
                }
                
                combined_results <- rbind(single_gene_cindex, relevant_pathways, comb_cindex)
                
                return(combined_results)
        }
        
        # Process Model 1
        model1_results <- process_model(
                model_name = "Model_1",
                gene_list = gene_list,
                single_gene_model = single_gene_model,
                multi_cox_model = multi_cox_model,
                surv_obj = surv_obj,
                cna_or_mut_col = "cna"
        )
        
        # Process Model 2
        model2_results <- process_model(
                model_name = "Model_2",
                gene_list = gene_list,
                single_gene_model = single_gene_model,
                multi_cox_model = multi_cox_model,
                surv_obj = surv_obj,
                cna_or_mut_col = "mut"
        )
        
        return(list(Model1 = model1_results, Model2 = model2_results))
}


# Function to analyze single genes and combinations for survival models
analyze_genes_for_combinations <- function(
                gene_list,
                processed_data,
                surv_obj,
                single_gene_model
) {
        # Helper function to process single genes and combinations
        process_model <- function(gene_list, processed_data, surv_obj,model, single_gene_model, cna_or_mut_col) {
                cat("Processing genes for", cna_or_mut_col, "\n")
                
                # Single gene C-index values
                single_gene_cindex <- single_gene_model[[model]]
                
                # Filter gene_list to include only genes present in single_gene_model
                gene_list <- gene_list[gene_list %in% single_gene_cindex$Gene]
                if (length(gene_list) == 0) {
                        cat("No valid genes for", cna_or_mut_col, ". Skipping...\n")
                        return(data.frame(Label = character(), C_Index = numeric(), Type = character()))
                }
                cat("Filtered gene list for", cna_or_mut_col, ":", paste(gene_list, collapse = ", "), "\n")
                
                single_gene_cindex <- single_gene_cindex[single_gene_cindex$Gene %in% gene_list, c("Gene", "Concordance_Index")]
                single_gene_cindex$Type <- "Single Gene"
                colnames(single_gene_cindex) <- c("Label", "C_Index", "Type")
                
                # Generate combinations of genes for 2 to length(gene_list)
                all_combinations <- list()
                for (k in 2:length(gene_list)) {
                        all_combinations <- c(all_combinations, combn(gene_list, k, simplify = FALSE))
                }
                
                # Penalized Cox model for combinations
                comb_cindex <- data.frame(Combination = character(), C_Index = numeric(), stringsAsFactors = FALSE)
                
                for (comb in all_combinations) {
                        # Check if CNA/Mutation columns are included for all genes in comb
                        # Also build model matrix step by step
                        model_matrix <- NULL
                        valid_comb <- TRUE
                        
                        # Create formula components dynamically
                        expr_terms <- character()
                        cna_mut_terms <- character()
                        interaction_terms <- character()
                        
                        for (gene in comb) {
                                expr_gene <- paste0("expr_", gene)
                                cna_mut_col_name <- paste0(cna_or_mut_col, "_", gene)
                                if (!(expr_gene %in% colnames(processed_data)) || !(cna_mut_col_name %in% colnames(processed_data))) {
                                        valid_comb <- FALSE
                                        break
                                }
                                
                                # Add terms
                                expr_terms <- c(expr_terms, paste0("processed_data[['", expr_gene, "']]"))
                                cna_mut_terms <- c(cna_mut_terms, paste0("factor(processed_data[['", cna_mut_col_name, "']])"))
                                interaction_terms <- c(interaction_terms, 
                                                       paste0("factor(processed_data[['", cna_mut_col_name, "']]) : factor(processed_data$condition)"))
                        }
                        
                        if (!valid_comb) next
                        
                        # Build formula string
                        # Base: expr terms + cna_mut terms + condition + interactions
                        full_formula_str <- paste("~", 
                                                  paste(c(expr_terms, cna_mut_terms, "factor(processed_data$condition)", interaction_terms), collapse = " + "))
                        
                        # Construct model matrix
                        model_matrix <- model.matrix(as.formula(full_formula_str))
                        
                        # Fit penalized Cox model
                        fit <- tryCatch({
                                set.seed(73)
                                cv.glmnet(model_matrix, surv_obj, family = "cox", alpha = 0)
                        }, error = function(e) {
                                NULL
                        })
                        
                        if (!is.null(fit)) {
                                best_lambda <- fit$lambda.min
                                if (!is.null(best_lambda)) {
                                        cox_fit <- glmnet(model_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                                        risk_score <- predict(cox_fit, newx = model_matrix, type = "link")
                                        cindex <- concordance.index(x = risk_score, surv.time = surv_obj[, 1], surv.event = surv_obj[, 2])$c.index
                                        
                                        comb_cindex <- rbind(comb_cindex, data.frame(Combination = paste(comb, collapse = " + "), C_Index = cindex))
                                }
                        }
                }
                
                if (nrow(comb_cindex) > 0) {
                        comb_cindex$Type <- "Combination"
                        colnames(comb_cindex) <- c("Label", "C_Index", "Type")
                }
                
                # Combine single gene and combination results
                combined_results <- rbind(single_gene_cindex, comb_cindex)
                
                return(combined_results)
        }
        
        # Process CNA-based models
        cna_results <- process_model(
                gene_list = gene_list,
                processed_data = processed_data,
                surv_obj = surv_obj,
                single_gene_model = single_gene_model,
                model='Model_1',
                cna_or_mut_col = "cna"
        )
        
        # Process Mutation-based models
        #mut_results <- process_model(
        #        gene_list = gene_list,
         #       processed_data = processed_data,
        #       surv_obj = surv_obj,
        #       single_gene_model = single_gene_model,
        #       model='Model_2',
        #       cna_or_mut_col = "mut"
        #  )
        
        return(list(CNA = cna_results))# , Mutation = mut_results))
}


# Function to plot C-Index results for separate analyses (CNA and Mutation)
plot_cindex_results_separate <- function(analyze_results, base_title = "C-Index Barplot") {
        library(ggplot2)
        library(RColorBrewer)
        
        cna_results <- analyze_results$CNA
        mut_results <- analyze_results$Mutation
        
        # Helper function to create a plot
        create_plot <- function(df, plot_title) {
                if (is.null(df) || nrow(df) == 0) {
                        return(NULL)
                }
                
                # Sort single genes and combinations
                sg <- df[df$Type == "Single Gene", , drop = FALSE]
                co <- df[df$Type == "Combination", , drop = FALSE]
                
                if (nrow(sg) > 0) {
                        sg <- sg[order(sg$C_Index, decreasing = TRUE), ]
                }
                if (nrow(co) > 0) {
                        co <- co[order(co$C_Index, decreasing = TRUE), ]
                }
                
                df <- rbind(sg, co)
                df$Label <- factor(df$Label, levels = df$Label)
                df$Type <- factor(df$Type, levels = c("Single Gene", "Combination"))
                
                # Assign colors: One palette for single genes, another for combinations
                sg_labels <- df$Label[df$Type == "Single Gene"]
                co_labels <- df$Label[df$Type == "Combination"]
                
                single_gene_colors <- if (length(sg_labels) > 0) {
                        setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(unique(sg_labels))),
                                 unique(sg_labels))
                } else NULL
                
                combination_colors <- if (length(co_labels) > 0) {
                        setNames(colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(co_labels))),
                                 unique(co_labels))
                } else NULL
                
                # Add ColorGroup column
                df$ColorGroup <- df$Label
                
                color_values <- c()
                if (!is.null(single_gene_colors)) {
                        color_values <- c(color_values, single_gene_colors)
                }
                if (!is.null(combination_colors)) {
                        color_values <- c(color_values, combination_colors)
                }
                
                # Create the plot
                p <- ggplot(df, aes(x = Label, y = C_Index, fill = ColorGroup)) +
                        geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.2) +
                        geom_text(aes(label = sprintf("%.2f", C_Index), y = C_Index + 0.03),
                                  vjust = 0, color = "black", size = 3.5, fontface = "bold") +
                        facet_grid(. ~ Type, scales = "free_x", space = "free_x") +
                        scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                        labs(
                                title = plot_title,
                                x = NULL,
                                y = "C-Index"
                        ) +
                        scale_fill_manual(values = color_values) +
                        theme_minimal(base_size = 14) +
                        theme(
                                plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "#34495e"),
                                strip.text = element_text(size = 13, face = "bold", color = "#34495e"),
                                strip.background = element_blank(),
                                axis.text.x = element_text(size = 10, face = "bold", color = "#2c3e50", angle = 30, hjust = 1),
                                axis.text.y = element_text(size = 10, face = "bold", color = "#2c3e50"),
                                panel.grid.major.y = element_line(color = "gray85", size = 0.5),
                                panel.grid.major.x = element_blank(),
                                panel.background = element_blank(),
                                axis.line = element_line(color = "gray50"),
                                axis.ticks = element_line(color = "gray50"),
                                panel.border = element_blank(),
                                legend.position = "none"
                        )
                
                return(p)
        }
        
        # Generate CNA plot
        cna_plot <- create_plot(cna_results, paste0(base_title, " - CNA"))
        # Generate Mutation plot
        mut_plot <- create_plot(mut_results, paste0(base_title, " - Mutation"))
        
        # Return both plots
        return(list(CNA = cna_plot, Mutation = mut_plot))
}


# Function to analyze genes and calculate hazard ratios (HR) for single and combined gene models
analyze_genes_with_hr <- function(gene_list, processed_data, surv_obj) {
        library(glmnet)
        library(survival)
        
        # Helper function to fit penalized Cox model and extract HRs
        fit_penalized_cox <- function(gene_combination, processed_data, surv_obj) {
                # Build formula components dynamically
                expr_terms <- character()
                cna_mut_terms <- character()
                interaction_terms <- character()
                
                for (gene in gene_combination) {
                        expr_gene <- paste0("expr_", gene)
                        cna_mut_col_name <- paste0("cna_", gene)
                        if (!(expr_gene %in% colnames(processed_data)) || !(cna_mut_col_name %in% colnames(processed_data))) {
                                stop(paste("Data for gene", gene, "is not available in processed_data"))
                        }
                        
                        expr_terms <- c(expr_terms, paste0("processed_data[['", expr_gene, "']]"))
                        cna_mut_terms <- c(cna_mut_terms, paste0("factor(processed_data[['", cna_mut_col_name, "']])"))
                        interaction_terms <- c(interaction_terms, 
                                               paste0("factor(processed_data[['", cna_mut_col_name, "']]) : factor(processed_data$condition)"))
                }
                
                # Build formula string
                full_formula_str <- paste("~", 
                                          paste(c(expr_terms, cna_mut_terms, "factor(processed_data$condition)", interaction_terms), collapse = " + "))
                
                # Construct model matrix
                model_matrix <- model.matrix(as.formula(full_formula_str))
                
                # Fit penalized Cox model
                fit <- tryCatch({
                        set.seed(73)
                        cv.glmnet(model_matrix, surv_obj, family = "cox", alpha = 0)
                }, error = function(e) {
                        NULL
                })
                
                if (is.null(fit)) {
                        stop("Penalized Cox model fitting failed")
                }
                
                best_lambda <- fit$lambda.min
                cox_fit <- glmnet(model_matrix, surv_obj, family = "cox", alpha = 0, lambda = best_lambda)
                
                # Extract coefficients and calculate hazard ratios
                coefficients <- as.matrix(cox_fit$beta)
                hr_values <- exp(coefficients)
                
                # Label terms appropriately
                term_labels <- rownames(coefficients)
                hr_data <- data.frame(
                        Term = term_labels,
                        HR = hr_values,
                        row.names = NULL
                )
                
                return(hr_data)
        }
        
        # Initialize results storage
        results <- list()
        
        # Process each gene individually
        for (gene in gene_list) {
                cat("Processing single gene model for:", gene, "\n")
                hr_data <- fit_penalized_cox(gene_combination = c(gene), processed_data = processed_data, surv_obj = surv_obj)
                hr_data$Model <- gene  # Label the model
                results[[gene]] <- hr_data
        }
        
        # Process combinations of genes
        if (length(gene_list) > 1) {
                for (k in 2:length(gene_list)) {
                        combinations <- combn(gene_list, k, simplify = FALSE)
                        for (comb in combinations) {
                                model_label <- paste(comb, collapse = " + ")
                                cat("Processing combination model for:", model_label, "\n")
                                hr_data <- fit_penalized_cox(gene_combination = comb, processed_data = processed_data, surv_obj = surv_obj)
                                hr_data$Model <- model_label  # Label the model
                                results[[model_label]] <- hr_data
                        }
                }
        }
        
        # Combine all results into a single data frame
        combined_results <- do.call(rbind, lapply(names(results), function(name) {
                data.frame(Model = name, results[[name]])
        }))
        
        return(combined_results)
}



format_terms <- function(term) {
        term <- gsub('processed_data\\[\\["expr_(.*?)"\\]\\]', "\\1 Expression", term)
        term <- gsub('factor\\(processed_data\\[\\["cna_(.*?)"\\]\\]\\)1', "\\1 CNA", term)
        term <- gsub('factor\\(processed_data\\$condition\\)resistant', "Response", term)
        term <- gsub("(.*?):", "\\1 x ", term) 
        term <- gsub('\\(Intercept\\)', "Intercept", term)
        return(term)
}



# Function to plot hazard ratios (HR) for multiple models
plot_hr_bar <- function(results_df, title_main = "Hazard Ratios for Models") {
        library(ggplot2)
        library(RColorBrewer)
        
        # Remove intercept terms
        results_df <- results_df[results_df$Term != "Intercept", ]
        
        # Keep the original term labels for coloring and ordering
        results_df$OriginalTerm <- results_df$Term
        
        # Adjust the display term for visualization purposes
        results_df$DisplayTerm <- with(results_df, ifelse(
                Model == "FYN", 
                gsub("FYN ", "", Term), # Remove "FYN " prefix for FYN model
                ifelse(Model == "TCF7L2", 
                       gsub("TCF7L2 ", "", Term), # Remove "TCF7L2 " prefix for TCF7L2 model
                       Term # Keep original terms for combined models
                )
        ))
        
        # Ensure models are in the desired order
        unique_models <- unique(results_df$Model)
        results_df$Model <- factor(results_df$Model, levels = unique_models)
        
        # Sort terms within each model in descending order of HR
        results_df <- do.call(rbind, lapply(split(results_df, results_df$Model), function(df) {
                df <- df[order(df$HR, decreasing = TRUE), ] # Sort by HR
                df$DisplayTerm <- factor(df$DisplayTerm, levels = df$DisplayTerm) # Factorize for display order
                df$OriginalTerm <- factor(df$OriginalTerm, levels = df$OriginalTerm) # Factorize for color
                return(df)
        }))
        
        # Assign unique colors to terms
        term_levels <- unique(results_df$OriginalTerm)
        term_colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(term_levels)), term_levels)
        
        # Create the plot
        p <- ggplot(results_df, aes(x = DisplayTerm, y = HR, fill = OriginalTerm)) +
                geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.2) +
                geom_text(aes(label = sprintf("%.2f", HR), y = HR + 0.03),
                          vjust = 0, color = "black", size = 2.5, fontface = "bold") +
                facet_grid(. ~ Model, scales = "free_x", space = "free_x") +
                scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                labs(
                        title = title_main,
                        x = NULL,
                        y = "Hazard Ratio (HR)"
                ) +
                scale_fill_manual(values = term_colors) +
                theme_minimal(base_size = 14) +
                theme(
                        plot.title = element_text(hjust = 0.5, size = 12, face = "bold", color = "#34495e"),
                        strip.text = element_text(size = 10, face = "bold", color = "#34495e"),
                        strip.background = element_blank(),
                        axis.text.x = element_text(size = 8, face = "bold", color = "#2c3e50", angle = 30, hjust = 1),
                        axis.text.y = element_text(size = 8, face = "bold", color = "#2c3e50"),
                        panel.grid.major.y = element_line(color = "gray85", size = 0.5),
                        panel.grid.major.x = element_blank(),
                        panel.background = element_blank(),
                        axis.line = element_line(color = "gray50"),
                        axis.ticks = element_line(color = "gray50"),
                        panel.border = element_blank(),
                        legend.position = "none"
                )
        
        return(p)
}



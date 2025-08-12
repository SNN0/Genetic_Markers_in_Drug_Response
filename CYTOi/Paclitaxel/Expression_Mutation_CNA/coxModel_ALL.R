

# Set the working directory to the location of CYTOi EMC data files
setwd('/path/to/your/CYTOi/Paclitaxel/Expression_Mutation_CNA')

#---------------------------------------------------------------------------------------------------------------
# Read sensitive CYTOi Paclitaxel EMC ALL sample and process expression, mutation, and CNA data

# Function to process TSV files for consistency
# process_tsv_files('TCGA - EC CYTOi Paclitaxel Resistant BRCA Samples CNA.tsv','TCGA - EC CYTOi Paclitaxel Resistant BRCA Samples Gene Expression.tsv')

# Define a conversion table for gene symbol changes
conversion_table <- data.frame(
        Old.Symbol = c("CCN1", "DSC3", "CELSR1", "DSP", "RIN2", "TCF4"),
        Converted.Symbol = c("PPP3CA", "DSC2", "ME2", "DSPP", "RASSF4", "TCF7L2")
)

# Process drug data for sensitive and resistant samples
CYTOi_Paclitaxel_ALL = process_drug_data(
        expr_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive Samples Gene Expression.tsv',
        cna_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive Samples CNA.tsv',
        expr_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant Samples Gene Expression.tsv',
        cna_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant Samples CNA.tsv',
        mut_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive Samples Mutation.tsv',
        mut_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant Samples Mutation.tsv',
        conversion_table = conversion_table
)

# Set condition factor levels for consistent modeling
CYTOi_Paclitaxel_ALL$condition = factor(CYTOi_Paclitaxel_ALL$condition, levels = c('sensitive', 'resistant'))

# Fetch clinical data for multiple TCGA projects
CYTOi_ALL_Clinical = get_multiple_clinical_data(tcga_project = c('TCGA-BRCA', 'TCGA-HNSC', 'TCGA-STAD', 'TCGA-LUAD', 'TCGA-UCS'))

# Check overlap between sample IDs in processed data and clinical data
sum(rownames(CYTOi_Paclitaxel_ALL) %in% CYTOi_ALL_Clinical$submitter_id) # 35 matches found

# Filter clinical data to retain only matching samples
CYTOi_ALL_Clinical_Paclitaxel = CYTOi_ALL_Clinical %>% filter(submitter_id %in% rownames(CYTOi_Paclitaxel_ALL))

# Create a survival object using clinical data
CYTOi_ALL_Paclitaxel_Surv = Surv(time = CYTOi_ALL_Clinical_Paclitaxel$overall_survival, event = CYTOi_ALL_Clinical_Paclitaxel$deceased)

#---------------------------------------------------------------------------------------------------------------
# Penalized Cox Regression Analysis

# Run penalized Cox regression models for the processed data
CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel = drug_small_sample_models_cond_penalized_with_pvalue_optimized(
        data = CYTOi_Paclitaxel_ALL,
        surv_obj = CYTOi_ALL_Paclitaxel_Surv
)

# Load enriched terms for pathway analysis
enrich_terms = read_excel('enriched_terms.xlsx')

# Extract unique genes from enriched terms
gen_list <- extract_unique_genes(enrich_terms)

# Process the penalized Cox models with the gene list
CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel <- process_models(CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel, gen_list)

# Save the processed Cox models
saveRDS(CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel, 'CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel.rds')

#---------------------------------------------------------------------------------------------------------------
# Multi-gene Cox Models

# Run multi-gene Cox models without combination terms
multiGeneCox_CYTOi_ALL = drug_multigene_models_cond_penalized(
        CYTOi_Paclitaxel_ALL,
        CYTOi_ALL_Paclitaxel_Surv,
        pathway_file = 'enriched_terms.xlsx',
        CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel,
        combination = F
)

# Save multi-gene Cox model results
saveRDS(multiGeneCox_CYTOi_ALL, 'CYTOi_EMC_ALL_Paclitaxel_MultiGeneCox.rds')



#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Survival Analysis

# Perform Kaplan-Meier analysis on filtered models
results <- kaplan_meier_analysis(
        model1_data = CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel$Filtered_Models$Model_1,
        model2_data = CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel$Filtered_Models$Model_2,
        surv_obj = CYTOi_ALL_Paclitaxel_Surv,
        processed_data = CYTOi_Paclitaxel_ALL,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results
saveRDS(results, 'CYTOi_EMC_ALL_Paclitaxel_KM.rds')

# Filter Kaplan-Meier results using enrichment table
results = filter_kaplan_meier_results(
        model1_results = results$Model1,
        model2_results = results$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results
saveRDS(results, 'CYTOi_EMC_ALL_Paclitaxel_KM_Filtered.rds')

#---------------------------------------------------------------------------------------------------------------
# Analyze Specific Genes and Pathways

# Analyze specific genes and pathways using a predefined list
results_FYN_MYL9_FLNB <- analyze_genes_and_pathways_for_models(
        gene_list = c("FYN", "MYL9", "FLNB"),
        excluded_pathways = c("Pathogenic Escherichia coli infection", 'Salmonella infection', 'Viral myocarditis', 'Shigellosis'),
        processed_data = CYTOi_Paclitaxel_ALL,
        surv_obj = CYTOi_ALL_Paclitaxel_Surv,
        single_gene_model = CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel,
        multi_cox_model = multiGeneCox_CYTOi_ALL,
        enrichment_table_path = "enriched_terms.xlsx"
)

# Generate barplot for analyzed models
model1_results = results_FYN_MYL9_FLNB$Model1
plot_model1_FYN_MYL9_FLNB <- cindex_barplot_vertical_v5_improved(
        model1_results,
        "Model_1",
        'CYTOi Paclitaxel All Samples'
)

# Save barplot
saveRDS(plot_model1_FYN_MYL9_FLNB, 'C-indexGraph/FYN_MYL9_FLNB_CYTOiPaclitaxel.rds')

#---------------------------------------------------------------------------------------------------------------
# Additional Analysis for Gene Combinations and HRs

# Analyze combinations of specific genes
results_FYN_LDHB = analyze_genes_for_combinations(
        gene_list = c('TCF7L2', 'LDHB'),
        processed_data = CYTOi_Paclitaxel_ALL,
        surv_obj = CYTOi_ALL_Paclitaxel_Surv,
        single_gene_model = CYTOi_EMC_ALL_Paclitaxel_Pen_CoxModel
)

# Plot results for gene combinations
results_FYN_LDHB_plot = plot_cindex_results_separate(
        results_FYN_LDHB,
        'Paclitaxel FYN-TCF7L2'
)

# Save CNA plot
saveRDS(results_FYN_LDHB_plot$CNA, 'C-indexGraph/FYN_LDHB_CYTOiPaclitaxel.rds')

#---------------------------------------------------------------------------------------------------------------
# Hazard Ratio Analysis

# Perform HR analysis for specific genes
A = analyze_genes_with_hr(
        gene_list = c('LDHB', 'TCF7L2'),
        processed_data = CYTOi_Paclitaxel_ALL,
        surv_obj = CYTOi_ALL_Paclitaxel_Surv
)

# Format terms for visualization
A$Term = format_terms(A$Term)

# Plot HR bar chart
LDHB_TCF_HR = plot_hr_bar(A)

# Save HR bar chart
saveRDS(LDHB_TCF_HR, 'C-indexGraph/LDHB_TCF_HR.rds')

#---------------------------------------------------------------------------------------------------------------
# Save Kaplan-Meier Plots for Specific Genes

TCF4_KM = CYTOi_EMC_ALL_Paclitaxel_KM_Filtered$Model1$Plots$TCF7L2_Interaction_cna_condition_HR
FLNB_KM = CYTOi_EMC_ALL_Paclitaxel_KM_Filtered$Model1$Plots$FLNB_Interaction_cna_condition_HR
LDHB_KM = CYTOi_EMC_ALL_Paclitaxel_KM_Filtered$Model1$Plots$LDHB_Interaction_cna_condition_HR
MMP7_KM = CYTOi_EMC_ALL_Paclitaxel_KM_Filtered$Model1$Plots$MMP7_Interaction_cna_condition_HR

saveRDS(TCF4_KM, 'C-indexGraph/TCF4_KM.rds')
saveRDS(FLNB_KM, 'C-indexGraph/FLNB_KM.rds')
saveRDS(LDHB_KM, 'C-indexGraph/LDHB_KM.rds')
saveRDS(MMP7_KM, 'C-indexGraph/MMP7_KM.rds')

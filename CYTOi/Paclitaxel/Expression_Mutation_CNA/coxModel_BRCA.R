
# Set the working directory to the location of CYTOi EMC data files for Paclitaxel BRCA analysis
setwd('/path/to/your/CYTOi/Paclitaxel/Expression_Mutation_CNA')

#---------------------------------------------------------------------------------------------------------------
# Read and process sensitive CYTOi Paclitaxel BRCA sample data (expression, mutation, CNA)


#process_tsv_files('TCGA - EMC CYTOi Paclitaxel Resistant BRCA Samples CNA.tsv','TCGA - EMC CYTOi Paclitaxel Resistant BRCA Samples Gene Expression.tsv')


# Process gene symbol conversion table
conversion_table <- data.frame(
        Old.Symbol = c("CCN1", "DSC3", "CELSR1", "DSP", "RIN2", "TCF4"),
        Converted.Symbol = c("PPP3CA", "DSC2", "ME2", "DSPP", "RASSF4", "TCF7L2")
)

# Process Paclitaxel BRCA data for sensitive and resistant samples
CYTOi_Paclitaxel_BRCA = process_drug_data(
        expr_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive BRCA Samples Gene Expression.tsv',
        cna_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive BRCA Samples CNA.tsv',
        expr_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant BRCA Samples Gene Expression.tsv',
        cna_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant BRCA Samples CNA.tsv',
        mut_sensitive_file = 'TCGA - EMC CYTOi Paclitaxel Sensitive BRCA Samples Mutation.tsv',
        mut_resistant_file = 'TCGA - EMC CYTOi Paclitaxel Resistant BRCA Samples Mutation.tsv',
        conversion_table = conversion_table
)

# Assign condition factor levels
CYTOi_Paclitaxel_BRCA$condition = factor(CYTOi_Paclitaxel_BRCA$condition, levels = c('sensitive', 'resistant'))

#---------------------------------------------------------------------------------------------------------------
# Clinical Data Processing

# Retrieve clinical data for TCGA-BRCA project
CYTOi_BRCA_Clinical = get_multiple_clinical_data(tcga_project = c('TCGA-BRCA'))

# Find matching samples between processed data and clinical data
sum(rownames(CYTOi_Paclitaxel_BRCA) %in% CYTOi_BRCA_Clinical$submitter_id) # Check the overlap (27 samples)

# Filter clinical data to match the processed samples
CYTOi_BRCA_Clinical_Paclitaxel = CYTOi_BRCA_Clinical %>% filter(submitter_id %in% rownames(CYTOi_Paclitaxel_BRCA))

# Create a survival object using clinical data
CYTOi_BRCA_Paclitaxel_Surv = Surv(
        time = CYTOi_BRCA_Clinical_Paclitaxel$overall_survival,
        event = CYTOi_BRCA_Clinical_Paclitaxel$deceased
)

#---------------------------------------------------------------------------------------------------------------
# Penalized Cox Regression Analysis

# Run penalized Cox regression models
CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel = drug_small_sample_models_cond_penalized_with_pvalue_optimized(
        data = CYTOi_Paclitaxel_BRCA,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv
)

# Load enriched terms and extract unique genes
enrich_terms = read_excel('enriched_terms.xlsx')
gen_list <- extract_unique_genes(enrich_terms)

# Process existing models with the gene list
CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel <- process_models(CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel, gen_list)

# Save the processed Cox models
saveRDS(CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel, 'CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel.rds')

#---------------------------------------------------------------------------------------------------------------
# Multi-Gene Cox Models

# Run multi-gene Cox models without combination terms
multiGeneCox_CYTOi_BRCA = drug_multigene_models_cond_penalized(
        CYTOi_Paclitaxel_BRCA,
        CYTOi_BRCA_Paclitaxel_Surv,
        pathway_file = 'enriched_terms.xlsx',
        CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel,
        combination = F
)

# Save multi-gene Cox model results
saveRDS(multiGeneCox_CYTOi_BRCA, 'CYTOi_EMC_BRCA_Paclitaxel_MultiGeneCox.rds')

# Run multi-gene Cox models with upregulated genes
multiGeneCox_CYTOi_BRCA_Upregulated = drug_upregulated_multigene_models_cond_penalized(
        CYTOi_Paclitaxel_BRCA,
        CYTOi_BRCA_Paclitaxel_Surv,
        pathway_file = 'enriched_terms.xlsx',
        CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel
)

# Save results for upregulated multi-gene Cox models
saveRDS(multiGeneCox_CYTOi_BRCA_Upregulated, 'CYTOi_EMC_BRCA_Paclitaxel_MultiGeneCox_Upregulated.rds')

#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Survival Analysis

# Perform Kaplan-Meier analysis for filtered models
results3 <- kaplan_meier_analysis(
        model1_data = CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel$Filtered_Models$Model_1,
        model2_data = CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel$Filtered_Models$Model_2,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv,
        processed_data = CYTOi_Paclitaxel_BRCA,
        hr_threshold = 0,
        plot = F
)


saveRDS(results3$Model1$Plots$CXCR4_Interaction_cna_condition_HR, 'yeniGraphlar/CXCR4_BRCA_KM.rds')


View(CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel$Filtered_Models$Model_1)

# Save Kaplan-Meier analysis results
saveRDS(results, 'CYTOi_EMC_BRCA_Paclitaxel_KM.rds')

# Filter Kaplan-Meier results using enrichment table
results <- filter_kaplan_meier_results(
        model1_results = results$Model1,
        model2_results = results$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results
saveRDS(results, 'CYTOi_EMC_BRCA_Paclitaxel_KM_Filtered.rds')

#---------------------------------------------------------------------------------------------------------------
# Plot Kaplan-Meier and Hazard Ratios for Genes and Pathways

# Analyze and visualize results for specific genes
results_FYN_MYL9_FLNB <- analyze_genes_and_pathways_for_models(
        gene_list = c("FYN", "TCF7L2"),
        excluded_pathways = c("Pathogenic Escherichia coli infection", 'Salmonella infection', 'Viral myocarditis', 'Shigellosis'),
        processed_data = CYTOi_Paclitaxel_BRCA,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv,
        single_gene_model = CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel,
        multi_cox_model = multiGeneCox_CYTOi_BRCA,
        enrichment_table_path = "enriched_terms.xlsx"
)

# Analyze and visualize results for specific genes
results_FYN_MYL9_FLNB <- analyze_genes_and_pathways_for_models(
        gene_list = c("FYN", "TCF7L2"),
        excluded_pathways = c('Salmonella infection','Melanogenesis','Arrhythmogenic right ventricular cardiomyopathy','Basal cell carcinoma','Thyroid cancer','Acute myeloid leukemia','Hepatocellular carcinoma','Breast cancer','Gastric cancer','Cushing syndrome','Prostate cancer','Colorectal cancer','Endometrial cancer','Human papillomavirus infection','Viral myocarditis',"Pathogenic Escherichia coli infection", 'Salmonella infection', 'Viral myocarditis', 'Shigellosis'),
        processed_data = CYTOi_Paclitaxel_BRCA,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv,
        single_gene_model = CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel,
        multi_cox_model = multiGeneCox_CYTOi_BRCA,
        enrichment_table_path = "enriched_terms.xlsx"
)

# Generate vertical barplot for the results
model1_results = results_FYN_MYL9_FLNB$Model1
plot_model1_FYN_MYL9_FLNB <- cindex_barplot_vertical_v5_improved(
        model1_results,
        "Model_1",
        'C-Index Barplot for CYTOi Paclitaxel:\nSingle Genes and Pathway Combinations'
)

# Save the plot
saveRDS(plot_model1_FYN_MYL9_FLNB, 'yeniGraphlar/FYN_TCF7L2_BRCA_Cindex.rds')


#---------------------------------------------------------------------------------------------------------------
# Additional Analysis for Gene Combinations and HRs

# Analyze combinations of specific genes
results_FYN_LDHB = analyze_genes_for_combinations(
        gene_list = c('FYN', 'TCF7L2'),
        processed_data = CYTOi_Paclitaxel_BRCA,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv,
        single_gene_model = CYTOi_EMC_BRCA_Paclitaxel_Pen_CoxModel
)

# Plot results for gene combinations
results_FYN_LDHB_plot = plot_cindex_results_separate(
        results_FYN_LDHB,
        'Paclitaxel FYN-TCF7L2'
)

# Save CNA plot
saveRDS(results_FYN_LDHB_plot$CNA, 'yeniGraphlar/FYN_TCF7L2_BRCA_Cindex.rds')

#---------------------------------------------------------------------------------------------------------------
# Hazard Ratio Analysis

# Perform HR analysis for specific genes
A = analyze_genes_with_hr(
        gene_list = c('FYN', 'CXCR4'),
        processed_data = CYTOi_Paclitaxel_BRCA,
        surv_obj = CYTOi_BRCA_Paclitaxel_Surv
)

# Format terms for visualization
A$Term = format_terms(A$Term)

# Plot HR bar chart
LDHB_TCF_HR = plot_hr_bar(A)

# Save HR bar chart
saveRDS(LDHB_TCF_HR, 'yeniGraphlar/CXCR4_FYN_HR_BRCA.rds')

#---------------------------------------------------------------------------------------------------------------



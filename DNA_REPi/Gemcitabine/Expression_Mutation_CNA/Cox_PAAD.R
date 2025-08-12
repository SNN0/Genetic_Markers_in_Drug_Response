
# Set the working directory for DNA_REPi Gemcitabine PAAD analysis
setwd('/path/to/your/DNA_REPi/Gemcitabine/Expression_Mutation_CNA')


#process_tsv_files('TCGA - EMC DNA_REPi Gemcitabine Resistant PAAD Samples CNA.tsv','TCGA - EMC DNA_REPi Gemcitabine Resistant PAAD Samples Gene Expression.tsv')


#---------------------------------------------------------------------------------------------------------------
# Read and process sensitive and resistant samples for DNA_REPi Gemcitabine PAAD data (expression, mutation, CNA)

# Process the input data files for sensitive and resistant samples
DNAi_Gemcitabine_PAAD = process_drug_data(
        expr_sensitive_file = 'TCGA - EMC DNA_REPi Gemcitabine Sensitive PAAD Samples Gene Expression.tsv',
        cna_sensitive_file = 'TCGA - EMC DNA_REPi Gemcitabine Sensitive PAAD Samples CNA.tsv',
        expr_resistant_file = 'TCGA - EMC DNA_REPi Gemcitabine Resistant PAAD Samples Gene Expression.tsv',
        cna_resistant_file = 'TCGA - EMC DNA_REPi Gemcitabine Resistant PAAD Samples CNA.tsv',
        mut_sensitive_file = 'TCGA - EMC DNA_REPi Gemcitabine Sensitive PAAD Samples Mutation.tsv',
        mut_resistant_file = 'TCGA - EMC DNA_REPi Gemcitabine Resistant PAAD Samples Mutation.tsv'
)

# Assign condition factor levels (sensitive and resistant)
DNAi_Gemcitabine_PAAD$condition = factor(DNAi_Gemcitabine_PAAD$condition, levels = c('sensitive', 'resistant'))

#---------------------------------------------------------------------------------------------------------------
# Clinical Data Processing

# Retrieve clinical data for TCGA-PAAD project
DNAi_PAAD_Clinical = get_multiple_clinical_data(tcga_project = c('TCGA-PAAD'))

# Check the overlap between processed data and clinical data
sum(rownames(DNAi_Gemcitabine_PAAD) %in% DNAi_PAAD_Clinical$submitter_id) # 

# Filter clinical data to retain only matching samples
DNAi_PAAD_Clinical_Gemcitabine = DNAi_PAAD_Clinical %>% filter(submitter_id %in% rownames(DNAi_Gemcitabine_PAAD))

# Create a survival object using clinical data
DNAi_PAAD_Gemcitabine_Surv = Surv(
        time = DNAi_PAAD_Clinical_Gemcitabine$overall_survival,
        event = DNAi_PAAD_Clinical_Gemcitabine$deceased
)

#---------------------------------------------------------------------------------------------------------------
# Penalized Cox Regression Analysis

# Run penalized Cox regression model for Gemcitabine
DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel = drug_small_sample_models_cond_penalized_with_pvalue_optimized(
        data = DNAi_Gemcitabine_PAAD,
        surv_obj = DNAi_PAAD_Gemcitabine_Surv
)

# Load enriched terms and extract unique genes for pathway analysis
enrich_terms = read_excel('enriched_terms.xlsx')
gen_list <- extract_unique_genes(enrich_terms)

# Process the existing models with the extracted gene list
DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel <- process_models(DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel, gen_list)

# Save the processed Cox models
saveRDS(DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel, 'DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel.rds')

#---------------------------------------------------------------------------------------------------------------
# Multi-Gene Cox Models

# Run multi-gene Cox models without combination terms
multiGeneCox_DNAi_PAAD = drug_multigene_models_cond_penalized(
        DNAi_Gemcitabine_PAAD,
        DNAi_PAAD_Gemcitabine_Surv,
        pathway_file = 'enriched_terms.xlsx',
        DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel,
        combination = F
)

# Save multi-gene Cox model results
saveRDS(multiGeneCox_DNAi_PAAD, 'DNAi_EMC_PAAD_Gemcitabine_MultiGeneCox.rds')

#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Survival Analysis

# Perform Kaplan-Meier analysis for filtered models
results <- kaplan_meier_analysis(
        model1_data = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel$Filtered_Models$Model_1,
        model2_data = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel$Filtered_Models$Model_2,
        surv_obj = DNAi_PAAD_Gemcitabine_Surv,
        processed_data = DNAi_Gemcitabine_PAAD,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results
saveRDS(results, 'DNAi_EMC_PAAD_Gemcitabine_KM.rds')

# Filter Kaplan-Meier results using enrichment table
results = filter_kaplan_meier_results(
        model1_results = results$Model1,
        model2_results = results$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results
saveRDS(results, 'DNAi_EMC_PAAD_Gemcitabine_KM_Filtered.rds')

# Perform Kaplan-Meier analysis for all genes
results_ALL <- kaplan_meier_analysis(
        model1_data = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel$Model_1,
        model2_data = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel$Model_2,
        surv_obj = DNAi_PAAD_Gemcitabine_Surv,
        processed_data = DNAi_Gemcitabine_PAAD,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results for all genes
saveRDS(results_ALL, 'DNAi_EMC_PAAD_Gemcitabine_KM_ALLGenes.rds')

# Filter Kaplan-Meier results for all genes using enrichment table
results_ALL = filter_kaplan_meier_results(
        model1_results = results_ALL$Model1,
        model2_results = results_ALL$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results for all genes
saveRDS(results_ALL, 'DNAi_EMC_PAAD_Gemcitabine_KM_ALLGenes_Filtered.rds')

#---------------------------------------------------------------------------------------------------------------
# Gene and Pathway Analysis for Specific Genes

results_PINK <- analyze_genes_and_pathways_for_models(
        gene_list = c("MGST2", "CDH1", "BMP4", "GPX8"),
        excluded_pathways = c("Amyotrophic lateral sclerosis", "Pathways of neurodegeneration - multiple diseases"),
        processed_data = DNAi_Gemcitabine_PAAD,
        surv_obj = DNAi_PAAD_Gemcitabine_Surv,
        single_gene_model = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel,
        multi_cox_model = multiGeneCox_DNAi_PAAD,
        enrichment_table_path = "enriched_terms.xlsx"
)

model1_results = results_PINK$Model1
model2_results = results_PINK$Model2

# Generate C-index barplots for Model_1 and Model_2
plot_model1_PINK <- cindex_barplot_vertical_v5_improved(model1_results, "Model_1", 'C-Index Barplot for DNAi Gemcitabine:\n Single Genes and Pathway Combinations')
plot_model2_PINK <- cindex_barplot_vertical_v5_improved(model2_results, "Model_2", 'C-Index Barplot for DNAi Gemcitabine:\n Single Genes and Pathway Combinations')

#---------------------------------------------------------------------------------------------------------------
# Gene Combination Analysis

results_MGST2_CDH1_BMP4_GPX8 = analyze_genes_for_combinations(
        gene_list = c("MGST2", "CDH1", "BMP4", "GPX8"),
        processed_data = DNAi_Gemcitabine_PAAD,
        surv_obj = DNAi_PAAD_Gemcitabine_Surv,
        single_gene_model = DNAi_EMC_PAAD_Gemcitabine_Pen_CoxModel
)

# Generate and save C-index barplots for gene combinations
results_MGST2_CDH1_BMP4_GPX8_plot = plot_cindex_results_separate(results_MGST2_CDH1_BMP4_GPX8, 'Gemcitabine PAAD Samples')
saveRDS(results_MGST2_CDH1_BMP4_GPX8_plot$CNA, 'Graph/MGST2_CDH1_BMP4_GPX8_Cindex.rds')

#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Plots for Specific Genes

# Generate and save Kaplan-Meier plots for specific genes
LSR_KM = DNAi_EMC_PAAD_Gemcitabine_KM_ALLGenes_Filtered$Model1$Plots$LSR_Interaction_cna_condition_HR
LSR_KM$plot$labels$title = "LSR - CNA x Response Interaction (HR = 4.15 )\nGemcitabine PAAD Samples"
LSR_KM$plot$theme$plot.title$size = 11
saveRDS(LSR_KM, 'Graph/LSR_KM_PAAD.rds')

TBC1D10C_KM = DNAi_EMC_PAAD_Gemcitabine_KM_ALLGenes$Model2$Plots$TBC1D10C_Interaction_mut_condition_HR
TBC1D10C_KM$plot$labels$title = "TBC1D10C - Mutation x Response (HR = 2.62 )\nGemcitabine PAAD Samples"
TBC1D10C_KM$plot$theme$plot.title$size = 11
saveRDS(TBC1D10C_KM, 'Graph/TBC1D10C_KM_PAAD.rds')



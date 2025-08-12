
# Set the working directory to the location of DNA_REPi Fluorouracil data files
setwd('/path/to/your/DNA_REPi/5-Fluorouracil/Expression_Mutation_CNA')

#---------------------------------------------------------------------------------------------------------------

#process_tsv_files('TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples CNA.tsv','TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples Gene Expression.tsv')


# Gene symbol conversion table for standardized processing
conversion_table <- data.frame(
        Old.Symbol = c("SH2D3C", "CCN1", "SH3D19", "RIN2", "LCP1", "SLA", "PMEPA1", "ENC1", "DSP", "HEBP1", "CNKSR3", "TCF4", "DLC1", "HOOK1", "SCG2"),
        Converted.Symbol = c("CHAT", "PPP3CA", "EBP", "RASSF4", "LPL", "SEPSECS", "STAG1", "CCL28", "DSPP", "NFE2L2", "MAGI1", "TCF7L2", "DYNLL1", "HK1", "MEN1")
)

# Process data for sensitive and resistant samples of Fluorouracil in STAD
DNAi_Fluoroucail_STAD = process_drug_data(
        expr_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive STAD Samples Gene Expression.tsv',
        cna_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive STAD Samples CNA.tsv',
        expr_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples Gene Expression.tsv',
        cna_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples CNA.tsv',
        mut_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive STAD Samples Mutation.tsv',
        mut_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples Mutation.tsv',
        conversion_table = conversion_table
)

# Assign factor levels for sensitive and resistant conditions
DNAi_Fluoroucail_STAD$condition = factor(DNAi_Fluoroucail_STAD$condition, levels = c('sensitive', 'resistant'))

#---------------------------------------------------------------------------------------------------------------
# Clinical Data Processing

# Retrieve clinical data for the TCGA-STAD project
DNAi_STAD_Clinical = get_multiple_clinical_data(tcga_project = c('TCGA-STAD'))

# Check the overlap between processed data and clinical data
sum(rownames(DNAi_Fluoroucail_STAD) %in% DNAi_STAD_Clinical$submitter_id) # 18 samples match

# Filter clinical data to retain only matching samples
DNAi_STAD_Clinical_Fluoroucail = DNAi_STAD_Clinical %>% filter(submitter_id %in% rownames(DNAi_Fluoroucail_STAD))

# Create a survival object using clinical data
DNAi_STAD_Fluoroucail_Surv = Surv(
        time = DNAi_STAD_Clinical_Fluoroucail$overall_survival,
        event = DNAi_STAD_Clinical_Fluoroucail$deceased
)

#---------------------------------------------------------------------------------------------------------------
# Penalized Cox Regression Analysis

# Run penalized Cox regression model for Fluorouracil in STAD
DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel = drug_small_sample_models_cond_penalized_with_pvalue_optimized(
        data = DNAi_Fluoroucail_STAD,
        surv_obj = DNAi_STAD_Fluoroucail_Surv
)

# Load enriched terms and extract unique genes for pathway analysis
enrich_terms = read_excel('enriched_terms.xlsx')
gen_list <- extract_unique_genes(enrich_terms)

# Process the existing models with the extracted gene list
DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel <- process_models(DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel, gen_list)

# Save the processed Cox models
saveRDS(DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel, 'DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel.rds')

#---------------------------------------------------------------------------------------------------------------
# Multi-Gene Cox Models

# Run multi-gene Cox models without combination terms
multiGeneCox_DNAi_STAD = drug_multigene_models_cond_penalized(
        DNAi_Fluoroucail_STAD,
        DNAi_STAD_Fluoroucail_Surv,
        pathway_file = 'enriched_terms.xlsx',
        DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel,
        combination = F
)

# Save multi-gene Cox model results
saveRDS(multiGeneCox_DNAi_STAD, 'DNAi_EMC_STAD_Fluorouracil_MultiGeneCox.rds')

#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Survival Analysis

# Perform Kaplan-Meier analysis for filtered models
results <- kaplan_meier_analysis(
        model1_data = DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel$Filtered_Models$Model_1,
        model2_data = DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel$Filtered_Models$Model_2,
        surv_obj = DNAi_STAD_Fluoroucail_Surv,
        processed_data = DNAi_Fluoroucail_STAD,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results
saveRDS(results, 'DNAi_EMC_STAD_Fluorouracil_KM.rds')

# Filter Kaplan-Meier results using enrichment table
results = filter_kaplan_meier_results(
        model1_results = results$Model1,
        model2_results = results$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results
saveRDS(results, 'DNAi_EMC_STAD_Fluorouracil_KM_Filtered.rds')

#---------------------------------------------------------------------------------------------------------------
# Single Gene and Pathway Analysis

# Analyze specific genes and pathways
results_PINK <- analyze_genes_and_pathways_for_models(
        gene_list = c("PINK1"),
        excluded_pathways = c("Amyotrophic lateral sclerosis", 'Pathways of neurodegeneration - multiple diseases'),
        processed_data = DNAi_Fluoroucail_STAD,
        surv_obj = DNAi_STAD_Fluoroucail_Surv,
        single_gene_model = DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel,
        multi_cox_model = multiGeneCox_DNAi_STAD,
        enrichment_table_path = "enriched_terms.xlsx"
)

# Generate and save barplot for specific gene results
model1_results = results_PINK$Model1
plot_model1_PINK <- cindex_barplot_horizontal_v5_improved(
        model1_results,
        "Model_1",
        'C-Index Barplot for DNAi 5-FU:\n Single Genes and Pathway Combinations'
)

# Save plot
saveRDS(plot_model1_PINK, 'Graph/PINK1_Cindex.rds')

#---------------------------------------------------------------------------------------------------------------
# KM Analysis and Combinations

# Analyze and visualize specific gene combinations
results_LYN_MYO1F_PINK1_TNFSF13B_PREX1 = analyze_genes_for_combinations(
        gene_list = c('LYN', 'MYO1F', 'PINK1', 'TNFSF13B', 'PREX1', 'MAGED1'),
        processed_data = DNAi_Fluoroucail_STAD,
        surv_obj = DNAi_STAD_Fluoroucail_Surv,
        single_gene_model = DNAi_EMC_STAD_Fluorouracil_Pen_CoxModel
)

# Plot the results for gene combinations
results_LYN_MYO1F_PINK1_TNFSF13B_PREX1_plot = plot_cindex_results_separate(
        results_LYN_MYO1F_PINK1_TNFSF13B_PREX1,
        '5-FU STAD Samples'
)

# Save plot for gene combinations
saveRDS(results_LYN_MYO1F_PINK1_TNFSF13B_PREX1_plot, 'Graph/LYN_MYO1F_PINK1_TNFSF13B_PREX1_Cindex.rds')

#---------------------------------------------------------------------------------------------------------------
# KM Plots for Specific Genes

# Generate and save KM plots for CDH1

CDH1_KM = DNAi_EMC_STAD_Fluorouracil_KM_ALLGenes_Filtered$Model2$Plots$CDH1_Mut_HR
CDH1_KM$plot$labels$title = "CDH1 - Mutation (HR = 1.18 )\n5-FU All Samples"
CDH1_KM$plot$theme$plot.title$size = 12
saveRDS(CDH1_KM, 'Graph/CDH1_KM_STADSamples.rds')

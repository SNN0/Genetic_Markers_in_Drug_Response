# Set the working directory to the location of DNA_REPi 5-Fluorouracil data files
setwd('/path/to/your/DNA_REPi/5-Fluorouracil/Expression_Mutation_CNA')



#process_tsv_files('TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples CNA.tsv','TCGA - EMC DNA_REPi Fluorouracil Resistant STAD Samples Gene Expression.tsv')

#---------------------------------------------------------------------------------------------------------------
# Read and process sensitive DNA_REPi 5-Fluorouracil data (expression, mutation, CNA)

# Load the gene expression, mutation, and CNA files for sensitive samples
expr_sensitive_file = read.table('TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Gene Expression.tsv', header = T)
mut_sensitive_file = read.table('TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Mutation.tsv', header = T)
cna_sensitive_file = read.table('TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples CNA.tsv', header = T)

# Remove unwanted samples from the data
olmayan = c('TCGA.AG.A01W', 'TCGA.AG.A01L', 'TCGA.AG.A01Y')
expr_sensitive_file = expr_sensitive_file %>% select(., -olmayan)
mut_sensitive_file = mut_sensitive_file %>% select(., -olmayan)
cna_sensitive_file = cna_sensitive_file %>% select(., -olmayan)

# Save the cleaned data to new files
write.table(expr_sensitive_file, 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Gene Expression-3.tsv', sep = '\t', row.names = F)
write.table(mut_sensitive_file, 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Mutation-3.tsv', sep = '\t', row.names = F)
write.table(cna_sensitive_file, 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples CNA-3.tsv', sep = '\t', row.names = F)

# Define a gene conversion table to standardize gene symbols
conversion_table <- data.frame(
        Old.Symbol = c("SH2D3C", "CCN1", "SH3D19", "RIN2", "LCP1", "SLA", "PMEPA1", "ENC1", "DSP", "HEBP1", "CNKSR3", "TCF4", "DLC1", "HOOK1", "SCG2"),
        Converted.Symbol = c("CHAT", "PPP3CA", "EBP", "RASSF4", "LPL", "SEPSECS", "STAG1", "CCL28", "DSPP", "NFE2L2", "MAGI1", "TCF7L2", "DYNLL1", "HK1", "MEN1")
)

# Process the data for both sensitive and resistant samples using the conversion table
DNAi_Fluoroucail_ALL = process_drug_data(
        expr_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Gene Expression-3.tsv',
        cna_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples CNA-3.tsv',
        expr_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant Samples Gene Expression.tsv',
        cna_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant Samples CNA.tsv',
        mut_sensitive_file = 'TCGA - EMC DNA_REPi Fluorouracil Sensitive Samples Mutation-3.tsv',
        mut_resistant_file = 'TCGA - EMC DNA_REPi Fluorouracil Resistant Samples Mutation.tsv',
        conversion_table = conversion_table
)

# Assign condition factor levels (sensitive and resistant)
DNAi_Fluoroucail_ALL$condition = factor(DNAi_Fluoroucail_ALL$condition, levels = c('sensitive', 'resistant'))

#---------------------------------------------------------------------------------------------------------------
# Clinical Data Processing

# Retrieve clinical data for multiple TCGA projects
DNAi_ALL_Clinical = get_multiple_clinical_data(tcga_project = c('TCGA-STAD', 'TCGA-PAAD', 'TCGA-ESCA', 'TCGA-READ'))

# Check the overlap between processed data and clinical data
sum(rownames(DNAi_Fluoroucail_ALL) %in% DNAi_ALL_Clinical$submitter_id) # 30 samples match

# Filter clinical data to retain only matching samples
DNAi_ALL_Clinical_Fluoroucail = DNAi_ALL_Clinical %>% filter(submitter_id %in% rownames(DNAi_Fluoroucail_ALL))

# Create a survival object using clinical data
DNAi_ALL_Fluoroucail_Surv = Surv(
        time = DNAi_ALL_Clinical_Fluoroucail$overall_survival,
        event = DNAi_ALL_Clinical_Fluoroucail$deceased
)

#---------------------------------------------------------------------------------------------------------------
# Penalized Cox Regression Analysis

# Run penalized Cox regression model for Fluorouracil
DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel = drug_small_sample_models_cond_penalized_with_pvalue_optimized(
        data = DNAi_Fluoroucail_ALL,
        surv_obj = DNAi_ALL_Fluoroucail_Surv
)

# Load enriched terms and extract unique genes for pathway analysis
enrich_terms = read_excel('enriched_terms.xlsx')
gen_list <- extract_unique_genes(enrich_terms)

# Process the existing models with the extracted gene list
DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel <- process_models(DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel, gen_list)

# Save the processed Cox models
saveRDS(DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel, 'DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel.rds')

#---------------------------------------------------------------------------------------------------------------
# Multi-Gene Cox Models

# Run multi-gene Cox models without combination terms
multiGeneCox_DNAi_ALL = drug_multigene_models_cond_penalized(
        DNAi_Fluoroucail_ALL,
        DNAi_ALL_Fluoroucail_Surv,
        pathway_file = 'enriched_terms.xlsx',
        DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel,
        combination = F
)

# Save multi-gene Cox model results
saveRDS(multiGeneCox_DNAi_ALL, 'DNAi_EMC_ALL_Fluorouracil_MultiGeneCox.rds')

#---------------------------------------------------------------------------------------------------------------
# Kaplan-Meier Survival Analysis

# Perform Kaplan-Meier analysis for filtered models
results <- kaplan_meier_analysis(
        model1_data = DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel$Filtered_Models$Model_1,
        model2_data = DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel$Filtered_Models$Model_2,
        surv_obj = DNAi_ALL_Fluoroucail_Surv,
        processed_data = DNAi_Fluoroucail_ALL,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results
saveRDS(results, 'DNAi_EMC_ALL_Fluorouracil_KM.rds')

# Filter Kaplan-Meier results using enrichment table
results = filter_kaplan_meier_results(
        model1_results = results$Model1,
        model2_results = results$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results
saveRDS(results, 'DNAi_EMC_ALL_Fluorouracil_KM_Filtered.rds')

# Perform Kaplan-Meier analysis for all genes
results_ALL <- kaplan_meier_analysis(
        model1_data = DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel$Model_1,
        model2_data = DNAi_EMC_ALL_Fluorouracil_Pen_CoxModel$Model_2,
        surv_obj = DNAi_ALL_Fluoroucail_Surv,
        processed_data = DNAi_Fluoroucail_ALL,
        hr_threshold = 1.1,
        plot = F
)

# Save Kaplan-Meier analysis results for all genes
saveRDS(results_ALL, 'DNAi_EMC_ALL_Fluorouracil_KM_ALLGenes.rds')

# Filter Kaplan-Meier results for all genes using enrichment table
results_ALL = filter_kaplan_meier_results(
        model1_results = results_ALL$Model1,
        model2_results = results_ALL$Model2,
        enrichment_table_path = 'enriched_terms.xlsx'
)

# Save filtered Kaplan-Meier results for all genes
saveRDS(results_ALL, 'DNAi_EMC_ALL_Fluorouracil_KM_ALLGenes_Filtered.rds')

#---------------------------------------------------------------------------------------------------------------
# Generate and save KM plots for BAMBI

BAMBI_KM = DNAi_EMC_ALL_Fluorouracil_KM_ALLGenes_Filtered$Model2$Plots$BAMBI_Mut_HR
BAMBI_KM$plot$labels$title = "BAMBI - Mutation (HR = 7.28 )\n5-FU All Samples"
BAMBI_KM$plot$theme$plot.title$size = 12
saveRDS(BAMBI_KM, 'Graph/BAMBI_KM_ALLSamples.rds')







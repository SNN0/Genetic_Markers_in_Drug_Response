# Analysis of Genetic Features' Impact on Drug Response and Survival



This project was developed to analyze the impact of genetic features (Copy Number Alteration - CNA, Mutation, Gene Expression) on patient survival in response to chemotherapeutic drugs (Paclitaxel, 5-Fluorouracil, Gemcitabine) in different cancer types (BRCA, STAD, PAAD). The analyses were performed using TCGA (The Cancer Genome Atlas) data, and the results are presented through statistical models and visualizations.



## Project Goal



The main goal of the project is to understand how specific genetic changes affect a patient group's sensitivity or resistance to a particular drug. These analyses aim to identify prognostic biomarkers for drug response and evaluate their impact on survival.



## Methodology



The project integrates genetic data (expression, CNA, mutation) and clinical data, and works with penalized Cox regression models and Kaplan-Meier analyses.



### Key Analysis Steps:



**Data Collection and Preparation:**

   - Gene expression, CNA, and mutation data are retrieved from TCGA.

   - Samples are segregated into **sensitive** and **resistant** groups based on drug response.

   - Gene symbols are standardized for consistency.

   - Patient survival data (overall survival time, vital status) is obtained from TCGA.

   - Genetic and clinical data are combined into a single dataset for each cancer/drug combination.

**Penalized Cox Regression Models:**

   * Penalized Cox regression models are created for each gene individually and for gene combinations using the `glmnet` R package.

   * Statistical metrics such as Hazard Ratio (HR), p-value, Confidence Interval (CI), Concordance Index (C-Index), and AIC are calculated for each model.

   * P-values are obtained using the bootstrapping method.

**Kaplan-Meier Survival Analysis:**

   * Kaplan-Meier curves are plotted to visualize survival rates for statistically significant genes or interactions.

**Visualization:**

   * Bar plots are created to compare C-Index values.

   * Bar plots are created to show Hazard Ratios.

   * All these plots are combined and arranged into a final figure format for publication.



## File Structure



* `functions.R`: This is the main library file containing all the core functions of the project. It includes functions for data processing, model building, and plotting. This file must be sourced before running the project.

* `coxModel_ALL.R`: Contains the main workflow for Paclitaxel analysis on a combined dataset of tumor samples (BRCA, HNSC, STAD, LUAD, UCS).

* `coxModel_BRCA.R`: Contains the main workflow for Paclitaxel analysis on BRCA (breast cancer) tumor samples.

* `Cox_ALL_5FU.R`: Contains the main workflow for the analysis of 5-Fluorouracil (5-FU) on a combined dataset of tumor samples (STAD, PAAD, ESCA, READ).

* `Cox_STAD.R`: Contains the main workflow for the analysis of 5-Fluorouracil (5-FU) on STAD (stomach cancer) tumor samples.

* `Cox_PAAD.R`: Contains the main workflow for Gemcitabine analysis on PAAD (pancreatic cancer) tumor samples.

* `Cox_ALL_Gem.R`: Contains the main workflow for Gemcitabine analysis on a combined dataset of tumor samples (PAAD, PCPG, LIHC, LUSC, SARC).

* `graphs.R`: Loads the `.rds` files of analysis results saved by the `coxModel\...R` files, adjusts the plot parameters, and combines them to create the final figures.



## Setup and Usage



**Prerequisites:** To run the project, you need a basic R programming environment and the libraries specified in the `functions.R` file (e.g., `TCGAbiolinks`, `glmnet`, `ggplot2`, `survminer`, `cowplot`, `readxl`). You can install the necessary libraries in the R console with the following command:

   ```R

   install.packages(c("TCGAbiolinks", "glmnet", "survminer", "ggplot2", "cowplot", "readxl", "gridExtra", "gplots", "ggpubr", "ggraph", "stringr", "boot", "dplyr", "tidyverse", "survcomp", "patchwork", "grid", "data.table"))

   ```

**File Paths:** At the beginning of each `coxModel\_...R` file, you need to update the data folder path (e.g., `setwd('/path/to/your/CYTOi/Paclitaxel/Expression_Mutation_CNA')`) to match your system.

**Execution:** To start the analyses, first source the `functions.R` file using the `source("functions.R")` command in R Studio or an R console. Then, run the relevant `coxModel_...R` file to start the analyses. Finally, run the `graphs.R` file to combine all the figures.

```R

# Load the functions.R file

source("functions.R")


# Run the relevant analysis file (e.g., for Paclitaxel)

source("coxModel\BRCA.R")

   
# Run the file that combines all figures

source("graphs.R")

```


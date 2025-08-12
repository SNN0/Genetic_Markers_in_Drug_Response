# Analysis of Genetic Features' Impact on Drug Response and Survival



This project was developed to analyze the impact of genetic features (Copy Number Alteration - CNA, Mutation, Gene Expression) on patient survival in response to chemotherapeutic drugs (Paclitaxel, 5-Fluorouracil, Gemcitabine) in different cancer types (BRCA, STAD, PAAD). The analyses were performed using TCGA (The Cancer Genome Atlas) data, and the results are presented through statistical models and visualizations.



## Project Goal



The main goal of the project is to understand how specific genetic changes affect a patient group's sensitivity or resistance to a particular drug. These analyses aim to identify prognostic biomarkers for drug response and evaluate their impact on survival.



## Methodology



The project integrates genetic data (expression, CNA, mutation) and clinical data, and works with penalized Cox regression models and Kaplan-Meier analyses.



### Key Analysis Steps:



**Data Collection and Preparation:**

&nbsp;   - Gene expression, CNA, and mutation data are retrieved from TCGA.

&nbsp;   - Samples are segregated into **sensitive** and **resistant** groups based on drug response.

&nbsp;   - Gene symbols are standardized for consistency.

&nbsp;   - Patient survival data (overall survival time, vital status) is obtained from TCGA.

&nbsp;   - Genetic and clinical data are combined into a single dataset for each cancer/drug combination.

2\.  \*\*Penalized Cox Regression Models:\*\*

&nbsp;   \* Penalized Cox regression models are created for each gene individually and for gene combinations using the `glmnet` R package.

&nbsp;   \* Statistical metrics such as Hazard Ratio (HR), p-value, Confidence Interval (CI), Concordance Index (C-Index), and AIC are calculated for each model.

&nbsp;   \* P-values are obtained using the bootstrapping method.

3\.  \*\*Kaplan-Meier Survival Analysis:\*\*

&nbsp;   \* Kaplan-Meier curves are plotted to visualize survival rates for statistically significant genes or interactions.

4\.  \*\*Visualization:\*\*

&nbsp;   \* Bar plots are created to compare C-Index values.

&nbsp;   \* Bar plots are created to show Hazard Ratios.

&nbsp;   \* All these plots are combined and arranged into a final figure format for publication.



\## File Structure



\* `functions.R`: This is the main library file containing all the core functions of the project. It includes functions for data processing, model building, and plotting. This file must be sourced before running the project.

\* `coxModel\_ALL.R`: Contains the main workflow for Paclitaxel analysis on a combined dataset of tumor samples (BRCA, HNSC, STAD, LUAD, UCS).

\* `coxModel\_BRCA.R`: Contains the main workflow for Paclitaxel analysis on BRCA (breast cancer) tumor samples.

\* `Cox\_ALL\_5FU.R`: Contains the main workflow for the analysis of 5-Fluorouracil (5-FU) on a combined dataset of tumor samples (STAD, PAAD, ESCA, READ).

\* `Cox\_STAD.R`: Contains the main workflow for the analysis of 5-Fluorouracil (5-FU) on STAD (stomach cancer) tumor samples.

\* `Cox\_PAAD.R`: Contains the main workflow for Gemcitabine analysis on PAAD (pancreatic cancer) tumor samples.

\* `Cox\_ALL\_Gem.R`: Contains the main workflow for Gemcitabine analysis on a combined dataset of tumor samples (PAAD, PCPG, LIHC, LUSC, SARC).

\* `graphs.R`: Loads the `.rds` files of analysis results saved by the `coxModel\_...R` files, adjusts the plot parameters, and combines them to create the final figures.



\## Setup and Usage



1\.  \*\*Prerequisites:\*\* To run the project, you need a basic R programming environment and the libraries specified in the `functions.R` file (e.g., `TCGAbiolinks`, `glmnet`, `ggplot2`, `survminer`, `cowplot`, `readxl`). You can install the necessary libraries in the R console with the following command:

&nbsp;   ```R

&nbsp;   install.packages(c("TCGAbiolinks", "glmnet", "survminer", "ggplot2", "cowplot", "readxl", "gridExtra", "gplots", "ggpubr", "ggraph", "stringr", "boot", "dplyr", "tidyverse", "survcomp", "patchwork", "grid", "data.table"))

&nbsp;   ```

2\.  \*\*File Paths:\*\* At the beginning of each `coxModel\_...R` file, you need to update the data folder path (e.g., `setwd('/path/to/your/CYTOi/Paclitaxel/Expression\_Mutation\_CNA')`) to match your system.

3\.  \*\*Execution:\*\* To start the analyses, first source the `functions.R` file using the `source("functions.R")` command in R Studio or an R console. Then, run the relevant `coxModel\_...R` file to start the analyses. Finally, run the `graphs.R` file to combine all the figures.

&nbsp;   ```R

&nbsp;   # Load the functions.R file

&nbsp;   source("functions.R")

&nbsp;   

&nbsp;   # Run the relevant analysis file (e.g., for Paclitaxel)

&nbsp;   source("coxModel\_BRCA.R")

&nbsp;   

&nbsp;   # Run the file that combines all figures

&nbsp;   source("graphs.R")

&nbsp;   ```


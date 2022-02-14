# Ropls-with-permutated-variable-selection
This script uses the Bioconductor package ropls by Etienne Thevenot to produce OPLS models. It performs variable selection using p(corr) and VIP. In addition to permutation post variable selection included in the ropls package also permutations pre variable selection with proceeding variable selection of every permutation resulting in p-values for R2 and Q2 including the variable selection procedure. It produces tables of all group comparisons including optional stratification by secondary ID.
Input is one matrix, one file with sampleID including the groups to be compared and one Configure file with desired settings.
The output is one html file per comparison of groups containing four selected models and one summary html file of all comparisons. Also tables with loadings are created.

# Installation
This package is tested on R 3.6.2, 4.0.0, 4.0.3, 4.0.5 and 4.1.2
Install needed packages BiocManager, Ropls, ggplot2, ggrepel, kableExtr, gridEstra, ggpubr, matrixStats, stringr, tryCatchLog, devtools and DescTools. These are installed if running the file dependencies.R using the following code
 
 ```
analysis_folder_name <- "Ropls-with-permutated-variable-selection"
source(paste(analysis_folder_name,"/dependencies.R", sep=""))
 ```

# Run example data
To run exampledata download folders and files in Ropls-with-permutated-variable-selection. Set workingdirectory to location of Ropls-with-permutated-variable-selection. To start the analysis run the file "Ropls_with_permutated_variable_selection_Run.R" using the following code

```
analysis_folder_name <- "Ropls-with-permutated-variable-selection"
source(paste(analysis_folder_name,"/Ropls_with_permutated_variable_selection_Run.R",sep=""))
```

# Run your own data
## Preparations:

1.	Prepare the datamatrix with subjectID in the first column and variable names in the first row. 
2.	All names has to be unique. Do not use the following symbols in subjectID or variable names; ?, $,%, ^, &, *, (, ),-,#, ?,,,<,>, /, |, \, [ ,] ,{, and };
3.	Missing values should be indicated by "", "NA" or "Inf"
4.	The datamatrix may contain numeric data, integers or categorical data as characters. Original categorical variable will be removed and replaced with dummy variables.
5. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. Also this file should not contain the symbols above and the samleID´s should agree with the subjectID´s in the datamatrix.
6. Both datamatrix and sampleID files should be saved as Tab delimited ("*.txt") files. 
7. Create a project folder for the analysis manually or use the Ropls-with-permutated-variable-selection-main folder. In this folder a subdirectory will be created called "outputR" where the results and the tables called model_table_to_analyse and the reordered_model_table_to_analyse describing which models are run will be stored. 
8. Create a folder and save your datamatrix and sampleID file here. Default name is "data_to_R_analysis" and is created by default in the project folder. You may also enter the path and foldername to your own input data folder.  
9. Save the Ropls_with_permutated_variable_selection_Configure_Get_Started.R and Ropls_with_permutated_variable_selection_Configure_Advanced.R in your project folder. The names of the files may be changed to your choice "Any_name_Advanced.R" and "Any_name_Get_Started.R".
10. Edit basic settings of parameters in file "Ropls_with_permutated_variable_selection_Configure_Get_Started.R" which contains the parameters that has to be entered including file names, folder names and column names discribing the groups to be compared. Default settings and advanced parameters may be altered in file "Ropls_with_permutated_variable_selection_Configure_Advanced.R".

## Run analysis:
When running Ropls_with_permutated_variable_selection_Run.R the data in Configure files are sent to file "Ropls_with_permutated_variable_selection_Models_of_each_comparison.Rmd" which renders one HTML file for each comparison containing score plots, loading plots, permutation pre and post varible selection plots and model statistics. When all comparisons have been performed "Ropls_with_permutated_variable_selection_Summary_of_models.Rmd" renders a summary HTML-file of all models containing tables of all selected models and all significant models.
Set working directory to the location of your project folder and run your analysis using the following code

```
setwd("enter path to where your project folder is stored")
directory_of_Ropls_with_permutated_variable_selection <- " enter path to where Ropls_with_permutated_variable_selection is stored"
analysis_folder_name <- "enter the name of your project folder" # 
source(paste(directory_of_Ropls_with_permutated_variable_selection,"/Ropls-with-permutated-variable-selection/Ropls_with_permutated_variable_selection_Run.R",sep=""))

```

## Features
The script generates four selected models: 
Model 1 uses p(corr) cutoff which is userset either to numeric value or corresponding to userset p[Pearson_pcorr_cutoff] and userset number of orthogonals with default 0. 
Model 2 uses p(corr) cutoff corresponding to userset p[Pearson_pcorr_cutoff] and number of orthogonals resulting in best performing model after variable selection. 
Model 3 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection. Only pcorr cutoff higher than p[Pearson_pcorr_cutoff] is used.
Model 4 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection while keeping the amount of variables down by only adding variables if a decrease in p(corr) cutoff less than userset Δpcorr cutoff results in an increase in Q2 of more than 1%.

### Best performing model during selection of amount of orthogonal variables and p(corr) cutoff
Best performing models are defined as models that after variable selection give high Q2, low difference between R2 and Q2 and also low p[Q2_perm_ post_vs] and low p[R2_perm_ post_vs]. The weight between low difference and low p[Q2_perm_ post_vs] is given by userset prefered_pR2_and_pQ2_permutated_post_vs with lower prefered_pR2_and_pQ2_permutated_post_vs giving more weight to p[Q2_perm_ post_vs] compared to high Q2 and low difference between R2 and Q2.

### Detailed descripting of method for selecting amount of orthogonal variables 
1)	Maximum number of orthogonals are userset by variable max_no_of_ortho with default setting 5 for Model2-4 and by setting no_of_ortho with default 0 for Model1.
2)	Selects models with p[R2_perm_ post_vs] and p[Q2_perm_ post_vs]< prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
3)	Selects models with max Q2 as long as adding an orthogonal increases Q2 more than 1%.
4)	If no model is found p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun.
5)	If there is more than one model with the same amount of orthogonals and the same Q2, the model with lowest amount of orthogonals after variable selection is chosen.

### Detailed descripting of method for selecting p(corr) cutoff
1)	From the start, all variables with |p(corr)|< p[Pearson_pcorr_cutoff] are removed.
2)	Amount of orthogonals selected when using pcorr cutoff set at p[Pearson_pcorr_cutoff] from Model2.
3)	Selects pcorr resulting in models with p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] < prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
4)	Model3: Selects pcorr resulting in models with max Q2 allowing for a decrease in Q2 by 1% if amount of variables are decreased.
5)	Model4: Selects pcorr resulting in models with max Q2 as long as adding variables by decreasing pcorr cutoff more than userset pcorr_diff increases Q2 more than 1%.
6)	If no model is found p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun. 

## Model statistics
### R2 and Q2
R2 is derived from the ropls package and is a measure of the fitness of the model and shows how much of the variation is explained by the model with a maximum value of 1. Q2 is a measure of the predictability of the model determined by 7 fold crossvalidation which is the default in the ropls package. The package gives R2 and Q2 pre and post variable selection.

### Permutations
In addition to the p[R2_perm_ post_vs] and p[Q2_perm_ post_vs] which is calculated by default by the ropls package this package also calculates p[R2_perm_ pre_vs] and p[Q2_perm_ pre_vs]. This includes randomisation of subject respons lables followed by variable selection and fitting of model post variable selection to obtain R2 and Q2 for permutated models. The same amount of orthogonals is used in the permutated models as in the unpermutated model under investigation. R2 and Q2 for the permutated models are compared to R2 and Q2 for the unpermutated models to generate p[R2_perm_ pre_vs] and p[Q2_perm_ pre_vs]. The number of permutations performed is userset by setting parameter no_permutations_post_vs, no_permutations_post_vs_selected_models and no_permutations_pre_vs with default 20.

## RMSE
RMSE is derived using ropls package and is a measure from the ropls package square root of the mean error between actual and predicted responses.


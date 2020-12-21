# Ropls-with-permutated-variable-selection
This script uses the package ropls by Etienne Thevenot to produce OPLS models. It performs variable selection using p(corr) and VIP. In addition to permutation after variable selection included in the ropls package also permutations before variable selection with proceeding variable selection of every permutation resulting in p-values for R2 and Q2 including the variable selection procedure. It produces tables of all group comparisons including optional stratification by secondary ID.
Input is one matrix, one file with sampleID including the groups to be compared and one Configure file with desired settings.
The output is one html file per comparison of groups and one summary html file of all comparisons. Also tables with loadings are created.

# Run example data
To run exampledata download folders and files in Ropls-with-permutated-variable-selection and run the file Ropls_with_permutated_variable_selection_Run.R

# Run your own data
## Preparations:
1.	Prepare the datamatrix with subjectID in the first column and variable names in the first row. 
2.	All names has to be unique. Do not use the following symbols in subjectID or variable names; ?, $,%, ^, &, *, (, ),-,#, ?,,,<,>, /, |, \, [ ,] ,{, and };
3.	Missing values should be indicated by “”, “NA” or “Inf”
4.	The datamatrix may contain numeric data, integers or categorical data as characters. Original categorical variable will be removed and replaced with dummy variables.
2. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. Also in this file there should be no symbols and the samleID´s should agree with the subjectID´s in the datamatrix.
3. Both datamatrix and sampleID files should be saved as Tab delimited (*.txt”) files. 
4. Create a project folder for the analysis manually. In this folder three subdirectories will be created automatically. One folder with default name data_to_R_analysis where your datamatrix and sampleID file should be saved, a second folder with default name Rscripts for the R-scripts and a third folder with default name outputR for the results and this is also where the model file called model_table_to_analyse and the reordered_model_table_to_analyse will be created. Model
5. Run the file with dependencies installing BiocManager, Ropls, ggplot2, ggrepel, kableExtr, gridEstra, ggpubr, matrixStats, stringr, tryCatchLog 

## Run analysis:
Enter desired settings of parameters in file “Ropls_with_permutated_variable_selection_Configure_Get_Started.R” which contains parameters that has to be entered. Default settings and advanced parameters may be altered in file Ropls_with_permutated_variable_selection_Configure_Advanced.R. When running Ropls_with_permutated_variable_selection_Run.R the data in Configure files are sent to file “Ropls_with_permutated_variable_selection_Models_of_each_comparison.Rmd” which renders one HTML file for each comparison. When all comparisons have been performed “Ropls_with_permutated_variable_selection_Summary_of_models.Rmd renders a summary HTML-file of all models containing tables of all selected models and all significant models.

## Features
The script generates three different models: 
Model 1 uses p(corr) cutoff which is userset either to numeric value or corresponding to userset p(pearson) cutoff and userset number of orthogonals with default 0. 
Model 2 uses p(corr) cutoff corresponding to userset p(pearson) cutoff and number of orthogonals resulting in best performing model after variable selection. 
Model 3 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection. Only pcorr cutoff higher that p(pearson) cutoff is used.
Model 4 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection while keeping the amount of variables down by only adding variables if a decrease in p(corr) cutoff less than userset Δpcorr cutoff results in an increase in Q2 of more than 1%.

### Choosing best performing model during selection fo amount of orthogonal variables and p(corr) cutoff
Best performing models are models that after variable selection give high Q2, low difference between R2 and Q2 and also low pQ2(perm_post_v.s.) and low pR2(perm after vs). The weight between low difference and low p(Q2_perm_after_vs) is given by userset prefered_pR2_and_pQ2_permutated_post_vs with lower prefered_pR2_and_pQ2_permutated_post_vs giving more weight to pQ2(perm_post_v.s.) compared to high Q2 and low difference between R2 and Q2.

### Detailed descripting of method for selecting amount of orthogonal variables 
1)	Maximum number of orthogonals are userset by variable max_no_of_ortho with default setting 5 for Model2-4 and by setting no_of_ortho with default 0 for Model1.
2)	Selects models with pR2(perm after vs) and pQ2(perm after vs)< prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
3)	Selects models with max Q2 as long as adding an orthogonal increases Q2 more than 1%.
4)	If no model is found pR2(perm after vs) and pQ2(perm after vs) limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun.
5)	If there is more than one model with the same amount of orthogonals and the same Q2, the model with lowest amount of orthogonals after variable selection is chosen.

### Detailed descripting of method for selecting p(corr) cutoff
1)	From the start, all variables with |p(corr)|< p[Pearson_pcorr_cutoff] are removed.
2)	Amount of orthogonals selected when using pcorr cutoff set at p[Pearson_pcorr_cutoff] from Model2.
3)	Selects pcorr resulting in models with pR2(perm after vs) and pQ2(perm after vs) < prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
4)	Model3: Selects pcorr resulting in models with max Q2 allowing for a decrease in Q2 by 1% if amount of variables are decreased.
5)	Model4: Selects pcorr resulting in models with max Q2 as long as adding variables by decreasing pcorr cutoff more than userset pcorr_diff increases Q2 more than 1%.
6)	If no model is found pR2(perm after vs) and pQ2(perm after vs) limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun. 

### Permutations



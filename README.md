# Ropls-with-permutated-variable-selection
This script uses the package ropls by Etienne Thevenot to produce OPLS models with variable selection using p(corr) and VIP. In addition to the standard permutation after variable selection also permutations before variable selection with proceeding variable selection of every permutation resulting in a p-value over the variable selection procedure. It produces tables of all group comparisons including optional stratification by secondary ID.
Input is one matrix, one file with sampleID including the groups to be compared and filling in Enter_data.R file with desired settings.
The output is one html file per comparison of groups and one summary html file of all comparisons. Also tables with loadings are created.

Preparations:
1.	Prepare the datamatrix with subjectID in the first column and variable names in the first row. 
2.	Be careful that all names are unique. Do not use the following symbols in subjectID or variable names; ?, $,%, ^, &, *, (, ),-,#, ?,,,<,>, /, |, \, [ ,] ,{, and };
3.	Missing values should be indicated by “”, “NA” or “Inf”
4.	The datamatrix may contain numeric data, integers or categorical data as characters. Original categorical variable will be removed and replaced with dummy variables.
2. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. Also in this file there should be no symbols and the samleID´s should agree with the subjectID´s in the datamatrix.
3. Both datamatrix and sampleID files should be saved as Tab delimited (*.txt”) files. 
4. Create a project folder for the analysis manually. In this folder three subdirectories will be created automatically. One folder called data_to_R_analysis where your datamatrix and sampleID file should be saved and this is also where the model file called model_table_to_analyse and the reordered_model_table_to_analyse is created, a second folder for the R-scripts and a third folder for the results.
5. Run the file with dependencies installing 

Performance:
Fill in the desired settings of parameters in section “#Enter data to send to Rmarkdown:” in file “Ropls_with_permutated_variable_selection_Enter_data.R” . Data that has to be edited from the default values is indicated by a “*”. Run the code acoording to instructions the notes in the file. Upload the datamatrix and the sampleID files when requested in section “#Upload files manually”. When running section “#Send data to Rmarkdown file to create html for each comparison” the data entered is sent to file “Ropls_with_permutated_variable_selection_Models_of_each_comparison.Rmd” which renders one HTML file for each comparison. When all comparisons have been performed “Ropls_with_permutated_variable_selection_Summary_of_models.Rmd renders a summary HTML-file of all models.

Features
The script generates three different models: 
Model 1 uses p(corr) cutoff which is userset either to numeric value or corresponding to userset p(pearson) cutoff and userset number of orthogonals with default 0. 
Model 2 uses p(corr) cutoff corresponding to userset p(pearson) cutoff and number of orthogonals resulting in best performing model after variable selection. 
Model 3 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection. Only pcorr cutoff higher that p(pearson) cutoff is used.
Model 4 uses both p(corr) cutoff and number of orthogonals resulting in best performing model after variable selection while keeping the amount of variables down by only adding variables if a decrease in p(corr) cutoff less than userset Δpcorr cutoff results in more than 1% higher Q2.

Best performing models are models that after variable selection give high Q2, low difference between R2 and Q2 and also low pQ2(perm_post_v.s.) and low pR2(perm after vs). The weight between low difference and low p(Q2_perm_after_vs) is given by userset prefered_pR2_and_pQ2_permutated_post_vs 

Method for selecting amount of orthogonal variables 
1)	Maximum number of orthogonals are userset by variable max_no_of_ortho.
2)	Selects models with pR2(perm after vs) and pQ2(perm after vs)< prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
3)	Selects models with max Q2 as long as adding an orthogonal increases Q2 more than 1%.
4)	If no model is found pR2(perm after vs) and pQ2(perm after vs) limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun.
5)	If models with the same amount of orthogonals and the same Q2 is found the model with lowest amount of orthogonals after variable selection is chosen.

Method for selecting p(corr) cutoff
1)	From the start, remove all variables with |p(corr)|< p[Pearson_pcorr_cutoff] 
2)	Amount of orthogonals selected when using pcorr cutoff set at p[Pearson_pcorr_cutoff] from Model2
3)	Selects pcorr resulting in models with pR2(perm after vs) and pQ2(perm after vs) < prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 < 0.2
4)	Selects pcorr resulting in models with max Q2 allowing for a decrease in Q2 by 1% if amount of variables are decreased. Applied for Model3.
5)	Selects pcorr resulting in models with max Q2 as long as adding variables by decreasing pcorr cutoff more than userset pcorr_diff increases Q2 more than 1%. Applied for Model4.
6)	If no model is found pR2(perm after vs) and pQ2(perm after vs) limit is increased by prefered_pR2_and_pQ2_permutated_post_vs and diff between R2 and Q2 is increased by 0.1 and selection is rerun. 

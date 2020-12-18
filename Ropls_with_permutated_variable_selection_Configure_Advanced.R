
# Set Advanced parameter settings
## Foldernames
foldername_Rmarkdownfiles <- "Rscript" #"foldername" Folder is automatically or manually created. This is the folder for the scripts .R and .rmd files.
foldername_of_input_matrix_and_sampleID <- "data_to_R_analysis" #"foldername" Folder is automatically created. This is the folder for datamatrix and sampleID files
foldername_output_reports <- "outputR" # "foldername"  Folder is automatically created below. In this folder one html file for each group comparison and a summary html file of all comparisons will be created.
foldername_function_file <- foldername_Rmarkdownfiles # foldername_Rmarkdownfiles or "foldername". This is the folder for the .R and .Rmd files.

## File names
filename_Rmarkdownfile_each_model <- "Ropls_with_permutated_variable_selection_Models_of_each_comparison.Rmd" # "filename.Rmd"
filename_Rmarkdownfile_summary <- "Ropls_with_permutated_variable_selection_Summary_of_models.Rmd" # "filename.Rmd"
filename_function_file <- "Ropls_with_permutated_variable_selection_Functions.R" # "filename.R" All functions used are saved in this file.

#Missing value tolerance
filter_percent_in_each_group <- 25 # Numeric. Missing value tolerance in each group which are compared.

##Variable selection
pcorr_cutoff_Model1_joint_models <- "according to p-value" #P(corr) cutoff in model1 for joint models. Either enter a value with p(corr)cutoff for all joint comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
pcorr_cutoff_Model1_stratified_models <- "according to p-value" #P(corr) cutoff in model1 for stratified models. Either enter a value with p(corr)cutoff for all stratified comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff

#Amount of orthogonals
no_of_ortho_pre_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for joint models. Either enter a value for nonstratified comparisons or a vector containing the number of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_pre_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for joint models. Either enter a value for joint comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
max_no_of_ortho_pre_vs_in_Model2 <- 5 # Max number of orthogonal variables in model pre variable selection of Model2 with pcorrcutoff according to pvalue
max_no_of_ortho_pre_vs_in_Model3_and_Model4 <- 5 # Max number of orthogonal variables in model pre variable selection of Model3 and Model4 with pcorr cutoff giving highest Q2
max_no_of_ortho_post_vs_in_Model2 <- 5 # Max number of orthogonal variables in model post variable selection of Model2 with pcorrcutoff according to pvalue
max_no_of_ortho_post_vs_in_Model3_and_Model4 <- 5 # Max number of orthogonal variables in model post variable selection of Model3 and Model4 with pcorr cutoff giving highest Q2

## Best performing model
prefered_pR2_and_pQ2_permutated_post_vs <- 0.05 # Prefered pR2 and pQ2 determined by permutation post variable selection during model optimization. The lower prefered pQ2 the larger weight will be given to pQ2 instead of Q2 as well as instead of diff between R2 and Q2. Minimum value is 1/no_permutations_post_vs.
pcorr_diff <- 0.02 

## Loading plots
variable_names_length <- "all" #Number of characters of variablenames shown in loadingplot or "all"
variable_names_position <- "all" # enter "beginning", "end" or "all" give selection of variablenames shown in loadingplot

## Running models
order_of_groups_numeric <- "correct" # numeric vector* containing correct order of groups or enter "correct". Deseased first and controls last in order to get direction of models with high score in diseased.
cluster <- "no" # "no" if analysis is run locally and "yes" if script is run using the bash scripts to run the R-scripts passing modelnumbers to run as arguments ("args").
name_intermediate_dir <- "default" # "name_of_directory" enter name of intermediate directory. On UPPMAX it is "SNIC_TMP" else enter "default"
each_model_or_summary <- "both" # enter "summary" if only summary should be run and "each" if only models of each comparison should be run. Otherwise enter "both". "Both" is not possible if cluster is "yes".






# Set Advanced parameter settings
## Foldernames
foldername_Rmarkdownfiles <- "Rscript" #"foldername". This is the folder for the scripts .R and .rmd files.
foldername_of_input_matrix_and_sampleID <- "data" #"foldername" Folder is automatically created if it does not excist. This is the folder for datamatrix and sampleID files
foldername_output_reports <- "outputR" # "foldername"  Folder is automatically created if it does not excist. In this folder one html file for each group comparison and a summary html file of all comparisons will be created.

directory_input_matrix_sampleID <-  paste(directory_of_analysis,"/",foldername_of_input_matrix_and_sampleID,"/",sep="")#may be altered to other than directory_of_analysis
directory_output_reports <- paste(directory_of_analysis,"/",foldername_output_reports,"/",sep="")  #may be altered to other than directory_of_analysis

## File names
filename_Rmarkdownfile_each_model <- "roplspvs_Models_of_each_comparison.Rmd" # "filename.Rmd"
filename_Rmarkdownfile_summary <- "roplspvs_Summary_of_models.Rmd" # "filename.Rmd"
filename_function_file <- "roplspvs_Functions.R" # "filename.R" All functions used are saved in this file.

##Variable selection
pcorr_cutoff_Model1_joint_models <- "according to p-value" #P(corr) cutoff in model1 for joint models. Either enter a value with p(corr)cutoff for all joint comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
pcorr_cutoff_Model1_stratified_models <- "according to p-value" #P(corr) cutoff in model1 for stratified models. Either enter a value with p(corr)cutoff for all stratified comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
variable_selection_using_VIP <- "no" # Enter "yes" if VIP is used during variable selection and "no" if only p(corr) is used.

#Amount of orthogonals
no_of_ortho_pre_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for joint models. Either enter a value for nonstratified comparisons or a vector containing the number of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_pre_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for joint models. Either enter a value for joint comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
max_no_of_ortho_pre_vs <- 5 # Max number of orthogonal variables in model pre variable selection of Model2 with pcorrcutoff according to pvalue
max_no_of_ortho_post_vs <- 5 # Max number of orthogonal variables in model post variable selection of Model2 with pcorrcutoff according to pvalue

## Best performing model
prefered_pR2_and_pQ2_permutated_post_vs <- 0.05 # Prefered pR2 and pQ2 determined by permutation post variable selection during model optimization. The lower prefered pQ2 the larger weight will be given to pQ2 instead of Q2 as well as instead of diff between R2 and Q2. Minimum value is 1/no_permutations_post_vs.
pcorr_diff <- 0.01

## Loading plots
variable_names_length <- "all" #Number of characters of variablenames shown in loadingplot or "all"
variable_names_position <- "all" # enter "beginning", "end" or "all" give selection of variablenames shown in loadingplot

## Running models
cluster <- "no" # "no" if analysis is run locally and "yes" if script is run using the bash scripts to run the R-scripts passing modelnumbers to run as arguments ("args").
name_intermediate_dir <- "default" # "name_of_directory" enter name of intermediate directory. On UPPMAX it is "SNIC_TMP" else enter "default"
each_model_or_summary <- "both" # enter "summary" if only summary should be run and "each" if only models of each comparison should be run. Otherwise enter "both". "Both" is not possible if cluster is "yes".
model_strategies_to_run <- c(1,2,3,4,5)

#Preparations
#1. Create a project folder for the analysis manually and change workingdirectory to this file.
#2. Download Rscript folder, Run.R, dependencies.R, Configure files

#Preparations to run example data
#3. Download data_to_R_analysis containing exampledata

# Set path where Ropls-with-permutated-variable-selection is stored
directory_of_Ropls_with_permutated_variable_selection <- getwd()

## Project settings
projectname <- "testproj MTBLS136" # "projectname" will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 201218 # "yymmdd" numeric date of analysis will appear in filenames.

## Data matrix
filename_matrix <- "s_MTBLS136_datamatrix.txt" #"filename.txt"
decimal_separator <- "dot" # "dot" or "comma"

## SampleID
filename_sampleID <- "sampleID.txt" #filename.txt"
colname_groupID <- "Factor_Value_CurrentPMH" # "column name" of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"Factor_Value_AgeAtBloodDraw" # "column name" column name of secondary id in sampleID file to stratify or write "joint" for no stratification

## Permutations
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_pre_vs <- 20 # Numeric. Number of permutations before variable selection in selected models.

## Variable selection
p_pearson_of_pcorr_cutoff <- 0.05 # P-value for p(corr) cutoff during variable selection

## Running models
setseedfirstmodel <- 200 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
order_of_groups_names <- c("E+P","E-only","Nonuser") # character vector containing correct order of group names or enter "correct". Deseased first and controls last in order to get direction of models with high score in diseased.
models_to_run <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,23) # numeric vector indicating models to run if all models are to be run enter "all"

## Paths of directories
directory_of_analysis <-  paste(location_of_analysis_folder,"/",analysis_folder_name,sep="")
directory_input_matrix_sampleID <-  paste(directory_of_analysis,"/",foldername_of_input_matrix_and_sampleID,"/",sep="")
directory_output_reports <- paste(directory_of_analysis,"/",foldername_output_reports,"/",sep="")


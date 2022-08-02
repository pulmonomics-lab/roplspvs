#Preparations
#1. Create a project folder for the analysis manually and change workingdirectory to this file.
#2. Download Rscript folder, Run.R, dependencies.R, Configure files

#Preparations to run example data
#3. Download data_to_R_analysis containing exampledata

## Project settings
projectname <- "testproj MTBLS136" # "projectname" will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 220614 # "yymmdd" numeric date of analysis will appear in filenames.

## Data matrix
filename_matrix <- "s_MTBLS136_datamatrix.txt" #"filename.txt"
decimal_separator <- "dot" # "dot" or "comma"
replace_0 <- "NA" # F or what to replace 0 values with. Replacing is performed before filtering.
filter_percent_in_each_group <- 25 # Numeric. Missing value tolerance in each group which are compared.
replace_NA <- "lld" # F or "lld" or value to replace with. If not replaced NAs will be imputed. Replaceing is performed after filtering. Llq is calculated by lld in dataset/3.
log_transform <- T # T or F

## SampleID
filename_sampleID <- "sampleID.txt" #filename.txt"
colname_groupID <- "Factor_Value_CurrentPMH" # "column name" of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"Factor_Value_AgeAtBloodDraw" # "column name" column name of secondary id in sampleID file to stratify or write "joint" for no stratification

## Permutations
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_over_vs <- 20 # Numeric. Number of permutations including variable selection in selected models.

## Variable selection
p_pearson_of_pcorr_cutoff <- 0.05 # P-value for p(corr) cutoff during variable selection

## Running models
setseedfirstmodel <- 200 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
order_of_groups <- c("E+P","E-only","Nonuser") # Character vector or numeric vector containing correct order of groups to compare or enter "correct" if order of levels in colname_groupID is already correct. Deseased first and controls last. This will define direction of scores as high in diseased.
models_to_run <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,23) # numeric vector indicating models to run if all models are to be run enter "all"


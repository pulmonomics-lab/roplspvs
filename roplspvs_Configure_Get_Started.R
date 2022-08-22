#Preparations
##1. Download roplspvs and see detailed description at https://github.com/MarikaStrom/roplspvs.git
##2. You are adviced to create a separate project folder for the analysis manually and change workingdirectory to this file. After editing parameters save this file in your project folder. Also roplspvs_Configure_Advanced.R may be edited and saved in project folder.

#Preparations of data
##1. Prepare the datamatrix with subjectID in the first column as row names and variable names in the first row as column names. Be careful that all names are unique.
##2. Prepare a sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with metadata (secondary ID for example gender) if stratification is desired.
##3. Save datamatrix and sampleID files as tabdelimited .txt files in a manually created folder in the project folder with default name "data" or enter path in advanced settings.

#Run analysis
##Run analysis by using the following code:
##source("<path to roplspvs>/dependencies.R")
##source("<path to roplspvs>/roplspvs_Run.R")

## Project settings
projectname <- "projectname" # "projectname" will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 220614 # yymmdd numeric date of analysis will appear in filenames.

## Data matrix
filename_matrix <- "filename.txt" #"filename.txt"
decimal_separator <- "dot" # "dot" or "comma" as decimal separator in filename_matrix
replace_0 <- F # F or "lld" or a value to replace 0 values with. Replacing is performed before filtering. Lld (lower limit of detection) is calculated by llq/3 with llq (lower limit of quantification) is the lowest value in the dataset.
filter_percent_in_each_group <- 25 # Numeric. Missing value tolerance in each group which are compared.
replace_NA <- F # F or "lld" or a value to replace NA with. If not replaced NAs will be imputed using NIPAL by ropls package. Replaceing is performed after filtering. Lld (lower limit of detection) is calculated by llq/3 with llq (lower limit of quantification) is the lowest value in the dataset.
log_transform <- F # T or F if datamatrix is transformed by natural logarithm

## SampleID
filename_sampleID <- "filename.txt" #"filename.txt"* containing metadata with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. All subjectIDs has to be unique and ordered in the same way as in datamatrix.
colname_groupID <- "column name" # "column name"* of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"column name" # "column name"* column name of secondary id in sampleID file to stratify or write "joint" for no stratification

## Permutations
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_over_vs <- 20 # Numeric. Number of permutations including variable selection in selected models.
no_permutations_sans_vs <- 20 # Numeric. Number of permutations sans variable selection

## Variable selection
p_pearson_of_pcorr_cutoff <- 0.05 # Numeric. P-value for p(corr) cutoff during variable selection in model 2 and minimum cutoff in model 3 and 4. Optionally in model 1.

## Running model settings
setseedfirstmodel <- 1 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
order_of_groups <- "correct" # # Character vector of group names or numeric vector of correct order of groups to compare or enter "correct" if order of levels in colname_groupID is already correct. Deseased first and controls last. This will define direction of scores as high in diseased.
comparisons_to_run <- "all" # numeric vector indicating which comparisons to run as numbered in the file "model_table_to_analyse". If all comparisons are to be run enter "all". model_table_to_analyse is created when starting the script. Edited model_table_to_analyse manually or indicate which comparisons to run.



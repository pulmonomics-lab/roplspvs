#Preparations
#1. Create a project folder for the analysis manually and change workingdirectory to this file.
#2. Download Rscript folder, Run.R, dependencies.R, Configure files

#Preparations of data
#1. Prepare the datamatrix with subjectID in the first column as rownames and variable names in the first row as columnnames. Be careful that all names are unique.
#2. Prepare a sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired.
#3. Save datamatrix and sampleID files as tabdelimited .txt files in a manually created folder in the project folder with default name "data_to_R_analysis" or enter path in advanced settings.

## Project settings
projectname <- "projectname" # "projectname" will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 201218 # yymmdd numeric date of analysis will appear in filenames.

## Data matrix
filename_matrix <- "filename.txt" #"filename.txt"
decimal_separator <- "dot" # "dot" or "comma" as decimal separator in filename_matrix

## SampleID
filename_sampleID <- "filename.txt" #"filename.txt"* containing metadata with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. All subjectIDs has to be unique and ordered in the same way as in datamatrix.
colname_groupID <- "column name" # "column name"* of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"column name" # "column name"* column name of secondary id in sampleID file to stratify or write "joint" for no stratification

## Permutations
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_over_vs <- 20 # Numeric. Number of permutations including variable selection in selected models.

## Variable selection
p_pearson_of_pcorr_cutoff <- 0.05 # Numeric. P-value for p(corr) cutoff during variable selection in model 2 and minimum cutoff in model 3 and 4. Optionally in model 1.

## Running model settings
setseedfirstmodel <- 1 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
order_of_groups <- "correct" # # Character vector or numeric vector containing correct order of groups to compare or enter "correct" if order of levels in colname_groupID is already correct. Deseased first and controls last. This will define direction of scores as high in diseased.
models_to_run <- "all" # numeric vector indicating models to run if all models are to be run enter "all"



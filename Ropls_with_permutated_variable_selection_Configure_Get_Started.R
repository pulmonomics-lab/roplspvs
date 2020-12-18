#Preparations 
#1. Prepare the datamatrix with subjectID in the first column as rownames and variable names in the first row as columnnames. Be careful that all names are unique.
#2. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. 
#3. Both datamatrix and sampleID files should be saved as tabdelimited .txt files. 
#4. Create a project folder for the analysis manually and change workingdirectory to this file

#Enter data to send to Rmarkdown:

## Project settings
directory_of_analysis <- getwd() # "path"* including name of project directory where subdirectories will be created.
projectname <- "testproj MTBLS136" # "projectname"* will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 201218 # "yymmdd"* numeric date of analysis will appear in filenames.

## Data matrix
filename_matrix <- "s_MTBLS136_datamatrix.txt" #"filename.txt"*
decimal_separator <- "dot" # "dot" or "comma"*

## SampleID
filename_sampleID <- "sampleID.txt" #filename.txt"*
colname_groupID <- "Factor_Value_CurrentPMH" # "column name"* of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"Factor_Value_AgeAtBloodDraw" # "column name"* column name of secondary id in sampleID file to stratify or write "joint" for no stratification

## Permutations
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_pre_vs <- 20 # Numeric. Number of permutations before variable selection in selected models.

setseedfirstmodel <- 200 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
p_pearson_of_pcorr_cutoff <- 0.05 # P-value for p(corr) cutoff during variable selection
order_of_groups_names <- c("E+P","E-only","Nonuser") # character vector* containing correct order of group names or enter "correct". Deseased first and controls last in order to get direction of models with high score in diseased.
models_to_run <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,23) # numeric vector indicating models to run if all models are to be run enter "all"

## Paths of directories
foldername_model_table_to_analyse <- foldername_output_reports #model_table_to_analyse will be saved in this folder.
directory_Rmarkdownfiles <- paste(directory_of_analysis,"/",foldername_Rmarkdownfiles,"/",sep="") # "path of directory/filename"
directory_input_matrix_sampleID <-  paste(directory_of_analysis,"/",foldername_of_input_matrix_and_sampleID,"/",sep="")
directory_output_reports <- paste(directory_of_analysis,"/",foldername_output_reports,"/",sep="")
directory_model_table_to_analyse <- paste(directory_of_analysis,"/",foldername_model_table_to_analyse,"/",sep="")
directory_function_file <- paste(directory_of_analysis,"/",foldername_function_file,"/", sep="")

#Upload files 
## 1.Datamatrix and sampleID files to foldername_of_input_matrix_and_sampleID
## 2. Scripts .R files to foldername_function_file and .rmd files to foldername_Rmarkdownfiles. These may be the same folder.

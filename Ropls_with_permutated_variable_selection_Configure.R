 #Preparations 
#1. Prepare the datamatrix with subjectID in the first column as rownames and variable names in the first row as columnnames. Be careful that all names are unique.
#2. A sampleID file with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. 
#3. Both datamatrix and sampleID files should be saved as tabdelimited .txt files. 
#4. Create a project folder for the analysis manually.


#Enter data to send to Rmarkdown:

directory_of_analysis <- "C:/Users/Marika/OneDrive - KI.SE/Documents/ropls_Rscript/testproj MTBLS136 201216/" # "path"* including name of project directory where subdirectories will be created.
projectname <- "testproj MTBLS136" # "projectname"* will appear in filenames and header of reports with underscores removed.
date_of_analysis <- 201216 # "yymmdd"* numeric date of analysis will appear in filenames.
foldername_Rmarkdownfiles <- "Rscript" #"foldername" Folder is automatically or manually created. This is the folder for the scripts .R and .rmd files.
filename_Rmarkdownfile_each_model <- "Ropls_with_permutated_variable_selection_vs10_Models_of_each_comparison.Rmd" # "filename.Rmd"
filename_Rmarkdownfile_summary <- "Ropls_with_permutated_variable_selection_vs10_Summary_of_models.Rmd" # "filename.Rmd"
foldername_of_input_matrix_and_sampleID <- "data_to_R_analysis" #"foldername" Folder is automatically created. This is the folder for datamatrix and sampleID files
filename_matrix <- "s_MTBLS136_datamatrix.txt" #"filename.txt"*
decimal_separator <- "dot" # "dot" or "comma"*
filename_sampleID <- "sampleID.txt" #filename.txt"*
colname_groupID <- "Factor_Value_CurrentPMH" # "column name"* of groups to compare in sampleID file
groupsnumeric <- "no" #write "yes" if you want tables in summary file to be sorted by numeric group belonging
colname_secID <-"Factor_Value_AgeAtBloodDraw" # "column name"* column name of secondary id in sampleID file to stratify or write "joint" for no stratification
foldername_output_reports <- "outputR" # "foldername"  Folder is automatically created below. In this folder one html file for each group comparison and a summary html file of all comparisons will be created.
foldername_function_file <- foldername_Rmarkdownfiles # foldername_Rmarkdownfiles or "foldername". This is the folder for the .R and .Rmd files.
filename_function_file <- "Ropls_with_permutated_variable_selection_vs10_Functions.R" # "filename.R" All functions used are saved in this file.
no_permutations_post_vs <- 20 # numeric. Number of permutations after variable selection during model selection.
no_permutations_post_vs_selected_models <- 20 # Numeric. Number of permutations after variable selection in selected models.
no_permutations_pre_vs <- 20 # Numeric. Number of permutations before variable selection in selected models.
filter_percent_in_each_group <- 25 # Numeric. Missing value tolerance in each group which are compared.
setseedfirstmodel <- 200 #Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
p_pearson_of_pcorr_cutoff <- 0.05 # P-value for p(corr) cutoff during variable selection for model 2 and optionally for model1
pcorr_cutoff_Model1_joint_models <- "according to p-value" #P(corr) cutoff in model1 for joint models. Either enter a value with p(corr)cutoff for all joint comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
pcorr_cutoff_Model1_stratified_models <- "according to p-value" #P(corr) cutoff in model1 for stratified models. Either enter a value with p(corr)cutoff for all stratified comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
no_of_ortho_pre_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for joint models. Either enter a value for nonstratified comparisons or a vector containing the number of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_pre_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model pre variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_joint_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for joint models. Either enter a value for joint comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
no_of_ortho_post_vs_Model1_stratified_models <- 0 #Number of orthogonal variables in model post variable selection of Model1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
max_no_of_ortho_pre_vs_in_Model2 <- 5 # Max number of orthogonal variables in model pre variable selection of Model2 with pcorrcutoff according to pvalue
max_no_of_ortho_pre_vs_in_Model3_and_Model4 <- 5 # Max number of orthogonal variables in model pre variable selection of Model3 and Model4 with pcorr cutoff giving highest Q2
max_no_of_ortho_post_vs_in_Model2 <- 5 # Max number of orthogonal variables in model post variable selection of Model2 with pcorrcutoff according to pvalue
max_no_of_ortho_post_vs_in_Model3_and_Model4 <- 5 # Max number of orthogonal variables in model post variable selection of Model3 and Model4 with pcorr cutoff giving highest Q2
filename_model_table_to_analyse <- "model_table_to_analyse.txt"
foldername_model_table_to_analyse <- "data_to_R_analysis" #model_table_to_analyse will be saved in this folder.
prefered_pR2_and_pQ2_permutated_post_vs <- 0.05 # Prefered pR2 and pQ2 determined by permutation post variable selection during model optimization. The lower prefered pQ2 the larger weight will be given to pQ2 instead of Q2 as well as instead of diff between R2 and Q2
pcorr_diff <- 0.02 
variable_names_length <- "all" #Number of characters of variablenames shown in loadingplot or "all"
variable_names_position <- "all" # enter "beginning", "end" or "all" give selection of variablenames shown in loadingplot
cluster <- "no" # "no" if analysis is run locally and "yes" if script is run using the bash scripts to run the R-scripts passing modelnumbers to run as arguments ("args").
each_model_or_summary <- "both" # enter "summary" if only summary should be run and "each" if only models of each comparison should be run. Otherwise enter "both". "Both" is not possible if cluster is "yes".
name_intermediate_dir <- "default" # "name_of_directory" enter name of intermediate directory. On UPPMAX it is "SNIC_TMP" else enter "default"

#Read paths of directories

directory_Rmarkdownfiles <- paste(directory_of_analysis,"/",foldername_Rmarkdownfiles,"/",sep="") # "path of directory/filename"
directory_input_matrix_sampleID <-  paste(directory_of_analysis,"/",foldername_of_input_matrix_and_sampleID,"/",sep="")
directory_output_reports <- paste(directory_of_analysis,"/",foldername_output_reports,"/",sep="")
directory_model_table_to_analyse <- paste(directory_of_analysis,"/",foldername_model_table_to_analyse,"/",sep="")
directory_and_filename_function_file <- paste(directory_of_analysis,"/",foldername_function_file,"/",filename_function_file,sep="")

#Load functions
source(directory_and_filename_function_file)

#Create directories.
makeproject(directory_of_analysis)

#Upload files manually
## 1.Datamatrix and sampleID files to foldername_of_input_matrix_and_sampleID
## 2. Scripts .R files to foldername_function_file and .rmd files to foldername_Rmarkdownfiles. These may be the same folder.

#Reorder groups to have deseased first if nessasary
setwd(directory_input_matrix_sampleID)
file_sampleID <- read.table(filename_sampleID,header=T, dec = ".", row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf",""), sep="\t")
levels_of_groups <- levels(as.factor(file_sampleID[,paste(colname_groupID)]))
levels_of_groups # Check if order of the groups are correct with deseased first and controls last
order_of_groups <- c(2,1,3) # numeric vector* containing correct order of groups or enter "correct"
reordered_levels_of_groups <- reorder_levels_of_groups() # Note that it has to be run even if order is "correct"
reordered_levels_of_groups #Check that reordered groups has deseased first and controls last.

#Create or load model table to analyse containing the comparisons to be modeled.

model_table_to_analyse <- create_or_load_Model_table_to_analyze() # Checks if model_table_to_analyse exist and loads or creats it in the folder foldername_model_table_to_analyse
model_table_to_analyse #check which models should be included in analysis
models_to_run <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,23) # numeric vector indicating models to run if all models are to be run enter "all"
model_table_to_analyse <- select_models_to_run(model_table_to_analyse, models_to_run)

#Send data to Rmarkdown file to create html files and lists of loadings

setwd(directory_output_reports)
if (each_model_or_summary=="each"| each_model_or_summary=="both") {
  #Send data to Rmarkdown file to create html for each comparison
  setwd(directory_output_reports)
  if (cluster=="yes") {args = commandArgs(trailingOnly=TRUE);model_number_to_run<-args;first_model_number_to_run<-model_number_to_run; last_model_number_to_run<-model_number_to_run} else {first_model_number_to_run<-1; last_model_number_to_run<-nrow(model_table_to_analyse)}
  for (i in first_model_number_to_run:last_model_number_to_run) {
    group1 <- paste(model_table_to_analyse$group1[i])
    group2 <- paste(model_table_to_analyse$group2[i])
    secID <- paste(model_table_to_analyse$secID[i])
    setseedno <- model_table_to_analyse$setseedno[i]
    pcorr_Model1 <- model_table_to_analyse$pcorr_Model1[i]
    ortho_pre_vs_Model1 <- model_table_to_analyse$ortho_pre_vs[i]
    ortho_post_vs_Model1 <-   model_table_to_analyse$ortho_post_vs[i]
    if (cluster=="yes" & (each_model_or_summary=="each" | each_model_or_summary=="both")) {
      rmarkdown::render(paste(directory_Rmarkdownfiles,"/",filename_Rmarkdownfile_each_model, sep =""),
                        output_file = paste(paste(directory_output_reports,"/",projectname, sep=""), date_of_analysis, group1, "vs", group2, secID, sep="_"),
                        intermediates_dir = Sys.getenv(paste(name_intermediate_dir)),
                        params=list(
                          colname_groupID = colname_groupID,
                          colname_secID = colname_secID,
                          group1 = group1,
                          group2 = group2,
                          secID = secID,
                          pcorr_Model1 = pcorr_Model1,
                          ortho_pre_vs_Model1 = ortho_pre_vs_Model1,
                          ortho_post_vs_Model1 = ortho_post_vs_Model1,
                          setseedno=setseedno,
                          p_pearson_of_pcorr_cutoff=p_pearson_of_pcorr_cutoff,
                          directory_input_matrix_sampleID=directory_input_matrix_sampleID,
                          filename_matrix=filename_matrix,
                          decimal_separator=decimal_separator,
                          variable_names_length=variable_names_length,
                          variable_names_position=variable_names_position,
                          filename_sampleID=filename_sampleID,
                          directory_output_reports=directory_output_reports,
                          projectname=projectname,
                          date_of_analysis=date_of_analysis,
                          no_permutations_post_vs=no_permutations_post_vs,
                          no_permutations_post_vs_selected_models=no_permutations_post_vs_selected_models,
                          no_permutations_pre_vs=no_permutations_pre_vs,
                          filter_percent_in_each_group=filter_percent_in_each_group,
                          directory_and_filename_function_file=directory_and_filename_function_file,
                          max_no_of_ortho_pre_vs_in_Model2=max_no_of_ortho_pre_vs_in_Model2,
                          max_no_of_ortho_pre_vs_in_Model3_and_Model4=max_no_of_ortho_pre_vs_in_Model3_and_Model4,
                          max_no_of_ortho_post_vs_in_Model2=max_no_of_ortho_post_vs_in_Model2,
                          max_no_of_ortho_post_vs_in_Model3_and_Model4=max_no_of_ortho_post_vs_in_Model3_and_Model4,
                          reordered_levels_of_groups=reordered_levels_of_groups,
                          prefered_pR2_and_pQ2_permutated_post_vs=prefered_pR2_and_pQ2_permutated_post_vs,
                          pcorr_diff=pcorr_diff
                        ))
    }
    
    if (cluster=="no" & (each_model_or_summary=="each" | each_model_or_summary=="both")) {
      rmarkdown::render(paste(directory_Rmarkdownfiles,"/",filename_Rmarkdownfile_each_model, sep =""),
                        output_file = paste(paste(directory_output_reports,"/",projectname, sep=""), date_of_analysis, group1, "vs", group2, secID, sep="_"),
                        params=list(
                          colname_groupID = colname_groupID,
                          colname_secID = colname_secID,
                          group1 = group1,
                          group2 = group2,
                          secID = secID,
                          pcorr_Model1 = pcorr_Model1,
                          ortho_pre_vs_Model1 = ortho_pre_vs_Model1,
                          ortho_post_vs_Model1 = ortho_post_vs_Model1,
                          setseedno=setseedno,
                          p_pearson_of_pcorr_cutoff=p_pearson_of_pcorr_cutoff,
                          directory_input_matrix_sampleID=directory_input_matrix_sampleID,
                          filename_matrix=filename_matrix,
                          decimal_separator=decimal_separator,
                          variable_names_length=variable_names_length,
                          variable_names_position=variable_names_position,
                          filename_sampleID=filename_sampleID,
                          directory_output_reports=directory_output_reports,
                          projectname=projectname,
                          date_of_analysis=date_of_analysis,
                          no_permutations_post_vs=no_permutations_post_vs,
                          no_permutations_post_vs_selected_models=no_permutations_post_vs_selected_models,
                          no_permutations_pre_vs=no_permutations_pre_vs,
                          filter_percent_in_each_group=filter_percent_in_each_group,
                          directory_and_filename_function_file=directory_and_filename_function_file,
                          max_no_of_ortho_pre_vs_in_Model2=max_no_of_ortho_pre_vs_in_Model2,
                          max_no_of_ortho_pre_vs_in_Model3_and_Model4=max_no_of_ortho_pre_vs_in_Model3_and_Model4,
                          max_no_of_ortho_post_vs_in_Model2=max_no_of_ortho_post_vs_in_Model2,
                          max_no_of_ortho_post_vs_in_Model3_and_Model4=max_no_of_ortho_post_vs_in_Model3_and_Model4,
                          reordered_levels_of_groups=reordered_levels_of_groups,
                          prefered_pR2_and_pQ2_permutated_post_vs=prefered_pR2_and_pQ2_permutated_post_vs,
                          pcorr_diff=pcorr_diff
                          
                        ))
      
      
    }
    
  }
}

if (each_model_or_summary=="summary"|(cluster=="no" & each_model_or_summary=="both")) {
  #Send data to Rmarkdown file to create html for summary of models
  rmarkdown::render(paste(directory_Rmarkdownfiles,"/",filename_Rmarkdownfile_summary, sep =""),
                    output_dir = paste(directory_output_reports,"/",sep=""),
                    output_file = paste(directory_output_reports,"/", "Summary_",projectname, "_",date_of_analysis, sep=""),
                    params=list(
                      directory_output_reports=directory_output_reports,
                      directory_and_filename_function_file=directory_and_filename_function_file,
                      projectname=projectname,
                      date_of_analysis=date_of_analysis,
                      directory_input_matrix_sampleID=directory_input_matrix_sampleID,
                      filename_matrix=filename_matrix,
                      filename_sampleID=filename_sampleID,
                      no_permutations_pre_vs=no_permutations_pre_vs,
                      no_permutations_post_vs_selected_models=no_permutations_post_vs_selected_models,
                      filter_percent_in_each_group=filter_percent_in_each_group,
                      groupsnumeric=groupsnumeric))
}



#Create SUSplots

#source(directory_and_filename_function_file)

#SUSplot1
 
#projectnamemodel2 <- "projectname"
#date_of_analysismodel2  <- 200621
#modelnamemodel2 <-  "resultModel2"
#group1model2<- "COPDex"
#group2model2 <- "COPDsm"
#secIDmodel2 <- "female"


#projectnamemodel1 <- "projectname"
#date_of_analysismodel1 <- 200621
#modelnamemodel1 <-  "resultModel2"
#group1model1 <- "COPDex"
#group2model1 <- "smoker"
#secIDmodel1 <- "female"

#plotSUSplotsameanalysis(directory_output_reports, projectnamemodel1, date_of_analysismodel1, modelnamemodel1, group1model1, group2model1, secIDmodel1, projectnamemodel2, date_of_analysismodel2, modelnamemodel2, group1model2, group2model2, secIDmodel2)


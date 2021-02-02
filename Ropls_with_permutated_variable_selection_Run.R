#Load configurations

source(paste(analysis_folder_name,"/Ropls_with_permutated_variable_selection_Configure_Advanced.R",sep=""))
if (file.exists(list.files(paste(getwd(),"/",analysis_folder_name,sep=""), pattern="\\Get_Started.R"))) {
  source(paste(analysis_folder_name,"/",list.files(paste(getwd(),"/",analysis_folder_name,sep=""), pattern="\\Get_Started.R"),sep=""))} else {
    source(paste(analysis_folder_name,"/",list.files(paste(getwd(),"/",analysis_folder_name,sep=""), pattern="\\Get_Started_example_data.R"),sep=""))}

#Read directories of Function files and Rmarkdown files

directory_Rmarkdownfiles <- paste(directory_of_Ropls_with_permutated_variable_selection,"/","Ropls-with-permutated-variable-selection","/",foldername_Rmarkdownfiles,"/",sep="") # "path of directory/filename"
directory_function_file <- paste(directory_of_Ropls_with_permutated_variable_selection,"/","Ropls-with-permutated-variable-selection","/",foldername_function_file,"/", sep="")
foldername_model_table_to_analyse <- foldername_output_reports #model_table_to_analyse will be saved in this folder.
directory_model_table_to_analyse <- paste(directory_of_analysis,"/",foldername_model_table_to_analyse,"/",sep="")

#Load functions

directory_and_filename_function_file <- paste(directory_function_file,"/",filename_function_file,sep="")
source(directory_and_filename_function_file)

#Create directories.
makeproject(directory_of_analysis)


#Reorder groups to have deseased first if nessasary
setwd(directory_input_matrix_sampleID)
file_sampleID <- read.table(filename_sampleID,header=T, dec = ".", row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf",""), sep="\t")
levels_of_groups <- levels(as.factor(file_sampleID[,paste(colname_groupID)]))
levels_of_groups # Check if order of the groups are correct with deseased first and controls last
reordered_levels_of_groups <- reorder_levels_of_groups() # Note that it has to be run even if order is "correct"
reordered_levels_of_groups #Check that reordered groups has deseased first and controls last.
order_of_groups_numeric <- match(levels_of_groups,reordered_levels_of_groups)

#Create or load model table to analyse containing the comparisons to be modeled.
filename_model_table_to_analyse <- "model_table_to_analyse.txt"
model_table_to_analyse <- create_or_load_Model_table_to_analyze() # Checks if model_table_to_analyse exist and loads or creats it in the folder foldername_model_table_to_analyse
model_table_to_analyse #check which models should be included in analysis
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
                          pcorr_diff=pcorr_diff,
                          variable_selection_using_VIP=variable_selection_using_VIP
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
                          pcorr_diff=pcorr_diff,
                          variable_selection_using_VIP=variable_selection_using_VIP
                          
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


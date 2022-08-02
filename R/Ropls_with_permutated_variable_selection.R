#' Opls with permutationed variable selection version 0.13.0
#'
#' Pipeline for OPLS models, variable selection and permutation
#'
#' Performs variable selection using p(corr) and optionally VIP. In addition to permutation post variable
#' selection included in the ropls package also permutations over variable selection with proceeding
#' variable selection of every permutation resulting in p-values for R2 and Q2 including the variable
#' selection procedure. It produces tables of all group comparisons including optional stratification
#' by secondary ID.
#'
#' @param directory_of_Ropls_with_permutated_variable_selection path to repository cloned from https://github.com/MarikaStrom/Ropls-with-permutated-variable-selection
#' @param directory_of_analysis path to folder where the outputR folder with results are created
#' @param projectname "projectname" will appear in filenames and header of reports with underscores removed.
#' @param date_of_analysis "yymmdd" Date of analysis will appear in filenames.
#' @param filename_matrix "filename.txt" file containing your data to be modeles. subjectID in the first column as rownames and variable names in the first row as columnnames. All names must be unique.
#' @param decimal_separator "dot" or "comma" as decimal separator in filename_matrix
#' @param filename_sampleID "filename.txt" file containing metadata with sampleID in the first column and containing one column with the groups to be compared and one column with secondaryID (for example gender) if stratification is desired. All subjectIDs has to be unique and ordered in the same way as in datamatrix.
#' @param colname_groupID "column name"* of groups to compare in sampleID file
#' @param groupsnumeric  "yes" if you want tables in summary file to be sorted by numeric group belonging. "no" if not.
#' @param colname_secID "column name"* column name of secondary id in sampleID file to stratify or write "joint" for no stratification
#' @param no_permutations_post_vs Numeric. Number of permutations post variable selection during model selection.
#' @param no_permutations_post_vs_selected_models Numeric. Number of permutations post variable selection in selected models.
#' @param no_permutations_over_vs Numeric. Number of permutations including variable selection in selected models.
#' @param p_pearson_of_pcorr_cutoff Numeric. P-value for p(corr) cutoff during variable selection in model 2 and minimum cutoff in model 3 and 4. Optionally in model 1.
#' @param setseedfirstmodel Numeric. Setseed of the first model. Second model will have setseedfirstmodel+1 etc.
#' @param order_of_groups Character vector or numeric vector containing correct order of groups or enter "correct" if order of levels in colname_groupID is already correct. Deseased first and controls last. This will define dirction of scores as high in diseased.
#' @param models_to_run Numeric vector indicating models to run if all models are to be run enter "all"
#' @param foldername_Rmarkdownfiles "foldername" Folder where .R and .rmd files required by the script are stored. It is located in the Ropls-with-permutated-variable-selection folder.
#' @param foldername_of_input_matrix_and_sampleID "foldername" Folder where datamatrix and sampleID files are stored.
#' @param foldername_output_reports "foldername"  Folder where one .html file and one .Rdata for each group comparison and a summary .html and a .Rdata file of all comparisons will be created. Folder is automatically created if it does not excist.
#' @param foldername_function_file "foldername" Folder where filename_function_file is stored.
#' @param directory_input_matrix_sampleID "path" Path to where "filename_matrix" and "filename_sampleID" files are stored.
#' @param directory_output_reports "path" Path to where "output_reports" are created in foldername_output_reports.
#' @param filename_Rmarkdownfile_each_model "filename.Rmd" rendered by function once per comparison to produce html report and .Rdata file for each comparison
#' @param filename_Rmarkdownfile_summary "filename.Rmd" rendered by function if "each_model_or_summary" is set to "both" or "summary" after all comparisons have been performed to produce a summary tables of all models.
#' @param filename_function_file "filename.R" All functions required are saved in this file.
#' @param filter_percent_in_each_group Numeric. Missing value tolerance in each group which are compared.
#' @param pcorr_cutoff_Model1_joint_models P(corr) cutoff in model1 for joint models. Either enter a value with p(corr)cutoff for all joint comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff.
#' @param pcorr_cutoff_Model1_stratified_models P(corr) cutoff in model1 for stratified models. Either enter a value with p(corr)cutoff for all stratified comparisons or vector containing p(corr) cutoff for each comparison in model_table_to_analyse or enter "according to p-value" to generate p(corr) corresponding to selected pvalue entered in p_pearson_of_pcorr_cutoff
#' @param variable_selection_using_VIP "yes" if VIP is used during variable selection and "no" if only p(corr) is used.
#' @param no_of_ortho_pre_vs_Model1_joint_models Numeric. Number of orthogonal variables in model pre variable selection of Model strategy 1 for joint models. Either enter a value for nonstratified comparisons or a vector containing the number of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
#' @param no_of_ortho_pre_vs_Model1_stratified_models Numeric. Number of orthogonal variables in model pre variable selection of Model strategy 1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model pre variable selection for each model in model_table_to_analyse.
#' @param no_of_ortho_post_vs_Model1_joint_models Numeric. Number of orthogonal variables in model post variable selection of Model strategy 1 for joint models. Either enter a value for joint comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
#' @param no_of_ortho_post_vs_Model1_stratified_models Numeric. Number of orthogonal variables in model post variable selection of Model strategy 1 for stratified models. Either enter a value for stratified comparisons or vector containing the amount of orthogonal variables investigated in model post variable selection for each model in model_table_to_analyse.
#' @param max_no_of_ortho_pre_vs Numeric. Max number of orthogonal variables in model pre variable selection
#' @param max_no_of_ortho_post_vs Numeric. Max number of orthogonal variables in model post variable selection
#' @param prefered_pR2_and_pQ2_permutated_post_vs Numeric. Prefered pR2 and pQ2 determined by permutation post variable selection during model optimization. The lower value the larger weight will be given to pQ2 instead of Q2 as well as instead of diff between R2 and Q2. THe value is applied during pcorr optimization in model strategy 3 and 4 and during selection of othogonals in model strategy 2 to 4 and while selecting best iteration model in model strategy 5. Minimum value is 1/no_permutations_post_vs.
#' @param pcorr_diff Numeric. A value by which pcorr cutoff is decreased to allow more variables in the model if Q2 is increased more than 1 percent. Used in model strategy 4.
#' @param variable_names_length Numeric. Number of characters of variablenames shown in loadingplot or "all"
#' @param variable_names_position "beginning", "end" or "all" give part of string of variablenames shown in loadingplot
#' @param cluster "no" if analysis is run locally and "yes" if script is run using the bash scripts to run the R-scripts passing modelnumbers to run as arguments ("args").
#' @param name_intermediate_dir "name_of_directory" enter name of intermediate directory. On UPPMAX it is "SNIC_TMP" else enter "default"
#' @param each_model_or_summary "summary" if only summary should be run and "each" if only models of each comparison should be run. "both" if first running each comparison followed by summary. "both" is not possible if cluster is "yes". "summary" requires that .Rdata files from "each" exists in outputR folder.
#' @param model_strategies_to_run Numeric vector indicating which model strategies to run. Model strategy 3 requires that 2 is run. Model strategy 4 requires that 3 and 2 is run.
#'
#' @return Outputs in a folder called OutputR in your analysis folder. one .html file and one .Rdata file
#' per group comparison  containing five models and one summary html file containing
#' tables of all models of all comparisons. Also tables with loadings are created.
#'
#' @examples
#' Ropls_with_permutated_variable_selection(directory_of_Ropls_with_permutated_variable_selection,directory_of_analysis,
#' projectname,date_of_analysis,filename_matrix,decimal_separator,filename_sampleID,colname_groupID,groupsnumeric,
#' colname_secID,no_permutations_post_vs,no_permutations_post_vs_selected_models,no_permutations_over_vs,
#' p_pearson_of_pcorr_cutoff,setseedfirstmodel,order_of_groups,models_to_run,model_strategies_to_run,
#' foldername_Rmarkdownfiles,foldername_of_input_matrix_and_sampleID,foldername_output_reports,foldername_function_file,
#' directory_input_matrix_sampleID,directory_output_reports,filename_Rmarkdownfile_each_model,filename_Rmarkdownfile_summary,
#' filename_function_file,filter_percent_in_each_group,pcorr_cutoff_Model1_joint_models,pcorr_cutoff_Model1_stratified_models,
#' variable_selection_using_VIP,no_of_ortho_pre_vs_Model1_joint_models,no_of_ortho_pre_vs_Model1_stratified_models,
#' no_of_ortho_post_vs_Model1_joint_models,no_of_ortho_post_vs_Model1_stratified_models, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs,
#' prefered_pR2_and_pQ2_permutated_post_vs,pcorr_diff,variable_names_length,variable_names_position,
#' cluster,name_intermediate_dir,each_model_or_summary,model_strategies_to_run)
#' @details In filename_sampleID and in filename_matrix do not use the following symbols in sampleID or variable names;
#' ?, $, %, ^, &, *, (, ), -, #, ?, ,, <, >, /, |, , ], {, } and [
#' Missing values should be indicated by "", "NA" or "Inf"
#' @export
Ropls_with_permutated_variable_selection <- function(directory_of_Ropls_with_permutated_variable_selection,directory_of_analysis,
                                                         projectname,date_of_analysis,filename_matrix,decimal_separator,filename_sampleID,colname_groupID,groupsnumeric,
                                                         colname_secID,no_permutations_post_vs,no_permutations_post_vs_selected_models,no_permutations_over_vs,
                                                         p_pearson_of_pcorr_cutoff,setseedfirstmodel,order_of_groups,models_to_run,
                                                         foldername_Rmarkdownfiles,foldername_of_input_matrix_and_sampleID,foldername_output_reports,foldername_function_file,
                                                         directory_input_matrix_sampleID,directory_output_reports,filename_Rmarkdownfile_each_model,filename_Rmarkdownfile_summary,
                                                         filename_function_file,replace_0, filter_percent_in_each_group, replace_NA, log_transform,pcorr_cutoff_Model1_joint_models,pcorr_cutoff_Model1_stratified_models,
                                                         variable_selection_using_VIP,no_of_ortho_pre_vs_Model1_joint_models,no_of_ortho_pre_vs_Model1_stratified_models,
                                                         no_of_ortho_post_vs_Model1_joint_models,no_of_ortho_post_vs_Model1_stratified_models,max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs,
                                                         prefered_pR2_and_pQ2_permutated_post_vs,pcorr_diff,variable_names_length,variable_names_position,
                                                         cluster,name_intermediate_dir,each_model_or_summary,model_strategies_to_run) {

#Read directories of Function files and Rmarkdown files

directory_Rmarkdownfiles <- paste(directory_of_Ropls_with_permutated_variable_selection,"/",foldername_Rmarkdownfiles,"/",sep="") # "path of directory/filename"
directory_function_file <- paste(directory_of_Ropls_with_permutated_variable_selection,"/",foldername_function_file,"/", sep="")
foldername_model_table_to_analyse <- foldername_output_reports #model_table_to_analyse will be saved in this folder.
directory_model_table_to_analyse <- paste(directory_of_analysis,"/",foldername_model_table_to_analyse,"/",sep="")

#Load functions

directory_and_filename_function_file <- paste(directory_function_file,"/",filename_function_file,sep="")
source(directory_and_filename_function_file)

#Create directories.
makeproject(directory_of_analysis, directory_output_reports)
if (length(dir(directory_input_matrix_sampleID)) ==0) {directory_input_matrix_sampleID <- paste(directory_of_Ropls_with_permutated_variable_selection,"/inst/extdata",sep="")}
#Reorder groups to have diseased first if necessary
setwd(directory_input_matrix_sampleID)
file_sampleID <- read.table(filename_sampleID,header=T, dec = ".", row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf",""), sep="\t")
levels_of_groups <- levels(as.factor(file_sampleID[,paste(colname_groupID)]))
levels_of_groups # Check if order of the groups are correct with diseased first and controls last
reordered_levels_of_groups <- reorder_levels_of_groups(order_of_groups,levels_of_groups) # Note that it has to be run even if order is "correct"
reordered_levels_of_groups #Check that reordered groups has diseased first and controls last.
order_of_groups_numeric <- match(levels_of_groups,reordered_levels_of_groups)

#Create or load model table to analyse containing the comparisons to be modeled.
filename_model_table_to_analyse <- "model_table_to_analyse.txt"
model_table_to_analyse <- create_or_load_Model_table_to_analyze(directory_input_matrix_sampleID, directory_model_table_to_analyse,filename_model_table_to_analyse,reordered_levels_of_groups, file_sampleID, filename_sampleID, colname_groupID, colname_secID, pcorr_cutoff_Model1_joint_models, no_of_ortho_pre_vs_Model1_joint_models, no_of_ortho_post_vs_Model1_joint_models, pcorr_cutoff_Model1_stratified_models, no_of_ortho_pre_vs_Model1_stratified_models, no_of_ortho_post_vs_Model1_stratified_models) # Checks if model_table_to_analyse exist and loads or creats it in the folder foldername_model_table_to_analyse

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
                          no_permutations_over_vs=no_permutations_over_vs,
                          replace_0=replace_0,
                          filter_percent_in_each_group=filter_percent_in_each_group,
                          replace_NA=replace_NA,
                          log_transform=log_transform,
                          directory_and_filename_function_file=directory_and_filename_function_file,
                          max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs,
                          max_no_of_ortho_post_vs=max_no_of_ortho_post_vs,
                          reordered_levels_of_groups=reordered_levels_of_groups,
                          prefered_pR2_and_pQ2_permutated_post_vs=prefered_pR2_and_pQ2_permutated_post_vs,
                          pcorr_diff=pcorr_diff,
                          variable_selection_using_VIP=variable_selection_using_VIP,
                          model_strategies_to_run=model_strategies_to_run
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
                          no_permutations_over_vs=no_permutations_over_vs,
                          replace_0=replace_0,
                          filter_percent_in_each_group=filter_percent_in_each_group,
                          replace_NA=replace_NA,
                          log_transform=log_transform,
                          directory_and_filename_function_file=directory_and_filename_function_file,
                          max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs,
                          max_no_of_ortho_post_vs=max_no_of_ortho_post_vs,
                          reordered_levels_of_groups=reordered_levels_of_groups,
                          prefered_pR2_and_pQ2_permutated_post_vs=prefered_pR2_and_pQ2_permutated_post_vs,
                          pcorr_diff=pcorr_diff,
                          variable_selection_using_VIP=variable_selection_using_VIP,
                          model_strategies_to_run=model_strategies_to_run

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
                      no_permutations_over_vs=no_permutations_over_vs,
                      no_permutations_post_vs_selected_models=no_permutations_post_vs_selected_models,
                      filter_percent_in_each_group=filter_percent_in_each_group,
                      groupsnumeric=groupsnumeric,
                      model_strategies_to_run=model_strategies_to_run,
                      max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs,
                      max_no_of_ortho_post_vs=max_no_of_ortho_post_vs,
                      prefered_pR2_and_pQ2_permutated_post_vs=prefered_pR2_and_pQ2_permutated_post_vs,
                      reordered_levels_of_groups=reordered_levels_of_groups,
                      pcorr_diff=pcorr_diff,
                      variable_selection_using_VIP=variable_selection_using_VIP))
}
}


#' Share and unique structures plot
#'
#' Compares two models from roplspvs package using p(corr) of the variables
#'
#' Compares selected models created by the function oplspvs.
#' A new model is fitted using the variables from both models. P(corr) from the original model is used for variables that were alo in the original models
#' and the refitted p(corr) is used for the variables that were uniquely in the other model. P(corr) is the correlation
#' between the scores of the model and the raw data.
#'
#'
#' @param directory_output_reports_modelX "path" Path to "outputR" where "output_reports" are created for model X .
#' @param projectname_modelX "projectname" as written in the Configure_Get_Strarted.R for model X
#' @param projectname_in_plot_modelX "projectname" as desired in the plot for model X
#' @param date_of_analysis_modelX date_of_analysis as given in Configure_Get_Started.R for model X
#' @param group1_modelX Diseased group of model X
#' @param group2_modelX Control group of model X
#' @param secID_modelX stratification group of model X. If no stratification write "joint"
#' @param result_modelX_name Name of model to be compared for model X. Either "result_Model1", "result_Model2", "result_Model3", "result_Model4" or "result_iterationmodel"
#' @param directory_output_reports_modelY "path" Path to "outputR" where "output_reports" are created for model Y .
#' @param projectname_modelY "projectname" as written in the Configure_Get_Started.R for model Y
#' @param projectname_in_plot_modelY "projectname" as desired in the plot for model Y
#' @param date_of_analysis_modelY date_of_analysis as given in Configure_Get_Strarted.R for model Y
#' @param group1_modelY Diseased group of model Y
#' @param group2_modelY Control group of model Y
#' @param secID_modelY stratification group of model Y. If no stratification write "joint"
#' @param result_modelY_name Name of model to be compared for model Y. Either "result_Model1", "result_Model2", "result_Model3", "result_Model4" or "result_iterationmodel"
#' @param variable_names_length Number of characters of variablenames shown in loadingplot or "all"
#' @param variable_names_position part of variablenames shown "beginning", "end" or "all"
#' @param title_length "long"  to write group comparison in title, "short" to write no comparison in title or FALSE to have no title
#' @param change_name_list list for changing names of variables in susplot in format c(c(change from,change to)), "edit" to change " " and "_" to "." or FALSE to not change
#' @param show_all_lables TRUE show all labels and FALSE excludes text labels that overlap too much
#'
#' @return Outputs a SUSplot figure
#'
#' @examples
#'
#' ## To compare two models from oplspvs:
#' SUSplot(directory_output_reports_modelX, projectname_modelX, projectname_in_plot_modelX, date_of_analysis_modelX, group1_modelX, group2_modelX, secID_modelX, result_modelX_name, directory_output_reports_modelY, projectname_modelY, projectname_in_plot_modelY, date_of_analysis_modelY, group1_modelY, group2_modelY, secID_modelY, result_modelY_name, variable_names_position, variable_names_length, title_length, change_name_list, show_all_lables)
#' @export
SUSplot <- function(directory_output_reports_modelX, projectname_modelX, projectname_in_plot_modelX, date_of_analysis_modelX, group1_modelX, group2_modelX, secID_modelX, result_modelX_name,
                    directory_output_reports_modelY, projectname_modelY, projectname_in_plot_modelY, date_of_analysis_modelY, group1_modelY, group2_modelY, secID_modelY, result_modelY_name, variable_names_position, variable_names_length, title_length,change_name_list, show_all_lables){
  library(ropls)
  library(ggplot2)
  library("ggrepel")
  modelX <- loadRData(directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename=result_modelX_name)
  variablelist1 <- modelX$pcorrlistaftervs
  modelY <- loadRData(directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename=result_modelY_name)
  variablelist2 <- modelY$pcorrlistaftervs
  subsetdatamatrix1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="subsetdatamatrix")
  class1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="classordered")
  subsetdatamatrix2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="subsetdatamatrix")
  class2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="classordered")




  #B is pcorr before variable selection variables from both models
  pcorrlist1B <- as.data.frame(as.data.frame(modelX$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
  rownames(pcorrlist1B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
  colnames(pcorrlist1B) <- "pcorrlist1B"
  choosecolumn <- rownames(pcorrlist1B) %in% intersect(row.names(pcorrlist1B), colnames(subsetdatamatrix1))
  pcorrlist1BE <- subset(pcorrlist1B, choosecolumn)# remove variables filtered away from subsetmatrix

  pcorrlist2B <- as.data.frame(as.data.frame(modelY$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
  rownames(pcorrlist2B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
  colnames(pcorrlist2B) <- "pcorrlist2B"
  choosecolumn <- rownames(pcorrlist2B) %in% intersect(row.names(pcorrlist2B), colnames(subsetdatamatrix2))
  pcorrlist2BE <- subset(pcorrlist2B, choosecolumn)# remove variables filtered away from subsetmatrix

  subsetdatamatrix1B <-subsetdatamatrix1[,row.names(pcorrlist1BE)]
  oplsda1C = opls(subsetdatamatrix1B, class1, predI = 1, orthoI = modelX$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

  subsetdatamatrix2B <-subsetdatamatrix2[,row.names(pcorrlist2BE)]
  oplsda2C = opls(subsetdatamatrix2B, class2, predI = 1, orthoI = modelY$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


  #C is fit model with variables from both models
  pcorrlist1C <-vector()
  for (i in 1:ncol(subsetdatamatrix1B)){
    pcorrtest <- cor.test(subsetdatamatrix1B[,i],oplsda1C@scoreMN,method="pearson")
    pcorrlist1C[i] <- pcorrtest$estimate
  }
  pcorrlist1C <- as.data.frame(pcorrlist1C)
  row.names(pcorrlist1C)<- colnames(subsetdatamatrix1B)
  pcorrlist2C<-vector()
  for (i in 1:ncol(subsetdatamatrix2B)){
    pcorrtest <- cor.test(subsetdatamatrix2B[,i],oplsda2C@scoreMN,method="pearson")
    pcorrlist2C[i] <- pcorrtest$estimate
  }
  pcorrlist2C <- as.data.frame(pcorrlist2C)
  row.names(pcorrlist2C)<- colnames(subsetdatamatrix2B)


  #A is fit model with variables from one model
  subsetdatamatrix1A <-subsetdatamatrix1[,row.names(variablelist1)]
  oplsda1A = opls(subsetdatamatrix1A, class1, predI = 1, orthoI = modelX$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

  subsetdatamatrix2A <-subsetdatamatrix2[,row.names(variablelist2)]
  oplsda2A = opls(subsetdatamatrix2A, class2, predI = 1, orthoI = modelY$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


  pcorrlist1A<-vector()
  for (i in 1:ncol(subsetdatamatrix1A)){
    pcorrtest1A <- cor.test(subsetdatamatrix1A[,i],oplsda1A@scoreMN,method="pearson")
    pcorrlist1A[i] <- pcorrtest1A$estimate
  }
  pcorrlist1A <- as.data.frame(pcorrlist1A)
  row.names(pcorrlist1A)<- colnames(subsetdatamatrix1A)

  pcorrlist2A<-vector()
  for (i in 1:ncol(subsetdatamatrix2A)){
    pcorrtest2A <- cor.test(subsetdatamatrix2A[,i],oplsda2A@scoreMN,method="pearson")
    pcorrlist2A[i] <- pcorrtest2A$estimate
  }
  pcorrlist2A <- as.data.frame(pcorrlist2A)
  row.names(pcorrlist2A)<- colnames(subsetdatamatrix2A)


  # Replace C pcorr in model with variables from both models with A pcorr form model with variables from one model
  for (i in 1:length(pcorrlist1A)){
    pcorrlist1C[rownames(pcorrlist1A),] <- pcorrlist1A[i]
  }

  for (i in 1:length(pcorrlist2A)){
    pcorrlist2C[rownames(pcorrlist2A),] <- pcorrlist2A[i]
  }

  pcorrlist1C <- as.data.frame(pcorrlist1C[unique(c(rownames(variablelist1),rownames(variablelist2))),])
  rownames(pcorrlist1C)<-unique(c(rownames(variablelist1),rownames(variablelist2)))
  pcorrlist2C <- as.data.frame(pcorrlist2C[unique(c(rownames(variablelist1),rownames(variablelist2))),])
  rownames(pcorrlist2C)<-unique(c(rownames(variablelist1),rownames(variablelist2)))
  colnames(pcorrlist1C)<-"pcorrlist1C"
  colnames(pcorrlist2C)<-"pcorrlist2C"
  SUSplot <- as.data.frame(cbind(pcorrlist1C,pcorrlist2C))
  SUSplot[is.na(SUSplot)] <- 0
  size <- 18
  colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(variablelist1),rownames(variablelist2)),rownames(SUSplot) %in% setdiff(rownames(variablelist2),rownames(variablelist1)),rownames(SUSplot) %in% intersect(rownames(variablelist1),rownames(variablelist2)))
  colorinplotvector <- vector()
  for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(group1_modelX," vs ",group2_modelX,secID_modelX)} else if(colorinplot[j,2]) {paste(group1_modelY," vs ",group2_modelY,secID_modelY)} else {"shared"}}
  susplotnames <- SUSplot
  rownames(susplotnames) <- gsub("_",".",rownames(susplotnames))
  rownames(susplotnames) <- gsub(" ",".",rownames(susplotnames))
  if (variable_names_position=="beginning") {
    rownames(susplotnames) <- substr(rownames(susplotnames), 1, variable_names_length)}
  if (variable_names_position=="end") {
    for (novariables in 1:nrow(susplotnames)) {
      rownames(susplotnames)[novariables] <- substr(rownames(susplotnames)[novariables], nchar(rownames(susplotnames)[novariables])-(variable_names_length-1), nchar(rownames(susplotnames)[novariables]))
      blankpos = StrPos(rownames(susplotnames)[novariables], '\\.', 1);
      if (!is.na(blankpos)&nchar(rownames(susplotnames)[novariables])>=variable_names_length) {
        rownames(susplotnames)[novariables] <- substr(rownames(susplotnames)[novariables], blankpos+1, nchar(rownames(susplotnames)[novariables]));
      }}}
if(is.logical(change_name_list)){if(change_name_list==F) susplotnames <- SUSplot}
  if(is.list(change_name_list)){
  change_name_list_df<-as.data.frame(change_name_list)
  names(change_name_list_df) <- NULL
  row.names(susplotnames)[match(change_name_list_df[1,],row.names(SUSplot))]<-change_name_list_df[2,]
  }


  colors_manual <- c("blue", "darkgreen","red")[c(sum(colorinplot[,1])!=0,sum(colorinplot[,2])!=0,sum(colorinplot[,3])!=0)]

  pC1 <- ggplot(SUSplot, aes(x=pcorrlist1C,y=pcorrlist2C, color=colorinplotvector))
  pC2 <- pC1 + geom_point()
  if (projectname_modelX==projectname_modelY) {
    if (title_length=="long") {
    pC3 <- pC2 + labs(y=paste(group1_modelY,"vs",group2_modelY, secID_modelY, "[p(corr)]"), x=paste(group1_modelX,"vs",group2_modelX, secID_modelX, "[p(corr)]"),
                     title=paste("SUS plot OPLS models" ,projectname_in_plot_modelX, group1_modelX,"vs", group2_modelX, secID_modelX, "\ncompared to ", projectname_in_plot_modelY, group1_modelY,"vs",group2_modelY, secID_modelY))
    } else {pC3 <- pC2 + labs(y=paste(group1_modelY,"vs",group2_modelY, secID_modelY, "[p(corr)]"), x=paste(group1_modelX,"vs",group2_modelX, secID_modelX, "[p(corr)]"),
                                title=paste("SUS plot OPLS models" ,projectname_in_plot_modelX))}
  } else if (title_length=="long") {
  pC3 <- pC2 + labs(y=paste(projectname_in_plot_modelY,group1_modelY,"vs",group2_modelY, secID_modelY, "[p(corr)]"), x=paste(projectname_in_plot_modelX, group1_modelX,"vs",group2_modelX, secID_modelX, "[p(corr)]"),
                    title=paste("SUS plot OPLS models" ,projectname_in_plot_modelX, group1_modelX,"vs", group2_modelX, secID_modelX, "\ncompared to ", projectname_in_plot_modelY, group1_modelY,"vs",group2_modelY, secID_modelY))
  } else {pC3 <- pC2 + labs(y=paste(projectname_in_plot_modelY,group1_modelY,"vs",group2_modelY, secID_modelY, "[p(corr)]"), x=paste(projectname_in_plot_modelX, group1_modelX,"vs",group2_modelX, secID_modelX, "[p(corr)]"),
                             title=paste("SUS plot OPLS models"))}
  if (title_length==FALSE) {pC4 <- pC3 + theme(plot.title = element_blank())
} else {pC4 <- pC3}
  pC5 <- pC4 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
  if (show_all_lables==TRUE) {
  pC6 <- pC5 + geom_text_repel(label=rownames(susplotnames), max.overlaps = Inf)
  } else {pC6 <- pC5 + geom_text_repel(label=rownames(susplotnames))}
  pC7 <- pC6 + theme(legend.position = "top",legend.title = element_blank())
  pC8 <- pC7 + scale_color_manual(values=colors_manual)
  pC9 <- pC8 + scale_y_continuous(breaks = c(seq(-1, 1, by=0.5)), minor_breaks = c(seq(-1, 1, by=0.1)),limits = c(-1, 1))
  pC10 <- pC9 + scale_x_continuous(breaks = c(seq(-1, 1, by=0.5)), minor_breaks = c(seq(-1, 1, by=0.1)),limits = c(-1, 1))
  pC10
}


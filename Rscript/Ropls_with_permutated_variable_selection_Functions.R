#Function file used in Ropls-with-permutated-variable-selection Version 0.15.0

subsetmatrixfunction <- function(sampleID,datamatrix,group1,group2,secID){
  if(secID=="joint"|secID=="no gender stratification"){
    ##Subset datamatrix
    subsetdatamatrix <- subset(datamatrix,sampleID[,paste(colname_groupID)]==group1 | sampleID[,paste(colname_groupID)]==group2)
  } else {
    ##Subset datamatrix
    subsetdatamatrix <- subset(datamatrix,sampleID[,paste(colname_groupID)]==group1 & sampleID[,paste(colname_secID)]==secID | sampleID[,paste(colname_groupID)]==group2 & sampleID[,paste(colname_secID)]==secID)
  }
subsetdatamatrix
}

subsetsampleIDfunction <- function(sampleID,group1,group2,secID){
  if(secID=="joint"|secID=="no gender stratification"){
    ## subset sampleID
    subsetsampleID <- subset(sampleID,sampleID[,paste(colname_groupID)]==group1 | sampleID[,paste(colname_groupID)]==group2)

  } else {
    ## subset sampleID
    subsetsampleID <- subset(sampleID,sampleID[,paste(colname_groupID)]==group1 & sampleID[,paste(colname_secID)]==secID | sampleID[,paste(colname_groupID)]==group2 & sampleID[,paste(colname_secID)]==secID)

  }
  subsetsampleID[,paste(colname_groupID)] <- as.factor(subsetsampleID[,paste(colname_groupID)])
  subsetsampleID
}

filterNAfunction <- function(subsetsampleID,subsetdatamatrix,group1,group2,secID,filter_percent_in_each_group){
#create matrix for each group
tsubsetdatamatrix <- t(subsetdatamatrix)
subsetdatamatrixgroup1 <- subset(subsetdatamatrix,subsetsampleID[,paste(colname_groupID)]==group1)
subsetdatamatrixgroup2 <- subset(subsetdatamatrix,subsetsampleID[,paste(colname_groupID)]==group2)
tsubsetdatamatrixgroup1 <- t(subsetdatamatrixgroup1)
tsubsetdatamatrixgroup2 <- t(subsetdatamatrixgroup2)

#count NAs for each peptide in each group
percentNAgroup1 <- rowCounts(tsubsetdatamatrixgroup1, value=NA, na.rm=FALSE)/ncol(tsubsetdatamatrixgroup1)*100
percentNAgroup2 <- rowCounts(tsubsetdatamatrixgroup2, value=NA, na.rm=FALSE)/ncol(tsubsetdatamatrixgroup2)*100

# keep all peptides that have values in at least percentage given in function of the subject in at least one group
percentNAs <- cbind(percentNAgroup1, percentNAgroup2)
minpercentNAs <- apply(percentNAs, 1, max)
minpercentNAsdf <- as.data.frame(minpercentNAs)
subsetdatamatrixNAfiltered <- subsetdatamatrix[,minpercentNAsdf<=filter_percent_in_each_group]
subsetdatamatrix <- subsetdatamatrixNAfiltered
subsetdatamatrix
}

   ##  ............................................................................
  ##  define function entering number of orthogonal variables and giving model statistics after selection default p(corr)>0.4   ####

  ### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
  ### oplsmodelwithvariableselection                                           ####


opls_model_with_variable_selection_trycatch <- function(subsetdatamatrix,ortho_pre_vs,ortho_post_vs,class, pcorr, printoptmodel="none",plotoptmodel="none", no_permutations_post_vs, variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){

  result <- tryCatch({
          result_sans_vs <- tryCatch({
        beforevsdata.oplsda <- opls(subsetdatamatrix, class, predI = 1, orthoI = ortho_pre_vs, scaleC="standard",info.txtC="none",fig.pdfC="none", permI=0)
        if (beforevsdata.oplsda@summaryDF$ort>max_no_of_ortho_pre_vs) {
          beforevsdata.oplsda <- opls(subsetdatamatrix, class, predI = 1, orthoI = max_no_of_ortho_pre_vs, scaleC="standard",info.txtC="none",fig.pdfC="none", permI=0)
        }
        beforevsdata.oplsda
      },
      error = function(e){
        beforevsdata.oplsda <- opls(subsetdatamatrix, class, predI = 1, orthoI = 0, scaleC="standard",info.txtC="none",fig.pdfC="none", permI=0)
        beforevsdata.oplsda
        }
      )
      beforevsdata.oplsda <- result_sans_vs
      vipropls <- getVipVn(beforevsdata.oplsda)
      summaryropls <- getSummaryDF(beforevsdata.oplsda)
      loadingropls <- getLoadingMN(beforevsdata.oplsda)
      Scoreofvariables <- getScoreMN(beforevsdata.oplsda)

      ### Calculate correlation between x and score to give p(corr)

      pcorrofvariables <-vector()
      for (i in 1:ncol(subsetdatamatrix)){
        pcorrofvariableshtest <- cor.test(subsetdatamatrix[,i],Scoreofvariables,method="pearson")
        pcorrofvariables[i] <- pcorrofvariableshtest$estimate
      }
      pcorrofvariables <- as.data.frame(pcorrofvariables)

      row.names(pcorrofvariables) <- colnames(subsetdatamatrix)
      pcorrofvariables <- na.omit(pcorrofvariables)
      ### pcorrlist of selected variables
      choosingsubset <- pcorrofvariables$pcorrofvariables>pcorr|pcorrofvariables$pcorrofvariables<(-(pcorr))
      choosingVIP <- vipropls>1
      choosingVIPandpcorr <- choosingsubset & choosingVIP
      subsetpcorrandVIPofvariables <- subset(pcorrofvariables,choosingVIPandpcorr)
      subsetpcorrofvariables <- subset(pcorrofvariables,choosingsubset)
      if (variable_selection_using_VIP=="yes") {
        pcorrlist_pre_vs_selected <- subsetpcorrandVIPofvariables} else {
          pcorrlist_pre_vs_selected <- subsetpcorrofvariables
      }

      if (nrow(pcorrlist_pre_vs_selected)>ortho_post_vs+1|(is.na(ortho_post_vs)&nrow(pcorrlist_pre_vs_selected)>1)) {

        ### model after variable selection
        aftervariableselectionsubsetdatamatrix <- subsetdatamatrix[,rownames(pcorrlist_pre_vs_selected)]


        result_post_vs <- tryCatch({
          aftervsdata.oplsda = opls(aftervariableselectionsubsetdatamatrix, class, predI = 1, orthoI = ortho_post_vs, scaleC="standard",info.txtC=printoptmodel,fig.pdfC=plotoptmodel,permI=0, .sinkC=NULL)

          if (aftervsdata.oplsda@summaryDF$ort>max_no_of_ortho_post_vs) {
            aftervsdata.oplsda = opls(aftervariableselectionsubsetdatamatrix, class, predI = 1, orthoI = max_no_of_ortho_post_vs, scaleC="standard",info.txtC=printoptmodel,fig.pdfC=plotoptmodel,permI=no_permutations_post_vs, .sinkC=NULL)

          } else {
          aftervsdata.oplsda = opls(aftervariableselectionsubsetdatamatrix, class, predI = 1, orthoI = ortho_post_vs, scaleC="standard",info.txtC=printoptmodel,fig.pdfC=plotoptmodel,permI=no_permutations_post_vs, .sinkC=NULL)
          }
          aftervsdata.oplsda
          },
        error = function(e){
          aftervsdata.oplsda = opls(aftervariableselectionsubsetdatamatrix, class, predI = 1, orthoI = 0, scaleC="standard",info.txtC=printoptmodel,fig.pdfC=plotoptmodel,permI=no_permutations_post_vs, .sinkC=NULL)
          aftervsdata.oplsda
          }
        )
        aftervsdata.oplsda <- result_post_vs
        viproplsaftervs <- getVipVn(aftervsdata.oplsda)
        summaryroplsaftervs <- getSummaryDF(aftervsdata.oplsda)
        loadingroplsaftervs <- getLoadingMN(aftervsdata.oplsda)
        scoreofvariablesaftervs <- getScoreMN(aftervsdata.oplsda)
        variables_no <-length(loadingroplsaftervs)
        resultaftervs <- cbind(pcorr,beforevsdata.oplsda@summaryDF$ort,variables_no,summaryroplsaftervs)

        pcorrlistaftervs <-data.frame(matrix(ncol=3))
        for (i in 1:ncol(aftervariableselectionsubsetdatamatrix)){
          pcorrlistaftervshtest <- cor.test(aftervariableselectionsubsetdatamatrix[,i],scoreofvariablesaftervs,method="pearson")
          pcorrlistaftervs[i,] <- c(pcorrlistaftervshtest$estimate,pcorrlistaftervshtest$conf.int[1],pcorrlistaftervshtest$conf.int[2])
        }
        colnames(pcorrlistaftervs) <- c("pcorrlistaftervs","conf.int.low","conf.int.high")
        row.names(pcorrlistaftervs) <- colnames(aftervariableselectionsubsetdatamatrix)

        if ("pR2Y" %in% colnames(resultaftervs)) {colnames(resultaftervs) <- c("pcorr cutoff","ortho pre v.s.","no. variables","R2X(cum)","R2Y(cum)","Q2(cum)","RMSEE","pred. post v.s.","ortho post v.s.","pR2Y permutated post v.s.","pQ2 permutated post v.s.")} else {
          colnames(resultaftervs) <- c("pcorr cutoff","ortho pre v.s.","no. variables","R2X(cum)","R2Y(cum)","Q2(cum)","RMSEE","pred. post v.s.","ortho post v.s.")
          }
        rownames(resultaftervs) <- "model"
        resultoplsmodelwithvariableselection <- list(beforevsdata.oplsda=beforevsdata.oplsda, aftervsdata.oplsda=aftervsdata.oplsda,resultaftervs=resultaftervs,variables_no=variables_no, scoreofvariablesaftervs=scoreofvariablesaftervs, loadingroplsaftervs=loadingroplsaftervs,pcorrlistaftervs=pcorrlistaftervs)
        resultoplsmodelwithvariableselection
      } else {
        resultaftervs <- data.frame(matrix(NA,nrow=1, ncol=11))
        colnames(resultaftervs)<-c("pcorr cutoff","ortho pre v.s.","no. variables","R2X(cum)","R2Y(cum)","Q2(cum)","RMSEE","pred. post v.s.","ortho post v.s.","pR2Y permutated post v.s.","pQ2 permutated post v.s.")
        row.names(resultaftervs)<-c("Total")
        resultoplsmodelwithvariableselection <- list(resultaftervs=resultaftervs)
        resultoplsmodelwithvariableselection
      }

    },
    error = function(e){
      resultaftervs <- data.frame(matrix(NA,nrow=1, ncol=11))
      colnames(resultaftervs)<-c("pcorr cutoff","ortho pre v.s.","no. variables","R2X(cum)","R2Y(cum)","Q2(cum)","RMSEE","pred. post v.s.","ortho post v.s.","pR2Y permutated post v.s.","pQ2 permutated post v.s.")
      row.names(resultaftervs)<-c("Total")
      resultoplsmodelwithvariableselection <- list(resultaftervs=resultaftervs)
      resultoplsmodelwithvariableselection
    }

)
	sinkout()
	result
}

sinkall <- function() {
  i <- sink.number(type = "message")
  while (i > 0) {
    sink(type = "message")
    i <- i - 1
  }
}

sinkout <- function() {
  i <- sink.number(type = "output")
  while (i > 1) {
    sink(type = "output")
    i <- i - 1
  }
}


  ##  ............................................................................
  ##  run original model with different amount of orthogal variables in original model
  model_pre_vs_with_different_amount_of_ortho_pre_vs_table <- function(subsetdatamatrix, no_of_orthogonal_in_model_pre_vs, class, no_permutations_post_vs,variable_selection_using_VIP){

    amountofvariablesinmodel <- data.frame()
  for (i in 0:no_of_orthogonal_in_model_pre_vs)
  {
    result <-  opls(subsetdatamatrix, class, predI = 1, orthoI = i,fig.pdfC="none", info.txtC="none", permI=no_permutations_post_vs)
    result <- getSummaryDF(result)
    amountofvariablesinmodel <- rbind(amountofvariablesinmodel,result)
  }
  amountofvariablesinmodelwithdiff <- amountofvariablesinmodel
  amountofvariablesinmodelwithdiff$"diff R2Y(cum)-Q2(cum)" <- amountofvariablesinmodel$`R2Y(cum)`- amountofvariablesinmodel$`Q2(cum)`
  colnames(amountofvariablesinmodelwithdiff["ort"]) <- "ortho pre v.s."
  row.names(amountofvariablesinmodelwithdiff) <- c(paste("model",1:nrow(amountofvariablesinmodelwithdiff)))
  amountofvariablesinmodelwithdiff
  }


  ## plot diff vs Q2 of original model

  plotof_model_pre_vs_with_different_amount_of_ortho_pre_vs <- function(table_different_amount_of_ortho_pre_vs){

    amountofvariablesinmodelwithdiff <- table_different_amount_of_ortho_pre_vs
    amountofvariablesinmodelwithdiffplot <- amountofvariablesinmodelwithdiff
    amountofvariablesinmodelwithdiffplot$Q2cum <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
    amountofvariablesinmodelwithdiffplot$diff <- amountofvariablesinmodelwithdiffplot$"diff R2Y(cum)-Q2(cum)"
     ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, diff)) +
       geom_point() +
     geom_text_repel(label=amountofvariablesinmodelwithdiffplot$ort)+
            labs(y="difference between R2Y(cum) and Q2(cum)")+
       labs(x="Q2(cum) pre variable selection")+
     theme(text=element_text(size=14), axis.text=element_text(size=14))

  }


  ##  ............................................................................
  ##  run optimized model with different amount of orthogal variables in original model                                                                       ####

  model_post_vs_with_different_amount_of_ortho_pre_vs_table <- function(subsetdatamatrix, no_of_orthogonal_in_model_pre_vs,ortho_post_vs,class, pcorr, no_permutations_post_vs, variable_selection_using_VIP,max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
    amountofvariablesinmodel <- data.frame()
    for (i in 0:no_of_orthogonal_in_model_pre_vs)
    {
        result <-  opls_model_with_variable_selection_trycatch(subsetdatamatrix, i,ortho_post_vs,class, pcorr, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)
        amountofvariablesinmodel <- rbind(amountofvariablesinmodel,result$resultaftervs)

    }
    amountofvariablesinmodelwithdiff <- amountofvariablesinmodel
    amountofvariablesinmodelwithdiff$"diff R2Y(cum)-Q2(cum)" <- amountofvariablesinmodel$`R2Y(cum)`- amountofvariablesinmodel$`Q2(cum)`
    row.names(amountofvariablesinmodelwithdiff) <- c(paste("model",1:nrow(amountofvariablesinmodelwithdiff)))
    amountofvariablesinmodelwithdiff
  }


  ##  ............................................................................
  ##  run optimized model with different amount of orthogal variables in original and optimized model                                                                       ####

  model_post_vs_table_with_different_amount_of_ortho_pre_and_post_vs <- function(subsetdatamatrix, no_of_orthogonal_in_model_pre_vs,no_of_orthogonal_in_model_post_vs,class, pcorr, no_permutations_post_vs,variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
  amountofvariablesinmodel <- data.frame()
  for (i in 0:no_of_orthogonal_in_model_pre_vs)
  {
    for (j in 0:no_of_orthogonal_in_model_post_vs)
    {
      result <-  opls_model_with_variable_selection_trycatch(subsetdatamatrix, i,j,class, pcorr, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)
      amountofvariablesinmodel <- rbind(amountofvariablesinmodel,result$resultaftervs)
    }
  }
  amountofvariablesinmodelwithdiff <- amountofvariablesinmodel
  amountofvariablesinmodelwithdiff$"diff R2Y(cum)-Q2(cum)" <- amountofvariablesinmodel$`R2Y(cum)`- amountofvariablesinmodel$`Q2(cum)`
  row.names(amountofvariablesinmodelwithdiff) <- c(paste("model",1:nrow(amountofvariablesinmodelwithdiff)))
  amountofvariablesinmodelwithdiff

  }

  ##plot of optimimized model table
  plotof_model_post_vs_table_with_different_amount_of_ortho_pre_and_post_vs <- function(
  model_post_vs_table, choosen_model){

    amountofvariablesinmodelwithdiffplot <- model_post_vs_table
    amountofvariablesinmodelwithdiffplot$`Q2cum` <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
    amountofvariablesinmodelwithdiffplot$diff <- amountofvariablesinmodelwithdiffplot$"diff R2Y(cum)-Q2(cum)"
    amountofvariablesinmodelwithdiffplot$"ortho pre v.s." <- as.character(amountofvariablesinmodelwithdiffplot$"ortho pre v.s.")
    ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, diff, color=c(amountofvariablesinmodelwithdiffplot$`ortho pre v.s.` == choosen_model$`ortho pre v.s.` & amountofvariablesinmodelwithdiffplot$`ortho post v.s.` == choosen_model$`ortho post v.s.`))) +
      geom_point() +
      geom_text_repel(label=paste(amountofvariablesinmodelwithdiffplot$"ortho pre v.s.",",",amountofvariablesinmodelwithdiffplot$"ortho post v.s.")) +
      labs(y="difference between R2Y(cum) and Q2(cum)", x="Q2(cum) post variable selection") +
      theme(text=element_text(size=14), axis.text=element_text(size=14)) +
      scale_color_manual(values=c("black","red"), name=NULL) +
      guides(fill=FALSE, color=FALSE)

  }
  ##plot pcorrtable
  plot_pcorrtable_diff_with_pcorr_labels <- function(
    pcorrtable,max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables){

    amountofvariablesinmodelwithdiffplot <- pcorrtable
    amountofvariablesinmodelwithdiffplot$`Q2cum` <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
    amountofvariablesinmodelwithdiffplot$diff <- amountofvariablesinmodelwithdiffplot$"diff R2Y(cum)-Q2(cum)"
    amountofvariablesinmodelwithdiffplot$"ortho pre v.s." <- as.character(amountofvariablesinmodelwithdiffplot$"ortho pre v.s.")
    ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, diff, color=c(amountofvariablesinmodelwithdiffplot$`no. variables`==max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables$`no. variables`))) +
      geom_point() +
      geom_text_repel(label=paste(amountofvariablesinmodelwithdiffplot$`pcorr cutoff`)) +
      labs(y="difference between R2Y(cum) and Q2(cum)", x="Q2(cum)") +
      theme(text=element_text(size=14), axis.text=element_text(size=14)) +
      scale_color_manual(values=c("black","red"), name=NULL) +
      guides(fill=FALSE, color=FALSE)
  }
  plot_pcorrtable_diff_with_amount_of_variables_labels <- function(
    pcorrtable,max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables){

    amountofvariablesinmodelwithdiffplot <- pcorrtable
    amountofvariablesinmodelwithdiffplot <- na.omit(amountofvariablesinmodelwithdiffplot)
    amountofvariablesinmodelwithdiffplot$`Q2cum` <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
    amountofvariablesinmodelwithdiffplot$diff <- amountofvariablesinmodelwithdiffplot$"diff R2Y(cum)-Q2(cum)"
    amountofvariablesinmodelwithdiffplot$"ortho pre v.s." <- as.character(amountofvariablesinmodelwithdiffplot$"ortho pre v.s.")
    ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, diff,color=c(amountofvariablesinmodelwithdiffplot$`no. variables`==max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables$`no. variables`))) +
      geom_point() +
      geom_text_repel(label=paste(amountofvariablesinmodelwithdiffplot$`no. variables`)) +
      labs(y="difference between R2Y(cum) and Q2(cum)", x="Q2(cum)") +
      theme(text=element_text(size=14), axis.text=element_text(size=14)) +
      scale_color_manual(values=c("black","red"), name=NULL) +
      guides(fill=FALSE, color=FALSE)
  }

  plot_pcorrtable_amount_of_variables_with_pcorr_labels <- function(
    pcorrtable,max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables){

    amountofvariablesinmodelwithdiffplot <- pcorrtable
    amountofvariablesinmodelwithdiffplot <- na.omit(amountofvariablesinmodelwithdiffplot)
    amountofvariablesinmodelwithdiffplot$`Q2cum` <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
     ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, `no. variables`,color=c(amountofvariablesinmodelwithdiffplot$`no. variables`==max_pcorrtable_with_max_Q2_and_diff_less_than_02_few_variables$`no. variables`))) +
      geom_point() +
      geom_text_repel(label=paste(amountofvariablesinmodelwithdiffplot$`pcorr cutoff`), max.overlaps=50) +
      labs(y="Number of variables", x="Q2(cum) post variable selection") +
      theme(text=element_text(size=14), axis.text=element_text(size=14)) +
      scale_color_manual(values=c("black","red"), name=NULL) +
      guides(fill=FALSE, color=FALSE)
  }

    ## plot diff vs Q2 optimized model

 plotof_model_post_vs_with_different_amount_of_ortho_pre_vs <- function(table_model_post_vs_different_amount_of_ort_in_model_pre_vs){

   amountofvariablesinmodelwithdiff <- table_model_post_vs_different_amount_of_ort_in_model_pre_vs
   amountofvariablesinmodelwithdiffplot <- amountofvariablesinmodelwithdiff
  amountofvariablesinmodelwithdiffplot$`Q2cum` <- amountofvariablesinmodelwithdiffplot$`Q2(cum)`
  amountofvariablesinmodelwithdiffplot$diff <- amountofvariablesinmodelwithdiffplot$"diff R2Y(cum)-Q2(cum)"
  amountofvariablesinmodelwithdiffplot$"ortho pre v.s." <- as.character(amountofvariablesinmodelwithdiffplot$"ortho pre v.s.")
  ggplot(amountofvariablesinmodelwithdiffplot, aes(Q2cum, diff)) +
    geom_point() +
    geom_text_repel(label=amountofvariablesinmodelwithdiffplot$"ortho pre v.s.") +
            labs(y="difference between R2Y(cum) and Q2(cum)", x="Q2(cum)") +
    theme(text=element_text(size=14), axis.text=element_text(size=14))

 }


  ##create table with different amount of pcorr during variable selection

  pcorrtableofmodelsaftervs_withdiff <- function(subsetdatamatrix, ortho_pre_vs,ortho_post_vs,class,pcorrvector, no_permutations_post_vs, variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
  pcorrtable <- data.frame()
  pcorrselections <- pcorrvector
  for (pcorr in pcorrselections){
    pcorrmodeli <- opls_model_with_variable_selection_trycatch(subsetdatamatrix, ortho_pre_vs=ortho_pre_vs,ortho_post_vs=ortho_post_vs,class=class, pcorr=pcorr, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)
    pcorrtable <- rbind(pcorrtable, pcorrmodeli$resultaftervs)
  }
  pcorrtablewithdiff <- pcorrtable
  pcorrtablewithdiff$"diff R2Y(cum)-Q2(cum)" <- pcorrtable$`R2Y(cum)`- pcorrtable$`Q2(cum)`

  row.names(pcorrtablewithdiff) <- c(paste("model",1:nrow(pcorrtablewithdiff)))
  pcorrtablewithdiff
  }


  plotloading <- function(model){
  loadingropls <- getLoadingMN(model)
  loadingplot <- as.data.frame(loadingropls)
  pB1 <- ggplot(loadingplot, aes(x=reorder(row.names(loadingplot),-loadingplot$p1),y=loadingplot$p1))
  pB2 <- pB1 +geom_col()
  pB3 <- pB2 + theme(axis.text.x = element_text(angle = 90))
  pB4 <- pB3 + labs(y="loadings (p)", x=element_blank(),title="Loadingplot")
  print(pB4)
  }

  plotpcorr <- function(model, variable_names_length, variable_names_position){
    pcorrplot <- model$pcorrlistaftervs
    rownames(pcorrplot) <- gsub("_", " ", rownames(pcorrplot))
    rownames(pcorrplot) <- gsub("[.]", " ", rownames(pcorrplot))
    if (variable_names_position=="beginning") {
      rownames(pcorrplot) <- substr(rownames(pcorrplot), 1, variable_names_length)}
    if (variable_names_position=="end") {
      for (novariables in 1:nrow(pcorrplot)) {
      rownames(pcorrplot)[novariables] <- substr(rownames(pcorrplot)[novariables], nchar(rownames(pcorrplot)[novariables])-(variable_names_length-1), nchar(rownames(pcorrplot)[novariables]))
      blankpos = StrPos(rownames(pcorrplot)[novariables], ' ', 1);
      if (!is.na(blankpos)&nchar(rownames(pcorrplot)[novariables])>=variable_names_length) {
        rownames(pcorrplot)[novariables] <- substr(rownames(pcorrplot)[novariables], blankpos+1, nchar(rownames(pcorrplot)[novariables]));
      }}}

    if (length(rownames(pcorrplot))>40) {
      if(max(nchar(rownames(pcorrplot[1])))>=80) {
        fontsize <- 4
        widthsize <- 80
        lineheight <- 0.7
        message("more than 40 variales, overlap")
      } else {
        fontsize <- 6
        widthsize <- 60
        lineheight <- 0.7
        message("more than 40 variables, no overlap")}} else {
          if(length(rownames(pcorrplot))>20) {
            if (max(nchar(rownames(pcorrplot[1])))>60) {
              fontsize <- 6
              widthsize <- 60
              lineheight <- 0.7
              message("20-40 variables, overlap")} else {
                fontsize <- 8
                widthsize <- 45
                lineheight <- 0.7
                message("20-40 variables, no overlap")}} else {
                  if(length(rownames(pcorrplot))>20) {
                    if (max(nchar(rownames(pcorrplot[1])))>45) {
                      fontsize <- 8
                      widthsize <- 45
                      lineheight <- 0.7
                      message("20-20 variables, overlap")} else {
                        fontsize <- 10
                        widthsize <- 30
                        lineheight <- 0.7
                        message("20-20 variables, no overlap")}} else {
                          if(length(rownames(pcorrplot))>10) {
                            if (max(nchar(rownames(pcorrplot[1])))>30) {
                              fontsize <- 10
                              widthsize <- 30
                              lineheight <- 0.7
                              message("10-20 variables, overlap")} else {
                                fontsize <- 12
                                widthsize <- 25
                                lineheight <- 0.7
                                message("10-20 variables, no overlap")
                              }} else {
                                if (max(nchar(rownames(pcorrplot[1])))>25) {
                                  fontsize <- 12
                                  widthsize <- 25
                                  lineheight <- 1
                                  message("<10 variables, overlap") } else {
                                    fontsize <- 15
                                    widthsize <- 20
                                    lineheight <- 1
                                    message("<10 variables, no overlap")
                                  }
                              }}}}
    pB1 <- ggplot(pcorrplot, aes(x=reorder(row.names(pcorrplot),-pcorrlistaftervs),y=pcorrlistaftervs))
    pB2 <- pB1 +geom_col()
    pB3 <- pB2 + theme(axis.text.x = element_text(angle = 90, size=fontsize, lineheight=lineheight,hjust=1,vjust=0.5), text=element_text(size=15), axis.text=element_text(size=15))
    pB4 <- pB3 + labs(y="p(corr)", x=element_blank(),title="P(corr) plot")
    pB5 <- pB4 + scale_x_discrete(labels = function(x) str_wrap(x, width = widthsize))
    pB6 <- pB5 + ylim(-1,1)
    pB7<-pB6 +
      geom_errorbar(
        aes(x=reorder(row.names(pcorrplot),-pcorrlistaftervs),
            ymin = conf.int.low,
            ymax = conf.int.high),
        color = "red"
      )
    if (nrow(pcorrplot)>50) {pB8 <- pB7 + theme(axis.text.x = element_blank())
    pB8} else {
      pB7}
  }

  plotpcorronly50variables <- function(model, variable_names_length, variable_names_position){
    pcorrplot <- model$pcorrlistaftervs
    rownames(pcorrplot) <- gsub("_", " ", rownames(pcorrplot))
    rownames(pcorrplot) <- gsub("[.]", " ", rownames(pcorrplot))

    pcorrplotdf <- tibble::rownames_to_column(pcorrplot, "VALUE")
    pcorrplotdfordered <- pcorrplotdf[order(abs(pcorrplotdf[,2]),decreasing = TRUE),]

    pcorrplot <- subset(pcorrplot, abs(pcorrplot$pcorrlistaftervs)>=abs(pcorrplotdfordered[50,2]))

    if (variable_names_position=="beginning") {
      rownames(pcorrplot) <- substr(rownames(pcorrplot), 1, variable_names_length)}
    if (variable_names_position=="end") {
      for (novariables in 1:nrow(pcorrplot)) {
        rownames(pcorrplot)[novariables] <- substr(rownames(pcorrplot)[novariables], nchar(rownames(pcorrplot)[novariables])-(variable_names_length-1), nchar(rownames(pcorrplot)[novariables]))
        blankpos = StrPos(rownames(pcorrplot)[novariables], ' ', 1);
        if (!is.na(blankpos)&nchar(rownames(pcorrplot)[novariables])>=variable_names_length) {
          rownames(pcorrplot)[novariables] <- substr(rownames(pcorrplot)[novariables], blankpos+1, nchar(rownames(pcorrplot)[novariables]));
        }}}
    if (length(rownames(pcorrplot))>40) {
      if(max(nchar(rownames(pcorrplot[1])))>=80) {
        fontsize <- 4
        widthsize <- 80
        lineheight <- 0.7
        message("more than 40 variales, overlap")
    } else {
   fontsize <- 6
    widthsize <- 60
   lineheight <- 0.7
      message("more than 40 variables, no overlap")}} else {
        if(length(rownames(pcorrplot))>30) {
          if (max(nchar(rownames(pcorrplot[1])))>60) {
            fontsize <- 6
            widthsize <- 60
            lineheight <- 0.7
            message("30-40 variables, overlap")} else {
              fontsize <- 8
              widthsize <- 45
              lineheight <- 0.7
              message("30-40 variables, no overlap")}} else {
        if(length(rownames(pcorrplot))>20) {
          if (max(nchar(rownames(pcorrplot[1])))>45) {
          fontsize <- 8
          widthsize <- 45
          lineheight <- 0.7
          message("20-30 variables, overlap")} else {
              fontsize <- 10
              widthsize <- 30
              lineheight <- 1
              message("20-30 variables, no overlap")}} else {
          if(length(rownames(pcorrplot))>10) {
            if (max(nchar(rownames(pcorrplot[1])))>30) {
              fontsize <- 10
              widthsize <- 30
              lineheight <- 1
              message("10-20 variables, overlap")} else {
  fontsize <- 12
  widthsize <- 25
  lineheight <- 1
  message("10-20 variables, no overlap")
              }} else {
                if (max(nchar(rownames(pcorrplot[1])))>25) {
                fontsize <- 12
                widthsize <- 25
                lineheight <- 1
                message("<10 variables, overlap") } else {
                  fontsize <- 15
                  widthsize <- 20
                  lineheight <- 1
                  message("<10 variables, no overlap")
                }
              }}}}



    pB1 <- ggplot(pcorrplot, aes(x=reorder(row.names(pcorrplot),-pcorrlistaftervs),y=pcorrlistaftervs))
    pB2 <- pB1 +geom_col()
    pB3 <- pB2 + theme(axis.text.x = element_text(angle = 90, size=fontsize, lineheight=lineheight,hjust=1,vjust=0.5), text=element_text(size=15), axis.text=element_text(size=15))
    pB4 <- pB3 + labs(y="p(corr)", x=element_blank(),title="P(corr) plot of 50 most contributing variables")
    pB5 <- pB4 + scale_x_discrete(labels = function(x) str_wrap(x, width = widthsize))
    pB6 <- pB5 + ylim(-1,1)
    pB7<-pB6 +
      geom_errorbar(
        aes(x=reorder(row.names(pcorrplot),-pcorrlistaftervs),
            ymin = conf.int.low,
            ymax = conf.int.high),
        color = "red")
    pB7
  }

  #library(bootstrap)
  #jackknifedloading <- jackknife(subsetdatamatrix[,row.names(loadingropls)[1]], getLoadingMN(aftervsdata.oplsda))

  plotscore <- function(model,classordered){
    scoresropls <- as.data.frame(model$scoreofvariablesaftervs)

  fontsize <- 15/4 + 15*(15/nrow(scoresropls)*3/4)
  if (fontsize>15) {fontsize <- 15}

  pC1 <-  ggplot(scoresropls, aes(x=row.names(scoresropls),y=scoresropls$p1, color=classordered))
  pC2 <- pC1 + geom_point()
  pC3 <- pC2 + theme(axis.text.x = element_text(angle = 90, size=fontsize))
  pC4 <- pC3 + labs(y="scores (t)", x="subject id",title="Predictive scores")
  pC5 <- pC4 + theme(legend.title = element_blank(), legend.justification = "top", text=element_text(size=15), axis.text=element_text(size=15))
  pC6 <- pC5 + scale_colour_manual(values=c("#0072B2","#D55E00"))
  pC7 <- pC6 + guides(col = guide_legend(reverse = TRUE))
  if (nrow(scoresropls)>60) {pC8 <- pC7 + theme(axis.text.x = element_blank())
  pC8} else {
    pC7}
  }

  plotrawdata <- function(variablename, Rdataname, dirRdata){
    load(paste(dirRdata,"/",Rdataname, sep=""))
    fontsize <- 15/4 + 15*(15/nrow(subsetdatamatrix)*3/4)
    if (fontsize>15) {fontsize <- 15}

    pC1 <-  ggplot(subsetdatamatrix, aes(x=row.names(subsetdatamatrix),y=subsetdatamatrix[,variablename], color=subsetsampleID[,paste(colname_groupID)]))
    pC2 <- pC1 + geom_point()
    pC3 <- pC2 + theme(axis.text.x = element_text(angle = 90, size=fontsize))
    pC4 <- pC3 + labs(y=variablename, x="subject id",title="Raw data plot")
    pC5 <- pC4 + theme(legend.title = element_blank(), legend.justification = "top", text=element_text(size=15), axis.text=element_text(size=15))
    pC6 <- pC5 + scale_colour_manual(values=c("#0072B2","#D55E00"))
    pC7 <- pC6 + guides(col = guide_legend(reverse = TRUE))
    pC7
  }


  plotboxplot <- function(model,subsetsampleID,colname_groupID,classordered){
    scoresropls <- as.data.frame(model$scoreofvariablesaftervs)
    df <- cbind((subsetsampleID[,paste(colname_groupID)]),as.data.frame(scoresropls$p1))
    colnames(df) <- c("sampleID","scores")
    df$sampleID <- as.character(df$sampleID)
    pwc <- df %>% t_test(scores ~ sampleID)
    pwc <- pwc %>% add_xy_position(x = "sampleID", step.increase = 0.2)
    pwc <- pwc %>% add_y_position(step.increase = 0.2)
    pwc$p<-formatC(signif(pwc$p,digits=1), digits=1,format="g")

    fontsize <- 15

    pC1 <- ggboxplot(df, x="sampleID" ,y="scores", fill="sampleID", order=levels(classordered))
    pC2 <- pC1 + theme(axis.text.x = element_text(size=fontsize))
    pC3 <- pC2 + labs(y="scores (t)", x=NULL,title="Predictive scores")
    pC4 <- pC3 + theme(legend.position = "none",  text=element_text(size=15), axis.text=element_text(size=15))
    pC5 <- pC4 + stat_pvalue_manual(pwc, label = "p = {p}", hide.ns = F, label.size = 4, bracket.size = 0.5)
    pC6 <- pC5 + geom_jitter(aes(fill=subsetsampleID[,paste(colname_groupID)]),width=0.2,shape=21, color="black",size=2, alpha=0.9)
    pC7 <- pC6 + scale_fill_manual(values=c("#0072B2","#D55E00"))
    pC7
  }

  ##  create iterationtable                                                   ####

  iterationtable <- function(subsetdatamatrix,ortho_pre_vs,ortho_post_vs,class, pcorrselectionsvector, no_permutations_post_vs, variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
  iterationtable <- data.frame()
  iterationpcorrlist <- list()
  iterationsummary <- list()
  iterationmatrix <- subsetdatamatrix
  iterationno <- 0
    for (i in pcorrselectionsvector){
      iterationno <- iterationno +1
      iterationi <- opls_model_with_variable_selection_trycatch(iterationmatrix, ortho_pre_vs=NA,ortho_post_vs=NA,class, pcorr=i, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)
      if (is.na(iterationi$resultaftervs$`no. variables`)){
        iterationi <- opls_model_with_variable_selection_trycatch(iterationmatrix, ortho_pre_vs=ortho_pre_vs,ortho_post_vs=NA,class, pcorr=i, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)}
      if (is.na(iterationi$resultaftervs$`no. variables`)){
        iterationi <- opls_model_with_variable_selection_trycatch(iterationmatrix, ortho_pre_vs=ortho_pre_vs,ortho_post_vs=ortho_post_vs,class, pcorr=i, no_permutations_post_vs=no_permutations_post_vs, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)}

      iterationmatrix <- iterationmatrix[,rownames(iterationi$loadingroplsaftervs)]

    iterationtable <- rbind(iterationtable, iterationi$resultaftervs)
    iterationpcorrlist[[iterationno]] <- iterationi$pcorrlistaftervs
    }
  row.names(iterationtable) <- c(paste("model",1:nrow(iterationtable)))
 iterationsummary[[1]] <- iterationtable
 iterationsummary[[2]] <- iterationpcorrlist
 iterationsummary
  }

  #   ____________________________________________________________________________
  #   randomize groupsettings of subsetsampleID                                                                         ####


  randomize_group <- function(subsetsampleID,colname_groupID){
    randomsubsetsampleID <- subsetsampleID
    nrd <- nrow(randomsubsetsampleID)
    randomsubsetsampleID[,paste(colname_groupID)] <- droplevels(randomsubsetsampleID[,paste(colname_groupID)])
    randomgroup <- sample(randomsubsetsampleID[,paste(colname_groupID)],nrd,replace=F)
    randomgroup
  }




  ##  ............................................................................
  ##  run permutated models                                                   ####

  permoplsmodelwithvariableselection <- function(subsetsampleID,colname_groupID,subsetdatamatrix,ortho_pre_vs,ortho_post_vs,class, pcorr, variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
    randomgroup <- randomize_group(subsetsampleID,colname_groupID)
    resultpermutation <-  opls_model_with_variable_selection_trycatch(subsetdatamatrix,ortho_pre_vs=ortho_pre_vs,ortho_post_vs=ortho_post_vs,class=randomgroup, pcorr,no_permutations_post_vs=0, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)

    corrcoff <- cor(as.numeric(as.factor(randomgroup)),as.numeric(as.factor(subsetsampleID[,paste(colname_groupID)])), method="pearson")
    resultpermutation <- cbind(resultpermutation$resultaftervs,corrcoff)
    if ("pR2Y permutated post v.s." %in% names(resultpermutation)) {resultpermutation$"pR2Y permutated post v.s." <- NULL}
    if ("pQ2 permutated post v.s." %in% names(resultpermutation)) {resultpermutation$"pQ2 permutated post v.s." <- NULL}
    resultpermutation

  }


  ##  ............................................................................
  ##  create table of randomized model                                        ####
  table_of_randomised_models_over_vs <- function(subsetsampleID,colname_groupID,subsetdatamatrix,ortho_pre_vs,ortho_post_vs,class, pcorr, no_permutations_over_vs, variable_selection_using_VIP,
                                                 max_no_of_ortho_pre_vs,
                                                 max_no_of_ortho_post_vs){
  permutated_models <- data.frame(matrix(NA,nrow=no_permutations_over_vs,ncol=10))
  for (i in 1:no_permutations_over_vs){
    permutated_models[i,] <- permoplsmodelwithvariableselection(subsetsampleID,colname_groupID,subsetdatamatrix,ortho_pre_vs=ortho_pre_vs,ortho_post_vs=ortho_post_vs,class, pcorr, variable_selection_using_VIP=variable_selection_using_VIP,
                                                                max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs,
                                                                max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)
  sinkout()
	}
  row.names(permutated_models) <- c(paste("model",1:nrow(permutated_models)))
  colnames(permutated_models) <- c("pcorr cutoff","ortho pre v.s.","no. variables","R2X(cum)","R2Y(cum)","Q2(cum)","RMSEE","pred. post v.s.","ortho post v.s.","corrcoff")
  permutated_models
  }


  ##  ............................................................................
  ##  calculate percent of R2 and Q2 in permuted model larger than original mo####

  percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated <- function(permutated_models, R2optimized, Q2optimized){


  percentageofpermmorethennotpermQ2 <- nrow(subset(permutated_models,permutated_models$"Q2(cum)">Q2optimized))/nrow(permutated_models)*100

  percentageofpermmorethennotpermR2 <- nrow(subset(permutated_models,permutated_models$"R2Y(cum)">R2optimized))/nrow(permutated_models)*100

  percentageofpermmorethennotperm <- rbind(percentageofpermmorethennotpermR2,percentageofpermmorethennotpermQ2)
  row.names(percentageofpermmorethennotperm) <- c("percentage of R2 in permutated models larger than R2 in unpermutated model=","percentage of Q2 in permutated models larger than Q2 in unpermutated model=")
  percentageofpermmorethennotperm

  }

  ##  calculate percent of R2 and Q2 in permuted model larger than original mo####

  pvalue_R2_and_Q2_in_perm <- function(permutated_models, R2optimized, Q2optimized){


    pvalueofpermQ2 <- nrow(subset(permutated_models,permutated_models$"Q2(cum)">Q2optimized))/nrow(permutated_models)
    if (pvalueofpermQ2==0) {pvalueofpermQ2<-1/nrow(permutated_models)}

    pvalueofpermR2 <- nrow(subset(permutated_models,permutated_models$"R2Y(cum)">R2optimized))/nrow(permutated_models)
    if (pvalueofpermR2==0) {pvalueofpermR2<-1/nrow(permutated_models)}

    pvalueofperm <- rbind(pvalueofpermR2,pvalueofpermQ2)
    row.names(pvalueofperm) <- c("gives pR2Y=","gives pQ2=")
    pvalueofperm

  }

  plotpermutations <- function(permutated_models){
  permutated_models$corrcoff <- abs(permutated_models$corrcoff)
  plot(permutated_models$corrcoff, permutated_models$`R2Y(cum)`, xlab="correlation coefficient between Y and permutated Y", ylab="Q2(cum) of permutated model")
  plot(permutated_models$corrcoff, permutated_models$`Q2(cum)`, xlab="correlation coefficient between Y and permutated Y", ylab="R2Y(cum) of permutated model")
  }

  calculatepforpermutation <- function(permutated_models, unpermutatedmodel){
    corrcoff <- 1
    resultunpermutatedmodel <- cbind(unpermutatedmodel$resultaftervs,corrcoff)
    resultunpermutatedmodel$"pR2Y permutated post v.s." <- NULL
    resultunpermutatedmodel$"pQ2 permutated post v.s." <- NULL
    permutated_modelsplusunpermutated <-  rbind(resultunpermutatedmodel , permutated_models)

    permutated_modelsplusunpermutated$corrcoff <- abs(permutated_modelsplusunpermutated$corrcoff)
    permutated_modelsplusunpermutated$unperm <- permutated_modelsplusunpermutated$corrcoff==1

    permutated_modelsplusunpermutated$`R2Ycum` <- permutated_modelsplusunpermutated$`R2Y(cum)`
    permutated_modelsplusunpermutated$`R2Ycum`[is.na(permutated_modelsplusunpermutated$`R2Ycum`)] <- 0

    permutated_modelsplusunpermutated$`Q2cum` <- permutated_modelsplusunpermutated$`Q2(cum)`
    permutated_modelsplusunpermutated$`Q2cum`[is.na(permutated_modelsplusunpermutated$`Q2cum`)] <- 0

    pcorrtestpermutationR2Y <- cor.test(permutated_modelsplusunpermutated$`R2Ycum`,permutated_modelsplusunpermutated$corrcoff,method="pearson")
    pcorrtestpermutationQ2 <- cor.test(permutated_modelsplusunpermutated$`Q2cum`,permutated_modelsplusunpermutated$corrcoff,method="pearson")
    pcorrtestpermutation <- list(pcorrtestpermutationR2Y,pcorrtestpermutationQ2)
    names(pcorrtestpermutation) <- c("R2Y","Q2")
    pforpermutationtable <- data.frame(nrow=5)
    reg_of_perm <- lm(permutated_modelsplusunpermutated$`Q2cum`~permutated_modelsplusunpermutated$corrcoff)
    pforpermutationtable[1] <- reg_of_perm$coefficients[1]
    pforpermutationtable[2] <- pcorrtestpermutation$R2Y$estimate
    pforpermutationtable[3] <- pcorrtestpermutation$Q2$estimate
    pforpermutationtable[4] <- pcorrtestpermutation$R2Y$p.value
    pforpermutationtable[5] <- pcorrtestpermutation$Q2$p.value


    names(pforpermutationtable) <-c("intercept of permutation","Correlation between R2Y for permutations over variable selection and correlation between permutated and unpermutated","Correlation between Q2 for permutations over variable selection. and correlation between permutated and unpermutated","P-value for correlation between R2Y for permutations over variable selection and correlation between permutated and unpermutated","P-value for correlation between Q2 for permutations over variable selection and correlation between permutated and unpermutated")
    t(pforpermutationtable)
    }

  plotpermutationswithoriginalmodelandreg <- function(permutated_models, subsetdatamatrix, ortho_pre_vs, ortho_post_vs,class,pcorr, variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
    corrcoff <- 1
    resultunpermutatedmodel <- cbind(opls_model_with_variable_selection_trycatch(subsetdatamatrix, ortho_pre_vs, ortho_post_vs,class,pcorr=pcorr, no_permutations_post_vs=0, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)$resultaftervs,corrcoff)
    colnames(permutated_models) <- colnames(resultunpermutatedmodel)
    permutated_modelsplusunpermutated <-  rbind(permutated_models, resultunpermutatedmodel)

    permutated_modelsplusunpermutated$corrcoff <- abs(permutated_modelsplusunpermutated$corrcoff)
    permutated_modelsplusunpermutated$unperm <- permutated_modelsplusunpermutated$corrcoff==1

    permutated_modelsplusunpermutated$`R2Ycum` <- permutated_modelsplusunpermutated$`R2Y(cum)`
    permutated_modelsplusunpermutated$`R2Ycum`[is.na(permutated_modelsplusunpermutated$`R2Ycum`)] <- 0

    pC1 <-  ggplot(permutated_modelsplusunpermutated, aes(x=corrcoff,y=`R2Ycum`, color=corrcoff==1))
    pC2 <- pC1 + geom_point()
    pC3 <- pC2 + labs(y="R2Y(cum) post variable selection", x=NULL,title="R2Y(cum)")
    pC4 <- pC3 + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
    pC5 <- pC4 + stat_cor(method = "pearson", color="black", size=5,label.y.npc="bottom",label.x.npc=0.5)
    pC6 <- pC5 + theme(legend.position = "top",text=element_text(size=15), axis.text=element_text(size=15))
    pC7 <- pC6 + scale_color_manual(values=c("red","blue"),labels = c("permutated","unpermutated"), name=NULL)

    permutated_modelsplusunpermutated$`Q2cum` <- permutated_modelsplusunpermutated$`Q2(cum)`
    permutated_modelsplusunpermutated$`Q2cum`[is.na(permutated_modelsplusunpermutated$`Q2cum`)] <- 0

    pC8 <-  ggplot(permutated_modelsplusunpermutated, aes(x=corrcoff,y=`Q2cum`, color=corrcoff==1))
    pC9 <- pC8 + geom_point()
    pC10 <- pC9 + labs(y="Q2(cum) post variable selection", x=NULL,title="Q2(cum)")
    pC11 <- pC10 + geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)
    pC12 <- pC11 + stat_cor(method = "pearson", color="black", size=5,label.y.npc="bottom",label.x.npc=0.5)
    pC13 <- pC12 + theme(legend.position = "top",text=element_text(size=15), axis.text=element_text(size=15))
    pC14 <- pC13 + scale_color_manual(values=c("red","blue"),labels = c("permutated","unpermutated"), name=NULL)

    grid.arrange(pC7, pC14, nrow = 1,top=text_grob("Permutation over variable selection", size = 20),bottom=text_grob(
      "correlation coefficient between Y and permutated Y", size=15))

}

  plotpermutationswithoriginalmodelandregusingggscatter <- function(permutated_models, subsetdatamatrix, ortho_pre_vs, orhoIoptimized,class,pcorr,variable_selection_using_VIP, max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs){
    corrcoff <- 1
    resultunpermutatedmodel <- cbind(opls_model_with_variable_selection_trycatch(subsetdatamatrix, ortho_pre_vs, ortho_post_vs,class,pcorr=pcorr,no_permutations_post_vs=0, variable_selection_using_VIP=variable_selection_using_VIP, max_no_of_ortho_pre_vs=max_no_of_ortho_pre_vs, max_no_of_ortho_post_vs=max_no_of_ortho_post_vs)$resultaftervs,corrcoff)
    permutated_modelsplusunpermutated <-  rbind(resultunpermutatedmodel , permutated_models)

    permutated_modelsplusunpermutated$corrcoff <- abs(permutated_modelsplusunpermutated$corrcoff)
    plot(permutated_modelsplusunpermutated$corrcoff, permutated_modelsplusunpermutated$`R2Y(cum)`, xlab="correlation coefficient between Y and permutated Y", ylab="Q2(cum) of permutated model",col = "#56B4E9"[permutated_modelsplusunpermutated$corrcoff==1])



permutated_modelsplusunpermutated$`R2Ycum` <- permutated_modelsplusunpermutated$`R2Y(cum)`
permutated_modelsplusunpermutated$`R2Ycum`[is.na(permutated_modelsplusunpermutated$`R2Ycum`)] <- 0
permutated_modelsplusunpermutated$unperm <- permutated_modelsplusunpermutated$corrcoff==1
ggscatter(as.data.frame(permutated_modelsplusunpermutated), x = "corrcoff", y = "R2Ycum",
          add = "reg.line", conf.int = F,
          xlab = "correlation coefficient between Y and permutated Y",
          ylab = "R2Y(cum) of permutated model",
          title = "Permutation plot", fill="unperm", color = "unperm", palette = c("blue", "Red"),
          all=TRUE, font.x=20, font.y=20, font.title=20, cor.coef.size = 8) + font("xy.text", size=20)+
  stat_cor(method = "pearson")  # Add correlation coefficient

plot(permutated_models$corrcoff, permutated_models$`Q2(cum)`, xlab="correlation coefficient between Y and permutated Y", ylab="R2Y(cum) of permutated model")
  }




  #Create summary table of models
  create_table_of_models <- function(resultmodelname, permutated_modelsname, percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name, directory_output_reports,projectname,date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    summarymodeltable <- data.frame()
    for (i in 1:length(selectedRows)) {
      load(selectedRows[i])
      rm(list=lsf.str())
resultmodel <- get(resultmodelname)
permutated_models <- get(permutated_modelsname)
percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated <- get(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name)
rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated==0] <-100/no_permutations_over_vs
pforpermutation <- calculatepforpermutation(permutated_models, unpermutatedmodel = resultmodel)

summarymodeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,resultmodel$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated)/100,t(pforpermutation))
      summarymodeltable <- rbind(summarymodeltable,summarymodeltablei)
    }
    colnames(summarymodeltable)[c(6,7,8,14,15,16)] <- c("pcorr cutoff","ortho pre v.s.","no. variables","ortho post v.s.","pR2Y permutated post v.s.","pQ2 permutated post v.s.")
    summarymodeltable
  }

  create_table_of_models_iterations <- function(directory_output_reports,projectname,date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    summarymodeltable <- data.frame()
    for (i in 1:length(selectedRows)) {
      load(selectedRows[i])
      rm(list=lsf.str())
     summarymodeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,best_performing_iteration_model_few_variables)
      summarymodeltable <- rbind(summarymodeltable,summarymodeltablei)
    }
    summarymodeltable$`diff R2Y(cum)-Q2(cum)` <-NULL
    summarymodeltable
  }


  load_model_result_with_perm <- function(directory_output_reports=directory_output_reports, projectname=projectname, date_of_analysis=date_of_analysis, firstgroup=firstgroup, secondgroup=secondgroup, secID=secID, resultmodelname=resultmodelname){

    load(paste(paste(directory_output_reports, "/",projectname, sep=""), date_of_analysis, firstgroup, "vs", secondgroup, secID,".Rdata", sep="_"))

    rm(list=lsf.str())
  resultmodel <- get(resultmodelname)
  if (resultmodelname=="result_Model3") {
    percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name<-"percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model3";
    permutated_modelsname<-"permutated_models_Model3"}
  if (resultmodelname=="result_Model2") {
    percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name<-"percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2";
    permutated_modelsname<-"permutated_models_Model2"}
  if (resultmodelname=="result_Model4") {
    percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name<-"percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model4";
    permutated_modelsname<-"permutated_models_Model4"}

  permutated_models <- get(permutated_modelsname)
  percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated <- get(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_name)
  rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
  percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated==0] <-100/no_permutations_over_vs
  pforpermutation <- calculatepforpermutation(permutated_models, unpermutatedmodel = resultmodel)

  modelresult_with_perm <- cbind(group1,ngroup1,group2,ngroup2,secID,resultmodel$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated)/100,t(pforpermutation))
   colnames(modelresult_with_perm)[c(6,7,8,14,15,16)] <- c("pcorr cutoff","ortho pre v.s.","no. variables","ortho post v.s.","pR2Y permutated post v.s.","pQ2 permutated post v.s.")
  modelresult_with_perm
  }

  #Create summary table SIMCA models
  create_table_SIMCA_model <- function(){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]

    summarySIMCAmodeltable <- data.frame()
    for (i in 1:length(selectedRows)) {
      load(selectedRows[i])
      summarySIMCAmodeltablei <- cbind(group1, " vs ",group2,secID,SIMCAmodelwithvs$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated(permutated_models,SIMCAmodelwithvs$resultaftervs$`R2Y(cum)` ,SIMCAmodelwithvs$resultaftervs$`Q2(cum)`)))

      summarySIMCAmodeltable <- rbind(summarySIMCAmodeltable,summarySIMCAmodeltablei)
    }
    summarySIMCAmodeltable
  }

 # create summary table for each comparison
  create_table_each_comparison <- function(directory_output_reports, projectname, date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    table_of_each_comparison <- list()
    for (i in 1:length(selectedRows)) {
      summaryeachcomparisontable <- data.frame()
      load(selectedRows[i])
      summarySIMCAmodeltablei <- cbind(group1, " vs ",group2,secID,SIMCAmodelwithvs$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_SIMCA))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summarySIMCAmodeltablei)

      summaryModel2tablei <- cbind(group1, " vs ",group2,secID,result_Model2$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryModel2tablei)

      summaryadditionallyinvestigatedmodeltablei <- cbind(group1, " vs ",group2,secID,resultmodeladditionallyinvestigated$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryadditionallyinvestigatedmodeltablei)
    table_of_each_comparison[[i]] <- list(summaryeachcomparisontable)
    }
    table_of_each_comparison
  }

  # create summary table for each comparison
  create_table_each_comparison_pcorr_by_pvalue_and_Q2max <- function(directory_output_reports, projectname, date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    table_of_each_comparison <- list()
    for (i in 1:length(selectedRows)) {
      summaryeachcomparisontable <- data.frame()
      load(selectedRows[i])

      summaryModel2tablei <- cbind(group1, " vs ",group2,secID,result_Model2$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryModel2tablei)

      summaryadditionallyinvestigatedmodeltablei <- cbind(group1, " vs ",group2,secID,resultmodeladditionallyinvestigated$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryadditionallyinvestigatedmodeltablei)
      table_of_each_comparison[[i]] <- list(summaryeachcomparisontable)
    }
    table_of_each_comparison
  }

  create_table_each_comparison_first_invstigated_pcorr_by_pvalue_and_Q2max <- function(directory_output_reports, projectname, date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    table_of_each_comparison <- list()
    for (j in 1:length(selectedRows)) {
      summaryeachcomparisontable <- data.frame()
      load(selectedRows[j])
      rm(list=lsf.str())
      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model1, unpermutatedmodel = result_Model1)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1==0] <- 100/no_permutations_over_vs
      summary_Model1_modeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model1$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summary_Model1_modeltablei)

      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model2 , unpermutatedmodel = result_Model2)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2==0] <-100/no_permutations_over_vs
      summaryModel2tablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model2$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryModel2tablei)

      pforpermutation <- calculatepforpermutation(permutated_models=permutated_modelsadditionallyinvestigated , unpermutatedmodel = resultmodeladditionallyinvestigated)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated==0] <-100/no_permutations_over_vs
      summaryadditionallyinvestigatedmodeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,resultmodeladditionallyinvestigated$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_additionallyinvestigated)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryadditionallyinvestigatedmodeltablei)
      table_of_each_comparison[[j]] <- list(summaryeachcomparisontable)

    }
    table_of_each_comparison
  }

  create_table_each_comparison_first_invstigated_pcorr_by_pvalue_and_Q2max_with_pvalue <- function(directory_output_reports, projectname, date_of_analysis){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    table_of_each_comparison <- list()
    for (j in 1:length(selectedRows)) {
      summaryeachcomparisontable <- data.frame()
      load(selectedRows[j])
      rm(list=lsf.str())
      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model1, unpermutatedmodel = result_Model1)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1==0] <- 100/no_permutations_over_vs
      summary_Model1_modeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model1$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model1)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summary_Model1_modeltablei)

      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model2 , unpermutatedmodel = result_Model2)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2==0] <-100/no_permutations_over_vs
      summaryModel2tablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model2$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model2)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summaryModel2tablei)

      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model3 , unpermutatedmodel = result_Model3)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model3) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model3[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model3==0] <-100/no_permutations_over_vs
      summary_Model3modeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model3$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model3)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summary_Model3modeltablei)
      table_of_each_comparison[[j]] <- list(summaryeachcomparisontable)

      pforpermutation <- calculatepforpermutation(permutated_models=permutated_models_Model4 , unpermutatedmodel = result_Model4)
      rownames(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model4) <- c("pR2Y permutated over v.s.", "pQ2 permutated over v.s.")
      percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model4[percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model4==0] <-100/no_permutations_over_vs
      summary_Model4modeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,result_Model4$resultaftervs,t(percent_R2_and_Q2_in_permutated_larger_than_in_unpermutated_Model4)/100, t(pforpermutation))
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summary_Model4modeltablei)
      table_of_each_comparison[[j]] <- list(summaryeachcomparisontable)

      best_performing_iteration_model_few_variables$`diff R2Y(cum)-Q2(cum)` <- NULL
      summary_Model5modeltablei <- cbind(group1,ngroup1,group2,ngroup2,secID,best_performing_iteration_model_few_variables,NA,NA,NA,NA,NA,NA,NA)
      colnames(summary_Model5modeltablei) <- colnames(summaryeachcomparisontable)
      summaryeachcomparisontable <- rbind(summaryeachcomparisontable,summary_Model5modeltablei)
      table_of_each_comparison[[j]] <- list(summaryeachcomparisontable)

    }
    table_of_each_comparison
  }

  #Create summary table of loadings
  create_table_of_loadings_with_pcorr <- function(resultmodelname, directory_output_reports, projectname, date_of_analysis, groupsnumeric,order_of_models){
    modelfilnames <- dir(directory_output_reports)
    selectedRowsRdata <- modelfilnames[modelfilnames %like% "%.Rdata"]
    selectedRows <- selectedRowsRdata[selectedRowsRdata %like% paste(paste(projectname,date_of_analysis,sep="_"),"%",sep="")]
    summarypcorrlist <- list()
    for (j in 1:length(selectedRows)) {
      load(selectedRows[j])
      resultmodel <- get(resultmodelname)
      pcorrlist <- resultmodel$pcorrlistaftervs
      if (is.null(pcorrlist)) {pcorrlist<-data.frame(matrix(NA))}
      pcorrlist <- tibble::rownames_to_column(pcorrlist, "VALUE")
      pcorrlist <- pcorrlist[order(abs(pcorrlist[,2]),decreasing = TRUE),]

        colnames(pcorrlist)[1] <- paste(group1, " vs ",group2,secID)
        colnames(pcorrlist)[2] <- paste(group1, " vs ",group2,secID,"p(corr)")
      summarypcorrlist[[j]] <- pcorrlist
    }
    max.rows <- max(unlist(lapply(summarypcorrlist, nrow), use.names = F))
  summarypcorrlistdf <- lapply(summarypcorrlist, function(x) {
    na.count <- max.rows - nrow(x)
    if (na.count > 0L) {
      na.dm <- matrix(NA, na.count, ncol(x))
      colnames(na.dm) <- colnames(x)
      rbind(x, na.dm)
    } else {
      x
    }
  })
  summarypcorrlistdf <- do.call(cbind, summarypcorrlistdf)
  if (groupsnumeric=="yes"){
    orderlist <- vector()
    for (i in 1:length(order_of_models)) {
    orderlist[c(i*2-1)] <- c(order_of_models[i]*2-1)
    orderlist[c(i*2)] <- order_of_models[i]*2
  }
    summarypcorrlistdf <- summarypcorrlistdf[,orderlist]
  }
  summarypcorrlistdf
  }

  #Calculate p(corr) from p-value
  #from formula derived from https://stats.stackexchange.com/questions/61026/can-p-values-for-pearsons-correlation-test-be-computed-just-from-correlation-co
  calculatepcorrfrompvalue <- function(selectpvalue,n){
    pcorrfrompvaluelist <- data.frame()
    for (i in 1:1000){
      r <- i/1000
      pcorrfrompvaluelist[i,1]=r
      z <- 0.5 * log((1+r)/(1-r))
      zse <- 1/sqrt(n-3)
      pcorrfrompvaluelist[i,2] <- min(pnorm(z, sd=zse), pnorm(z, lower.tail=F, sd=zse))*2
    }
    colnames(pcorrfrompvaluelist) <- c("pcorr","pvalue")
    pcorrfrompvaluelist$pvaluediff <- pcorrfrompvaluelist$pvalue - selectpvalue
    pcorrfrompvalue <- subset(pcorrfrompvaluelist,abs(pcorrfrompvaluelist$pvaluediff) == min(abs(pcorrfrompvaluelist$pvaluediff)))
    pcorrfrompvalue$pcorr
  }

  selecting_model_with_diff_less_than_02_and_max_Q2 <- function(modeltable, prefered_pR2_and_pQ2_permutated_post_vs) {


    nrowmodeltable <- 0
    pvalueR2andQ2 <- 0
    diffR2Q2 <- 0.1
    while (nrowmodeltable==0 & (diffR2Q2<1 | pvalueR2andQ2<1)) {
      pvalueR2andQ2 <- pvalueR2andQ2+prefered_pR2_and_pQ2_permutated_post_vs
      diffR2Q2 <- diffR2Q2+0.1
      modeltable_with_diff_less_than_02_and_pvalue_less_than_given_value <-
        subset(
          modeltable,
          modeltable$"pR2Y permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"pQ2 permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"diff R2Y(cum)-Q2(cum)"<=diffR2Q2)
      nrowmodeltable <- nrow(modeltable_with_diff_less_than_02_and_pvalue_less_than_given_value)
    }
      maxQ2 <- max(modeltable_with_diff_less_than_02_and_pvalue_less_than_given_value$`Q2(cum)`)

      model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2 <-
        subset(
          modeltable_with_diff_less_than_02_and_pvalue_less_than_given_value,
          modeltable_with_diff_less_than_02_and_pvalue_less_than_given_value$`Q2(cum)`>(maxQ2-0.01*maxQ2)
        )

      model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2 <-
        subset(
          model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2,
          model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2$`ortho post v.s.`==min(model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2$`ortho post v.s.`)
        )
      model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2 <- model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2[1,]


    model_with_diff_less_than_02_and_pvalue_less_than_given_value_and_max_Q2
  }

  selecting_model_with_diff_less_than_02_max_Q2_low_pperm_and_high_pcorr_few_variables <- function(modeltable, prefered_pR2_and_pQ2_permutated_post_vs, pcorr_diff) {


    nrowmodeltable <- 0
    pvalueR2andQ2 <- 0
    diffR2Q2 <- 0.1
    while (nrowmodeltable==0 & (diffR2Q2<1 | pvalueR2andQ2<1)) {
      pvalueR2andQ2 <- pvalueR2andQ2+prefered_pR2_and_pQ2_permutated_post_vs
      diffR2Q2 <- diffR2Q2+0.1
      modeltable_selected <-
        subset(
          modeltable,
          modeltable$"pR2Y permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"pQ2 permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"diff R2Y(cum)-Q2(cum)"<=diffR2Q2)
      nrowmodeltable <- nrow(modeltable_selected)
    }

      maxQ2 <- max(modeltable_selected$`Q2(cum)`)
      modeltable_selected_reduced <- modeltable_selected

      model_selected_with_max_Q2 <-
        subset(
          modeltable_selected_reduced,
          modeltable_selected_reduced$`Q2(cum)`==max(modeltable_selected_reduced$`Q2(cum)`)
        )
      model_selected_with_max_Q2 <- model_selected_with_max_Q2[1,]



      for (i in seq(1,pcorr_diff , by=-0.001)) {
      maxQ2_at_min_pcorr <- 0

      while (maxQ2_at_min_pcorr!=maxQ2) {

      maxQ2 <- max(model_selected_with_max_Q2$`Q2(cum)`)

      maxQ2diff <- (modeltable_selected$`pcorr cutoff`-model_selected_with_max_Q2$`pcorr cutoff`)/i*0.01*maxQ2

      model_selected_with_max_Q2 <-
        subset(
          modeltable_selected,
          modeltable_selected$`Q2(cum)`>=(maxQ2-maxQ2diff)
        )
      model_selected_with_max_Q2 <- model_selected_with_max_Q2[1,]
      maxQ2_at_min_pcorr <- model_selected_with_max_Q2$`Q2(cum)`
      }
      }


    model_selected_with_max_Q2
  }

  selecting_model_with_diff_less_than_02_max_Q2_low_pperm_and_high_pcorr_many_variables <- function(modeltable, prefered_pR2_and_pQ2_permutated_post_vs) {


    nrowmodeltable <- 0
    pvalueR2andQ2 <- 0
    diffR2Q2 <- 0.1
    while (nrowmodeltable==0 & (diffR2Q2<1 | pvalueR2andQ2<1)) {
      pvalueR2andQ2 <- pvalueR2andQ2+prefered_pR2_and_pQ2_permutated_post_vs
      diffR2Q2 <- diffR2Q2+0.1
      modeltable_selected <-
        subset(
          modeltable,
          modeltable$"pR2Y permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"pQ2 permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"diff R2Y(cum)-Q2(cum)"<=diffR2Q2)

      maxQ2 <- max(modeltable_selected$`Q2(cum)`)
      model_selected_with_max_Q2 <- modeltable_selected
      maxQ2_at_min_pcorr <- 0
      while (maxQ2_at_min_pcorr!=maxQ2) {
        maxQ2 <- max(model_selected_with_max_Q2$`Q2(cum)`)
        model_selected_with_max_Q2 <-
          subset(
            modeltable_selected,
            modeltable_selected$`Q2(cum)`==max(modeltable_selected$`Q2(cum)`)
          )

        model_selected_with_max_Q2 <- model_selected_with_max_Q2[1,]
        maxQ2_at_min_pcorr <- model_selected_with_max_Q2$`Q2(cum)`
        maxQ2diff <- (modeltable_selected$`pcorr cutoff`-model_selected_with_max_Q2$`pcorr cutoff`)/0.05*0.01

        model_selected_with_max_Q2 <-
          subset(
            modeltable_selected,
            modeltable_selected$`Q2(cum)`>=(maxQ2_at_min_pcorr-maxQ2diff)
          )
        model_selected_with_max_Q2 <- model_selected_with_max_Q2[1,]
        maxQ2_at_min_pcorr <- max(model_selected_with_max_Q2$`Q2(cum)`)
      }
      nrowmodeltable <- nrow(model_selected_with_max_Q2)
    }
    model_selected_with_max_Q2
  }



  selecting_model_with_diff_less_than_02_max_Q2_low_pperm_and_few_orthogonals_old <- function(modeltable, prefered_pR2_and_pQ2_permutated_post_vs) {


  nrowmodeltable <- 0
  pvalueR2andQ2 <- 0
  diffR2Q2 <- 0.1
  while (nrowmodeltable==0 & (diffR2Q2<1 | pvalueR2andQ2<1)) {
    pvalueR2andQ2 <- pvalueR2andQ2+prefered_pR2_and_pQ2_permutated_post_vs
    diffR2Q2 <- diffR2Q2+0.1
  modeltable_selected <-
    subset(
      modeltable,
      modeltable$"pR2Y permutated post v.s."<=pvalueR2andQ2 &
        modeltable$"pQ2 permutated post v.s."<=pvalueR2andQ2 &
        modeltable$"diff R2Y(cum)-Q2(cum)"<=diffR2Q2)

  maxQ2 <- max(modeltable_selected$`Q2(cum)`)
  model_selected_and_max_Q2 <- modeltable_selected
  maxQ2_at_min_tot_no_ortho <- 0
  while (maxQ2_at_min_tot_no_ortho!=maxQ2) {
    maxQ2 <- max(model_selected_and_max_Q2$`Q2(cum)`)
  model_selected_and_max_Q2 <-
    subset(
      modeltable_selected,
      modeltable_selected$`Q2(cum)`>(maxQ2-0.01)
    )
  totalnoortho_in_modeltable_selected <- model_selected_and_max_Q2$`ortho pre v.s.`+ model_selected_and_max_Q2$`ortho post v.s.`
  model_selected_and_max_Q2 <-
    subset(
      model_selected_and_max_Q2,
      totalnoortho_in_modeltable_selected==min(totalnoortho_in_modeltable_selected)
    )

  model_selected_and_max_Q2 <-
    subset(
      model_selected_and_max_Q2,
      model_selected_and_max_Q2$`ortho post v.s.`==min(model_selected_and_max_Q2$`ortho post v.s.`)
    )
  maxQ2_at_min_tot_no_ortho <- max(model_selected_and_max_Q2$`Q2(cum)`)
  }


  totalnoortho <- modeltable_selected$`ortho pre v.s.`+
    modeltable_selected$`ortho post v.s.`

  maxQ2diff <- (min(totalnoortho_in_modeltable_selected) - totalnoortho)*0.01

    model_selected_and_max_Q2 <-
    subset(
      modeltable_selected,
      modeltable_selected$`Q2(cum)`>=(maxQ2-maxQ2diff)
  )
   totalnoortho <- model_selected_and_max_Q2$`ortho pre v.s.`+
      model_selected_and_max_Q2$`ortho post v.s.`
    model_selected_and_max_Q2 <-
      subset(
        model_selected_and_max_Q2,
        totalnoortho==min(totalnoortho)
      )

    model_selected_and_max_Q2 <-
    subset(
      model_selected_and_max_Q2,
      model_selected_and_max_Q2$`ortho post v.s.`==min(model_selected_and_max_Q2$`ortho post v.s.`)
    )
    model_selected_and_max_Q2 <- model_selected_and_max_Q2[1,]
  nrowmodeltable <- nrow(model_selected_and_max_Q2)
  }
  model_selected_and_max_Q2
  }

  selecting_model_with_diff_less_than_02_max_Q2_low_pperm_and_few_orthogonals <- function(modeltable, prefered_pR2_and_pQ2_permutated_post_vs) {

    nrowmodeltable <- 0
    pvalueR2andQ2 <- 0
    diffR2Q2 <- 0.1
    while (nrowmodeltable==0 & (diffR2Q2<1 | pvalueR2andQ2<1)) {
      pvalueR2andQ2 <- pvalueR2andQ2+prefered_pR2_and_pQ2_permutated_post_vs
      diffR2Q2 <- diffR2Q2+0.1
      modeltable_selected <-
        subset(
          modeltable,
          modeltable$"pR2Y permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"pQ2 permutated post v.s."<=pvalueR2andQ2 &
            modeltable$"diff R2Y(cum)-Q2(cum)"<=diffR2Q2)
      nrowmodeltable <- nrow(modeltable_selected)
    }

      maxQ2 <- max(modeltable_selected$`Q2(cum)`)
      model_selected_and_max_Q2 <- subset(modeltable_selected, modeltable_selected$`Q2(cum)`== maxQ2)
      model_selected_and_max_Q2 <- model_selected_and_max_Q2[1,]

      for (i in seq(0.001,0.01, by=0.001)) {
      maxQ2_at_min_tot_no_ortho <- 0
      while (maxQ2_at_min_tot_no_ortho!=maxQ2) {
        maxQ2 <- max(model_selected_and_max_Q2$`Q2(cum)`)

      totalnoortho <- modeltable_selected$`ortho pre v.s.`+
        modeltable_selected$`ortho post v.s.`
      totalnoortho_in_modeltable_selected <- model_selected_and_max_Q2$`ortho pre v.s.`+ model_selected_and_max_Q2$`ortho post v.s.`

      maxQ2diff <- (min(totalnoortho_in_modeltable_selected) - totalnoortho)*i*maxQ2

      model_selected_and_max_Q2 <-
        subset(
          modeltable_selected,
          modeltable_selected$`Q2(cum)`>=(maxQ2-maxQ2diff)
        )
      totalnoortho <- model_selected_and_max_Q2$`ortho pre v.s.`+
        model_selected_and_max_Q2$`ortho post v.s.`
      model_selected_and_max_Q2 <-
        subset(
          model_selected_and_max_Q2,
          totalnoortho==min(totalnoortho)
        )

      model_selected_and_max_Q2 <-
        subset(
          model_selected_and_max_Q2,
          model_selected_and_max_Q2$`ortho post v.s.`==min(model_selected_and_max_Q2$`ortho post v.s.`)
        )
      model_selected_and_max_Q2 <- model_selected_and_max_Q2[1,]
      maxQ2_at_min_tot_no_ortho <- max(model_selected_and_max_Q2$`Q2(cum)`)
      }
      }

    model_selected_and_max_Q2
  }


 exist_in_rdata <- function(variablexist, filename_rdata, directory_output_reports, directory_and_filename_rdata) {
   if (filename_rdata %in% dir(directory_output_reports)) {
     test <- load(directory_and_filename_rdata)
     if(variablexist %in% test) {
       variablexistout <- get(variablexist)
       variablexistout
     }
   }
 }

 plotSUSplot <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secIDmodel1,secIDmodel2){
   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secIDmodel1, variablename=variablenamemodel1)
   pcorrlist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secIDmodel2, variablename=variablenamemodel2)
   pcorrlist2 <- model2$pcorrlistaftervs
   pcorrlist1B <- as.data.frame(as.data.frame(model1$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   pcorrlist2B <- as.data.frame(as.data.frame(model2$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   for (i in 1:length(pcorrlist1)){
     pcorrlist1B[rownames(pcorrlist1),] <- pcorrlist1[i]
   }

   for (i in 1:length(pcorrlist2)){
     pcorrlist2B[rownames(pcorrlist2),] <- pcorrlist2[i]
   }

   SUSplot <- as.data.frame(cbind(pcorrlist1B,pcorrlist2B))
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(pcorrlist1),rownames(pcorrlist2)),rownames(SUSplot) %in% setdiff(rownames(pcorrlist2),rownames(pcorrlist1)),rownames(SUSplot) %in% intersect(rownames(pcorrlist1),rownames(pcorrlist2)))
   colorinplotvector <- vector()
   for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)} else if(colorinplot[j,2]) {paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)} else {"shared by models"}}
   pC1 <- ggplot(SUSplot, aes(x=pcorrlist1B,y=pcorrlist2B, color=colorinplotvector))
   pC2 <- pC1 + geom_point()
   pC3 <- pC2 + labs(y=paste("Cluster ", clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2), x=paste("Cluster", clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1),title=paste("SUS plot ROPLS-model cluster ", clustermodel1,":",group1model1," vs ", clustermodel1,":",group2model1,"\ncompared to cluster ", clustermodel2, ":", group1model2,"vs",clustermodel2, ":",group2model2))
   pC4 <- pC3 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
   pC5 <- pC4 + geom_text_repel(label=rownames(SUSplot))
   pC6 <- pC5 + theme(legend.title = element_blank())
   pC7 <- pC6 + scale_color_manual(values=c("blue", "darkgreen","red"))
   pC8 <- pC7 + xlim(-0.8,0.8) + ylim(-0.8,0.8)
   pC8
 }

 plotSUSplot_no_gender_stratification <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secID, variable_names_position, variable_names_length){
   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename=variablenamemodel1)
   pcorrlist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename=variablenamemodel2)
   pcorrlist2 <- model2$pcorrlistaftervs
   pcorrlist1B <- as.data.frame(as.data.frame(model1$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   pcorrlist2B <- as.data.frame(as.data.frame(model2$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   for (i in 1:length(pcorrlist1)){
     pcorrlist1B[rownames(pcorrlist1),] <- pcorrlist1[i]
   }

   for (i in 1:length(pcorrlist2)){
     pcorrlist2B[rownames(pcorrlist2),] <- pcorrlist2[i]
   }

   SUSplot <- as.data.frame(cbind(pcorrlist1B,pcorrlist2B))
   SUSplot[is.na(SUSplot)] <- 0
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(pcorrlist1),rownames(pcorrlist2)),rownames(SUSplot) %in% setdiff(rownames(pcorrlist2),rownames(pcorrlist1)),rownames(SUSplot) %in% intersect(rownames(pcorrlist1),rownames(pcorrlist2)))
   colorinplotvector <- vector()
   for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)} else if(colorinplot[j,2]) {paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)} else {"shared by models"}}
   susplotnames <- SUSplot
   if (variable_names_position=="beginning") {
     rownames(susplotnames) <- substr(rownames(susplotnames), 1, variable_names_length)}
   if (variable_names_position=="end") {
     rownames(susplotnames) <- substr(rownames(susplotnames), nchar(rownames(susplotnames))-(variable_names_length-1), nchar(rownames(susplotnames)))}

   pC1 <- ggplot(SUSplot, aes(x=pcorrlist1B,y=pcorrlist2B, color=colorinplotvector))
   pC2 <- pC1 + geom_point()
   pC3 <- pC2 + labs(y=paste("Cluster ", clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2), x=paste("Cluster", clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1),title=paste("SUS plot ROPLS-model cluster ", clustermodel1,":",group1model1," vs ", clustermodel1,":",group2model1,"\ncompared to cluster ", clustermodel2, ":", group1model2,"vs",clustermodel2, ":",group2model2))
   pC4 <- pC3 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
   pC5 <- pC4 + geom_text_repel(label=rownames(susplotnames))
   pC6 <- pC5 + theme(legend.title = element_blank())
   pC7 <- pC6 + scale_color_manual(values=c("blue", "darkgreen","red"))
   pC8 <- pC7 + scale_x_continuous(breaks = c(seq(-0.8, 0.8, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9 <- pC8 + scale_y_continuous(breaks = c(seq(-0.8, 0.8, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9
 }

 SUSplot_table_no_gender_stratification <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secID){

   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename=variablenamemodel1)
   pcorrlist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename=variablenamemodel2)
   pcorrlist2 <- model2$pcorrlistaftervs
   pcorrlist1B <- as.data.frame(as.data.frame(model1$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
  rownames(pcorrlist1B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
  colnames(pcorrlist1B) <- "pcorrlist1B"
   pcorrlist2B <- as.data.frame(as.data.frame(model2$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
  rownames(pcorrlist2B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   for (i in 1:length(pcorrlist1)){
   pcorrlist1B[rownames(pcorrlist1),] <- pcorrlist1[i]
   }

   for (i in 1:length(pcorrlist2)){
   pcorrlist2B[rownames(pcorrlist2),] <- pcorrlist2[i]
    }
 colnames(pcorrlist1B)<-paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)
 colnames(pcorrlist2B)<-paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)
  SUSplot <- as.data.frame(cbind(pcorrlist1B,pcorrlist2B))
  size <- 18
  colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(pcorrlist1),rownames(pcorrlist2)),rownames(SUSplot) %in% setdiff(rownames(pcorrlist2),rownames(pcorrlist1)),rownames(SUSplot) %in% intersect(rownames(pcorrlist1),rownames(pcorrlist2)))

   SUSplotordered <- SUSplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
  colorinplotordered<- colorinplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
  setcolor <- vector()
  setcolor[colorinplotordered[,1]] <- "lightblue"
  setcolor[colorinplotordered[,2]] <- "lightgreen"
  setcolor[colorinplotordered[,3]] <- "orange"
  kbl(SUSplotordered) %>%
    kable_styling(bootstrap_options = "condensed", full_width = T) %>%
    column_spec(1, background= setcolor)
 }

 SUSplot_table_no_gender_stratification_new_model <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secID){

   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename=variablenamemodel1)
   variablelist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename=variablenamemodel2)
   variablelist2 <- model2$pcorrlistaftervs
   subsetdatamatrix1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename="subsetdatamatrix")
   class1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename="class")
   subsetdatamatrix2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename="subsetdatamatrix")
   class2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename="class")


   #B is pcorr before variable selection variables from both models
   pcorrlist1B <- as.data.frame(as.data.frame(model1$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   choosecolumn <- rownames(pcorrlist1B) %in% intersect(row.names(pcorrlist1B), colnames(subsetdatamatrix1))
   pcorrlist1BE <- subset(pcorrlist1B, choosecolumn)# remove variables filtered away from subsetmatrix

   pcorrlist2B <- as.data.frame(as.data.frame(model2$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   choosecolumn <- rownames(pcorrlist2B) %in% intersect(row.names(pcorrlist2B), colnames(subsetdatamatrix2))
   pcorrlist2BE <- subset(pcorrlist2B, choosecolumn)# remove variables filtered away from subsetmatrix

   subsetdatamatrix1B <-subsetdatamatrix1[,row.names(pcorrlist1BE)]
   oplsda1C = opls(subsetdatamatrix1B, class1, predI = 1, orthoI = model1$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

   subsetdatamatrix2B <-subsetdatamatrix2[,row.names(pcorrlist2BE)]
   oplsda2C = opls(subsetdatamatrix2B, class2, predI = 1, orthoI = model2$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


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
   oplsda1A = opls(subsetdatamatrix1A, class1, predI = 1, orthoI = model1$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

   subsetdatamatrix2A <-subsetdatamatrix2[,row.names(variablelist2)]
   oplsda2A = opls(subsetdatamatrix2A, class2, predI = 1, orthoI = model2$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


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
   colnames(pcorrlist1C)<-paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)
   colnames(pcorrlist2C)<-paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)
   SUSplot <- as.data.frame(cbind(pcorrlist1C,pcorrlist2C))
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(variablelist1),rownames(variablelist2)),rownames(SUSplot) %in% setdiff(rownames(variablelist2),rownames(variablelist1)),rownames(SUSplot) %in% intersect(rownames(variablelist1),rownames(variablelist2)))

   SUSplotordered <- SUSplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   colorinplotordered<- colorinplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   setcolor <- vector()
   setcolor[colorinplotordered[,1]] <- "lightblue"
   setcolor[colorinplotordered[,2]] <- "lightgreen"
   setcolor[colorinplotordered[,3]] <- "orange"

   printSUSplot_table <- kbl(SUSplotordered) %>%
     kable_styling(bootstrap_options = "condensed", full_width = T) %>%
     column_spec(1, background= setcolor)
   printSUSplot_table
   }

 plotSUSplot_no_gender_stratification_new_model <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secID, variable_names_position, variable_names_length){
   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename=variablenamemodel1)
   variablelist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename=variablenamemodel2)
   variablelist2 <- model2$pcorrlistaftervs
   subsetdatamatrix1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename="subsetdatamatrix")
   class1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename="class")
   subsetdatamatrix2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename="subsetdatamatrix")
   class2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename="class")


   #B is pcorr before variable selection variables from both models
   pcorrlist1B <- as.data.frame(as.data.frame(model1$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   choosecolumn <- rownames(pcorrlist1B) %in% intersect(row.names(pcorrlist1B), colnames(subsetdatamatrix1))
   pcorrlist1BE <- subset(pcorrlist1B, choosecolumn)# remove variables filtered away from subsetmatrix

   pcorrlist2B <- as.data.frame(as.data.frame(model2$pcorrlistaftervs)[unique(c(rownames(variablelist1),rownames(variablelist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(variablelist1),rownames(variablelist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   choosecolumn <- rownames(pcorrlist2B) %in% intersect(row.names(pcorrlist2B), colnames(subsetdatamatrix2))
   pcorrlist2BE <- subset(pcorrlist2B, choosecolumn)# remove variables filtered away from subsetmatrix

   subsetdatamatrix1B <-subsetdatamatrix1[,row.names(pcorrlist1BE)]
   oplsda1C = opls(subsetdatamatrix1B, class1, predI = 1, orthoI = model1$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

   subsetdatamatrix2B <-subsetdatamatrix2[,row.names(pcorrlist2BE)]
   oplsda2C = opls(subsetdatamatrix2B, class2, predI = 1, orthoI = model2$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


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
   oplsda1A = opls(subsetdatamatrix1A, class1, predI = 1, orthoI = model1$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)

   subsetdatamatrix2A <-subsetdatamatrix2[,row.names(variablelist2)]
   oplsda2A = opls(subsetdatamatrix2A, class2, predI = 1, orthoI = model2$resultaftervs$`ortho post v.s.`, scaleC="standard",info.txtC="none",fig.pdfC="none",permI=0)


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
   for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)} else if(colorinplot[j,2]) {paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)} else {"shared by models"}}
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

   pC1 <- ggplot(SUSplot, aes(x=pcorrlist1C,y=pcorrlist2C, color=colorinplotvector))
   pC2 <- pC1 + geom_point()
   pC3 <- pC2 + labs(y=paste("Cluster ", clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2), x=paste("Cluster", clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1),title=paste("SUS plot ROPLS-model cluster ", clustermodel1,":",group1model1," vs ", clustermodel1,":",group2model1,"\ncompared to cluster ", clustermodel2, ":", group1model2,"vs",clustermodel2, ":",group2model2))
   pC4 <- pC3 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
   pC5 <- pC4 + geom_text_repel(label=rownames(susplotnames))
   pC6 <- pC5 + theme(legend.title = element_blank())
   pC7 <- pC6 + scale_color_manual(values=c("blue", "darkgreen","red"))
   pC8 <- pC7 + scale_x_continuous(breaks = c(seq(-1, 1, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9 <- pC8 + scale_y_continuous(breaks = c(seq(-1, 1, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9
 }


 SUSplot_table_no_gender_stratification <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, variablenamemodel1, group1model1, group2model1, projectnamemodel2, date_of_analysismodel2, variablenamemodel2, group1model2, group2model2, secID){

   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secID, variablename=variablenamemodel1)
   pcorrlist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secID, variablename=variablenamemodel2)
   pcorrlist2 <- model2$pcorrlistaftervs
   pcorrlist1B <- as.data.frame(as.data.frame(model1$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   pcorrlist2B <- as.data.frame(as.data.frame(model2$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   for (i in 1:length(pcorrlist1)){
     pcorrlist1B[rownames(pcorrlist1),] <- pcorrlist1[i]
   }

   for (i in 1:length(pcorrlist2)){
     pcorrlist2B[rownames(pcorrlist2),] <- pcorrlist2[i]
   }
   colnames(pcorrlist1B)<-paste(clustermodel1,":",group1model1," vs ",clustermodel1,":",group2model1)
   colnames(pcorrlist2B)<-paste(clustermodel2,":",group1model2," vs ",clustermodel2,":",group2model2)
   SUSplot <- as.data.frame(cbind(pcorrlist1B,pcorrlist2B))
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(pcorrlist1),rownames(pcorrlist2)),rownames(SUSplot) %in% setdiff(rownames(pcorrlist2),rownames(pcorrlist1)),rownames(SUSplot) %in% intersect(rownames(pcorrlist1),rownames(pcorrlist2)))

   SUSplotordered <- SUSplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   colorinplotordered<- colorinplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   setcolor <- vector()
   setcolor[colorinplotordered[,1]] <- "lightblue"
   setcolor[colorinplotordered[,2]] <- "lightgreen"
   setcolor[colorinplotordered[,3]] <- "orange"
   kbl(SUSplotordered) %>%
     kable_styling(bootstrap_options = "condensed", full_width = T) %>%
     column_spec(1, background= setcolor)
 }

 SUSplot_table_new_model <- function(directory_output_reports_modelX, projectname_modelX, projectname_in_plot_modelX, date_of_analysis_modelX, variablename_modelX, group1_modelX, group2_modelX, secID_modelX, result_modelX_name, directory_output_reports_modelY, projectname_modelY, projectname_in_plot_modelY, date_of_analysis_modelY, group1_modelY, group2_modelY, secID_modelY, result_modelY_name){

   modelX <- loadRData(directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename=result_modelX_name)
   variablelist1 <- modelX$pcorrlistaftervs
   modelY <- loadRData(directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename=result_modelY_name)
   variablelist2 <- modelY$pcorrlistaftervs
   subsetdatamatrix1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="subsetdatamatrix")
   class1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="class")
   subsetdatamatrix2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="subsetdatamatrix")
   class2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="class")

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
   colnames(pcorrlist1C)<-paste(projectname_in_plot_modelX," ",group1_modelX," vs ",group2_modelX,secID_modelX)
   colnames(pcorrlist2C)<-paste(projectname_in_plot_modelY," ",group1_modelY," vs ",group2_modelY,secID_modelY)
   SUSplot <- as.data.frame(cbind(pcorrlist1C,pcorrlist2C))
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(variablelist1),rownames(variablelist2)),rownames(SUSplot) %in% setdiff(rownames(variablelist2),rownames(variablelist1)),rownames(SUSplot) %in% intersect(rownames(variablelist1),rownames(variablelist2)))

   SUSplotordered <- SUSplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   colorinplotordered<- colorinplot[order(abs(SUSplot[,1]),decreasing = TRUE),]
   setcolor <- vector()
   setcolor[colorinplotordered[,1]] <- "lightblue"
   setcolor[colorinplotordered[,2]] <- "lightgreen"
   setcolor[colorinplotordered[,3]] <- "orange"

   printSUSplot_table <- kbl(SUSplotordered) %>%
     kable_styling(bootstrap_options = "condensed", full_width = T) %>%
     column_spec(1, background= setcolor)
   printSUSplot_table
 }

 plotSUSplot_new_model <- function(directory_output_reports_modelX, projectname_modelX, projectname_in_plot_modelX, date_of_analysis_modelX, variablename_modelX, group1_modelX, group2_modelX, secID_modelX, result_modelX_name, directory_output_reports_modelY, projectname_modelY, projectname_in_plot_modelY, date_of_analysis_modelY, group1_modelY, group2_modelY, secID_modelY, result_modelY_name){
   modelX <- loadRData(directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename=result_modelX_name)
   variablelist1 <- modelX$pcorrlistaftervs
   modelY <- loadRData(directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename=result_modelY_name)
   variablelist2 <- modelY$pcorrlistaftervs
   subsetdatamatrix1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="subsetdatamatrix")
   class1 <- loadRData(directory_output_reports=directory_output_reports_modelX, projectname=projectname_modelX, date_of_analysis=date_of_analysis_modelX, firstgroup=group1_modelX, secondgroup=group2_modelX, secID=secID_modelX, variablename="class")
   subsetdatamatrix2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="subsetdatamatrix")
   class2 <- loadRData(directory_output_reports=directory_output_reports_modelY, projectname=projectname_modelY, date_of_analysis=date_of_analysis_modelY, firstgroup=group1_modelY, secondgroup=group2_modelY, secID=secID_modelY, variablename="class")




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
   for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(group1_modelX," vs ",group2_modelX,secID_modelX)} else if(colorinplot[j,2]) {paste(group1_modelY," vs ",group2_modelY,secID_modelY)} else {"shared by models"}}
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

   pC1 <- ggplot(SUSplot, aes(x=pcorrlist1C,y=pcorrlist2C, color=colorinplotvector))
   pC2 <- pC1 + geom_point()
   pC3 <- pC2 + labs(y=paste(projectname_in_plot_modelY,group1_modelY,"vs",group2_modelY, secID_modelY), x=paste(projectname_in_plot_modelX, group1_modelX,"vs",group2_modelX, secID_modelX),title=paste("SUS plot ROPLS models" ,projectname_in_plot_modelX, group1_modelX,"vs", group2_modelX, secID_modelX, "\ncompared to ", projectname_in_plot_modelY, group1_modelY,"vs",group2_modelY, secID_modelY))
   pC4 <- pC3 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
   pC5 <- pC4 + geom_text_repel(label=rownames(susplotnames))
   pC6 <- pC5 + theme(legend.title = element_blank())
   pC7 <- pC6 + scale_color_manual(values=c("blue", "darkgreen","red"))
   pC8 <- pC7 + scale_x_continuous(breaks = c(seq(-1, 1, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9 <- pC8 + scale_y_continuous(breaks = c(seq(-1, 1, by=0.4)), minor_breaks = c(seq(-1, 1, by=0.2)))
   pC9
 }



 loadRData <- function(directory_output_reports, projectname, date_of_analysis, firstgroup, secondgroup, secID, variablename) {
   #loads an RData file, and returns it
   load(paste(paste(directory_output_reports, projectname, sep=""), date_of_analysis, firstgroup, "vs", secondgroup, secID,".Rdata", sep="_"))
   get(variablename)
 }

 plotSUSplotsameanalysis <- function(directory_output_reports, projectnamemodel1, date_of_analysismodel1, modelnamemodel1, group1model1, group2model1, secIDmodel1, projectnamemodel2, date_of_analysismodel2, modelnamemodel2, group1model2, group2model2, secIDmodel2){
   model1 <- loadRData(directory_output_reports, projectname=projectnamemodel1, date_of_analysis=date_of_analysismodel1, firstgroup=group1model1, secondgroup=group2model1, secID=secIDmodel1, variablename=modelnamemodel1)
   pcorrlist1 <- model1$pcorrlistaftervs
   model2 <- loadRData(directory_output_reports, projectname=projectnamemodel2, date_of_analysis=date_of_analysismodel2, firstgroup=group1model2, secondgroup=group2model2, secID=secIDmodel2, variablename=modelnamemodel2)
   pcorrlist2 <- model2$pcorrlistaftervs
   pcorrlist1B <- as.data.frame(as.data.frame(model1$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist1B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist1B) <- "pcorrlist1B"
   pcorrlist2B <- as.data.frame(as.data.frame(model2$beforevsdata.oplsda@loadingMN)[unique(c(rownames(pcorrlist1),rownames(pcorrlist2))),])
   rownames(pcorrlist2B) <- unique(c(rownames(pcorrlist1),rownames(pcorrlist2)))
   colnames(pcorrlist2B) <- "pcorrlist2B"
   for (i in 1:length(pcorrlist1)){
     pcorrlist1B[rownames(pcorrlist1),] <- pcorrlist1[i]
   }

   for (i in 1:length(pcorrlist2)){
     pcorrlist2B[rownames(pcorrlist2),] <- pcorrlist2[i]
   }

   SUSplot <- as.data.frame(cbind(pcorrlist1B,pcorrlist2B))
   size <- 18
   colorinplot <- cbind(rownames(SUSplot) %in% setdiff(rownames(pcorrlist1),rownames(pcorrlist2)),rownames(SUSplot) %in% setdiff(rownames(pcorrlist2),rownames(pcorrlist1)),rownames(SUSplot) %in% intersect(rownames(pcorrlist1),rownames(pcorrlist2)))
   colorinplotvector <- vector()
   for (j in 1:nrow(colorinplot)) {colorinplotvector[j] <- if (colorinplot[j,1]) {paste(group1model1," vs ",group2model1, "\nin ", secIDmodel1)} else if(colorinplot[j,2]) {paste(group1model2," vs ",group2model2, "\nin ", secIDmodel2)} else {"shared by models"}}
   pC1 <- ggplot(SUSplot, aes(x=pcorrlist1B,y=pcorrlist2B, color=colorinplotvector))
   pC2 <- pC1 + geom_point()
   pC3 <- pC2 + labs(y=paste(group1model2," vs ",group2model2, " in ", secIDmodel2), x=paste(group1model1," vs ",group2model1, " in ", secIDmodel1),title=paste("SUS plot ROPLS-model ", group1model1," vs ", group2model1, " in ", secIDmodel1, "\ncompared to ", group1model2,"vs",group2model2, " in ", secIDmodel2))
   pC4 <- pC3 + theme(text=element_text(size=size), axis.text=element_text(size=size), title = element_text(size=size))
   pC5 <- pC4 + geom_text_repel(label=rownames(SUSplot))
   pC6 <- pC5 + theme(legend.position = "right",legend.title = element_blank())
   pC7 <- pC6 + scale_color_manual(values=c("blue", "darkgreen","red"))
   pC8 <- pC7 + xlim(-0.8,0.8) + ylim(-0.8,0.8)
   pC8
 }

 loadRData <- function(directory_output_reports, projectname, date_of_analysis, firstgroup, secondgroup, secID, variablename) {
   #loads an RData file, and returns it
   load(paste(paste(directory_output_reports, projectname, sep=""), date_of_analysis, firstgroup, "vs", secondgroup, secID,".Rdata", sep="_"))
   get(variablename)
 }
 make_combinations <- function(x) {

   l <- length(x)
   mylist <- lapply(2:l, function(y) {
     combn(x, y, simplify = FALSE)
   })
   mylist

 }

 create_table_to_analyze <- function(file_sampleID,colname_groupID,colname_secID,pcorr_cutoff_Model1_joint_models,no_of_ortho_pre_vs_Model1_joint_models,no_of_ortho_post_vs_Model1_joint_models,pcorr_cutoff_Model1_stratified_models,no_of_ortho_pre_vs_Model1_stratified_models,no_of_ortho_post_vs_Model1_stratified_models){

   #without stratification by secondary ID
   twowaycomparison <- make_combinations(levels(as.factor(file_sampleID[,paste(colname_groupID)])))
   twowaycomparisondf <- as.data.frame(twowaycomparison[[1]])
   nrow(t(twowaycomparisondf))
   model_table_to_analyse <- data.frame(matrix(NA,nrow=nrow(t(twowaycomparisondf))))
   model_table_to_analyse$group1 <-  t(twowaycomparisondf)[,1]
   model_table_to_analyse$group2 <-  t(twowaycomparisondf)[,2]
   model_table_to_analyse$secID <- "joint"
   model_table_to_analyse$pcorr_Model1 <- pcorr_cutoff_Model1_joint_models #either enter a choosen value or vector or enter "according to p-value"
   model_table_to_analyse$ortho_pre_vs <- no_of_ortho_pre_vs_Model1_joint_models #either enter a choosen value or vector
   model_table_to_analyse$ortho_post_vs <- no_of_ortho_post_vs_Model1_joint_models #either enter a choosen value or vector

   if (colname_secID!="joint"){
     # with stratification by seconary ID
     for (i in 1:length(levels(as.factor(file_sampleID[,paste(colname_secID)])))){
       twowaycomparison <- make_combinations(levels(as.factor(file_sampleID[,paste(colname_groupID)])))
       twowaycomparisondf <- as.data.frame(twowaycomparison[[1]])
       nrow(t(twowaycomparisondf))
       model_table_to_analysei <- data.frame(matrix(NA,nrow=nrow(t(twowaycomparisondf))))
       model_table_to_analysei$group1 <-  t(twowaycomparisondf)[,1]
       model_table_to_analysei$group2 <-  t(twowaycomparisondf)[,2]
       model_table_to_analysei$secID <- levels(as.factor(file_sampleID[,paste(colname_secID)]))[i]
       model_table_to_analysei$pcorr_Model1 <- pcorr_cutoff_Model1_stratified_models #either enter a choosen value or vector or enter "according to p-value"
       model_table_to_analysei$ortho_pre_vs <- no_of_ortho_pre_vs_Model1_stratified_models #either enter a choosen value or vector
       model_table_to_analysei$ortho_post_vs <- no_of_ortho_post_vs_Model1_stratified_models #either enter a choosen value or vector
       model_table_to_analyse <- rbind(model_table_to_analyse,model_table_to_analysei)
     }
   }
   model_table_to_analyse[,1] <- NULL
   model_table_to_analyse$setseedno <- c(setseedfirstmodel+1):c(setseedfirstmodel+nrow(model_table_to_analyse))
   model_table_to_analyse
 }

 create_table_to_analyze_reordered <- function(reordered_levels_of_groups,file_sampleID,colname_groupID,colname_secID,pcorr_cutoff_Model1_joint_models,no_of_ortho_pre_vs_Model1_joint_models,no_of_ortho_post_vs_Model1_joint_models,pcorr_cutoff_Model1_stratified_models,no_of_ortho_pre_vs_Model1_stratified_models,no_of_ortho_post_vs_Model1_stratified_models){

   #without stratification by secondary ID
   twowaycomparison <- make_combinations(reordered_levels_of_groups)
   twowaycomparisondf <- as.data.frame(twowaycomparison[[1]])
   model_table_to_analyse <- data.frame(matrix(NA,nrow=nrow(t(twowaycomparisondf))))
   model_table_to_analyse$group1 <-  t(twowaycomparisondf)[,1]
   model_table_to_analyse$group2 <-  t(twowaycomparisondf)[,2]
   model_table_to_analyse$secID <- "joint"
   model_table_to_analyse$pcorr_Model1 <- pcorr_cutoff_Model1_joint_models #either enter a choosen value or vector or enter "according to p-value"
   model_table_to_analyse$ortho_pre_vs <- no_of_ortho_pre_vs_Model1_joint_models #either enter a choosen value or vector
   model_table_to_analyse$ortho_post_vs <- no_of_ortho_post_vs_Model1_joint_models #either enter a choosen value or vector

   if (colname_secID!="joint"){
     # with stratification by seconary ID
     for (i in 1:length(levels(as.factor(file_sampleID[,paste(colname_secID)])))){
       twowaycomparison <- make_combinations(reordered_levels_of_groups)
       twowaycomparisondf <- as.data.frame(twowaycomparison[[1]])
       model_table_to_analysei <- data.frame(matrix(NA,nrow=nrow(t(twowaycomparisondf))))
       model_table_to_analysei$group1 <-  t(twowaycomparisondf)[,1]
       model_table_to_analysei$group2 <-  t(twowaycomparisondf)[,2]
       model_table_to_analysei$secID <- levels(as.factor(file_sampleID[,paste(colname_secID)]))[i]
       model_table_to_analysei$pcorr_Model1 <- pcorr_cutoff_Model1_stratified_models #either enter a choosen value or vector or enter "according to p-value"
       model_table_to_analysei$ortho_pre_vs <- no_of_ortho_pre_vs_Model1_stratified_models #either enter a choosen value or vector
       model_table_to_analysei$ortho_post_vs <- no_of_ortho_post_vs_Model1_stratified_models #either enter a choosen value or vector
       model_table_to_analyse <- rbind(model_table_to_analyse,model_table_to_analysei)
     }
   }
   model_table_to_analyse[,1] <- NULL
   model_table_to_analyse$setseedno <- c(setseedfirstmodel+1):c(setseedfirstmodel+nrow(model_table_to_analyse))
   model_table_to_analyse
 }

 create_or_load_Model_table_to_analyze <- function(directory_input_matrix_sampleID, directory_model_table_to_analyse,filename_model_table_to_analyse,reordered_levels_of_groups, file_sampleID, filename_sampleID, colname_groupID, colname_secID, pcorr_cutoff_Model1_joint_models, no_of_ortho_pre_vs_Model1_joint_models, no_of_ortho_post_vs_Model1_joint_models, pcorr_cutoff_Model1_stratified_models, no_of_ortho_pre_vs_Model1_stratified_models, no_of_ortho_post_vs_Model1_stratified_models){
   setwd(directory_input_matrix_sampleID)
   file_sampleID <- read.table(filename_sampleID,header=T, dec = ".", row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t")
   setwd(directory_model_table_to_analyse)
   if (filename_model_table_to_analyse %in% dir(directory_model_table_to_analyse )) {
     model_table_to_analyse <- as.data.frame(read.table(filename_model_table_to_analyse, header=T, dec = ".", row.names=1, check.names = FALSE, na.strings=c("", "NA", "Inf"), sep="\t"))
   } else {

     model_table_to_analyse <- create_table_to_analyze(file_sampleID, colname_groupID, colname_secID, pcorr_cutoff_Model1_joint_models, no_of_ortho_pre_vs_Model1_joint_models, no_of_ortho_post_vs_Model1_joint_models, pcorr_cutoff_Model1_stratified_models, no_of_ortho_pre_vs_Model1_stratified_models, no_of_ortho_post_vs_Model1_stratified_models)
     write.table(model_table_to_analyse, paste("first_created",filename_model_table_to_analyse,sep=""), row.names=T,quote=F,sep="\t")

     model_table_to_analyse_reordered <- create_table_to_analyze_reordered(reordered_levels_of_groups, file_sampleID, colname_groupID, colname_secID, pcorr_cutoff_Model1_joint_models, no_of_ortho_pre_vs_Model1_joint_models, no_of_ortho_post_vs_Model1_joint_models, pcorr_cutoff_Model1_stratified_models, no_of_ortho_pre_vs_Model1_stratified_models, no_of_ortho_post_vs_Model1_stratified_models)
     write.table(model_table_to_analyse_reordered, paste("reordered",filename_model_table_to_analyse,sep="_"), row.names=T,quote=F,sep="\t")
     model_table_to_analyse<-model_table_to_analyse_reordered

   if (models_to_run!="all") {
   model_table_to_analyse <- select_models_to_run(model_table_to_analyse, models_to_run)
   write.table(model_table_to_analyse, paste("selected",filename_model_table_to_analyse,sep="_"), row.names=T,quote=F,sep="\t")
   }
 write.table(model_table_to_analyse, filename_model_table_to_analyse, row.names=T,quote=F,sep="\t")
   }
   model_table_to_analyse
 }


 makeproject <- function(directory_of_analysis, directory_output_reports) {

   if (dir.exists(directory_output_reports)) {
     message("Not performed. The output report folder already exists.")
   } else {
     # main dir

     setwd(directory_of_analysis)
     # sub dir
     l.subdir <- c(directory_output_reports)
     sapply(l.subdir,dir.create)
   }
 }

plotScorepredO1 <- function(resultmodel) {
  toplot <-cbind(resultmodel$aftervsdata.oplsda@scoreMN,resultmodel$aftervsdata.oplsda@orthoScoreMN[,1])
  toplot <- as.data.frame(toplot)
  ggplot(toplot, aes(p1, V2)) +
    geom_point() +
    geom_text_repel(label=rownames(toplot))
}
plotloadingpredO1 <- function(resultmodel) {
  toplot <-cbind(resultmodel$aftervsdata.oplsda@loadingMN,resultmodel$aftervsdata.oplsda@orthoLoadingMN[,1])
  toplot <- as.data.frame(toplot)
  ggplot(toplot, aes(p1, V2)) +
    geom_point() +
    geom_text_repel(label=rownames(toplot))
}



reorder_levels_of_groups <- function(order_of_groups,levels_of_groups) {
  if (order_of_groups=="correct") {
    reordered_levels_of_groups<-levels_of_groups} else {
  if (is.numeric(order_of_groups)) {
    reordered_levels_of_groups <- levels_of_groups[order_of_groups]} else {
    reordered_levels_of_groups <- order_of_groups}}
  reordered_levels_of_groups
}

select_models_to_run <- function(model_table_to_analyse,models_to_run) {
  if (models_to_run!="all") {
  model_table_to_analyse <- model_table_to_analyse[paste(models_to_run),]
  }
  model_table_to_analyse
}

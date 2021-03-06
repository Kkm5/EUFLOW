#require(IdMappingAnalysis)
#require(mvbutils)
options(stringsAsFactors = FALSE)
options(digits=3)
#
# #############example

#' WorkflowPathIdTag
#'
#' places tag on to an identifier
#' @param x Dataframe which represent sthe Evaluation dataset or the Reference Dataset
#' @param status tag to determine the evaluation or reference status
#' @param data.type Workflow category
#' @param version Workflow path number
#' @return the tagged identifier for further processing
#' @export
#'
WorkflowPathIdTag<-function(x,status,data.type,version){
  if (status == "workflow_path")
    status_tag<-paste(row.names(x),"_path","_",as.character(data.type),"_",as.character(version),sep="")
  if (status == "reference")
    status_tag<-paste(row.names(x),"_reference",sep="")
  return(status_tag)
}


#' WorkflowEvaluationData
#'
#' @param EvaluationExperimentSet         Merged set of Workflow options on the same samples.
#' @param ReferenceSet                    The Reference data set which is related to the Evaluation set as determined by the Model quality
#' @return An object of the class list which the first item is the reference dataset and all following list items are workflow options on the same samples
#' @export


WorkflowPathData<-function(EvaluationExperimentSet,ReferenceSet){

    names(EvaluationExperimentSet)[1]<- "Symbol"
    names(ReferenceSet)[1]<- "Symbol"
    row.names(EvaluationExperimentSet)<-EvaluationExperimentSet$Symbol
    row.names(ReferenceSet)<-as.character(ReferenceSet$Symbol)
    WorkflowList<-strsplit(row.names(EvaluationExperimentSet),"_")
    WorkflowNameVector<-sapply(WorkflowList, "[", 1)
    WorkflowOptionVector<-sapply(WorkflowList, "[", 2)
    #WorkflowNumber<-substr(WorkflowOptionVector,nchar(WorkflowOptionVector),nchar(WorkflowOptionVector))
    return(list(ReferenceSet,(split(EvaluationExperimentSet,WorkflowOptionVector))))

}


#' BuildEvaluationStructure
#'
#' @param Workflow.Data A list of length 2. The first item is the reference dataset, and the second item is a dataframe combining all the evaluation datasets.
#' @param ReferenceTag A string to use as a suffix for ID's from the reference dataset.
#' @param EvaluationTag A string to use as a suffix for ID's from the evaluation dataset.
#' @return Merged.options A dataframe with renamed labels
#' @export
BuildEvaluationStructure<-function(Workflow.Data,ReferenceTag="Protein",EvaluationTag="RNASeq"){
    EvaluationList<-Workflow.Data[[2]]
    Merged.options<-Workflow.Data[[1]]
    row.names(Merged.options)<-WorkflowPathIdTag(Merged.options,status="reference",ReferenceTag,1)
    for(o in c(1:length(EvaluationList))) {
        Evaluation_dataframe<-as.data.frame(EvaluationList[[o]])
        SymbolList<-strsplit(row.names(Evaluation_dataframe),"_")
        row.names(Evaluation_dataframe)<-sapply(SymbolList, "[", 1)
        row.names(Evaluation_dataframe)<-WorkflowPathIdTag(Evaluation_dataframe,status="workflow_path",EvaluationTag[o],o)
        Merged.options=rbind(Merged.options,Evaluation_dataframe)
    }
    return(Merged.options)
}




#' WorkflowPathMap
#'
#' Make a workflow map
#'
#' @param Merged.options A concatenation of the evaluation and reference data set features.
#' @return A data frame with two columns: "drivers" (the IDs for the features in the reference data set), and "workflow_options_merged":
#' @details  The "drivers" are the ID
#' @export
#'
WorkflowPathMap <- function(Merged.options){
    reference<-row.names(Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "reference",])
    workflow_path_data<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "path",]
    count.options<-function(x) {
        paste(row.names(workflow_path_data[sapply(strsplit(row.names(workflow_path_data),"_"), "[", 4) == x,]))
    }
    imax<-max(unique(sapply(strsplit(row.names(workflow_path_data),"_"), "[", 4)),na.rm = TRUE)
    workflow_path_matrix<-sapply(1:imax,count.options)
    if(class(workflow_path_matrix)=="matrix"){
        workflow_paths_combined<-sapply(1:dim(workflow_path_matrix)[1],function(i){paste0(as.character(workflow_path_matrix[i,]),collapse=",")})
        WorkflowMap<-data.frame(reference,workflow_paths_combined)
    }
    if(class(workflow_path_matrix)!="matrix"){
        ref_ids<-sapply(strsplit(row.names(Merged.options),"_"), "[", 1)[1:length(reference)]
        ref_ids_df<-as.data.frame(ref_ids)
        ref_ids_df$positions<-row.names(ref_ids_df)
        workflow_path_df<-as.data.frame(unlist(workflow_path_matrix))
        names(workflow_path_df)<-"path_ids"
        Symbol<-sapply(strsplit(as.character(workflow_path_df$path_ids),"_"), "[", 1)
        workflow_path_df$positions<-as.character(match(Symbol,ref_ids))
        workflow_paths_combined_df<-merge(workflow_path_df,ref_ids_df,by="positions")
        workflow_paths_combined<-sapply(ref_ids,function(i){paste0(workflow_paths_combined_df[workflow_paths_combined_df$ref_ids ==as.character(i),]$path_ids,collapse=",")})
        workflow_paths_combined<-as.data.frame(workflow_paths_combined)
        WorkflowMap<-data.frame(reference,workflow_paths_combined$workflow_paths_combined)
        names(WorkflowMap)<-c("reference","workflow_paths_combined")
    }
    return(WorkflowMap)
}


#' WorkflowPathModelQuality
#'
#' Produces a 'model quality object'
#'
#' @param Merged.options See other doc.
#' @return A "model quality object".
#' @export
#'
#'
WorkflowPathModelQuality<-function(Merged.options){
    WorkflowMap.object<-WorkflowPathMap(Merged.options)
    IdMap.example<-IdMap(DF=WorkflowMap.object,name="Workflowmap.object", primaryKey="reference",secondaryKey="workflow_paths_combined")
    WorkflowMap.object$workflow_paths_combined = as.character(WorkflowMap.object$workflow_paths_combined)
    secondaryIDs<-unlist(strsplit(WorkflowMap.object$workflow_paths_combined,","))
    uniquePairs_workflow <- as.UniquePairs.IdMap(IdMap.example,secondaryIDs)
    reference<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"),"[",2) == "reference",]
    names(reference)[names(reference)=="Symbol"] <- "reference"
    reference$reference<-row.names(reference)
    evaluation<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"),"[",2) != "reference",]
    names(evaluation)[names(evaluation)=="Symbol"] <- "workflow_paths_combined"
    evaluation$workflow_paths_combined<-row.names(evaluation)
    Model.quality.object<-CorrData(uniquePairs_workflow,reference,evaluation)
    return(Model.quality.object)
}

#' ModelQualityPairs
#'
#' User specifies the model quality used to evaluate the workflow
#' @param Model.quality.object The pairs (what kind of object?) with data to be used to estimate model quality
#' @param method Not sure what this parameter is.
#' @return "Values of model quality per pair"
#' @export
#'
ModelQualityPairs<-function(Model.quality.object, method=method){
  Model.quality<- Corr(Model.quality.object,method=method,verbose=FALSE)
  return(Model.quality)
}



#' ModelQualityFitToClusters
#'
#' "Performs the EM algorithm over the bootstrap values"
#' @param Multiple parameters - to be completed later
#' @return "a dataframe of posterior odds and posterior probability variance"
#' @export
#'
ModelQualityFitToClusters<-function(Y, Ysigsq,
         bootModel,
         piStart = c(0.5, 0.5),
         VStart = c(0.1,0.1),
         psiStart = c(0,0.1),
         NinnerLoop = 1,
         nReps=500,
         psi0Constraint,
         V0Constraint,
         sameV=FALSE,
         estimatesOnly=TRUE,
         printMe = TRUE,
         plotMe = TRUE,
         testMe=FALSE,
         Ntest = 5000,
         seed) {
  ### EM algorithm for 2 clusters,
  ### with constraints on the cluster means and variances, and known data variances
  if(testMe) {
    if(missing(seed)) .Random.seed <<- Random.seed.save
    else if(!is.na(seed)) .Random.seed <<- seed
    #  NA ==>  a new dataset.
    simPsi = c(0, 0.4)  ##
    simPi = c(2/3, 1/3)
    simData = data.frame(G = 1+rbinom(Ntest, 1, simPi[2]))
    simV = c(0.05^2, 0.05^2)
    simData$Ysigsq = rgamma(Ntest, 10, 400)
    simData$sd = sqrt(simV[simData$G] +simData$Ysigsq)
    simData = within(simData, Y <- simPsi[G] + rnorm(Ntest)*sqrt(simV[G]) + rnorm(Ntest)*sqrt(Ysigsq))
    print(summary(simData$Y))
    Y = simData$Y
    Ysigsq = simData$Ysigsq
  }
  ###############  Begining of EM algorithm  ##############
  piStar = piStart
  VStar = VStart
  psiStar = psiStart
  stopMe = FALSE
  iRep = 0
  while(1) {
    iRep = iRep + 1
    #    catn(", ", missing(V0Constraint))
    if(!missing(V0Constraint))
      VStar[1] = V0Constraint
    if(!missing(psi0Constraint))
      psiStar[1] = psi0Constraint
    #    print(psiStar)
    piStarOdds = piStar[2]/piStar[1]
    piStarOddsGK = piStarOdds *
      dnorm(Y, psiStar[2], sqrt(VStar[2] + Ysigsq)) /
      dnorm(Y, psiStar[1], sqrt(VStar[1] + Ysigsq))
    piStarGK = cbind(1/(1+piStarOddsGK), piStarOddsGK/(1+piStarOddsGK))
    EstarN = apply(piStarGK, 2, sum)
    piStar = apply(piStarGK, 2, mean)
    psiHat = psiStar
    VHat = VStar
    for(iRepInner in 1:NinnerLoop) {
      varHatTotal = colSums(outer(Y, psiHat, "-")^2 * piStarGK)
      sigsqTotal  = Ysigsq %*% piStarGK
      VHat = pmax(0, varHatTotal - sigsqTotal) / EstarN
      if(!missing(V0Constraint))
        VHat[1] = V0Constraint
      psiHat = colSums( Y %*% (piStarGK
                               / outer(Ysigsq, VHat, "+"))) /
        colSums( piStarGK
                 / outer(Ysigsq, VHat, "+"))
      if(!missing(psi0Constraint))
        psiHat[1] = psi0Constraint
      if(sameV)
        VHat[1] = VHat[2] = mean(VHat[1], VHat[2])
      if(max(abs(psiHat-psiStar), abs(VHat-VStar)) < 1e-7)
        stopMe = TRUE;
      psiHat -> psiStar
      VHat -> VStar
    }
    if(iRep >= nReps) stopMe = TRUE
    if(stopMe) break
  }
  cat(iRep, ifelse(iRep==nReps, ".  Loop exhausted.", ".  Converged."), "\n")

  if(plotMe) {
    options(echo=F)
    plot(col="blue", type = "l",main ="Mixture density", Ytemp<-seq(-1,1,length=100),
         xlab= as.character(names(bootModel[3])), cex.lab = 1.5, ylab="Density", sub = as.character(paste(names(bootModel[1]),"vs",names(bootModel[2]))),
         piStar[1]*dnorm(Ytemp, psiStar[1], sqrt(VStar[1] + mean(Ysigsq)))
         +
           piStar[2]*dnorm(Ytemp, psiStar[2], sqrt(VStar[2]  + mean(Ysigsq)))
    )
    for(g in 1:2) lines(col="blue", lty=2, Ytemp<-seq(-1,1,length=100),
                        piStar[g]*dnorm(Ytemp, psiStar[g], sqrt(VStar[g]  + mean(Ysigsq))))
    lines(density(Y), lwd=2, col="black")
    abline(v = 0)
    #abline(v = 0.38352)
    ###  Should we make a better choice than the means of the Ysigsq?
    if(testMe) lines(col="red", Ytemp<-seq(-1,1,length=100),
                     piStar[1]*dnorm(Ytemp, simPsi[1], sqrt(simV[1] + mean(simData$Ysigsq)))
                     +
                       piStar[2]*dnorm(Ytemp, simPsi[2], sqrt(simV[2]  + mean(simData$Ysigsq)))
    )
    legend(x=par("usr")[1], y=par("usr")[4],
           legend=c(ifelse(testMe, "truth", ""),
                    "data smooth", "estimate", " x or 0 component", " + component"),
           col=c("red", "black", "blue", "blue", "blue"),
           lty=c(ifelse(testMe, 1,0),1,1,2,2),
           lwd=c(1,2,1,1,1)
    )
    options(echo=T)
  }
  estimates = c(pi1=piStar[2], psi0=psiHat[1],
                psi1=psiHat[2], Var0=VStar[1], Var1=VStar[2])

  posteriorOdds =
    piStar[2]*dnorm(Y, psiHat[2], sqrt(VStar[2] + Ysigsq)) /
    piStar[1]/dnorm(Y, psiHat[1], sqrt(VStar[1] + Ysigsq))
  postProb = posteriorOdds/(1+posteriorOdds)
  postProbVar = Ysigsq * (postProb*(1-postProb))^2 *
    ((Y-psiHat[1])/(VStar[1]+Ysigsq) - (Y-psiHat[2])/(VStar[2]+Ysigsq))^2

  if(testMe) {
    simTruth = c(pi1=simPi[2], psi0=simPsi[1],
                 psi1=simPsi[2], Var0=simV[1], Var1=simV[2])
    estimates = data.frame(row.names=c("true","estimated"),
                           rbind(simTruth, estimates))
  }
  if(estimatesOnly) return(estimates)
  else {
    attr(x=posteriorOdds, which="estimates") = estimates
    return(data.frame(posteriorOdds,postProbVar))
  }


}

#' EstimatePosteriorProbability
#'
#' Calls the fit2clusters.workflow function and creates a dataframe with bootstap values, posterior odds, posterior probability variance, for each identifier pair
#' @param Model quality data and Model Quality
#' @return a evaulatioin ready dataframe of the posterior estimates
#' @export
#'
EstimatePosteriorProbability<-function(Model.quality.object,Model.Quality,postProb=NULL,postProbVar=NULL){
  bootstrap<-Bootstrap(Model.quality.object,Fisher=TRUE,verbose=FALSE)
  bootModel<-as.data.frame(bootstrap)
  bootModel<-bootModel[complete.cases(bootModel),]
  pairs<-bootModel[,1:2]
  #bootModel<-bootModel[c(1:96,99:134),]
  EMtest<-ModelQualityFitToClusters(bootModel$corr, bootModel$sd^2,bootModel,psi0Constraint=0, sameV=T,estimatesOnly=F,seed=Random.seed.save)
  #again part of the EMtest
  postProbs<-as.vector(EMtest[[1]]/(1+EMtest[[1]]))  #not needed at the output is a datafrmae from fit2clusters
  postProbVar <-as.vector(EMtest[[2]])
  bootMergedWithPairs = merge(data.frame(postProbs=postProbs, postProbVar=postProbVar,bootModel), pairs)
  return(bootMergedWithPairs)
}

#' expectedUtility
#'
#' Given the Posterior estimates evaluate the expected utility for a users determined Utp,Lfp,deltaPlus,guarantee
#' @param dataset the resulting dataframe from Workflow.posteriorestimate
#' @param a label for the workflow
#' @param Loss of a false positive
#' @param Utility of a true positive
#' @param deltaPlus parameter to estimate error
#' @param guarantee starting point
#' @return the expected utilty and associated estimates
#' @export
#'
WorkflowPathExpectedUtility<-function(dataset, label="", Lfp=1,Utp=1,deltaPlus=1,guarantee=1e-5)
    {
    postProbVar = pmax(dataset$postProbVar, guarantee)
    PrPlus = sum(dataset$postProbs/postProbVar)/
        sum(1/postProbVar)
#     result = data.frame(label=label,
#                         Utp=Utp, Lfp=Lfp, deltaPlus=deltaPlus,
#                         nPairs=nrow(dataset),
#                         PrPlus= PrPlus,
#                         PrTrue= PrTrue<-PrPlus / deltaPlus,
#                         PrFalse= PrFalse<-1 - PrTrue,
#                         Utrue=  Utrue<-PrTrue * Utp,
#                         Lfalse= Lfalse<-PrFalse * Lfp,
#                         Eutility1= Utrue-Lfalse,
#                         Eutility= nrow(dataset)*(Utrue-Lfalse))
    result = data.frame(nPairs=nrow(dataset),PrPlus= PrPlus,PrTrue= PrTrue<-PrPlus / deltaPlus,PrFalse= PrFalse<-1 - PrTrue,Utrue=  Utrue<-PrTrue * Utp,Lfalse= Lfalse<-PrFalse * Lfp,Eutility1= Utrue-Lfalse,Eutility= nrow(dataset)*(Utrue-Lfalse))
    rownames(result) = label
    return(result)
}


#' WorkflowEvaluationTable
#'
#' WorkflowEvaluationTable:  Produces an expected utility table for guidance and evaluation
#'
#' @param Posterior.dataframe Produced by Workflow.posteriorestimate().
#' @return Nicely formatted table of posterior probabilities, Pr(+) and Pr(-), standard deviations, model quality scores, and biases.
#' @export
WorkflowEvaluationTable<-function(Posterior.dataframe,Lfp=1,Utp=1,deltaPlus=1,guarantee=1e-5){
    WorkflowStats<-data.frame(sapply(strsplit(Posterior.dataframe$workflow_paths_combined,"_"),"[",1),sapply(strsplit(Posterior.dataframe$workflow_paths_combined,"_"),"[",4),Posterior.dataframe)
    colnames(WorkflowStats)[1]<-"Marker"
    colnames(WorkflowStats)[2]<-"WorkflowID"
    WorkflowLabelDF<-WorkflowStats[WorkflowStats$WorkflowID==1,]
    WorkflowLabel<-as.data.frame(strsplit(WorkflowLabelDF$workflow_paths_combined[1],"_"))[3,]
    types<-2:max(WorkflowStats$WorkflowID)
    Evaluation.table<-WorkflowPathExpectedUtility(dataset=WorkflowStats[WorkflowStats$WorkflowID==1,],Lfp=Lfp,Utp=Utp,deltaPlus=deltaPlus,guarantee=guarantee,label=WorkflowLabel)
    for(i in types){
        WorkflowLabelDF<-WorkflowStats[WorkflowStats$WorkflowID==i,]
        WorkflowLabel<-as.data.frame(strsplit(WorkflowLabelDF$workflow_paths_combined[i],"_"))[3,]
        Evaluation.table<-rbind(Evaluation.table,WorkflowPathExpectedUtility(dataset=WorkflowStats[WorkflowStats$WorkflowID==i,],Lfp=Lfp,Utp=Utp,deltaPlus=deltaPlus,guarantee=guarantee,label=WorkflowLabel))
    }
    return(Evaluation.table)
}


#' ExpectedUtilityPlot
#'
#' ExpectedUtilityPlot:  Produces a plot for the expected utility of workflow path
#'
#' @param Evaluation Table
#' @return Plot of Expected utility and number of filters applied
#' @export
ExpectedUtilityPlot<-function(Evaluation.table,EvaluationOrder="TEU"){
    i<-as.numeric(dim(Evaluation.table)[1])-1
    if(EvaluationOrder=="TEU"){
        plot(NA,xlim=c(-0.2,5.2),ylim=c(-100,300),main=paste(EvaluationOrder," vs Number of Filters Applied"),xlab="Number of Filters Applied",ylab=as.character(EvaluationOrder)) # make an empty plot
        points(c(0:i),Evaluation.table$Eutility,type="b",pch=1,lwd=2)
        text(c(0:i),Evaluation.table$Eutility,labels=row.names(Evaluation.table),pos=c(3,3,3))
    }
    if(EvaluationOrder=="MEU"){
        plot(NA,xlim=c(-0.2,5.2),ylim=c(-1,3),main=paste(EvaluationOrder," vs Number of Filters Applied"),xlab="Number of Filters Applied",ylab=as.character(EvaluationOrder)) # make an empty plot
        points(c(0:i),Evaluation.table$Eutility1,type="b",pch=1,lwd=2)
        text(c(0:i),Evaluation.table$Eutility1,labels=row.names(Evaluation.table),pos=c(3,3,3))
    }
}





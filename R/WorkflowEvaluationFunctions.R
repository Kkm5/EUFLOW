require(IdMappingAnalysis)
# #require(mvbutils)
options(stringsAsFactors = FALSE)
#
# #############example

.assign.status<-function(x,status,data.type,version){
  if (status == "workflow_option")
    status_tag<-paste(row.names(x),"_WFO","_",as.character(data.type),"_",as.character(version),sep="")
  if (status == "driver")
    status_tag<-paste(row.names(x),"_DRIVER",sep="")
  return(status_tag)
}


#' WorkflowEvaluationData
#'
#' @param EvaluationExperimentSet         Merged set of Workflow options on the same samples.
#' @param ReferenceSet                    The Reference data set which is related to the Evaluation set as determined by the Model quality
#' @return An object of the class list which the first item is the reference dataset and all following list items are workflow options on the same samples
#' @export


WorkflowEvaluationData<-function(EvaluationExperimentSet,ReferenceSet){

    names(EvaluationExperimentSet)[1]<- "Symbol"
    names(ReferenceSet)[1]<- "Symbol"
    row.names(EvaluationExperimentSet)<-EvaluationExperimentSet$Symbol
    row.names(ReferenceSet)<-ReferenceSet$Symbol
    WorkflowList<-strsplit(row.names(EvaluationExperimentSet),"_")
    WorkflowNameVector<-sapply(WorkflowList, "[", 1)
    WorkflowOptionVector<-sapply(WorkflowList, "[", 2)
    #WorkflowNumber<-substr(WorkflowOptionVector,nchar(WorkflowOptionVector),nchar(WorkflowOptionVector))
    return(list(ReferenceSet,(split(EvaluationExperimentSet,WorkflowOptionVector))))

}


#' Merge tag options
#'
#' @param Workflow.Data object where reference is the first item an the evaluation data is all further items in list
#' @param ReferenceTag which sets a particular symbol for the reference set tag
#' @param EvaluationTag which sets a particular symbol for the evaluation set tag
#' @return Merged.options a dataframe with renamed labels
#' @export

merge.tag.options<-function(Workflow.Data,ReferenceTag="P",EvaluationTag="RS"){
    EvaluationList<-Workflow.Data[[2]]
    Merged.options<-Workflow.Data[[1]]
    row.names(Merged.options)<-.assign.status(Merged.options,status="driver",ReferenceTag,1)
    for(o in c(1:length(EvaluationList))) {
        Evaluation_dataframe<-as.data.frame(EvaluationList[[o]])
        SymbolList<-strsplit(row.names(Evaluation_dataframe),"_")
        row.names(Evaluation_dataframe)<-sapply(SymbolList, "[", 1)
        row.names(Evaluation_dataframe)<-.assign.status(Evaluation_dataframe,status="workflow_option",EvaluationTag,o)
        Merged.options=rbind(Merged.options,Evaluation_dataframe)
    }
    return(Merged.options)
    }

# Merged.options<-merge.tag.options(Workflow.Data)

#' make.workflow.map
#'
#' Make a workflow map
#'
#' @param Merged options
#' @return WorkflowMap.object a dataframe of reference tags mapped to the respective evaluation tags
#' @export

make.workflow.map <- function(Merged.options){
  drivers<-row.names(Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "DRIVER",])
  workflow_options_data<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"), "[", 2) == "WFO",]
  workflow_options<-row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 2) == "WFO",])
  imax<-max(unique(sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4)),na.rm = TRUE)
  imin<-min(unique(sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4)),na.rm = TRUE)
  workflow_options_merged<-paste(row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4) == imin,]),row.names(workflow_options_data[sapply(strsplit(row.names(workflow_options_data),"_"), "[", 4) == imax,]),sep=",")
  WorkflowMap<-data.frame(drivers,workflow_options_merged)
  ###when your ready class(Merged.options)<- "WorkflowMap"
  return(WorkflowMap)
}


Workflow.Criterion<-function(Model.quality.object){
  Model.quality<- Corr(Model.quality.object,method="spearman",verbose=TRUE)
  return(Model.quality)
}


fit2clusters.workflow<-function(Y, Ysigsq,
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
    abline(v = 0.38352)
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

Workflow.posteriorestimate<-function(Model.quality.object,Model.Quality,postProb=NULL,postProbVar=NULL){
  bootstrap<-Bootstrap(Model.quality.object,Fisher=TRUE,verbose=TRUE)
  bootModel<-as.data.frame(bootstrap)
  bootModel<-bootModel[complete.cases(bootModel),]
  pairs<-bootModel[,1:2]
  #bootModel<-bootModel[c(1:96,99:134),]
  EMtest<-fit2clusters.workflow(bootModel$corr, bootModel$sd^2,bootModel,psi0Constraint=0, sameV=T,estimatesOnly=F,seed=Random.seed.save)
  #again part of the EMtest
  postProbs<-as.vector(EMtest[[1]]/(1+EMtest[[1]]))  #not needed at the output is a datafrmae from fit2clusters
  postProbVar <-as.vector(EMtest[[2]])
  bootMergedWithPairs = merge(data.frame(postProbs=postProbs, postProbVar=postProbVar,bootModel), pairs)
  return(bootMergedWithPairs)
}


expectedUtility<-function(dataset, label="", Lfp=1,Utp=1,deltaPlus=1,guarantee=1e-5)
    {
    postProbVar = pmax(dataset$postProbVar, guarantee)
    PrPlus = sum(dataset$postProbs/postProbVar)/
        sum(1/postProbVar)
    result = data.frame(label=label,
                        Utp=Utp, Lfp=Lfp, deltaPlus=deltaPlus,
                        nPairs=nrow(dataset),
                        PrPlus= PrPlus,
                        PrTrue= PrTrue<-PrPlus / deltaPlus,
                        PrFalse= PrFalse<-1 - PrTrue,
                        Utrue=  Utrue<-PrTrue * Utp,
                        Lfalse= Lfalse<-PrFalse * Lfp,
                        Eutility1= Utrue-Lfalse,
                        Eutility= nrow(dataset)*(Utrue-Lfalse))
    rownames(result) = label
    return(result)
}


#' Workflow.Evaluation.table
#'
#' Workflow.Evaluation.table:  Produces an expected utility table for guidance and evaluation
#'
#' @param Posterior.dataframe Produced by Workflow.posteriorestimate().
#' @return Nicely formatted table of posterior probabilities, Pr(+) and Pr(-), standard deviations, model quality scores, and biases.


Workflow.Evaluation.table<-function(Posterior.dataframe,Lfp=1,Utp=1,deltaPlus=1,guarantee=1e-5){
    WorkflowStats<-data.frame(sapply(strsplit(Posterior.dataframe$workflow_options_merged,"_"),"[",1),sapply(strsplit(Posterior.dataframe$workflow_options_merged,"_"),"[",4),Posterior.dataframe)
    colnames(WorkflowStats)[1]<-"Marker"
    colnames(WorkflowStats)[2]<-"WorkflowID"
    og<-expectedUtility(label="Use All", dataset=WorkflowStats,Lfp=Lfp,Utp=Utp,deltaPlus=deltaPlus,guarantee=guarantee)
    for(p in unique(WorkflowStats$WorkflowID)){
        set<-expectedUtility(label=as.character(WorkflowStats$WorkflowID[as.numeric(p)]), dataset=WorkflowStats[WorkflowStats$WorkflowID==as.numeric(p),],Lfp=Lfp,Utp=Utp,deltaPlus=deltaPlus,guarantee=guarantee)
        og=rbind(og,set)
    }
    return(og)
}






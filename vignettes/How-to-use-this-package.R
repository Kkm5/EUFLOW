## ----setup, include=FALSE------------------------------------------------
library(knitr)
knitr::opts_chunk$set(echo = TRUE)
options(stringsAsFactors = FALSE)
system("touch I_ran_Rmd")

## ------------------------------------------------------------------------
search() 


## ------------------------------------------------------------------------
ls(pos=2)

## ------------------------------------------------------------------------
data(RNASEQDATA)
RNASEQDATA[1:9,1:5]
RNASEQDATA[68:77,1:5]

## ------------------------------------------------------------------------
data(RPPADATA.original)
RPPADATA<-RPPADATA.original
RPPADATA[1:9,1:5]

## ------------------------------------------------------------------------

EvaluationExperimentSet<-RNASEQDATA
ReferenceSet<-RPPADATA

## ----eval=FALSE----------------------------------------------------------
#  
#  Workflow.Data<-WorkflowEvaluationData(EvaluationExperimentSet,ReferenceSet)
#  Merged.options<-merge_tag_options(Workflow.Data)
#  Model.quality.object<-Model.quality.list(Merged.options)
#  Model.Quality<-Workflow.Criterion(Model.quality.object)
#  Posterior.dataframe<-Workflow.posteriorestimate(Model.quality.object,Model.Quality)
#  Workflow.Evaluation.table(Posterior.dataframe)
#  Workflow.Evaluation.table(Posterior.dataframe,deltaPlus = 2)
#  Workflow.Evaluation.table(Posterior.dataframe,Utp=3)
#  Workflow.Evaluation.table(Posterior.dataframe,Utp=1)


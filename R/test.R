
setwd("C:/Users/Mcdade/Desktop/EUFLOW")
RNASEQDATA<-read.csv(file="RNASEQDATA.csv",header=TRUE)
RPPADATA<-read.csv(file="RPPADATA.original.csv",header=TRUE)
EvaluationExperimentSet<-RNASEQDATA
#EvaluationExperimentSet<-rbind(RNASEQDATA,RANDOMSET)
ReferenceSet<-RPPADATA
Workflow.Data<-WorkflowEvaluationData(EvaluationExperimentSet,ReferenceSet)
Merged.options<-merge.tag.options(Workflow.Data)
Model.quality.object<-Model.quality.list(Merged.options)
Model.Quality<-Workflow.Criterion(Model.quality.object)
Posterior.dataframe<-Workflow.posteriorestimate(Model.quality.object,Model.Quality)
Workflow.Evaluation.table(Posterior.dataframe)
Workflow.Evaluation.table(Posterior.dataframe,deltaPlus = 2)
Workflow.Evaluation.table(Posterior.dataframe,Utp=3)
Workflow.Evaluation.table(Posterior.dataframe,Utp=1)

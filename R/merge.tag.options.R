

RNASEQDATA<-read.csv(file="RNASEQDATA.csv",header=TRUE)
RPPADATA<-read.csv(file="RPPADATA.original.csv",header=TRUE)
EvaluationExperimentSet<-RNASEQDATA
ReferenceSet<-RPPADATA
Workflow.Data<-WorkflowEvaluationData(EvaluationExperimentSet,ReferenceSet)

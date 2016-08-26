

RNASEQDATA<-read.csv(file="RNASEQDATA.csv",header=TRUE)
RPPADATA<-read.csv(file="RPPADATA.original.csv",header=TRUE)
#EvaluationExperimentSet<-RNASEQDATA
EvaluationExperimentSet<-rbind(RNASEQDATA,RANDOMSET)
ReferenceSet<-RPPADATA
Workflow.Data<-WorkflowEvaluationData(EvaluationExperimentSet,ReferenceSet)
Merged.options<-merge.tag.options(Workflow.Data)
WorkflowMap.object<-make.workflow.map(Merged.options)

IdMap.example<-IdMap(DF=WorkflowMap.object,name="Workflowmap.object", primaryKey="drivers",secondaryKey="workflow_options_merged")

secondaryIDs<-unlist(strsplit(WorkflowMap.object$workflow_options_merged,","))
uniquePairs_workflow <- as.UniquePairs.IdMap(IdMap.example,secondaryIDs)


P<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"),"[",2) == "DRIVER",]
names(P)[names(P)=="Symbol"] <- "drivers"
P$drivers<-row.names(P)
RSALL<-Merged.options[sapply(strsplit(row.names(Merged.options),"_"),"[",2) != "DRIVER",]
names(RSALL)[names(RSALL)=="Symbol"] <- "workflow_options_merged"
RSALL$workflow_options_merged<-row.names(RSALL)
Model.quality.object<-CorrData(uniquePairs_workflow,P,RSALL)
Model.Quality<-Workflow.Criterion(Model.quality.object)
Posterior.dataframe<-Workflow.posteriorestimate(Model.quality.object,Model.Quality)
Workflow.Evaluation.table(Posterior.dataframe)
Workflow.Evaluation.table(Posterior.dataframe,deltaPlus = 2)
Workflow.Evaluation.table(Posterior.dataframe,Utp=3)
Workflow.Evaluation.table(Posterior.dataframe,Utp=1)

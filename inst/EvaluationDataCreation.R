fn <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62944/suppl/GSE62944%5F06%5F01%5F15%5FTCGA%5F24%5FNormal%5FCancerType%5FSamples.txt.gz"
download.file(fn,destfile="tmp.txt.gz")


fn <- "http://s.wordpress.org/resources/survey/wp2011-survey.tar.gz"
download.file(fn,destfile="tmp.tar.gz")
untar("GSE62944_RAW.tar",list=TRUE)  ## check contents
untar("GSE62944_RAW.tar")
## or, if you just want to extract the target file:
untar("tmp.tar.gz",files="wp2011-survey/anon-data.csv")
X <- read.csv("wp2011-survey/anon-data.csv")

untar("GSE62944_RAW.tar",files="GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz")




NormalFeatureCounts<-read.table(file="GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt.gz", header=TRUE)
CancerFeatureCounts<-read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt.gz", header=TRUE)
CancerFPKM<-read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FPKM.txt.gz", header=TRUE)
CancerTPM<-read.table(file="GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_TPM.txt.gz", header=TRUE)

#################Cancer Type file formatting
CancerTCGAtype<-read.table(file="GSE62944_06_01_15_TCGA_24_CancerType_Samples.txt.gz",header=TRUE)
CancerTCGAtype<-rbind(names(CancerTCGAtype),CancerTCGAtype)
names(CancerTCGAtype)<-c("TCGA","Cancer")
CancerTCGAtype$TCGA<-gsub("[-]",".",CancerTCGAtype$TCGA) 


######Cancer type selection and file creation
BRCA.names<-CancerTCGAtype[CancerTCGAtype$Cancer=="BRCA",]

OV.names<-CancerTCGAtype[CancerTCGAtype$Cancer=="OV",]

BRCACancerFeatureCounts<-CancerFeatureCounts[,names(CancerFeatureCounts) %in% BRCA.names$TCGA]
OVCancerFeatureCounts<-CancerFeatureCounts[,names(CancerFeatureCounts) %in% OV.names$TCGA]


COMPLETERPPADATA<-RPPADATA$X
OVCancerFeatureCountsRPPAOverlap<-OVCancerFeatureCounts[row.names(OVCancerFeatureCounts) %in% COMPLETERPPADATA,]
OV.name.split<-strsplit(names(OVCancerFeatureCountsRPPAOverlap),"[.]")
OV.name.matrix<-apply(as.matrix(OV.name.split.df),MARGIN=1,paste,sep=".")
just3<-as.data.frame(OV.name.matrix)
names(just3)<-c("first","second","third")
OV.tcga.names<-paste(as.character(just3$first),as.character(just3$second),as.character(just3$third),sep=".")

OVCancerFeatureCountsDATA<-OVCancerFeatureCountsRPPAOverlap
names(OVCancerFeatureCountsDATA)<-OV.tcga.names

RNASEQDATA<-read.csv(file="RNASEQDATA.csv",head=TRUE)
OVMERGE<-log2(OVCancerFeatureCountsDATA)
OVMERGE<-cbind(row.names(OVMERGE),OVMERGE)
names(OVMERGE)[names(OVMERGE)=="row.names(OVMERGE)"] <- "X"
OVMERGE$X<-paste(OVMERGE$X,"_v3",sep="")
row.names(OVMERGE)<-NULL
OVMERGE<-OVMERGE[,names(OVMERGE) %in% names(RNASEQDATA)]

write.csv(OVCancerFeatureCountsDATA,file="OVCancerFeatureCountsDATA.csv")





row.names(BRCARPPADATA)<-BRCARPPADATA$COMMON
BRCARPPADATA.original<-BRCARPPADATA

rppa_names_split<-as.data.frame(strsplit(names(BRCARPPADATA),"[.]"))
rppa_names_split_df<-as.data.frame(rppa_names_split)
rppa_names_df<-as.data.frame(apply(as.matrix(rppa_names_split_df),MARGIN=1,paste,sep="."))
names(rppa_names_df)<-c("first","second","third")
rppa_names_reformed<-paste(as.character(rppa_names_df$first),as.character(rppa_names_df$second),as.character(rppa_names_df$third),sep=".")
names(BRCARPPADATA)<-rppa_names_reformed

BRCARPPADATAFINAL<-BRCARPPADATA[, colSums(is.na(BRCARPPADATA)) != nrow(BRCARPPADATA)]
X<-row.names(BRCARPPADATAFINAL)
BRCARPPADATAFINAL<-cbind(X,BRCARPPADATAFINAL)
row.names(BRCARPPADATAFINAL)<-NULL


BRCASALMONDATAFINAL<-BRCASALMONDATA[,names(BRCASALMONDATA) %in% names(BRCARPPADATAFINAL)]
BRCARPPADATAFINAL<-BRCARPPADATAFINAL[,names(BRCARPPADATAFINAL) %in% names(BRCASALMONDATAFINAL)]

COMPLETERPPADATA<-RPPADATA$X[as.character(RPPADATA$X) %!in% c("BIRC2","DIABLO","NCOA3","RB1","YAP1")]
BRCARPPADATAFINAL<-BRCARPPADATAFINAL[c(names(BRCASALMONDATAFINAL))]

BRCARPPADATAFINAL<-BRCARPPADATAFINAL[BRCARPPADATAFINAL$X %in% COMPLETERPPADATA,]

MissingBRCASALMONDATA<-c(paste(c("BIRC2","DIABLO","NCOA3","RB1","YAP1"),c("_v1"),sep=""),paste(c("BIRC2","DIABLO","NCOA3","RB1","YAP1"),c("_v2"),sep=""),paste(c("BIRC2","DIABLO","NCOA3","RB1","YAP1"),c("_v3"),sep=""),paste(c("BIRC2","DIABLO","NCOA3","RB1","YAP1"),c("_v4"),sep=""),paste(c("BIRC2","DIABLO","NCOA3","RB1","YAP1"),c("_v5"),sep=""))

COMPLETEBRCASALMONDATA<-BRCASALMONDATA$X[as.character(BRCASALMONDATA$X) %!in% MissingBRCASALMONDATA]

BRCASALMONDATAFINAL<-BRCASALMONDATAFINAL[BRCASALMONDATAFINAL$X %in% COMPLETEBRCASALMONDATA,]



BRCASALMONTHRESHOLDDATA<-BRCASALMONDATA[1:67,]
BRCASALMONOVER1000<-as.data.frame(apply(BRCASALMONTHRESHOLDDATA[,2:407],1,mean)>1000)
names(BRCASALMONOVER1000)<-"OVER1000"
BRCASALMONOVER5000<-as.data.frame(apply(BRCASALMONTHRESHOLDDATA[,2:407],1,mean)>5000)
names(BRCASALMONOVER5000)<-"OVER5000"
BRCASALMONOVER10000<-as.data.frame(apply(BRCASALMONTHRESHOLDDATA[,2:407],1,mean)>10000)
names(BRCASALMONOVER10000)<-"OVER10000"

BRCASALMONTHRESHOLDDATA$X<-sapply(strsplit(as.character(BRCASALMONTHRESHOLDDATA$X),"_"),"[",1)

BRCASALMONTHRESHOLD<-data.frame(BRCASALMONTHRESHOLDDATA$X,BRCASALMONOVER1000$OVER1000,BRCASALMONOVER5000$OVER5000,BRCASALMONOVER10000$OVER10000)
names(BRCASALMONTHRESHOLD)<-c("X","OVER1000","OVER5000","OVER10000")


MissingBRCASymbol<-c("BIRC2","DIABLO","NCOA3","RB1","YAP1")
COMPLETETHRESHOLDIDS<-BRCASALMONTHRESHOLDDATA$X[as.character(BRCASALMONTHRESHOLDDATA$X) %!in% MissingBRCASymbol]

BRCASALMONDATAFINAL<-BRCASALMONDATA[,names(BRCASALMONDATA) %in% names(BRCARPPADATAFINAL)]
BRCASALMONTHRESHOLDDATAFINAL<-BRCASALMONTHRESHOLDDATA[BRCASALMONTHRESHOLDDATA$X %in% COMPLETETHRESHOLDIDS,]
BRCASALMONTHRESHOLDDATAFINAL<-BRCASALMONTHRESHOLDDATAFINAL[,names(BRCASALMONTHRESHOLDDATAFINAL) %in% names(BRCARPPADATAFINAL)]


OVER1000<-BRCASALMONTHRESHOLD[BRCASALMONTHRESHOLD$OVER1000==TRUE,]$X
OVER5000<-BRCASALMONTHRESHOLD[BRCASALMONTHRESHOLD$OVER5000==TRUE,]$X
OVER10000<-BRCASALMONTHRESHOLD[BRCASALMONTHRESHOLD$OVER10000==TRUE,]$X

ALLGENESDATA<-BRCASALMONTHRESHOLDDATAFINAL
OVER1000DATA<-BRCASALMONTHRESHOLDDATAFINAL[BRCASALMONTHRESHOLDDATAFINAL$X %in% OVER1000,]
OVER5000DATA<-BRCASALMONTHRESHOLDDATAFINAL[BRCASALMONTHRESHOLDDATAFINAL$X %in% OVER5000,]
OVER10000DATA<-BRCASALMONTHRESHOLDDATAFINAL[BRCASALMONTHRESHOLDDATAFINAL$X %in% OVER10000,]


ALLGENESDATA$X<-paste(ALLGENESDATA$X,"_v1",sep="")
OVER1000DATA$X<-paste(OVER1000DATA$X,"_v2",sep="")
OVER5000DATA$X<-paste(OVER5000DATA$X,"_v3",sep="")
OVER10000DATA$X<-paste(OVER10000DATA$X,"_v4",sep="")


BRCATHRESHOLD<-rbind(ALLGENESDATA,OVER1000DATA,OVER5000DATA,OVER10000DATA)
save(BRCATHRESHOLD,file="BRCATHRESHOLD.rda")


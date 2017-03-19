data(brca.tcga.salmon.gene.counts)

library(biomaRt)

variation = useEnsembl(biomart="ensembl", dataset="hsapiens_snp")

listFilters(mart)

listAttributes(mart)

mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)

Describe_hgnc<-getBM(attributes = c("hgnc_symbol","description"),
                filters    = "with_hgnc",
                values     = TRUE, 
                mart       = mart)

Describe_hgnc_map<-Describe_hgnc[Describe_hgnc$hgnc_symbol %in% row.names(BRCARPPADATA),]
row.names(Describe_hgnc_map)<-NULL
write.csv(Describe_hgnc_map, file="Describe_hgnc_map.csv")
save(Describe_hgnc_map, file="Describe_hgnc_map.rda")

Phenotype_hgnc<-getBM(attributes = c("hgnc_symbol","phenotype_description"),
                     filters    = "with_hgnc",
                     values     = TRUE, 
                     mart       = mart)

Phenotype_hgnc_map<-Phenotype_hgnc[Phenotype_hgnc$hgnc_symbol %in% row.names(BRCARPPADATA),]
row.names(Phenotype_hgnc_map)<-NULL
write.csv(Phenotype_hgnc_map,file="Phenotype_hgnc_map.csv")
save(Phenotype_hgnc_map,file="Phenotype_hgnc_map.rda")

All_hgnc<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters    = "with_hgnc",
      values     = TRUE, 
      mart       = mart)

All_hgnc_map<-All_hgnc[All_hgnc$hgnc_symbol %in% BRCARPPADATA$COMMON,]
row.names(All_hgnc_map)<-NULL

################v1
BRCA_SALMON_ALL_HGNC<-brca.tcga.salmon.gene.counts[row.names(brca.tcga.salmon.gene.counts) %in% All_hgnc_map$ensembl_gene_id,]

Biotype_proteincoding<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters    = "biotype",
                values     = "protein_coding", 
                mart       = mart)

Biotype_proteincoding_map<-Biotype_proteincoding[Biotype_proteincoding$hgnc_symbol %in% BRCARPPADATA$COMMMON]
row.names(Biotype_proteincoding_map)<-NULL


Biotype_pseudogene<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                             filters    = "biotype",
                             values     = "transcribed_processed_pseudogene", 
                             mart       = mart)

Biotype_pseudogene_map<-Biotype_pseudogene[Biotype_pseudogene$hgnc_symbol %in% BRCARPPADATA$COMMON,]
row.names(Biotype_pseudogene_map)<-NULL

Somatic_variation<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                             filters    = "somatic_variation_source",
                             values     = TRUE, 
                             mart       = mart)

SNP_valid<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters    = "with_validated_snp",
                         values     = TRUE, 
                         mart       = mart)


filterType("biotype",mart)

filterOptions("biotype",mart)


Transmembrane<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters    = "with_tmhmm",
                          values     = TRUE, 
                          mart       = mart)

Transmembrane_map<-Transmembrane[Transmembrane$hgnc_symbol %in% BRCARPPADATA$COMMON,]
No_Transmembrane_map<-All_hgnc_map[All_hgnc_map$ensembl_gene_id %!in% Transmembrane_map$ensembl_gene_id,]


###########v2
BRCA_SALMON_TRANSMEMBRANE<-brca.tcga.salmon.gene.counts[row.names(brca.tcga.salmon.gene.counts) %in% Transmembrane_map$ensembl_gene_id,]

###########v3
BRCA_SALMON_NO_TRANSMEMBRANE<-brca.tcga.salmon.gene.counts[row.names(brca.tcga.salmon.gene.counts) %in% No_Transmembrane_map$ensembl_gene_id,]

write.csv(Transmembrane_map,file="Transmembrane_map.csv")
write.csv(No_Transmembrane_map,file="No_Transmembrane_map.csv")

FormatDataFile<-function(data,map,version){
    names_split<-as.data.frame(strsplit(names(data),"[.]"))
    names_split_df<-as.data.frame(names_split)
    names_df<-as.data.frame(apply(as.matrix(names_split_df),MARGIN=1,paste,sep="."))
    names(names_df)<-c("first","second","third")
    names_reformed<-paste(as.character(names_df$first),as.character(names_df$second),as.character(names_df$third),sep=".")
    names(data)<-names_reformed
    X<-row.names(data)
    formatdata<-cbind(X,data)
    formatdata<-merge(map,formatdata,by.x="ensembl_gene_id",by.y="X")
    formatdata$ensembl_gene_id<-NULL
    colnames(formatdata)[1]<-"X"
    formatdata$X<-paste(formatdata$X,version,sep="_")
    row.names(formatdata)<-NULL
    #return()
    return(formatdata)
}

BRCARPPADATA$COMMON<-NULL
BRCARPPADATA$GENE_ID<-NULL



Tigrfam<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters    = "with_tigrfam",
                     values     = TRUE, 
                     mart       = mart)

Tigrfam_map<-Tigrfam[Tigrfam$hgnc_symbol %in% BRCARPPADATA$COMMON,]


Pfam<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters    = "with_protein_feature_pfam",
                     values     = TRUE, 
                     mart       = mart)

Pfam_map<-Pfam[Pfam$hgnc_symbol %in% BRCARPPADATA$COMMON,]


Panther<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
            filters    = "with_hmmpanther",
            values     = TRUE, 
            mart       = mart)

Panther_map<-Panther[Panther$hgnc_symbol %in% BRCARPPADATA$COMMON,]


Low_complexity<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
            filters    = "with_low_complexity",
            values     = TRUE, 
            mart       = mart)

Low_complexity_map<-Low_complexity[Low_complexity$hgnc_symbol %in% BRCARPPADATA$COMMON,]
No_Low_Complexity_map<-All_hgnc_map[All_hgnc_map$ensembl_gene_id %!in% Low_complexity_map$ensembl_gene_id,]

write.csv(Low_complexity_map,file="Low_complexity_map.csv")
write.csv(No_Low_Complexity_map,file="No_Low_Complexity_map.csv")

############v4
BRCA_SALMON_LOW_COMPLEXITY<-brca.tcga.salmon.gene.counts[row.names(brca.tcga.salmon.gene.counts) %in% Low_complexity_map$ensembl_gene_id,]

############v5
BRCA_SALMON_NO_LOW_COMPLEXITY<-brca.tcga.salmon.gene.counts[row.names(brca.tcga.salmon.gene.counts) %in% No_Low_Complexity_map$ensembl_gene_id,]


v1<-FormatDataFile(BRCA_SALMON_ALL_HGNC,All_hgnc_map,"v1")
v2<-FormatDataFile(BRCA_SALMON_TRANSMEMBRANE,Transmembrane_map,"v2")
v3<-FormatDataFile(BRCA_SALMON_NO_TRANSMEMBRANE,No_Transmembrane_map,"v3")
v4<-FormatDataFile(BRCA_SALMON_LOW_COMPLEXITY,Low_complexity_map,"v4")
v5<-FormatDataFile(BRCA_SALMON_NO_LOW_COMPLEXITY,No_Low_Complexity_map,"v5")

BRCASALMONDATA<-rbind(v1,v2,v3,v4,v5)





BiocManager::install("TCGAbiolinks")
BiocManager::install("survival")
BiocManager::install("survminer")
BiocManager::install("maftools")

library("TCGAbiolinks")
library("survival")
library("survminer")

library(maftools)
library(ggplot2)
library(pheatmap)

library(dplyr)
library(stringr)

library(purrr)

data_mut <- read.delim("data_mutations.txt", header = TRUE, sep = '\t')
clinical_data <- read.table("data_clinical_patient.txt", header = TRUE, sep = '\t')
rna_seq <- read.csv("RNAseq_BRCA.csv")
datatest <- read.table("data_mutations.txt", header = TRUE, sep = '\t')


###

#isolating common patient data, finding unique patients
unique_pat_mut <- data_mut
unique_pat_mut$Tumor_Sample_Barcode = substr(unique_pat_mut$Tumor_Sample_Barcode,1,
                                             nchar(unique_pat_mut$Tumor_Sample_Barcode)-3)
unique_pat_mut <- unique(unique_pat_mut$Tumor_Sample_Barcode)
unique_pat_clin <- unique(clinical_data$PATIENT_ID)
#finding common patients between data_clinical_patient and data_mutations
K <- intersect(unique_pat_clin,unique_pat_mut)
#editing and pulling out patient ID from RNAseq_BRCA
unique_rna_data <- colnames(rna_seq)
unique_rna_data <- unique_rna_data[-1]
unique_rna_data = substr(unique_rna_data,1,nchar(unique_rna_data)-16)
for (i in 1:length(unique_rna_data)) {
  unique_rna_data[i] <- chartr('.', '-', unique_rna_data[i])
}
unique_rna_data <- unique(unique_rna_data)
#array of common patients with clinical, mutation and expression data.
common_pat <- intersect(K, unique_rna_data)


common_patient_mut <- subset(data_mut, grepl(paste(common_pat, collapse='|'), data_mut$Tumor_Sample_Barcode))


#mutation consequence (specific) frequency
mutation.consequence_freq <- as.data.frame(table(common_patient_mut$Consequence))
colnames(mutation.consequence_freq)[1] <- "mutation_consequence"


#variant classification frequency
variant_freq <- as.data.frame(table(common_patient_mut$Variant_Classification))
colnames(variant_freq)[1] <- "variant_class"

ggplot(data=variant_freq, aes(x=variant_class, y=Freq)) + 
  geom_col() + theme(axis.text.x = element_text(angle = 45,hjust=1))


#impact level frequency
impact_freq <- as.data.frame(table(common_patient_mut$IMPACT))
colnames(impact_freq)[1] <- "impact"


#mutation type (general) frequency
mutation.type_freq <- as.data.frame(table(common_patient_mut$VARIANT_CLASS))
colnames(mutation.type_freq)[1] <- "mutation_type"


ggplot(data=mutation.type_freq, aes(x=mutation_type, y=Freq))+
  geom_col(aes(fill=mutation_type))

#variant type frequency
var.type_freq <- as.data.frame(table(common_patient_mut$Variant_Type))
colnames(var.type_freq)[1] <- "variant_type"

ggplot(data=var.type_freq, aes(x=variant_type, y=Freq))+
  geom_col( aes(fill=variant_type))


#gene (HUGO) name frequency
#sample.name <- as.data.frame(table(common_patient_mut$Tumor_Sample_Barcode))
removed_silent_and_missense <- subset(common_patient_mut, common_patient_mut$Variant_Classification != "Missense_Mutation" & common_patient_mut$Variant_Classification != "Silent")
unique_rows <- !duplicated(removed_silent_and_missense[,c("Variant_Classification", "Hugo_Symbol", "Tumor_Sample_Barcode")])
removed_silent_and_missense_unique <- removed_silent_and_missense[unique_rows,]

hugo <- as.data.frame(table(removed_silent_and_missense_unique$Hugo_Symbol))
hugo.ordered <- hugo[order(-hugo$Freq),]

ggplot(data=hugo.ordered[1:15,], aes(x=Var1, y=Freq))+
  geom_col()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_x_discrete(limits = hugo.ordered[1:15,]$Var1)


##Heatmap from tutorial 4

#generate oncoplot matrix with edits  



removed <- subset(common_patient_mut, 
                  (common_patient_mut$IMPACT == "HIGH" |
                    common_patient_mut$IMPACT == "MODERATE") &
                    !(common_patient_mut$Variant_Type == "INS" & (common_patient_mut$Tumor_Seq_Allele1 == "-" & common_patient_mut$Tumor_Seq_Allele2 == "-")))

cnv_events = unique(removed$Variant_Classification)

oncomat = reshape2::dcast(
  data = removed,
  formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
  fun.aggregate = function(x, cnv = cnv_events) {
    x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit
    xad = x[x %in% cnv]
    xvc = x[!x %in% cnv]
    
    if (length(xvc) > 0) {
      xvc = ifelse(test = length(xvc) > 1,
                   yes = 'Multi_Hit',
                   no = xvc)
    }
    
    
    x = ifelse(
      test = length(xad) > 0,
      yes = paste(xad, xvc, sep = ';'),
      no = xvc
    )
    x = gsub(pattern = ';$',
             replacement = '',
             x = x)
    x = gsub(pattern = '^;',
             replacement = '',
             x = x)
    return(x)
  },
  value.var = 'Variant_Classification',
  fill = '',
  drop = FALSE
)
rownames(oncomat) = oncomat$Hugo_Symbol
oncomat <- oncomat[,-1]

# Get the number of non-empty string columns for each row
num_non_empty_cols <- data.frame(row.names(oncomat), integer(nrow(oncomat)))
colnames(num_non_empty_cols) <- c("HUGO", "num_nonzero_col")


num_non_empty_cols[, 2] <- apply(oncomat != "", 1, sum)
# Order the rows based on the number of non-empty string columns
ordered_rows <- num_non_empty_cols[order(-num_non_empty_cols$num_nonzero_col),]

oncomat.ordered <- oncomat[ordered_rows$HUGO,]

#oncomat.ordered <- oncomat[order(-hugo$Freq),] #this orders according to hugo, but hugo is different here



#transform into binary matrix

mat <- oncomat.ordered
mat[mat!=""]=1
mat[mat==""]=0
mat <- apply(mat, 2 ,as.numeric)
mat <- as.matrix(mat)
rownames(mat)  <-  row.names(oncomat.ordered)

reduce.mat <- mat[1:10,]
res <- pheatmap(reduce.mat,
                cluster_rows = F,
                show_colnames=FALSE)

cluster <-  as.data.frame(cutree(res$tree_col, k = 4))
colnames(cluster)[1] <- 'clusternum'
row.names(cluster) <- substr(row.names(cluster), 1, nchar(row.names(cluster)) - 3)


clin_patient_w_data <- subset(clinical_data, grepl(paste(common_pat, collapse='|'), clinical_data$PATIENT_ID))
clin_patient_w_data$clusternum <- cluster$clusternum
#remove patients without DSS
pat_w_dss <- subset(clin_patient_w_data, clin_patient_w_data$DSS_STATUS != "")
#converts to 0,1
pat_w_dss$DSS_STATUS <- as.numeric(substr(pat_w_dss$DSS_STATUS, 1, 1))


#subset so only include data that has subtype (is not empty)
survival_disease_subtype <- subset(pat_w_dss, pat_w_dss$SUBTYPE != "")
  

#remove BRCA_Normal (not enough data) norm gene?
survival_disease_subtype_no_norm <- subset(survival_disease_subtype, survival_disease_subtype$SUBTYPE != "BRCA_Normal")




#plotting KM plot
#fit, do not change the OS_months or DSS status. just change the "subtype" part to the name of the column with the mutation data that you added
fit = survfit(Surv(OS_MONTHS, DSS_STATUS) ~ clusternum, data=survival_disease_subtype_no_norm)
print(fit)

#plot 
ggsurvplot(fit, data=survival_disease_subtype_no_norm, pval=T, risk.table=T, risk.table.col="strata",
           risk.table.height=0.35)



cars <- mtcars











#maftools

laml = read.maf(maf = "testdata.txt")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
#oncoplot(maf = lapl, top = 10, sepwd_genes = 10000000000000, sepwd_samples = 5000000000)
#wayy too much data to use oncoplot with, ask TA


#maftools removes all silent mutations


PIK3CA <- subset(common_patient_mut, common_patient_mut$Hugo_Symbol == "PIK3CA")
unique_rows <- !duplicated(PIK3CA[,c("Variant_Classification", "Tumor_Sample_Barcode")])
unique_PIK3CA <- PIK3CA[unique_rows,]
sum(unique_PIK3CA$Variant_Classification != c("Silent"))


test <- subset(clin_patient_w_data, clin_patient_w_data$DSS_STATUS == "0:ALIVE OR DEAD TUMOR FREE" | clin_patient_w_data$DSS_STATUS == "1:DEAD WITH TUMOR")

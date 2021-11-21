####Score our signature

# patient<-"..." #this usually gets read in from the DE script that I have in the RNA pipeline so you will need to modify this to fit your code to read in the patient id
# tpm<-read.delim("<..._exp.genes.results>",sep="\t",header=T,stringsAsFactors=F) ### read in *_exp.genes.result file

##########################
# Modified by Rachel

library(dplyr)
library(tibble)

# Fetch command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript patientIPASS.R [patient_id] [input_tsv] [gene_col] [expr_col]", call.=FALSE)
}

patient <- args[1]
input_filename <- args[2]
gene_col <- args[3]
expr_col <- args[4]

# If input file doesn't exist, abort with message
if (!(file.exists(input_filename))) {
  stop("Input file does not exist", call.=FALSE)
}

tpm <- read.delim(input_filename, sep="\t", header=T, stringsAsFactors=F)
ipass<-read.delim("/app/IPASS_Genes.txt", sep="\t", header=T, stringsAsFactors=F)

##########################
tpm$gene = tpm[,gene_col]
# rownames(tpm)<-tpm$gene
tpm <- tpm %>% remove_rownames %>% column_to_rownames(var=gene_col)


index<-which(tpm$gene %in% ipass$gene)

subTpm<-tpm[index,]
expTpm<-as.matrix(subTpm[,expr_col])
colnames(expTpm)<-patient
rownames(expTpm)<-rownames(subTpm)
expTpm[expTpm==0]<-0.0001

#calculate score
logexpScore<-log(expTpm,base=10)
Score<-colMeans(logexpScore)
expScore<-rbind(expTpm,Score)
texpScore<-as.data.frame(t(expScore))

#normalise PaedSigScore between -1 and 1

minScore<--1.8 #from the cohort that developed the IPASS
maxScore<-2 #from the chohort that developed the IPASS
scaledScore<-(-1)+((texpScore$Score-minScore)*(1-(-1)))/(maxScore-minScore)
texpScore$IPASS<-scaledScore

filename<-paste0(patient,"_IPASS.txt")
write.table(texpScore,filename,sep="\t",row.names=F,quote=F)

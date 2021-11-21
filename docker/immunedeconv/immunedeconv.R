#################################################################
## Script to run immunedeconv on RSEM gene expression data     ##
##                                                             ##
## Inputs:                                                     ##
##   - RSEM expression TSV (gene level)                        ##
## Outputs:                                                    ##
##   - Immunedeconv summary                                    ##
#################################################################

library(tidyverse)
library(immunedeconv)

# Fetch command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript immunedeconv.R [patient_id] [input_tsv] [gene_col] [expr_col]", call.=FALSE)
}

# Get input file name from arguments
patient_id <- args[1]
input_file <- args[2]
gene_col <- args[3]
expr_col <- args[4]

# If input file doesn't exist, abort with message
if (!(file.exists(input_file))) {
  stop("Input file does not exist", call.=FALSE)
}

# Expression RSEM file to dataframe
expr_df <- read.csv(input_file, header=TRUE, sep='\t')
print(head(expr_df))
edited_df <- select(expr_df, gene_col, expr_col)

edited_df <- edited_df %>% remove_rownames %>% column_to_rownames(var=gene_col)
colnames(edited_df) <- c('patient')


res_quantiseq <- immunedeconv::deconvolute(edited_df, "quantiseq", tumor = T, scale_mrna = T)
colnames(res_quantiseq) <- c("cell_type", "cell_fraction")
# res <- res_quantiseq %>% immunedeconv::map_result_to_celltypes(c("T cell CD4+"), "quantiseq")
knitr::kable(res_quantiseq, digits=2)
# knitr::kable(list(res, res_quantiseq), digits=2)


# res_mcp <- immunedeconv::deconvolute(edited_df, "mcp_counter")
# knitr::kable(res_mcp, digits=2)

# Requires at least 2 samples
# res_timer <- immunedeconv::deconvolute(edited_df, "timer", indications=c('sarc'))
# knitr::kable(res_timer, digits=2)

# Should scaling be enabled?
res_epic <- immunedeconv::deconvolute(edited_df, "epic", tumor = T, scale_mrna = T)
colnames(res_epic) <- c("cell_type", "cell_fraction")
knitr::kable(res_epic, digits=2)


# Write results to separate TSV files
write.table(res_quantiseq, file=paste(patient_id, '_quantiseq.tsv', sep=""), row.names=FALSE, sep="\t", quote=FALSE)
write.table(res_epic, file=paste(patient_id, '_epic.tsv', sep=""), row.names=FALSE, sep="\t", quote=FALSE)
# write.table(res_mcp, file=paste(patient_id, '_mcp-counter.tsv', sep=""), row.names=FALSE, sep="\t", quote=FALSE)


# Mapping cell types (is this necessary??? Loss of info)
# ctm <- immunedeconv::cell_type_map
# immunedeconv::map_result_to_celltypes(res_quantiseq, c('T cell CD8+', 'T cell CD4+', 'B cell', 'Monocyte', 'Myeloid dendritic cell', 'NK cell', 'Neutrophil', 'Mast cell', 'Eosinophil', 'Basophil', 'other cell', 'cytotoxicity score', 'immune score', 'microenvironment score', 'stroma score'))
# immunedeconv::map_result_to_celltypes(res_epic, c('T cell CD8+', 'T cell CD4+', 'B cell', 'Monocyte', 'Myeloid dendritic cell', 'NK cell', 'Neutrophil', 'Mast cell', 'Eosinophil', 'Basophil', 'other cell', 'cytotoxicity score', 'immune score', 'microenvironment score', 'stroma score'))
# immunedeconv::map_result_to_celltypes(res_mcp, c('T cell CD8+', 'T cell CD4+', 'B cell', 'Monocyte', 'Myeloid dendritic cell', 'NK cell', 'Neutrophil', 'Mast cell', 'Eosinophil', 'Basophil', 'other cell', 'cytotoxicity score', 'immune score', 'microenvironment score', 'stroma score'))


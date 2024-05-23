rm(list=ls())
library(ontologyIndex)
library(ontologySimilarity)
args = commandArgs(trailingOnly=TRUE)

id=args[1]
code=args[2]
refg=args[3]
#dat <- read.csv('/run/data_dependencies/omim_annotate/HGMD_feature_20190510.csv')
#dat <- dat[!is.na(dat$hpo_id),]
dat <- read.csv(paste0('/run/data_dependencies/omim_annotate/',refg,'/HGMD_phen.tsv'),
                sep='\t')

library(dplyr)
get_HPO_list <- function(df1){
  df2 <- df1 %>%
    dplyr::group_by(acc_num, phen_id, gene_sym) %>%
    dplyr::summarise(HPO = paste0(hpo_id, collapse = " ")) %>%
    dplyr::ungroup() %>% as.data.frame()
  df2$HPO_list <- lapply(df2$HPO, function(i) unlist(strsplit(i, " ")))
  #df2$Gene <- lapply(df2$Gene, function(i) unique(strsplit(i, ",")[[1]]))
  #rownames(df2) <- ifelse(is.na(df2$OMIM_ID), df2$Disease_Name, df2$OMIM_ID)
  return(df2)
}


# Load HPO_obo
HPO_obo <- get_OBO("/run/data_dependencies/omim_annotate/hp.obo",propagate_relationships = c('is_a','part_of'), extract_tags = "minimal")

# set simi_thresh
simi_thresh <- 0
# Set pwd:
out_pwd <- "./"


# In public release, there might be empty HGMD phenotype file
if (dim(dat)[1] == 0) {
  col_names <- c('acc_num', 'phen_id', 'gene_sym', 'HPO', 'HPO_list', 'Similarity_Score')
  dat2 <- data.frame(matrix(nrow = 0, ncol = length(col_names)))
  colnames(dat2)<- col_names

} else {
  dat <- dat[!is.na(dat$hpo_id),]
  dat2_ori <- get_HPO_list(dat)
  #dat2_ori <- dat2_ori[with(dat2_ori,order(hgvs)),]
  #dat2_ori$Gene <- unlist(lapply(dat2_ori$Gene, function(x) x[[1]]))
  dat2 <- dat2_ori
  
  # Load patient HPO
  HPO.file <- file.path(id)
  HPO.orig <- read.table(HPO.file, sep = "\t", fill=T, header=F, stringsAsFactors=F)
  HPO <- HPO.orig$V1

  # remove terms without a HPO ID
  HPO <- HPO[grepl("HP:", HPO)]
  HPO <- list(HPO)

  sim_mat <- get_asym_sim_grid(HPO, dat2$HPO_list, ontology=HPO_obo)
  dat2$Similarity_Score <- as.vector(sim_mat)
  dat2$HPO_list <- unlist(lapply(dat2$HPO_list, function(x) paste0(unlist(x), collapse = "|")))
  dat2 <- dat2[order(dat2$Similarity_Score, decreasing = T),]
}

output_file_name2 <- paste0("/out/",code,"-cz")
write.table(dat2, output_file_name2, sep = "\t", quote = F, row.names = F)




# Load OMIM gene-disease data
genemap2 <-readRDS("/run/data_dependencies/omim_annotate/hg19/genemap2_v2022.rds")

## ---- Load OMIM Phenotype  ----
HPO_orig <- read.table("/run/data_dependencies/omim_annotate/hg19/HPO_OMIM.tsv", sep="\t", header= T, stringsAsFactors = FALSE, comment.char = "", fill = TRUE, quote = "\"")
OMIM_HPO <- HPO_orig[,c('OMIM_ID', 'DiseaseName', 'HPO_ID')]
# rename colnames
colnames(OMIM_HPO) = c("OMIM_ID", "Disease_Name", "HPO_ID")
OMIM_HPO_cl <- unique(OMIM_HPO)


# Function to get_HPO_list 
get_HPO_list <- function(df1){
  df2 <- df1 %>% 
    dplyr::group_by(OMIM_ID, Disease_Name) %>% 
    dplyr::summarise(HPO = paste0(HPO_ID, collapse = " ")) %>% 
    dplyr::ungroup() %>% as.data.frame()
  df2$HPO_list <- lapply(df2$HPO, function(i) unlist(strsplit(i, " ")))
  #rownames(df2) <- ifelse(is.na(df2$OMIM_ID), df2$Disease_Name, df2$OMIM_ID)
  return(df2)
}
# Prepare OMIM HPO list for ontology similarity comparision
OMIM_HPO_all <- get_HPO_list(OMIM_HPO_cl[, c("OMIM_ID", "Disease_Name","HPO_ID")])

HPO.file <- file.path(id)
#HPO.file <- '/mnt/atlas_local/sasidhar/data/1127895/P-031-101_HPOTerms.txt'
HPO.orig <- read.table(HPO.file, sep = "\t", fill=T, header=F, stringsAsFactors=F)
HPO <- HPO.orig$V1
# remove terms without a HPO ID
HPO <- HPO[grepl("HP:", HPO)]
HPO <- list(HPO)

 
# Generate HPO list
#Unknown_Disease <- data.frame(Disease_Name = "Unknown",
#                              HPO_ID = HPO,
#                              OMIM_ID = 'NA')

#Unknown_Disease_all <- get_HPO_list(Unknown_Disease)

# compare disease similarity
#sim_mat <- get_asym_sim_grid(Unknown_Disease_all$HPO_list, OMIM_HPO_all$HPO_list, ontology=HPO_obo)
sim_mat <- get_asym_sim_grid(HPO, OMIM_HPO_all$HPO_list, ontology=HPO_obo)


# OMIM_HPO_all$HPO_term <- HPO_obo$name[OMIM_HPO_all$HPO]
OMIM_HPO_all$Similarity_Score <- as.vector(sim_mat)
# convert HPO ID to HPO term
OMIM_HPO_all$HPO_term <- unlist(lapply(OMIM_HPO_all$HPO_list, function(x) paste0(HPO_obo$name[unlist(x)], collapse = "|")))

## Add gene - disease relationship
OMIM_HPO_all_wGene <- merge(unique(genemap2[,c("Pheno_ID","Approved_Gene_Symbol","Ensembl_Gene_ID","Entrez_Gene_ID")]),
                            OMIM_HPO_all[,c("OMIM_ID","Disease_Name","Similarity_Score","HPO_term")],  by.y = "OMIM_ID", by.x = "Pheno_ID")
colnames(OMIM_HPO_all_wGene)[2] <- 'Gene_Symbol'
OMIM_HPO_all_order <-OMIM_HPO_all_wGene[order(OMIM_HPO_all_wGene$Similarity_Score, decreasing = T),]
# write.table(OMIM_HPO_all_order, paste0(unknown_disease_path,"HPOsimi_all_",output_file_name), sep = "\t", quote = F, row.names = F)

# OMIM_HPO_all_filt <- head(OMIM_HPO_all_order, n = No_candidate)
OMIM_HPO_all_filt <- OMIM_HPO_all_order[OMIM_HPO_all_order$Similarity_Score >= simi_thresh,]
output_file_name2 <- paste0("/out/",code,"-dx")
write.table(OMIM_HPO_all_filt, output_file_name2, sep = "\t", quote = F, row.names = F)




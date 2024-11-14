# This script is used to expand the variants in the HGMD vcf file to gDNA level. 
# (Originall written by @Dongxue Mao, reorganized by @Zhijian Yu)



# Install and load libraries ---------------------------------------------
required_packages <- c("vcfR", "stringr", "stringi", "reshape2", 
                       "tidyr", "plyr", "tidyr", "data.table", "tictoc")

missing_packages <- required_packages[!(required_packages %in% 
                                       rownames(installed.packages()))]

if(length(missing_packages)) {
  install.packages(missing_packages)
}

library(vcfR)
library(stringr)
library(stringi)
library(reshape2)
library(tidyr)
library(plyr) # could potentially be replaced by dplyr, needs testing
library("data.table")
library(tictoc)


# Prepare for data processing ---------------------------------------------

# arg1: hgmd input vcf.gz
# arg2: out file directory
# arg3: reference file: AA_abbreviations.txt
# arg4: genome version

## parse the args
args = commandArgs(trailingOnly = T)

print("Printing out arguments;")
print(paste("args1:", args[1], sep = " "))

print(paste("args2:", args[2], sep = " "))

print(paste("args3:", args[3], sep = " "))

print(paste("args4:", args[4], sep = " "))

print("")

genome_version <- args[4]

#set working directory
setwd(args[2])

#create output folders
system("mkdir -p {Expand_Result,Expand_Result_final,TransVarInput,TransVarOut}")



# Step 1 - ParseInfo_mc ----------------------------------------------------------------------------------------------------

## step 1  parse the infor in the hgmd vcf save as RDS
## step 2: dcast infor column and combine with other columns

## parse the infor in the hgmd vcf save as RDS ####
print("Runnning Step 1 --- ParseInfo_mc")
tic("step1")

# load the vcf
vcf <- read.vcfR(args[1], verbose = F)
vcf.df <- vcf@fix
vcf.df <- as.data.frame(vcf.df)


# seperate the info column to key, value
# rbind to long df
vcf_dt <- as.data.table(vcf.df) # Convert to data.table for efficiency

# Split the INFO column into fields, handling ";" within quotes
vcf_dt[, fields := stri_split_regex(INFO, ';(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)', simplify = FALSE)]

# Unnest the fields into long format
result_dt <- vcf_dt[, .(field = unlist(fields)), by = ID]

# Split each field into key and value based on the first "=" using sub
result_dt[, key := sub("=.*", "", field)]
result_dt[, value := sub("^[^=]*=", "", field)]

# Select and reorder the desired columns, including 'field'
results_df <- result_dt[, .(field, key, value, id = ID)]

# Convert to data.frame
results_df <- as.data.frame(results_df)

print("saving result")
saveRDS(results_df, file.path(args[2], "hgmd_mc.rds"))
toc()

# Step 2 - Combine info ----------------------------------------------------------------------------------------------------

print("Runnning Step 2 --- CombineInfo")
tic("step2")

## dcast infor column and combine with other columns ####
info.l <- readRDS(file.path(args[2], "hgmd_mc.rds"))
info.w <- reshape2::dcast(info.l, id ~ key, value.var = "value")

# combine with other columns
vcf <- read.vcfR(args[1], verbose = FALSE)
vcf.df <- vcf@fix
vcf.df <- as.data.frame(vcf.df)

vcf.all <- merge(vcf.df[, 1:7],
                 info.w,
                 by.x = "ID",
                 by.y = "id",
                 all = T)

saveRDS(vcf.all, file.path(args[2], "hgmd_all.rds"))
toc()


# Step 3&4 - pro2genome ----------------------------------------------------------------------------------------------------

tic("step 3 and 4")
print("Runnning Step 3&4 --- pro2genome")

# parse pro2g result from TransVar
vcf.all <- readRDS(file.path(args[2], "hgmd_all.rds"))

# keep CLASS == "DM" and "DM?"
vcf.all <- vcf.all[vcf.all$CLASS %in% c("DM", "DM?"), ]
write.table(
  vcf.all,
  file.path(args[2], "hgmd_DM.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(vcf.all, file.path(args[2], "hgmd_DM.rds"))

# load var in HGMD
hgmd.orig <- readRDS(file.path(args[2], "hgmd_DM.rds"))
hgmd <- hgmd.orig
# parse PROT column
hgmd$PROT2 <-  gsub(".[0-9]+%3A", ":", hgmd$PROT)
hgmd$PROT2 <- gsub("%3B", ";", hgmd$PROT2)
hgmd$PROT2 <- gsub("%3D", "=", hgmd$PROT2)
hgmd$PROT2 <- gsub("\\(|\\)", "", hgmd$PROT2)

## pre-filter the HGMD entries ####
# parse PROT column
hgmd.part <- hgmd[, c("ID", "GENE", "PROT", "PROT2")]

# drop fs var (no expansion for fs)
hgmd.filt <- hgmd.part[!grepl("fs", hgmd.part$PROT), ]

# drop var with multiple protein change (no expansion for them) (there are 48)
hgmd.filt <- hgmd.filt[!grepl(";", hgmd.filt$PROT2), ]

# drop synonymous var (will parse as non-coding)
hgmd.filt <- hgmd.filt[!grepl("=", hgmd.filt$PROT2), ]

# drop stop loss or stop gain (with * in the PROT)
hgmd.filt <- hgmd.filt[!grepl("\\*", hgmd.filt$PROT2), ]

# drop inframe insertions and delins (ins and delins)
hgmd.filt <- hgmd.filt[!grepl("ins", hgmd.filt$PROT2), ]
hgmd.filt <- hgmd.filt[!grepl("dup", hgmd.filt$PROT2), ]

# drop inframe deletions that affect more than 1 aa
del <- hgmd.filt[grepl("del", hgmd.filt$PROT2), ]
del <- separate(
  data = del,
  col = "PROT2",
  into = c("NP_ID", "pChange"),
  remove = F,
  sep = ":"
)
del_multi_AA <- del[grepl("_", del$pChange), ]
hgmd.filt <- hgmd.filt[!(hgmd.filt$ID %in% del_multi_AA$ID), ]

hgmd.non.coding <- hgmd[!(hgmd$ID %in% hgmd.filt$ID), ]

# check No of del with 2 AA del
del_multi_AA$pChange <- gsub("del", "", del_multi_AA$pChange)
del_multi_AA$pChange <- gsub("p\\.", "", del_multi_AA$pChange)
del_multi_AA <- separate(
  data = del_multi_AA,
  col = "pChange",
  into = c("p_start", "p_end"),
  remove = F,
  sep = "_"
)
del_multi_AA$p_start <- str_extract(del_multi_AA$p_start, "[0-9]+$")
del_multi_AA$p_end <- str_extract(del_multi_AA$p_end, "[0-9]+$")
del_multi_AA$del_plength <- as.numeric(del_multi_AA$p_end)  - as.numeric(del_multi_AA$p_start) + 1

# drop ones that are NA in PROT
hgmd.filt <- hgmd.filt[!is.na(hgmd.filt$PROT), ]

# save hgmd.filt
write.table(
  hgmd.filt,
  file.path(args[2], "Expand_Result", "hgmd.filt.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(hgmd.filt, file.path(args[2], "Expand_Result", "hgmd.filt.rds"))

## parse the single AA del ####
# Add alt of the AA (all possible change of change, except the ref)
del.1aa <- hgmd.filt[grepl("del", hgmd.filt$PROT2), ]
del.1aa <- separate(
  data = del.1aa,
  col = "PROT2",
  into = c("NP_ID", "pChange"),
  remove = F,
  sep = ":"
)
del.1aa$pChange <- gsub("del", "", del.1aa$pChange)
del.1aa$pChange <- gsub("p\\.", "", del.1aa$pChange)
del.1aa$pLoc <- str_extract(del.1aa$pChange, "[0-9]+$")
del.1aa$pRef <- str_extract(del.1aa$pChange, "^[A-z]+")

# load AA table
AA.table <- read.delim2(args[3], sep = "\t", header = T)
AA.filt <- AA.table[!grepl("\\/", AA.table$Full_Name), ]
AA.3letter <- AA.filt$ThreeLetter
AA.1letter <- AA.filt$OneLetter
# c('Ile','Phe','Lys','Tyr','Val','Glu','Met','Thr','Leu','Asp','Gly','Gln','Pro','Ser','His','Asn','Cys','Arg','Ala','Trp','Ter')

del.1aa2 <- del.1aa[rep(seq_len(nrow(del.1aa)), each = 20), ]
del.1aa2$pAlt <- AA.3letter

# remove Synonymous
del.1aa3 <- del.1aa2[del.1aa2$pRef != del.1aa2$pAlt, ]

# format the expanded AA to PROT2 format
del.1aa3$PROT2_exp <- paste0(del.1aa3$NP_ID,
                             ":p.",
                             del.1aa3$pRef,
                             del.1aa3$pLoc,
                             del.1aa3$pAlt)

##TransVar with del.1aa ####

var4transVar.Del <-  unique(del.1aa3$PROT2_exp)
write.table(
  var4transVar.Del,
  file.path(args[2], "TransVarInput", "del_var.txt"),
  quote = F,
  col.names = F,
  row.names = F
)

# Warning from TransVar [_annotate_snv_protein] warning: unknown alternative: TER, ignore alternative.


## run TransVar step 1 ####
system(
  paste(
    "conda run -n hgmd transvar panno --refseq -l TransVarInput/del_var.txt --idmap protein_id --refversion",
    genome_version,
    "> TransVarOut/out_del_1aa.txt"
  ),
  intern = FALSE
)

# load result
del2g.orig <- read.delim2(
  file.path(args[2], "TransVarOut", "out_del_1aa.txt"),
  header = T,
  stringsAsFactors = F
) # 59567-2022version 80678-2024version
del2g <- del2g.orig[del2g.orig$transcript != ".", ]

# del2g <- del2g[grepl("candidate_codons=",del2g$info),] # 50711
var <- unique(del2g$input)
failed.del2g <- del.1aa3[!(del.1aa3$PROT2_exp %in% var) &
                           del.1aa3$pAlt != "Ter", ]

length(unique(failed.del2g$ID))

# add the HGMD info
# 1. add the HGMD ID
del2g.wID <- merge(del2g,
                   del.1aa3[, c("PROT2_exp", "ID", "PROT2")],
                   by.x = "input",
                   by.y = "PROT2_exp",
                   all.x = T)

# save the result
write.table(
  del2g.wID,
  file.path(args[2], "Expand_Result", "del2g.wID.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(del2g.wID, file.path(args[2], "Expand_Result", "del2g.wID.rds"))


##TransVar with NP ####

# select var to run TransVar: non del
NP_df <- hgmd.filt[!(hgmd.filt$PROT2 %in% del.1aa3$PROT2), ]

var4transVar.p <-  unique(NP_df$PROT2)
table(grepl("\\?", var4transVar.p))

write.table(
  var4transVar.p[!grepl("\\?", var4transVar.p)],
  file.path(args[2], "TransVarInput", "NP_var.txt"),
  quote = F,
  col.names = F,
  row.names = F
)

NP_df_exp <- NP_df
NP_df_exp <- separate(
  data = NP_df_exp,
  col = "PROT2",
  into = c("NP_ID", "pChange"),
  remove = F,
  sep = ":"
)
NP_df_exp$pChange <- gsub("p\\.", "", NP_df_exp$pChange)
NP_df_exp$pAlt <- str_extract(NP_df_exp$pChange, "([A-z]|\\?)+$")
NP_df_exp$pRef <- str_extract(NP_df_exp$pChange, "^[A-z]+")
NP_df_exp$pLoc <- str_extract(NP_df_exp$pChange, "[0-9]+")

# AA in different format 1 letter or 3 letter
# parse them seperately

NP_df_exp_3aa <-  NP_df_exp[NP_df_exp$pAlt %in% AA.3letter, ]
NP_df_exp_1aa <- NP_df_exp[NP_df_exp$pAlt %in% LETTERS, ]
NP_df_exp_qMark <- NP_df_exp[NP_df_exp$pAlt == "?", ]

# 3 letters
NP_df_exp_3aa <- NP_df_exp_3aa[rep(seq_len(nrow(NP_df_exp_3aa)), each = 20), ]
NP_df_exp_3aa$pAlt <- AA.3letter
NP_df_exp_3aa <- NP_df_exp_3aa[NP_df_exp_3aa$pRef != NP_df_exp_3aa$pAlt, ]
NP_df_exp_3aa$PROT2_exp <- paste0(
  NP_df_exp_3aa$NP_ID,
  ":p.",
  NP_df_exp_3aa$pRef,
  NP_df_exp_3aa$pLoc,
  NP_df_exp_3aa$pAlt
)

# 1 letter
NP_df_exp_1aa <- NP_df_exp_1aa[rep(seq_len(nrow(NP_df_exp_1aa)), each = 20), ]
NP_df_exp_1aa$pAlt <- AA.1letter
NP_df_exp_1aa <- NP_df_exp_1aa[NP_df_exp_1aa$pRef != NP_df_exp_1aa$pAlt, ]
NP_df_exp_1aa$PROT2_exp <- paste0(
  NP_df_exp_1aa$NP_ID,
  ":p.",
  NP_df_exp_1aa$pRef,
  NP_df_exp_1aa$pLoc,
  NP_df_exp_1aa$pAlt
)

pRef <- as.data.frame(table(NP_df_exp_qMark$pRef))
pRef <- pRef$Var1
# check if all are start loss
length(setdiff(pRef, c("M", "Met"))) == 0

# rename pRef
NP_df_exp_qMark$pRef = "Met"
NP_df_exp_qMark <- NP_df_exp_qMark[rep(seq_len(nrow(NP_df_exp_qMark)), each = 19), ]
AA.nonStart = AA.3letter[AA.3letter != "Met"]

NP_df_exp_qMark$pAlt <- AA.nonStart
NP_df_exp_qMark$PROT2_exp <- paste0(
  NP_df_exp_qMark$NP_ID,
  ":p.",
  NP_df_exp_qMark$pRef,
  NP_df_exp_qMark$pLoc,
  NP_df_exp_qMark$pAlt
)

NP_df_exp_final <- rbind(NP_df_exp_3aa, NP_df_exp_1aa, NP_df_exp_qMark)

var4transVar.missense <-  unique(c(
  NP_df_exp_3aa$PROT2_exp,
  NP_df_exp_1aa$PROT2_exp,
  NP_df_exp_qMark$PROT2_exp
)) # 2186825
write.table(
  var4transVar.missense,
  file.path(args[2], "TransVarInput", "NP_missense_var.txt"),
  quote = F,
  col.names = F,
  row.names = F
)

##run TransVar step 2 ####
system(
  paste(
    "conda run -n hgmd transvar panno -l TransVarInput/NP_var.txt --refseq --idmap protein_id --refversion",
    genome_version,
    "> TransVarOut/out_NP.txt"
  ),
  intern = FALSE
)

system(
  paste(
    "conda run -n hgmd transvar panno -l TransVarInput/NP_missense_var.txt --refseq --idmap protein_id --refversion",
    genome_version,
    "> TransVarOut/out_NP_missense.txt"
  ),
  intern = FALSE
)

## load transVar NP result ####
# NP
NP2g.orig <- read.delim2(
  file.path(args[2], "TransVarOut", "out_NP.txt"),
  header = T,
  stringsAsFactors = F
) # 166544
# write.table(NP2g.orig, file.path(args[2],"TransVarOut","out_NP.txt"), col.names = T, row.names = F, quote = F, sep = "\t")

NP2g <- NP2g.orig[NP2g.orig$transcript != ".", ]
NP2g <- NP2g[grepl("candidate_codons=", NP2g$info), ] # all failed are M1?, start loss, can be rescued in NP2g_exp

# check failed
var <- unique(NP2g$input) # 164045
# failed.NP.var <- var4transVar.p[!(var4transVar.p %in% var)] # 2303

# NP exp
NP2g_exp.orig <- fread(
  file.path(args[2], "TransVarOut", "out_NP_missense.txt"),
  header = T,
  stringsAsFactors = F,
  sep = "\t"
)

NP2g_exp <- NP2g_exp.orig[NP2g_exp.orig$transcript != ".", ]
# NP2g_exp <- NP2g_exp[grepl("candidate_codons=",NP2g_exp$info),] # 2167881

# check failed
var_exp <- unique(NP2g_exp$input)

NP2g_exp.wID <- merge(
  NP2g_exp,
  NP_df_exp_final[, c("PROT2_exp", "ID", "PROT2")],
  by.x = "input",
  by.y = "PROT2_exp",
  all.x = T
)
write.table(
  NP2g_exp.wID,
  file.path(args[2], "Expand_Result", "NP2g_exp.wID.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(NP2g_exp.wID,
        file.path(args[2], "Expand_Result", "NP2g_exp.wID.rds"))


## ADD hgmd ID to the NP2g result ####
NP2g.wID <- merge(NP2g,
                  hgmd.filt,
                  by.x = "input",
                  by.y = "PROT2",
                  all.x = T)
NP2g.wID$input_type <- "hgmd.NP_ID"
NP2g.wID$PROT2 <- NP2g.wID$input

# save the result
write.table(
  NP2g.wID,
  file.path(args[2], "Expand_Result", "NP2g.wID.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(NP2g.wID, file.path(args[2], "Expand_Result", "NP2g.wID.rds"))

## with gDNA ####

## step 1: gDNA to protein change. ####
# most failed due to change in the Gene name/ transcript ID / protein location
# for those, run TransVar with gDNA -> pro change
# then run TransVar with pro change -> get the gDNA change

var <- unique(c(NP2g.wID$input, del2g.wID$PROT2, NP2g_exp.wID$PROT2))
failed.NP2g <- hgmd.filt[!(hgmd.filt$PROT2 %in% var), ]
failed.ID <- unique(failed.NP2g$ID)

# get the gDNA for the failed.NP2g (failed NP and Gene)
gDNA.failed.NP2g <- hgmd[hgmd$ID %in% failed.ID, c("ID", "CHROM", "POS", "REF", "ALT", "PROT2")]
gDNA.failed.NP2g$gDNA <- paste0(
  "chr",
  gDNA.failed.NP2g$CHROM,
  ":g.",
  gDNA.failed.NP2g$POS,
  gDNA.failed.NP2g$REF,
  ">",
  gDNA.failed.NP2g$ALT
)

write.table(
  unique(gDNA.failed.NP2g$gDNA),
  file.path(args[2], "TransVarInput", "gDNA_var.txt"),
  quote = F,
  col.names = F,
  row.names = F
)

## run TransVar step 3 ####
system(
  paste(
    "conda run -n hgmd transvar ganno -l TransVarInput/gDNA_var.txt --refseq --idmap protein_id --refversion",
    genome_version,
    "> TransVarOut/out_gDNA.txt"
  ),
  intern = FALSE
)

# extract the ref and alt protein for future mapping
gDNA.failed.NP2g.part <- gDNA.failed.NP2g[, c("ID", "gDNA", "PROT2")]
gDNA.failed.NP2g.part <- separate(
  data = gDNA.failed.NP2g.part,
  col = "PROT2",
  into = c("NP_ID", "pChange"),
  remove = F,
  sep = ":"
)
gDNA.failed.NP2g.part$pChange <- gsub("p\\.", "", gDNA.failed.NP2g.part$pChange)
gDNA.failed.NP2g.part$pAlt <- str_extract(gDNA.failed.NP2g.part$pChange, "([A-z]|\\?)+$")
gDNA.failed.NP2g.part$pRef <- str_extract(gDNA.failed.NP2g.part$pChange, "^[A-z]+")
gDNA.failed.NP2g.part$pLoc <- str_extract(gDNA.failed.NP2g.part$pChange, "[0-9]+")

# replace 3 letter AA to 1 letter AA (TransVar using 1 letter AA)
gDNA.failed.NP2g.part$pAlt <- mapvalues(gDNA.failed.NP2g.part$pAlt,
                                        from = AA.table$ThreeLetter,
                                        to = AA.table$OneLetter)
gDNA.failed.NP2g.part$pRef <- mapvalues(gDNA.failed.NP2g.part$pRef,
                                        from = AA.table$ThreeLetter,
                                        to = AA.table$OneLetter)
# consider stop: U and Ter
gDNA.failed.NP2g.part$pRef <- mapvalues(gDNA.failed.NP2g.part$pRef,
                                        from = c("U", "Ter"),
                                        to = c("*", "*"))

# load result
pDNA2pro.orig <- read.delim2(
  file.path(args[2], "TransVarOut", "out_gDNA.txt"),
  header = T,
  stringsAsFactors = F
)
pDNA2pro <- pDNA2pro.orig[grepl("p.", pDNA2pro.orig$coordinates.gDNA.cDNA.protein.), c("input", "gene", "coordinates.gDNA.cDNA.protein.", "info")]

# parse to info and extract the aliases(NP ID)
# vectorized transformation
library(dplyr)
info.df <- pDNA2pro  |>
  mutate(original_info = as.character(info)) |>
  separate_rows(info, sep = ";") |>
  rename(fields = info) |>
  separate(
    fields,
    into = c("key", "value"),
    sep = "=",
    fill = "right",
    extra = "merge",
    remove = FALSE
  ) |>
  rename(coordinates = coordinates.gDNA.cDNA.protein.)  |>
  select(fields, key, value, original_info, coordinates, input) |>
  rename(info = original_info)
library(plyr)

info.w.orig <- reshape2::dcast(
  info.df,
  formula = input + coordinates + info ~ key,
  value.var = "value",
  fill = NULL
)
info.w <- info.w.orig[, c("input", "coordinates", "aliases")]
info.w <- unique(info.w)
info.w <- info.w %>%
  separate(
    col = "coordinates",
    c("gDNA", "cDNA", "protein"),
    sep = "/",
    remove = T
  )
info.w <- info.w[, c("input", "aliases", "protein")]
colnames(info.w) <- c("input", "TransVar_NP_ID", "TransVar_pChange")
# extract TransVar_pRef, TransVar_pAlt
info.w$TransVar_pChange <- gsub("p\\.", "", info.w$TransVar_pChange)
info.w$TransVar_pAlt <- str_extract(info.w$TransVar_pChange, "([A-z]|\\?|\\*)+$")
info.w$TransVar_pRef <- str_extract(info.w$TransVar_pChange, "^[A-z]+|\\*")
info.w$TransVar_pLoc <- str_extract(info.w$TransVar_pChange, "[0-9]+")

# compare with original HGMD pRef, pAlt
pDNA2pro.wID <- merge(
  gDNA.failed.NP2g.part,
  info.w,
  by.x = "gDNA",
  by.y = "input",
  all.y = T
)
pDNA2pro.wID.all <- merge(
  gDNA.failed.NP2g.part,
  info.w,
  by.x = "gDNA",
  by.y = "input",
  all = T
)
ID.failed <- gDNA.failed.NP2g.part$ID[!(gDNA.failed.NP2g.part$gDNA %in% info.w$input)]  # 150
# keep the entries where both pRef and pAlt matches

# Some of the cases are difficult
# 1. inframe del --- parse separately
# 2. start loss --- parse separately
# 3. reversed p ref/alt from HGMD (most likely error from hgmd)  handled by the for loop bellow

pDNA2pro.wID.selected <- data.frame()

# Vectorized Transformation with dplyr
library(dplyr)

pDNA2pro.wID.selected <- pDNA2pro.wID %>%
  mutate(
    priority = case_when(
      pRef == TransVar_pRef & pAlt == TransVar_pAlt ~ 1,
      pAlt == TransVar_pRef & pRef == TransVar_pAlt ~ 2,
      pRef == TransVar_pRef &
        paste0(pAlt, pRef) == TransVar_pAlt ~ 3,
      pRef == TransVar_pRef & pRef == "M" & pLoc == "1" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(priority)) %>%
  group_by(ID) %>%
  slice(1) %>%
  select(-priority) |>
  ungroup()

# reload plyr to mask functions for compatability
library(plyr)

# pDNA2pro.failed <- gDNA.failed.NP2g.part[gDNA.failed.NP2g.part$ID %in% ID.failed,]
pDNA2pro.failed <- pDNA2pro.wID.all[pDNA2pro.wID.all$ID %in% ID.failed, ]
length(unique(pDNA2pro.failed$gDNA))
# "CD169908" "CD198588", TransVar failed to map to a pChange, therefore, no expansion can be done for this two cases.

# expand the pDNA2pro.wID
pDNA2pro.wID_exp <- pDNA2pro.wID[rep(seq_len(nrow(pDNA2pro.wID)), each = 20), ]
pDNA2pro.wID_exp$exp.pAlt <- AA.filt$OneLetter

pDNA2pro.wID_exp <- pDNA2pro.wID_exp[pDNA2pro.wID_exp$TransVar_pRef != pDNA2pro.wID_exp$exp.pAlt, ]
pDNA2pro.wID_exp$PROT2_exp <- paste0(
  pDNA2pro.wID_exp$TransVar_NP_ID,
  ":p.",
  pDNA2pro.wID_exp$TransVar_pRef,
  pDNA2pro.wID_exp$TransVar_pLoc,
  pDNA2pro.wID_exp$exp.pAlt
)

write.table(
  unique(pDNA2pro.wID_exp$PROT2_exp),
  file.path(args[2], "TransVarInput", "gDNA2pro_exp_var.txt"),
  quote = F,
  col.names = F,
  row.names = F
)

## run transvar step 4 ####
system(
  paste(
    "conda run -n hgmd transvar panno -l TransVarInput/gDNA2pro_exp_var.txt --refseq --idmap protein_id --refversion",
    genome_version,
    "> TransVarOut/out_gDNA2pro_exp.txt"
  ),
  intern = FALSE
)

## step 2: then map back to gDNA ####

# load result
gDNA2pro2g.orig <- read.delim2(
  file.path(args[2], "TransVarOut", "out_gDNA2pro_exp.txt"),
  header = T,
  stringsAsFactors = F
)
gDNA2pro2g <- gDNA2pro2g.orig[gDNA2pro2g.orig$transcript != ".", ]
gDNA2pro2g <- gDNA2pro2g[grepl("candidate_codons=", gDNA2pro2g$info), ] # same nrow as gDNA2pro2g.orig
gDNA2pro2g$idx <- rownames(gDNA2pro2g)
g2p2g.info <- gDNA2pro2g[, c("idx", "info")]

colnames(gDNA2pro2g)[1] <- "PROT2_exp"
gDNA2pro2g.wID <- merge(
  gDNA2pro2g,
  pDNA2pro.wID_exp,
  by.x = "PROT2_exp",
  by.y = "PROT2_exp",
  all = T
)

## save the gDNA2pro2g result ####
write.table(
  gDNA2pro2g.wID,
  file.path(args[2], "Expand_Result", "gDNA2pro2g.wID.txt"),
  sep = "\t",
  col.names = T,
  row.names = F,
  quote = F
)
saveRDS(gDNA2pro2g.wID,
        file.path(args[2], "Expand_Result", "gDNA2pro2g.wID.rds"))
toc()


# Step 5 - gDNA2VCF ----------------------------------------------------------------------------------------------------
# step5 extract the gDNA chagne and convert to VCF formats

print("Runnning Step 5 --- gDNA2VCF")
tic("step 5")

g.built <- genome_version
del.orig <- readRDS(file.path(args[2], "Expand_Result", "del2g.wID.rds"))
NP_exp.orig <- readRDS(file.path(args[2], "Expand_Result", "NP2g_exp.wID.rds"))
gDNA.orig <- readRDS(file.path(args[2], "Expand_Result", "gDNA2pro2g.wID.rds"))

# hgmd.orig <- readRDS("hgmd.filt.rds")
hgmd.orig <- readRDS(file.path(args[2], "hgmd_DM.rds"))
hgmd <- hgmd.orig[, c("ID", "RANKSCORE", "CLASS")]

# function to parse the infor column
# parse to info and extract the all the possible gDNA changes
library(dplyr) #use dplyr instead of plyr here
parse_info <- function(input.df) {
  input.df$idx <- as.numeric(rownames(input.df))
  infor.wIdx <- input.df[, c("idx", "info")]
  
  # Vectorized Transformation
  info.df <- infor.wIdx %>%
    separate_rows(info, sep = ";") %>%
    separate(info,
             into = c("key", "value"),
             sep = "=",
             fill = "right") %>%
    select(idx, key, value)
  
  info.df$idx <- as.numeric(info.df$idx)
  info.df.w <- reshape2::dcast(
    info.df,
    formula = idx ~ key,
    value.var = "value",
    fill = NULL
  )
  
  if (sum(info.df.w$idx !=  infor.wIdx$idx) == 0) {
    output.df <- cbind(input.df, info.df.w)
    output.df$idx <- NULL
    output.df$coordinates.gDNA <- gsub("\\/.*", "", output.df$coordinates.gDNA.cDNA.protein.)
    
    return(output.df)
  } else {
    print("error")
  }
  
}

# parse the info
del.parse.info <- parse_info(del.orig)
saveRDS(del.parse.info,
        file.path(args[2], "Expand_Result", "Step5_del.parse.info.rds"))


NP_exp.orig <- as.data.table(NP_exp.orig)
colnames(NP_exp.orig)[5] <- "coordinates.gDNA.cDNA.protein."
NP_exp.parse.info <- parse_info(NP_exp.orig)
saveRDS(
  NP_exp.parse.info,
  file.path(args[2], "Expand_Result", "Step5_NP_exp.parse.info.rds")
)

gDNA.parse.info <- parse_info(gDNA.orig)
saveRDS(
  gDNA.parse.info, 
  file.path(args[2], "Expand_Result","Step5_gDNA.parse.info.rds")
)

del.parse.info <- readRDS(file.path(args[2], "Expand_Result", "Step5_del.parse.info.rds"))
NP_exp.parse.info <- readRDS(file.path(args[2], "Expand_Result", "Step5_NP_exp.parse.info.rds"))
gDNA.parse.info <- readRDS(file.path(args[2], "Expand_Result", "Step5_gDNA.parse.info.rds"))


library(plyr) #go back with plyr
if (T) {
  del.failed.idx <- is.na(del.parse.info$candidate_mnv_variants) &
    is.na(del.parse.info$candidate_snv_variants) &
    is.na(del.parse.info$coordinates.gDNA)
  del.failed <- del.parse.info[del.failed.idx, ]
  del.filt <- del.parse.info[!del.failed.idx, ]
  
  del.gDNA <- del.filt[, c(
    "ID",
    "PROT2",
    "input",
    "candidate_mnv_variants",
    "candidate_snv_variants",
    "coordinates.gDNA"
  )]
  
  # parse the candidate variants
  
  del.gDNA.mnv <- del.filt[!is.na(del.filt$candidate_mnv_variants), c("ID", "PROT2", "input", "candidate_mnv_variants")]
  del.gDNA.snv <- del.filt[!is.na(del.filt$candidate_snv_variants), c("ID", "PROT2", "input", "candidate_snv_variants")]
  del.gDNA.orig <- del.filt[!is.na(del.filt$coordinates.gDNA), c("ID", "PROT2", "input", "coordinates.gDNA")]
  # colnames(del.gDNA.mnv) <- colnames(del.gDNA.snv) <- colnames(del.gDNA.orig) <- c("ID","PROT2_exp","input","candidate_variants")
  colnames(del.gDNA.mnv) <- colnames(del.gDNA.snv) <- colnames(del.gDNA.orig) <- c("ID", "PROT2", "PROT2_exp", "candidate_variants") # edits Jan 24, 2023
  del.gDNA <- rbind(del.gDNA.snv, del.gDNA.orig)
  
  del.gDNA <- del.gDNA %>%
    mutate(candidate_variants = strsplit(candidate_variants, ",")) %>%
    unnest(candidate_variants) %>%
    as.data.frame()
  
  NP_exp.failed.idx <- is.na(NP_exp.parse.info$candidate_mnv_variants) &
    is.na(NP_exp.parse.info$candidate_snv_variants) &
    is.na(NP_exp.parse.info$coordinates.gDNA)
  NP_exp.failed <- NP_exp.parse.info[NP_exp.failed.idx, ]
  NP_exp.filt <- NP_exp.parse.info[!NP_exp.failed.idx, ]
  
  NP_exp.gDNA <- NP_exp.filt[, c(
    "ID",
    "PROT2",
    "input",
    "candidate_mnv_variants",
    "candidate_snv_variants",
    "coordinates.gDNA"
  )]
  
  # parse the candidate variants
  
  NP_exp.gDNA.mnv <- NP_exp.filt[!is.na(NP_exp.filt$candidate_mnv_variants), c("ID", "PROT2", "input", "candidate_mnv_variants")]
  NP_exp.gDNA.snv <- NP_exp.filt[!is.na(NP_exp.filt$candidate_snv_variants), c("ID", "PROT2", "input", "candidate_snv_variants")]
  NP_exp.gDNA.orig <- NP_exp.filt[!is.na(NP_exp.filt$coordinates.gDNA), c("ID", "PROT2", "input", "coordinates.gDNA")]
  colnames(NP_exp.gDNA.mnv) <- colnames(NP_exp.gDNA.snv) <- colnames(NP_exp.gDNA.orig) <- c("ID", "PROT2", "PROT2_exp", "candidate_variants")
  NP_exp.gDNA <- rbind(NP_exp.gDNA.snv, NP_exp.gDNA.mnv, NP_exp.gDNA.orig)
  
  NP_exp.gDNA <- NP_exp.gDNA %>%
    mutate(candidate_variants = strsplit(candidate_variants, ",")) %>%
    unnest(candidate_variants) %>%
    as.data.frame()
  
  ## gDNA ####
  
  gDNA.failed.idx <- is.na(gDNA.parse.info$candidate_mnv_variants) &
    is.na(gDNA.parse.info$candidate_snv_variants) &
    is.na(gDNA.parse.info$coordinates.gDNA)
  gDNA.failed <- gDNA.parse.info[gDNA.failed.idx, ]
  gDNA.filt <- gDNA.parse.info[!gDNA.failed.idx, ]
  
  gDNA.gDNA <- gDNA.filt[, c(
    "ID",
    "PROT2",
    "PROT2_exp",
    "candidate_mnv_variants",
    "candidate_snv_variants",
    "coordinates.gDNA"
  )]
  # parse the candidate variants
  
  gDNA.gDNA.mnv <- gDNA.filt[!is.na(gDNA.filt$candidate_mnv_variants), c("ID", "PROT2", "PROT2_exp", "candidate_mnv_variants")]
  gDNA.gDNA.snv <- gDNA.filt[!is.na(gDNA.filt$candidate_snv_variants), c("ID", "PROT2", "PROT2_exp", "candidate_snv_variants")]
  gDNA.gDNA.orig <- gDNA.filt[!is.na(gDNA.filt$coordinates.gDNA), c("ID", "PROT2", "PROT2_exp", "coordinates.gDNA")]
  colnames(gDNA.gDNA.mnv) <- colnames(gDNA.gDNA.snv) <- colnames(gDNA.gDNA.orig) <- c("ID", "PROT2", "PROT2_exp", "candidate_variants")
  gDNA.gDNA <- rbind(gDNA.gDNA.snv, gDNA.gDNA.mnv, gDNA.gDNA.orig)
  
  gDNA.gDNA <- gDNA.gDNA %>%
    mutate(candidate_variants = strsplit(candidate_variants, ",")) %>%
    unnest(candidate_variants) %>%
    as.data.frame()
  
  ## Try convert gDNA to VCF format with code ####
  # from del
  del.gDNA$gChr <- str_extract(del.gDNA$candidate_variants,
                               "^chr([:digit:]+)|^(chr[A-Z])")
  del.gDNA$gChr <- gsub("chr", "", del.gDNA$gChr)
  del.gDNA$gLoc <- str_extract(del.gDNA$candidate_variants, "g.([:digit:]+)")
  del.gDNA$gLoc <- gsub("g.", "", del.gDNA$gLoc)
  # del.gDNA$gRef <- str_extract(del.gDNA$candidate_variants, "g.(([:digit:]|\\_)+|[A-Z]+)")
  del.gDNA$gRef_tmp <- str_extract(del.gDNA$candidate_variants, "(([A-z]|\\>)+)$")
  del.gDNA$gRef <- str_extract(del.gDNA$gRef_tmp, "([A-Z]+)")
  del.gDNA$gRef_tmp <- NULL
  del.gDNA$gAlt <- str_extract(del.gDNA$candidate_variants, "([A-Z])+$")
  del.gDNA$HGMD_exp_type <- "Del_to_Missense"
  
  # from gDNA
  gDNA.gDNA$gChr <- str_extract(gDNA.gDNA$candidate_variants,
                                "^chr([:digit:]+)|^(chr[A-Z])")
  gDNA.gDNA$gChr <- gsub("chr", "", gDNA.gDNA$gChr)
  gDNA.gDNA$gLoc <- str_extract(gDNA.gDNA$candidate_variants, "g.([:digit:]+)")
  gDNA.gDNA$gLoc <- gsub("g.", "", gDNA.gDNA$gLoc)
  # gDNA.gDNA$gRef <- str_extract(gDNA.gDNA$candidate_variants, "g.(([:digit:]|\\_)+|[A-Z]+)")
  gDNA.gDNA$gRef_tmp <- str_extract(gDNA.gDNA$candidate_variants, "(([A-z]|\\>)+)$")
  gDNA.gDNA$gRef <- str_extract(gDNA.gDNA$gRef_tmp, "([A-Z]+)")
  gDNA.gDNA$gRef_tmp <- NULL
  gDNA.gDNA$gAlt <- str_extract(gDNA.gDNA$candidate_variants, "([A-Z])+$")
  
  # add the HGMD_exp_type
  PROT2_Alt <- str_extract(gDNA.gDNA$PROT2, "([A-z])+$")
  PROT2_exp_Alt <- str_extract(gDNA.gDNA$PROT2_exp, "([A-z])+$")
  gDNA.gDNA$HGMD_exp_type[PROT2_Alt == PROT2_exp_Alt] <- "Same_pChange"
  gDNA.gDNA$HGMD_exp_type[PROT2_Alt != PROT2_exp_Alt] <- "Different_pChange"
  gDNA.gDNA$HGMD_exp_type[grepl("\\?", gDNA.gDNA$PROT2)] <- "Start_Loss"
  gDNA.gDNA$HGMD_exp_type[grepl("\\*", gDNA.gDNA$PROT2_exp)] <- "Stop_Loss"
  
  table(gDNA.gDNA$HGMD_exp_type)
  
  # from NP_exp.gDNA
  NP_exp.gDNA$gChr <- str_extract(NP_exp.gDNA$candidate_variants,
                                  "^chr([:digit:]+)|^(chr[A-Z])")
  NP_exp.gDNA$gChr <- gsub("chr", "", NP_exp.gDNA$gChr)
  NP_exp.gDNA$gLoc <- str_extract(NP_exp.gDNA$candidate_variants, "g.([:digit:]+)")
  NP_exp.gDNA$gLoc <- gsub("g.", "", NP_exp.gDNA$gLoc)
  # NP_exp.gDNA$gRef <- str_extract(NP_exp.gDNA$candidate_variants, "g.(([:digit:]|\\_)+|[A-Z]+)")
  NP_exp.gDNA$gRef_tmp <- str_extract(NP_exp.gDNA$candidate_variants, "(([A-z]|\\>)+)$")
  NP_exp.gDNA$gRef <- str_extract(NP_exp.gDNA$gRef_tmp, "([A-Z]+)")
  NP_exp.gDNA$gRef_tmp <- NULL
  NP_exp.gDNA$gAlt <- str_extract(NP_exp.gDNA$candidate_variants, "([A-Z])+$")
  
  # add the HGMD_exp_type
  PROT2_Alt <- str_extract(NP_exp.gDNA$PROT2, "([A-z])+$")
  PROT2_exp_Alt <- str_extract(NP_exp.gDNA$PROT2_exp, "([A-z])+$")
  # check if the length of Alt, exp_Alt are the same
  PROT2_Alt_nchar <- nchar(PROT2_Alt)
  PROT2_exp_Alt_nchar <- nchar(PROT2_exp_Alt)
  table(PROT2_Alt_nchar == PROT2_exp_Alt_nchar) # all T
  
  NP_exp.gDNA$HGMD_exp_type[PROT2_Alt == PROT2_exp_Alt] <- "Same_pChange"
  NP_exp.gDNA$HGMD_exp_type[PROT2_Alt != PROT2_exp_Alt] <- "Different_pChange"
  NP_exp.gDNA$HGMD_exp_type[grepl("\\?", NP_exp.gDNA$PROT2)] <- "Start_Loss"
  NP_exp.gDNA$HGMD_exp_type[grepl("\\*", NP_exp.gDNA$PROT2_exp)] <- "Stop_Loss"
  
  table(is.na(NP_exp.gDNA$HGMD_exp_type))
  
  ## add the RankScore ####
  HGMD_exp.all <- rbind(del.gDNA, NP_exp.gDNA, gDNA.gDNA)
  
  HGMD_exp.all <- merge(HGMD_exp.all,
                        hgmd,
                        by.x = "ID",
                        by.y = "ID",
                        all.x = T)
  saveRDS(HGMD_exp.all,
          file.path(args[2], "Expand_Result_final", "hgmd.Coding.rds"))
  write.table(
    HGMD_exp.all,
    file.path(args[2], "Expand_Result_final", "hgmd.Coding.txt"),
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
  )
  
}

toc()

# Step 6 - nonCoding ----------------------------------------------------------------------------------------------------
print("Runnning Step 6 --- nonCoding")
tic("step 6")

hgmd.orig <- readRDS(file.path(args[2], "hgmd_DM.rds"))
hgmd.filt <- readRDS(file.path(args[2], "Expand_Result", "hgmd.filt.rds"))

hgmd <- hgmd.orig
hgmd$PROT2 <-  gsub(".[0-9]+%3A", ":", hgmd$PROT)
hgmd$PROT2 <- gsub("%3B", ";", hgmd$PROT2)
hgmd$PROT2 <- gsub("%3D", "=", hgmd$PROT2)
hgmd$PROT2 <- gsub("\\(|\\)", "", hgmd$PROT2)

## pre-filter the HGMD entries ####

# parse PROT column
hgmd <- hgmd[, c("ID", "CHROM", "POS", "REF", "ALT", "GENE", "PROT", "PROT2")]

## var to drop ####
# frame shift (drop)
hgmd$fs <- grepl("fs", hgmd$PROT2)

# drop var with multiple protein change (no expansion for them) (there are 48)
hgmd$multi_pChange <- grepl(";", hgmd$PROT2)

# drop stop loss or stop gain (with * in the PROT)
hgmd$StopLossGain <- grepl("\\*", hgmd$PROT2)

# drop inframe insertions and delins (ins and delins), keep del
hgmd$ins <- grepl("ins", hgmd$PROT2)
hgmd$dup <- grepl("dup", hgmd$PROT2)

# drop inframe deletions that affect more than 1 aa
# del2 <- hgmd.filt[grepl("del", hgmd.filt$PROT2),]
del <- hgmd[grepl("del", hgmd$PROT2), ]
del <- separate(
  data = del,
  col = "PROT2",
  into = c("NP_ID", "pChange"),
  remove = F,
  sep = ":"
)
del_multi_AA <- del[grepl("_", del$pChange), ] # dim 2319
hgmd$del_multi_AA <- hgmd$ID %in% del_multi_AA$ID

# drop ones that are NA in PROT
hgmd$NA_PROT2 <- is.na(hgmd$PROT2)

#  synonymous var (will parse as non-coding)
hgmd$synonymous <- grepl("=", hgmd$PROT2)

#  hgmd drop
drop_for_pExp <- hgmd$fs + hgmd$multi_pChange + hgmd$StopLossGain + hgmd$ins + hgmd$dup + hgmd$del_multi_AA + hgmd$synonymous + hgmd$NA_PROT2
table(drop_for_pExp)
# same number as hgmd.filt, confirmed

# for non-coding
# first, not in hgmd.filt
# then drop hgmd$fs + hgmd$multi_pChange + hgmd$StopLossGain + hgmd$ins + hgmd$dup + hgmd$del_multi_AA
drop_for_noncoding <- hgmd$fs + hgmd$multi_pChange + hgmd$StopLossGain + hgmd$ins + hgmd$dup + hgmd$del_multi_AA + (hgmd$ID %in% hgmd.filt$ID)
table(drop_for_noncoding) # 32165 remaining for non-coding

hgmd.non.coding <- hgmd[drop_for_noncoding == 0, ]

hgmd.non.coding$RefLength <- nchar(hgmd.non.coding$REF)
hgmd.non.coding$AltLength <- nchar(hgmd.non.coding$ALT)
hgmd.non.coding$delLength <- hgmd.non.coding$RefLength - hgmd.non.coding$AltLength

hgmd.non.coding$ExpPos_start <- as.numeric(hgmd.non.coding$POS) - 2
hgmd.non.coding$ExpPos_stop <- as.numeric(hgmd.non.coding$POS) + hgmd.non.coding$RefLength + 2
hgmd.non.coding.part <- hgmd.non.coding[, c(
  "ID",
  "CHROM",
  "POS",
  "REF",
  "ALT",
  "GENE",
  "PROT2",
  "RefLength",
  "AltLength",
  "delLength",
  "ExpPos_start",
  "ExpPos_stop"
)]

Ref.length.df <- as.data.frame(table(hgmd.non.coding.part$RefLength[hgmd.non.coding.part$RefLength <= 10]))
Ref.length.df$Var1 <- as.character(Ref.length.df$Var1)
Ref.length.df <- rbind(Ref.length.df, c(">10", sum(hgmd.non.coding.part$RefLength > 10)))

del.length.df <- as.data.frame(table(hgmd.non.coding.part$delLength[hgmd.non.coding.part$delLength <= 10 &
                                                                      hgmd.non.coding.part$delLength >= 0]))
del.length.df$Var1 <- as.character(del.length.df$Var1)
del.length.df <- rbind(del.length.df, c(">10", sum(hgmd.non.coding.part$delLength > 10)))
del.length.df <- rbind(del.length.df, c("<0", sum(hgmd.non.coding.part$delLength < 0)))

hgmd.non.coding.select <- hgmd.non.coding.part[hgmd.non.coding.part$RefLength <= 10 &
                                                 hgmd.non.coding.part$delLength >= 0, ]

hgmd.non.coding.all <- merge(
  hgmd.non.coding.select,
  hgmd.orig[, c("ID", "RANKSCORE")],
  by.x = "ID",
  by.y = "ID",
  all.x = T
)
hgmd.non.coding.all$HGMD_Exp <- "nonCoding"
hgmd.non.coding.final <- hgmd.non.coding.all[, c(
  "CHROM",
  "ExpPos_start",
  "ExpPos_stop",
  "RANKSCORE",
  "HGMD_Exp",
  "RefLength",
  "delLength"
)]

saveRDS(
  hgmd.non.coding.final,
  file.path(args[2], "Expand_Result_final", "hgmd.nonCoding.rds")
)
write.table(
  hgmd.non.coding.final,
  file.path(args[2], "Expand_Result_final", "hgmd.nonCoding.txt"),
  sep = "\t",
  col.names = T,
  quote = F,
  row.names = F
)

saveRDS(hgmd.non.coding.all,
        file.path(args[2], "Expand_Result", "hgmd.nonCoding.all.rds"))
write.table(
  hgmd.non.coding.all,
  file.path(args[2], "Expand_Result", "hgmd.nonCoding.all.txt"),
  sep = "\t",
  col.names = T,
  quote = F,
  row.names = F
)


# Step 7 dropped
# Step 8 - reformat4AI ----------------------------------------------------------------------------------------------------
print("Runnning Step 8 (skipping step7) --- reformat4AI")
tic("step 8")
# reformat expansion
Coding <- fread(file.path(args[2], "Expand_Result_final", "hgmd.Coding.txt"))
nonCoding <- fread(file.path(args[2], "Expand_Result_final", "hgmd.nonCoding.txt"))

Coding <- Coding[, c("ID",
                     "HGMD_exp_type",
                     "gAlt",
                     "gRef",
                     "RANKSCORE",
                     "CLASS",
                     "gLoc",
                     "gChr")]
colnames(Coding) <- c("hgmdID",
                      "c_HGMD_Exp",
                      "alt",
                      "ref",
                      'c_RANKSCORE',
                      'CLASS',
                      'new_pos',
                      "new_chr")

Coding$new_chr[Coding$new_chr == "X"] <- "23"
Coding$new_chr[Coding$new_chr == "Y"] <- "24"



## One potential issue!!!!!!: it keeps col.names = T when writing the output, which could cause the issue in reading this tsv file in R; However, it doesn't affect results when testing with pd in Python
write.table(
  Coding,
  file.path(args[2], "Expand_Result_final", "hgmd_c.tsv"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)

nonCoding <- nonCoding[, c(
  "RANKSCORE",
  "HGMD_Exp",
  "RefLength",
  "delLength",
  "ExpPos_start",
  "ExpPos_stop",
  "CHROM"
)]
colnames(nonCoding) <- c(
  "nc_RANKSCORE",
  'nc_HGMD_Exp',
  'hgmd_refLength',
  'hgmd_delLength',
  "new_start",
  'new_stop',
  "new_chr"
)

nonCoding$new_chr[nonCoding$new_chr == "X"] <- "23"
nonCoding$new_chr[nonCoding$new_chr == "Y"] <- "24"

## One potential issue!!!!!!: it keeps col.names = T when writing the output, which could cause the issue in reading this tsv file in R; However, it doesn't affect results when testing with pd in Python
write.table(
  nonCoding,
  file.path(args[2], "Expand_Result_final", "hgmd_nc.tsv"),
  sep = "\t",
  col.names = T,
  row.names = T,
  quote = F
)
toc()
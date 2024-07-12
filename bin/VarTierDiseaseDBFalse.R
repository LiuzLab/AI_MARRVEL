#!/usr/bin/env Rscript

library(tidyr)
library(readr)
library(tibble)
library(dplyr)
library(data.table)
# set adj.k value to adjust Tier
adj.k <- 1.05

# set if diseases DB will be included

args <- commandArgs(trailingOnly = TRUE)
refg <- args[1]

cat("STEP0. load database: OMIM inheritance ####\n")
omim_inheritance <- file.path(paste0("var_tier/", refg, "/genemap2.Inh.F.txt"))
in.f <- "scores.csv"
out.f <- "Tier.v2.tsv"

genemap2.Inh.F <- readr::read_tsv(omim_inheritance)

# load the input ####
anno.columns <- c(
  "varId_dash",
  "zyg",
  "geneSymbol",
  "geneEnsId",
  "gnomadAF",
  "gnomadAFg",
  "omimSymptomSimScore",
  "hgmdSymptomSimScore",
  "IMPACT",
  "Consequence",
  "hgmdVarFound",
  "clinvarSignDesc",
  "spliceAImax"
)

anno.orig <- readr::read_csv(in.f) 
anno <- anno.orig[, anno.columns]

# rename the col
colnames(anno) <- c(
  "Uploaded_variation",
  "GT",
  "SYMBOL",
  "Gene",
  "gnomadAF",
  "gnomadAFg",
  "omimSymptomSimScore",
  "hgmdSymptomSimScore",
  "IMPACT",
  "Consequence",
  "hgmdVarFound",
  "clinvarSignDesc",
  "spliceAImax"
)

anno_noGeneID <- anno |> 
  mutate(
    IMPACT.max = case_when(
      IMPACT == "MODIFIER" ~ 1,
      IMPACT == "LOW" ~ 2,
      IMPACT == "MODERATE" ~ 3,
      T ~ 4
    )
  ) |> 
  mutate(
    IMPACT = NULL,
    TierAD = 5 - IMPACT.max,
    TierAR = ifelse(IMPACT.max >= 3, 3, 4),
    TierAR.adj = ifelse(IMPACT.max >= 3, 3, 4),
    No.Var.HM = as.numeric(IMPACT.max >= 3),
    No.Var.H = as.numeric(IMPACT.max == 4),
    No.Var.M = as.numeric(IMPACT.max == 3),
    No.Var.L = as.numeric(IMPACT.max <= 2)
  ) 
  

# For rows with Gene ID ->  Tier analysis
anno <- anno |> filter(Gene != "-") |> distinct() |> 
  mutate(
    IMPACT = ifelse(spliceAImax >= 0.8, "HIGH", IMPACT),
    IMPACT.no = case_when(
      IMPACT == "MODIFIER" ~ 1,
      IMPACT == "LOW" ~ 2,
      IMPACT == "MODERATE" ~ 3, 
      T ~ 4)
    ) 
  


# Implement the rule based Tier classification ####

# take IMPACT HIGH, MODERATE, LOW, or MODIFIER as H,M,L,L
# Tier at Gene level

cat("STEP 1. classify T4: genes with only L variants, nonCodingVar2keep\n")

# STEP 1.1 Keep Genes with IMPACT H or M variants
Gene.all <- unique(anno$Gene)
Gene.IMPACT <- unique(anno[, c("Gene", "IMPACT")])

Gene.IMPACT <- anno |> select(Gene, IMPACT) |> distinct()
Gene.IMPACT.HM <- Gene.IMPACT |> filter(IMPACT %in% c("HIGH", "MODERATE"))

# GeneTier4 will be removed for remaining Tier analysis
GeneTier4 <- Gene.all[!(Gene.all %in% Gene.IMPACT.HM$Gene)]

anno <- anno |> 
  mutate(TierAD = ifelse(Gene %in% GeneTier4, 4, NA)) |> 
  mutate(TierAR = ifelse(Gene %in% GeneTier4, 4, NA)) 


# STEP 1.2 Keep Var with max IMPACT == MODIFIER/LOW and ClinVarClinSig == conflict
# note: many variants with “Conflicting” are non-coding and coding for different transcripts, only keep the var with max impact as MODIFIER

# find the non-Coding Var to keep
nonCodingVar2keep <- c()

# split var to two groups:  ####

# Run Full Tier analysis:
# non - GeneTier4 and nonCodingVar2keep
idx <- is.na(anno$TierAD) 
Var.RunTier <- unique(anno$Uploaded_variation[idx])

# Only Count No.Var.L：
# Others:
# check if there are variants without Tier
ParseVar.NoTier <- F
if (F %in% idx) {
  ParseVar.NoTier <- T
  Var.NoTier <- unique(anno$Uploaded_variation[!idx])
}

anno.f1 <- anno[idx, ]

cat("STEP 2: Get the most severe IMPACT for each variant ####\n")
VEP.f1 <- anno.f1

# it won't be the same order but will work
VEP.f1.maxIMAPCT <- VEP.f1 |> 
  filter(is.na(TierAD)) |> 
  group_by(Uploaded_variation, Gene) |> 
  mutate(
    IMPACT.max = max(IMPACT.no)
  ) |> 
  ungroup()

# find all intronic Var with AF == 0, will adjust the Tier if in trans with a high impact variant. ####

rareVar <- VEP.f1.maxIMAPCT$gnomadAF == "-" &
  VEP.f1.maxIMAPCT$gnomadAFg == "-"
intronicVar <- grepl("intron_variant", VEP.f1.maxIMAPCT$Consequence) &
  VEP.f1.maxIMAPCT$IMPACT.max == 1
rareIntronicVar <- VEP.f1.maxIMAPCT[rareVar &
                                      intronicVar, c("Uploaded_variation", "Gene")]
rareIntronicVar <- unique(rareIntronicVar)

# for each Variant keep the most severe IMPACT only
VEP.f2 <- distinct(
  VEP.f1.maxIMAPCT[, c("Uploaded_variation",
  "SYMBOL",
  "Gene",
  "GT",
  "IMPACT.max",
  "omimSymptomSimScore",
  "hgmdSymptomSimScore")])

cat("STEP 3: assign Tiers\n")

# VEP.f2$GT <- ifelse(VEP.f2$GT=="HET", "0/1","1/1")
# though Tier classification is Gene based, Tier should be assigned to each Variant

VEP.Tier.part1 <- data.frame(matrix(ncol = ncol(VEP.f2) + 6, nrow = 0))
colnames(VEP.Tier.part1) <- c(
  colnames(VEP.f2),
  "TierAR",
  "TierAD",
  "No.Var.HM",
  "No.Var.H",
  "No.Var.M",
  "No.Var.L"
)

VEP.Tier.part1 <- data.frame(matrix(ncol = ncol(VEP.f2) + 6, nrow = 0))
colnames(VEP.Tier.part1) <- c(
  colnames(VEP.f2),
  "TierAR",
  "TierAD",
  "No.Var.HM",
  "No.Var.H",
  "No.Var.M",
  "No.Var.L"
)

cat("STEP 4: Find all genes with No.Var.H == 1, adjust the Tier for rareIntronicVar\n")

VEP.Tier.part1 <- VEP.f2 |> 
  group_by(Gene) |>
  group_modify(~{
    g.df <- .x
    
    # AD regardless of the GT
    g.df <- g.df |> mutate(TierAD = if_else(IMPACT.max < 3, 4, 5 - IMPACT.max))
    
    # AR
    g.df2 <- g.df
    
    # duplicate rows with HOM var
    if ("HOM" %in% g.df$GT) {
      g.df2 <- bind_rows(g.df2, filter(g.df, GT == "1/1"))
    }
    
    IMPACT.max <- g.df2$IMPACT.max
    # Calculate TierAR
    if (sum(IMPACT.max == 4) >= 2) {
      g.df$TierAR <- 5 - g.df$IMPACT.max
      g.df$TierAR[g.df$IMPACT.max == 3] <- 1.5
    } else if (sum(IMPACT.max == 4) == 1 & sum(IMPACT.max == 3) >= 1) {
      g.df$TierAR[g.df$IMPACT.max >= 3] <- 1.5
      g.df$TierAR[g.df$IMPACT.max < 3] <- 5 - g.df$IMPACT.max[g.df$IMPACT.max < 3]
    } else if (sum(IMPACT.max == 4) == 1) {
      g.df$TierAR[g.df$IMPACT.max == 4] <- 3
      g.df$TierAR[g.df$IMPACT.max <= 3] <- 5 - g.df$IMPACT.max[g.df$IMPACT.max <= 3]
    } else if (sum(IMPACT.max == 3) >= 2) {
      g.df$TierAR[g.df$IMPACT.max == 3] <- 2
      g.df$TierAR[g.df$IMPACT.max < 3] <- 5 - g.df$IMPACT.max[g.df$IMPACT.max < 3]
    } else if (sum(IMPACT.max == 3) == 1) {
      g.df$TierAR[g.df$IMPACT.max == 3] <- 3
      g.df$TierAR[g.df$IMPACT.max < 3] <- 5 - g.df$IMPACT.max[g.df$IMPACT.max < 3]
    } else {
      g.df$TierAR <- 5 - g.df$IMPACT.max
    }
    # Calculate the count columns
    g.df$No.Var.HM <- sum(IMPACT.max >= 3)
    g.df$No.Var.H <- sum(IMPACT.max == 4)
    g.df$No.Var.M <- sum(IMPACT.max == 3)
    g.df$No.Var.L <- sum(IMPACT.max == 2)
    
    g.df

  }) |>
  ungroup()

# If TierAR.adj is needed and not calculated earlier, you can add it here
# For example, if it's the same as TierAR:
 VEP.Tier.part1 <- VEP.Tier.part1 |>
  mutate(TierAR.adj = TierAR)

VEP.Tier.part1$TierAR.adj <- VEP.Tier.part1$TierAR


# add the TierAR.adj to the VEP.Tier df

VEP.Tier.part1 <- VEP.Tier.part1[, c(
  "Uploaded_variation",
  "Gene",
  "GT",
  "IMPACT.max",
  "TierAD",
  "TierAR",
  "TierAR.adj",
  "No.Var.HM",
  "No.Var.H",
  "No.Var.M",
  "No.Var.L"
)]

cat("STEP 5: add the No.Var for Var.NoTier\n")
if (ParseVar.NoTier) {
  VEP.Tier.part2 <- anno[anno$Uploaded_variation %in% Var.NoTier, c("Uploaded_variation", "Gene", "GT")]
  VEP.Tier.part2$IMPACT.max <- 1
  VEP.Tier.part2[, c("TierAD", "TierAR", "TierAR.adj")] <- 4
  VEP.Tier.part2[, c("No.Var.HM", "No.Var.H", "No.Var.M", "No.Var.L")] <- NA
  
  VEP.Tier <- rbind(VEP.Tier.part1, VEP.Tier.part2)
} else {
  VEP.Tier <- VEP.Tier.part1
}


VEP.Tier.wGene <- VEP.Tier |>
  group_by(Gene) |>
  group_modify(~{
    g.df <- .x
    
    # Duplicate rows with HOM var
    g.df2 <- if ("HOM" %in% g.df$GT) {
      bind_rows(g.df, filter(g.df, GT == "1/1"))
    } else {
      g.df
    }
    
    IMPACT.max <- g.df2$IMPACT.max
    
    g.df |>
      mutate(
        No.Var.HM = sum(IMPACT.max >= 3),
        No.Var.H = sum(IMPACT.max == 4),
        No.Var.M = sum(IMPACT.max == 3),
        No.Var.L = sum(IMPACT.max == 2)
      )
  }) |>
  ungroup() |>
  select(-GT)  # Remove the GT column

# If you need to ensure the same column order as the original VEP.Tier
VEP.Tier.wGene <- VEP.Tier.wGene |>
  select(names(VEP.Tier)[names(VEP.Tier) != "GT"])

VEP.Tier.wGene <- as_tibble(VEP.Tier.wGene)

cat("STEP 6: add the Tiers from anno_noGeneID\n")
if ("-" %in% anno$Gene) {
  VEP.Tier.final <- rbind(VEP.Tier.wGene, anno_noGeneID)
}
VEP.Tier.final <- VEP.Tier.wGene

colnames(VEP.Tier.final)[colnames(VEP.Tier.final) == "IMPACT.max"] <- "IMPACT.from.Tier"
VEP.Tier.final$TierAR.adj[is.na(VEP.Tier.final$TierAR.adj)] <- VEP.Tier.final$TierAR[is.na(VEP.Tier.final$TierAR.adj)]

# add inheritance information
colnames(genemap2.Inh.F)[1] <- "Gene"
VEP.Tier.wInh <- merge(
  VEP.Tier.final,
  genemap2.Inh.F,
  by.x = "Gene",
  by.y = "Gene",
  all.x = T
)
VEP.Tier.wInh$dominant[is.na(VEP.Tier.wInh$dominant)] <- 0
VEP.Tier.wInh$recessive[is.na(VEP.Tier.wInh$recessive)] <- 0
VEP.Tier.wInh$AD.matched <- ifelse(VEP.Tier.wInh$TierAD <= 2 &
                                     (VEP.Tier.wInh$dominant == 1), 1, 0)
VEP.Tier.wInh$AR.matched <- ifelse(VEP.Tier.wInh$TierAR <= 2 &
                                     (VEP.Tier.wInh$recessive == 1), 1, 0)

write.table(
  VEP.Tier.wInh,
  out.f,
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

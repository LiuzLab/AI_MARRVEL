# Run Tier

# Updates history ####
# Update on 2021.Sep13
#' Updates:
#' 1. corrected chromosome number 23->X, 24-> Y
#' 2. add interaction term of phenotype and Tier. 
#' 3. for non-coding var, if ClinVarClinSig == conflict ->  push the Tier higher (if IMPACT == MODIFIER/LOW and ClinVarClinSig == conflict, Tier = Tier + 1)
#' 4. Keep all variant that in-trans (sudo) with a strong Var in the same gene, adjust the Tier to 3.5?
# Updates on 2022.Apr5
#' 1. remove the pre-processing of OMIM inheritance to VarTier_preprocessOMIM.R
#' 2. Add the number of variants (H/M/L impact) per gene: "No.Var.HM","No.Var.H","No.Var.M","No.Var.L"
#' 3. Add the option if IncludeClinVarHGMD 
#' 
#' Updates on 2022.Sep29
#' ClinVar clinvarSignDesc is in dif format after optimizing the pipeline, update the code to be compatible for that 
#' gnomadAF, gnomadAFg is "-" not "NA", update the code for that
#' add IncludeClinVarHGMD for line 101


#' change genemap2.Inh into 2 features, AD, AR, save to tsv file
#' keep GeneID == "-" assinge Tier value
#' input  output /out/explain/ "ID" + "_Tier"
#' 

# parameter settings ####

# set adj.k value to adjust Tier
adj.k <- 1.05

# set if diseases DB will be included
IncludeClinVarHGMD <- FALSE

# load database: OMIM inheritance ####
# genemap2.Inh.F <- readRDS(file.path("/run/genemap2.Inh.F.rds"))
#genemap2.Inh.F <- read.delim2(file.path("/run/data_dependencies/var_tier/hg19/genemap2.Inh.F.txt"))
# genemap2.Inh.F <- read.delim2("/houston_10t/dongxue/MARRVEL_AI/resource/OMIM/out/genemap2.Inh.F.txt")
args = commandArgs(trailingOnly=TRUE)
refg <- args[2]
# load database: OMIM inheritance ####
# genemap2.Inh.F <- readRDS(file.path("/run/genemap2.Inh.F.rds"))
genemap2.Inh.F <- read.delim2(file.path(paste0('/run/data_dependencies/var_tier/',refg,'/genemap2.Inh.F.txt')))

# set the input output file / directory ####



# in.f.path <- "/out/input/"
in.f.path <- "/out/rami-test/"
out.f.path <- "/out/tier-test-false/"

# When test the code without docker
#in.f.path <- out.f.path <- "/out/"

# in.f <- paste0(args[1],"_Tier")
# in.f <- paste0(args[1],"_scores.csv")
in.f <- paste0(args[1],"_scores.csv")
out.f <- paste0(args[1],"_Tier.v2.tsv")

# load libraries ####
library("data.table")
library("dplyr")



# load the input ####
anno.columns <- c('varId_dash','zyg','geneSymbol','geneEnsId','gnomadAF','gnomadAFg','omimSymptomSimScore','hgmdSymptomSimScore',
                  'IMPACT','Consequence','hgmdVarFound','clinvarSignDesc','spliceAImax')

anno.orig <- read.table(file.path(in.f.path,in.f), stringsAsFactors = F, sep = ",", header = T)
# To test:
# score.path <- "/houston_20t/chaozhong/MARRVEL_AI_2022/explain/906246_raw.csv"
# score.path <- "/houston_20t/chaozhong/MARRVEL_AI_2022/opt_09302022/output/rami-test/906246_scores.csv"
# score.path <- "/houston_30t/chaozhong/NDV_model/fake_VCF_test/out/rami-test/NDVG_causal_scores.csv"
# anno.orig <- read.csv(score.path)
anno <- anno.orig[,anno.columns]

# rename the col
colnames(anno) <- c("Uploaded_variation","GT","SYMBOL","Gene","gnomadAF","gnomadAFg",
                      "omimSymptomSimScore", "hgmdSymptomSimScore","IMPACT",
                      "Consequence","hgmdVarFound","clinvarSignDesc","spliceAImax")

if ("-" %in% anno$Gene) {
  # for rows without gene ID -> Assign Tier now:
  anno_noGeneID <- anno[anno$Gene == "-",c("Uploaded_variation","Gene","IMPACT")]
  # convert IMPACT to numeric value
  anno_noGeneID$IMPACT.max <- ifelse(anno_noGeneID$IMPACT == "MODIFIER", 1, ifelse(anno_noGeneID$IMPACT == "LOW", 2, ifelse(anno_noGeneID$IMPACT == "MODERATE",3,4)) )
  anno_noGeneID$IMPACT <- NULL
  # "TierAD","TierAR","TierAR.adj"
  anno_noGeneID$TierAD <- 5 - anno_noGeneID$IMPACT.max
  anno_noGeneID[,c("TierAR","TierAR.adj")] <- ifelse(anno_noGeneID$IMPACT.max >= 3, 3,4)
  anno_noGeneID$No.Var.HM <- as.numeric(anno_noGeneID$IMPACT.max >= 3)
  anno_noGeneID$No.Var.H <- as.numeric(anno_noGeneID$IMPACT.max == 4)
  anno_noGeneID$No.Var.M <- as.numeric(anno_noGeneID$IMPACT.max == 3)
  anno_noGeneID$No.Var.L <- as.numeric(anno_noGeneID$IMPACT.max <=2)
}



# For rows with Gene ID ->  Tier analysis
anno <- anno[anno$Gene != "-",]
anno <- unique(anno)
  

# Tier ####


# STEP 0. Adjust variant IMPACT for ClinVar, HGMD, SpliceAImax >= 0.8  ####

if (IncludeClinVarHGMD){
  anno$IMPACT <- ifelse(grepl("Pathogenic|Likely_pathogenic",anno$clinvarSignDesc) | anno$hgmdVarFound == 1 , "HIGH",anno$IMPACT)
}

# label Splice AI > 0.8 as HIGH
anno$IMPACT <- ifelse(anno$spliceAImax >=0.8, "HIGH",anno$IMPACT)
# convert IMPACT to numeric value
anno$IMPACT.no <- ifelse(anno$IMPACT == "MODIFIER", 1, ifelse(anno$IMPACT == "LOW", 2, ifelse(anno$IMPACT == "MODERATE",3,4)) )
  

# Implement the rule based Tier classification ####
  
# take IMPACT HIGH, MODERATE, LOW, or MODIFIER as H,M,L,L
# Tier at Gene level

# STEP 1. classify T4: genes with only L variants, nonCodingVar2keep ####

# STEP 1.1 Keep Genes with IMPACT H or M variants
Gene.all <- unique(anno$Gene)
Gene.IMPACT <- unique(anno[,c("Gene","IMPACT")])
Gene.IMPACT.HM <- Gene.IMPACT[Gene.IMPACT$IMPACT %in% c("HIGH","MODERATE"),]
  
# GeneTier4 will be removed for remaining Tier analysis
GeneTier4 <- Gene.all[!(Gene.all %in% Gene.IMPACT.HM$Gene)]
anno$TierAD[anno$Gene %in% GeneTier4] = anno$TierAR[anno$Gene %in% GeneTier4] <- 4

# STEP 1.2 Keep Var with max IMPACT == MODIFIER/LOW and ClinVarClinSig == conflict 
# note: many variants with “Conflicting” are non-coding and coding for different transcripts, only keep the var with max impact as MODIFIER

# find the non-Coding Var to keep
nonCodingVar2keep <- c()

if(IncludeClinVarHGMD){
  nonCodingModifier.idx <- which(grepl("Conflicting_interpretations_of_pathogenicity", anno$clinvarSignDesc) & (anno$IMPACT %in% c("MODIFIER")))
  nonCodingModifier <- anno$Uploaded_variation[nonCodingModifier.idx]
  nonCodingModifier.df <- anno[anno$Uploaded_variation %in% nonCodingModifier,]
  
  for(var in unique(nonCodingModifier.df$Uploaded_variation)){
    var.df <- nonCodingModifier.df[nonCodingModifier.df$Uploaded_variation == var,]
    var.impact <- unique(var.df$IMPACT)
    if(length(var.impact) == 1 ){
      if(var.impact == "MODIFIER"){
        nonCodingVar2keep <- c(nonCodingVar2keep,var)
      }
    }
  }
}


# split var to two groups:  ####

# Run Full Tier analysis:
# non - GeneTier4 and nonCodingVar2keep
idx <- is.na(anno$TierAD) | (anno$Uploaded_variation %in% nonCodingVar2keep)
Var.RunTier <- unique(anno$Uploaded_variation[idx])

# Only Count No.Var.L：
# Others:
# check if there are variants without Tier
ParseVar.NoTier <- F
if(F %in% idx){
  ParseVar.NoTier <- T
  Var.NoTier <- unique(anno$Uploaded_variation[!idx])
}

anno.f1 <- anno[idx,]

# set nonCodingVar2keep IMPACT as 3
anno.f1$IMPACT.no[anno.f1$Uploaded_variation %in% nonCodingVar2keep] <- 3


# STEP 2: Get the most severe IMPACT for each variant ####

VEP.f1 <- anno.f1

# for Gene-Var pair with multiple IMPACT, take the MAX
VEP.f1.maxIMAPCT <- data.frame(matrix(ncol = ncol(VEP.f1) + 1, nrow = 0))
colnames(VEP.f1.maxIMAPCT) <- c(colnames(VEP.f1),"IMPACT.max")
for (var in Var.RunTier){
  var.df <- VEP.f1[VEP.f1$Uploaded_variation == var,]
  # check no. of gene for this variant
  if(length(unique(var.df$Gene)) == 1){
    var.df$IMPACT.max <- max(var.df$IMPACT.no)
    VEP.f1.maxIMAPCT <- rbind(VEP.f1.maxIMAPCT,var.df)
  } else {
    for (gene in unique(var.df$Gene)){
      var.g.df <- var.df[var.df$Gene == gene,]
      var.g.df$IMPACT.max <- max(var.g.df$IMPACT.no)
      VEP.f1.maxIMAPCT <- rbind(VEP.f1.maxIMAPCT,var.g.df)
    }
  }
}
  
# find all intronic Var with AF == 0, will adjust the Tier if in trans with a high impact variant. ####

rareVar <- VEP.f1.maxIMAPCT$gnomadAF == "-" & VEP.f1.maxIMAPCT$gnomadAFg == "-"
intronicVar <- grepl("intron_variant",VEP.f1.maxIMAPCT$Consequence)  & VEP.f1.maxIMAPCT$IMPACT.max == 1
rareIntronicVar <- VEP.f1.maxIMAPCT[rareVar & intronicVar, c("Uploaded_variation","Gene")]
rareIntronicVar <- unique(rareIntronicVar)

# for each Variant keep the most severe IMPACT only
VEP.f2 <- unique(VEP.f1.maxIMAPCT[,c("Uploaded_variation", "SYMBOL","Gene","GT","IMPACT.max","omimSymptomSimScore","hgmdSymptomSimScore")])

# STEP 3: assign Tiers ####
  
# VEP.f2$GT <- ifelse(VEP.f2$GT=="HET", "0/1","1/1")

# though Tier classification is Gene based, Tier should be assigned to each Variant 
  
VEP.Tier.part1 <- data.frame(matrix(ncol = ncol(VEP.f2) + 6, nrow = 0))
colnames(VEP.Tier.part1) <- c(colnames(VEP.f2),"TierAR","TierAD","No.Var.HM","No.Var.H","No.Var.M","No.Var.L")

for (g in unique(VEP.f2$Gene)){
  
  g.df <- VEP.f2[VEP.f2$Gene == g,]
  
  # X linked, seudo autosomal region -> treat as autosomal
  # if (sum(grepl("X",g.df$Uploaded_variation)) >0) {
  #   g.df$TierAD = g.df$TierAR = ifelse(g.df$IMPACT.max <3, 4, 5-g.df$IMPACT.max)
  # }
  
  # AD regardless of the GT
  g.df$TierAD <- ifelse(g.df$IMPACT.max <3, 4, 5-g.df$IMPACT.max)
  
  # AR
  g.df2 <- g.df
  
  # duplicate rows with HOM var
  if ("HOM" %in% g.df$GT){
    g.df2 <- rbind(g.df, g.df[g.df$GT == "1/1",])
  }
  IMPACT.max <- g.df2$IMPACT.max
  if(sum(IMPACT.max==4) >=2){
    # >=2 strong alleles
    g.df$TierAR <- 5-g.df$IMPACT.max
    # update 202221005: H+H: M-> change Tie from 2 to 1.5
    g.df$TierAR[g.df$IMPACT.max ==3] <- 1.5
  } else if (sum(IMPACT.max==4) ==1 & sum(IMPACT.max==3) >=1) {
    # H + M
    g.df$TierAR[g.df$IMPACT.max >=3] <- 1.5 # (need more thought, 2 or 1.5 ???)
    g.df$TierAR[g.df$IMPACT.max <3] <- 5-g.df$IMPACT.max[g.df$IMPACT.max <3]
  } else if (sum(IMPACT.max==4) ==1 ) {
    # H only
    g.df$TierAR[g.df$IMPACT.max ==4] <- 3 # (need more thought, 2 or 3 ???)
    g.df$TierAR[g.df$IMPACT.max <=3] <- 5-g.df$IMPACT.max[g.df$IMPACT.max <=3]
  } else if (sum(IMPACT.max==3) >=2 ) {
    # M + M
    g.df$TierAR[g.df$IMPACT.max ==3] <- 2 # (need more thought, 2 or 3 ???)
    g.df$TierAR[g.df$IMPACT.max < 3] <- 5-g.df$IMPACT.max[g.df$IMPACT.max < 3]
  }  else if (sum(IMPACT.max==3) == 1 ) {
    # M only
    g.df$TierAR[g.df$IMPACT.max ==3] <- 3 # (need more thought, 2 or 3 ???)
    g.df$TierAR[g.df$IMPACT.max < 3] <- 5-g.df$IMPACT.max[g.df$IMPACT.max < 3]
  }
  
  # count the number of high IMPACT var per gene. 
  No.Var.HM <- sum(IMPACT.max >=3)
  No.Var.H <- sum(IMPACT.max == 4)
  No.Var.M <- sum(IMPACT.max ==3)
  No.Var.L <- sum(IMPACT.max == 2)
  g.df$No.Var.HM <- No.Var.HM
  g.df$No.Var.H <- No.Var.H
  g.df$No.Var.M <- No.Var.M
  g.df$No.Var.L <- No.Var.L
  
  VEP.Tier.part1 <- rbind(VEP.Tier.part1,g.df)
}
  
# STEP 4: Find all genes with No.Var.H == 1, adjust the Tier for rareIntronicVar ####

if (IncludeClinVarHGMD){
  gene.w1HighVar <- VEP.Tier.part1$Gene[VEP.Tier.part1$No.Var.H == 1]
  rareIntronicVar.2adj <- rareIntronicVar[rareIntronicVar$Gene %in% gene.w1HighVar,]
  
  # adjust the Tier for rareIntronicVar based on similarity
  phenoSimi <- unique(VEP.f1.maxIMAPCT[,c("Gene","omimSymptomSimScore","hgmdSymptomSimScore")])
  
  # check for non numeric values and convert to 0
  omim_NA_idx <- which(is.na(as.numeric(phenoSimi$omimSymptomSimScore)))
  hgmd_NA_idx <- which(is.na(as.numeric(phenoSimi$hgmdSymptomSimScore)))
  # print(hgmd_NA_idx)
  phenoSimi$omimSymptomSimScore[omim_NA_idx] <- 0
  phenoSimi$hgmdSymptomSimScore[hgmd_NA_idx] <- 0
  
  phenoSimi$omimSymptomSimScore <- as.numeric(phenoSimi$omimSymptomSimScore)
  phenoSimi$hgmdSymptomSimScore <- as.numeric(phenoSimi$hgmdSymptomSimScore)
  
  phenoSimi <- phenoSimi[phenoSimi$omimSymptomSimScore != 0 | phenoSimi$hgmdSymptomSimScore != 0,]
  
  rareIntronicVar.2adj.wPhenoSimi <- merge(rareIntronicVar.2adj, phenoSimi, by.x = "Gene", by.y = "Gene")
  rareIntronicVar.2adj.wPhenoSimi$simi <- apply(rareIntronicVar.2adj.wPhenoSimi[, 3:4], 1, max)
  rareIntronicVar.2adj.wPhenoSimi$TierAR.adj <- 5-(5-4)*exp(adj.k*rareIntronicVar.2adj.wPhenoSimi$simi)
  
  rareIntronicVar.df <- rareIntronicVar.2adj.wPhenoSimi[,c("Uploaded_variation","Gene","TierAR.adj")]
  VEP.Tier.part1 <- merge(VEP.Tier.part1,rareIntronicVar.df,
                          by.x <- c("Uploaded_variation","Gene"),
                          by.y = c("Uploaded_variation","Gene"),
                          all.x = T)
} else {
  VEP.Tier.part1$TierAR.adj <- VEP.Tier.part1$TierAR
}

# add the TierAR.adj to the VEP.Tier df

VEP.Tier.part1 <- VEP.Tier.part1[,c("Uploaded_variation","Gene","GT","IMPACT.max","TierAD","TierAR","TierAR.adj","No.Var.HM","No.Var.H","No.Var.M","No.Var.L")]

# STEP 5: add the No.Var for Var.NoTier  ####
if (ParseVar.NoTier){
  VEP.Tier.part2 <- anno[anno$Uploaded_variation %in% Var.NoTier,c("Uploaded_variation","Gene","GT")]
  VEP.Tier.part2$IMPACT.max <- 1
  VEP.Tier.part2[,c("TierAD","TierAR","TierAR.adj")] <- 4
  VEP.Tier.part2[,c("No.Var.HM","No.Var.H","No.Var.M","No.Var.L")] <- NA
  
  VEP.Tier <- rbind(VEP.Tier.part1, VEP.Tier.part2)
} else {
  VEP.Tier <- VEP.Tier.part1
}

VEP.Tier.wGene <- data.frame(matrix(ncol = ncol(VEP.Tier), nrow = 0))
for (g in unique(VEP.Tier$Gene)){

  g.df <- VEP.Tier[VEP.Tier$Gene == g,]
  # AR
  g.df2 <- g.df
  
  # duplicate rows with HOM var
  if ("HOM" %in% g.df$GT){
    g.df2 <- rbind(g.df, g.df[g.df$GT == "1/1",])
  }
  IMPACT.max <- g.df2$IMPACT.max
  
  No.Var.HM <- sum(IMPACT.max >=3)
  No.Var.H <- sum(IMPACT.max == 4)
  No.Var.M <- sum(IMPACT.max ==3)
  No.Var.L <- sum(IMPACT.max == 2)
  g.df$No.Var.HM <- No.Var.HM
  g.df$No.Var.H <- No.Var.H
  g.df$No.Var.M <- No.Var.M
  g.df$No.Var.L <- No.Var.L
  
  VEP.Tier.wGene <- rbind(VEP.Tier.wGene,g.df)
}
# remove the GT column
VEP.Tier.wGene$GT <- NULL


# STEP 6: add the Tiers from anno_noGeneID
if ("-" %in% anno$Gene) {
  VEP.Tier.final <- rbind(VEP.Tier.wGene, anno_noGeneID)
}
VEP.Tier.final <- VEP.Tier.wGene

colnames(VEP.Tier.final)[colnames(VEP.Tier.final) == "IMPACT.max"] <- "IMPACT.from.Tier"
VEP.Tier.final$TierAR.adj[is.na(VEP.Tier.final$TierAR.adj)] <- VEP.Tier.final$TierAR[is.na(VEP.Tier.final$TierAR.adj)]
  

# add inheritance information
colnames(genemap2.Inh.F)[1] <- "Gene"
VEP.Tier.wInh <- merge(VEP.Tier.final, genemap2.Inh.F, by.x = "Gene", by.y = "Gene", all.x = T)
VEP.Tier.wInh$dominant[is.na(VEP.Tier.wInh$dominant)] <- 0
VEP.Tier.wInh$recessive[is.na(VEP.Tier.wInh$recessive)] <- 0
VEP.Tier.wInh$AD.matched <- ifelse(VEP.Tier.wInh$TierAD <=2 & (VEP.Tier.wInh$dominant == 1),1,0)
VEP.Tier.wInh$AR.matched <- ifelse(VEP.Tier.wInh$TierAR <=2 & (VEP.Tier.wInh$recessive == 1),1,0)

write.table(VEP.Tier.wInh,file.path(out.f.path,out.f), sep = "\t", quote = F, col.names = T, row.names = F)

# write.table(VEP.Tier.wInh,"/houston_20t/dongxue/MARRVEL_AI/V2/FeatureEng/TestOut/906246_Tier.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
# write.table(VEP.Tier.wInh,"/houston_30t/dongxue/AI/UDN_Sample_2023/ModifiedAIMCode/NDVG_causal_Tier.v2.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
### Plot for Tier adjustment ####




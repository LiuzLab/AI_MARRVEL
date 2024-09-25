#NOTE: Replace the download links with newest ones and then parse the files

library(tidyverse)
library(httr)

#download new hpo-related data files
hp_obo_file <- "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/hp.obo"
hpo_omim_file <- "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/phenotype.hpoa"
genemap2_file <- "https://data.omim.org/downloads/E8eFWaP3SOu67gXVTVKiGA/genemap2.txt"

GET(hp_obo_file, write_disk("hp.obo", overwrite = TRUE))
GET(hpo_omim_file, write_disk("phenotype.hpoa", overwrite = TRUE))
GET(genemap2_file, write_disk("genemap2.txt", overwrite = TRUE))


#read the new data file
df_new <- read_delim("phenotype.hpoa", delim = '\t', skip = 4)

old_colnames <- c("OMIM_ID",
                  "HPO_ID",
                  "DiseaseName",
                  "Onset",
                  "Frequency",
                  "Sex",
                  "Aspect")

#make the format of new data file same as the old one
df_new_omim <- df_new %>%
  separate(
    database_id,
    into = c("database", "omim"),
    sep = ":",
    remove = FALSE
  ) %>%
  filter(database == "OMIM")

df_new_omim_clean <- df_new_omim %>%
  select(omim, hpo_id, disease_name, onset, frequency, sex, aspect)

colnames(df_new_omim_clean) <- old_colnames
df_new_omim_clean$OMIM_ID <- as.numeric(df_new_omim_clean$OMIM_ID)

#save the new HPO_OMIM.tsv
write_tsv(df_new_omim_clean, 'HPO_OMIM.tsv')

#parse the genemap2 table
system("cat genemap2.txt | ./parseGeneMap2_output.py > genemap2_parsed.txt",
       intern = FALSE)

#organize the table
old_colnames <- c(
  "MIM_Number",
  "Pheno_ID",
  "Phenotypes",
  "Approved_Gene_Symbol",
  "Gene_Symbols",
  "Ensembl_Gene_ID",
  "Entrez_Gene_ID",
  "inheritance",
  "Pheno_mapKey"
)

gm2_new_parsed <- read_delim('genemap2_parsed.txt')
gm2_new_parsed_selected <- gm2_new_parsed |>
  select(
    MIM_Number,
    Phenotype_MIM_Number,
    Phenotype_Name,
    Approved_Gene_Symbol,
    Gene_Symbols,
    Ensembl_Gene_ID,
    Entrez_Gene_ID,
    Inheritance,
    Mapping_Key
  )
colnames(gm2_new_parsed_selected) <- old_colnames

#save genemap2 data as rds
#the file name here is not changed, even though the newest data is from 2024
saveRDS(gm2_new_parsed_selected, file = "genemap2_v2022.rds")

#remove unnecessary file
file.remove("phenotype.hpoa", "genemap2.txt", "genemap2_parsed.txt")

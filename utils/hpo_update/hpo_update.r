library(tidyverse)
library(httr)

#download new hpo-related data files
hp_obo_file <- "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/hp.obo"
hpo_omim_file <- "https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2024-08-13/phenotype.hpoa"
genemap2_file <- "https://data.omim.org/downloads/E8eFWaP3SOu67gXVTVKiGA/genemap2.txt"

GET(hp_obo_file, write_disk("hp.obo", overwrite = TRUE))
GET(hpo_omim_file, write_disk("phenotype.hpoa", overwrite = TRUE))
GET(genemap2_file, write_disk("genemap2.txt", overwrite = TRUE))


#use old dataframe for reference; read the new data file
old_colnames <- c("OMIM_ID",
                  "HPO_ID",
                  "DiseaseName",
                  "Onset",
                  "Frequency",
                  "Sex",
                  "Aspect")
df_new <- read_delim("phenotype.hpoa", delim = '\t', skip = 4)

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

#remove unnecessary file
file.remove("phenotype.hpoa")

#save the new HPO_OMIM.tsv
write_tsv(df_new_omim_clean, 'HPO_OMIM.tsv')

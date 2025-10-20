setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# Load packages
library(duckplyr)
library(tidyverse)

# import metadata
metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/UC_MP_Emperor_Map.csv")
metadata_summary_2 <- metadata %>% 
  ungroup() %>% 
  group_by(Disease, `current_5ASA`, ASA_exposure) %>% 
  summarize(count = n(), .groups = "drop")

# import feature table
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/MZmine_outout_2_join_aligner_02min/fbmn_02min_join_aligner_iimn_fbmn_quant.csv")
feature_table <- feature_table %>%
  rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)

# import annotations
annotations <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/MSV000082094_2/MSV000082094_IBD_multiplex_libraries.tsv")
annotations$`#Scan#` <- as.character(annotations$`#Scan#`)
annotations <- annotations |> 
  dplyr::rename(Feature = `#Scan#`)

# transpose data
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzXML")) |> 
  t() |>  
  as.data.frame() |>  
  rownames_to_column("filename") |> 
  dplyr::mutate(id = sub("_.*", "", filename)) |> 
  dplyr::select(filename, id, everything())

data_5ASA <- data_transpose %>%
  dplyr::select(id, c(`108`, `964`, `2381`, `3589`,`3614`))

keep_ids <- c("108","964","2381","3589","3614")

# OLD -> NEW names
name_map <- c(
  "108"  = "5-aminosalicylic_acid",
  "964"  = "5-Aminosalicylic acid_propanoic-acid",
  "2381" = "5-aminosalicylic acid_phenylpropionic acid",
  "3589" = "5-Aminosalicylic acid_palmitic-acid",
  "3614" = "5-Aminosalicylic acid_oleic-acid")

# Rename only present columns (dplyr::rename uses new = old)
old_ids_present <- intersect(keep_ids, names(data_transpose))

data_5ASA <- data_transpose %>%
  dplyr::select(id, tidyselect::any_of(old_ids_present)) %>%
  dplyr::rename(!!! setNames(old_ids_present, unname(name_map[old_ids_present])))

merged_df_UC <- data_5ASA %>%
  dplyr::left_join(metadata, by = c("id" = "vial_name")) %>%
  dplyr::select(id, tidyselect::any_of(unname(name_map)), Disease, current_5ASA, ASA_exposure, everything())

#write_csv(merged_df_UC, "MSV000082094_40_patients_merged_df_UC.csv")


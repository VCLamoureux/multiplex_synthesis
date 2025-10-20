setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# Load packages
library(tidyverse)

# import table
metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/metadata (1).csv")

# import feature table
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/MZmine_output/fbmn_pediatric_iimn_fbmn_quant.csv")
feature_table <- feature_table %>%
  rename_all(~gsub(" Peak area", "", .))
colnames(feature_table)

# import annotation table
annotations <- fread("OneDrive-UniversityofCalifornia,SanDiegoHealth/Download_public_data/Pediatric_cohort_Shirley/99cebd129c6748e3b67a8834589b9faa-merged_results_with_gnps.tsv")
annotations$`#Scan#` <- as.character(annotations$`#Scan#`)
annotations <- annotations |> 
  dplyr::rename(Feature = `#Scan#`)

# get info feature
info_feature <- feature_table[,1:3]
colnames(info_feature) <- c("Feature", "mz", "RT")
info_feature$Feature <- as.character(info_feature$Feature)
info_feature_name <- left_join(info_feature, annotations, by = "Feature")

# transpose data
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzML")) |> 
  t() |>  
  as.data.frame() |>  
  rownames_to_column("filename")

# Merge with metadata using the matching file name columns
merged_df <- data_transpose |> 
  dplyr::mutate(Status = case_when(
    str_detect(filename, "HC") ~ "healthy",
    str_detect(filename, "PT") ~ "IBD", TRUE ~ NA_character_)) |> 
  dplyr::select(Status, filename, everything())
merged_df$filename <- gsub(".mzML", "", merged_df$filename)

merged_df_metadata <- merged_df |> 
  left_join(metadata, by = c("filename" = "sample_ID")) |> 
  dplyr::filter(!is.na(Status))

merged_df_metadata_synthesis <- merged_df_metadata |> 
  dplyr::select(`5239`,`35389`,`31962`,`64841`,`127282`, `123913`, `125288`, `127323`, filename, Status, Medications, ALT, AST, Alkaline_Phosphatase, Cohort)

keep_ids <- c("5239","35389","31962","64841","127282", "123913", "125288", "127323")
name_map <- c(
  "5239" = "5-aminosalicylic_acid",
  "35389" = "5-Aminosalicylic acid_propanoic-acid",
  "31962" = "5-Aminosalicylic acid_Succinic acid",
  "64841" = "5-aminosalicylic acid_phenylpropionic acid",
  "127282" = "5-Aminosalicylic acid_oleic-acid", 
  "123913" = "5-Aminosalicylic acid_ricinoleic-acid_1", 
  "125288" = "5-Aminosalicylic acid_ricinoleic-acid_2", 
  "127323" = "5-Aminosalicylic acid_12-hydroxy-stearic-acid")

# Select the columns, then rename using the map (only those present are changed)
merged_df_metadata_synthesis_renamed <- merged_df_metadata_synthesis %>%
  dplyr::select(filename, any_of(keep_ids), Medications, Cohort) %>%
  dplyr::rename_with(~ unname(ifelse(.x %in% names(name_map), name_map[.x], .x)), .cols = -filename)

# export table
#write_csv(merged_df_metadata_synthesis_renamed, "Pediatrics_5ASA_derivatives.csv")


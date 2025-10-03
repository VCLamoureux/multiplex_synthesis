# Set your working directory
#setwd("")
setwd("/Users/vincentlamoureux/")

# These packages are part of CRAN repository
#install.packages("data.table", dependencies = TRUE)
#install.packages("tidyverse", dependencies = TRUE)
#install.packages("pheatmap", dependencies = TRUE)

## Load the packages required for the analysis
library(data.table)
library(tidyverse)
library(pheatmap)
library(circlize)

# Specify the folder path - it should be the folder inside the working directory 
folder_path <- "Desktop/output_fasst_5ASA_distinct_V3/"

## Download/Import the ReDU metadata file - it should be in the working directory folder and NOT be in the sub-folder with the csv files from the Fast Search

# Define the filename for the ReDU metadata 
processed_redu_metadata <- "OneDrive - University of California, San Diego Health/Abubaker/all_sampleinformation_september2025.tsv"

# Check if the pre-processed metadata file exists in the working directory
if (!file.exists(file.path(getwd(), processed_redu_metadata))) {
  redu_url <- "https://redu.gnps2.org/dump"
  options(timeout = 600) # If 10 min is not enough, add more time 
  download.file(redu_url, file.path(getwd(), processed_redu_metadata), mode = "wb")
  redu_metadata <- data.table::fread(processed_redu_metadata)
} else {
  redu_metadata <- data.table::fread(processed_redu_metadata)
}

keep_cols <- c("compound_name", "lib_usi", "USI", "Charge", "Cosine", "Matching Peaks", "Dataset")

read_one <- function(path) {
  df <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()))
  df$Compound <- tools::file_path_sans_ext(basename(path))
  missing <- setdiff(keep_cols, names(df))
  if (length(missing)) df[missing] <- NA_character_
  dplyr::select(df, dplyr::any_of(c("Compound", keep_cols)))
}

file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

molecules_interest <- purrr::map_dfr(file_list, read_one)
molecules_interest <- molecules_interest %>%
  dplyr::mutate(
    USI = str_squish(as.character(USI)),
    USI = na_if(USI, ""),
    USI = na_if(USI, "NA"),
    USI = na_if(USI, "na")) %>%
  dplyr::filter(!is.na(USI))


molecules_interest_filtered <- molecules_interest |> 
  dplyr::filter(Cosine >= 0.7 & `Matching Peaks` >= 4 )

MassiveID_filename <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the USI column in the molecules_interest
molecules_interest_filtered$USI <- vapply(molecules_interest_filtered$USI, MassiveID_filename, FUN.VALUE = character(1))

# Prepare the ReDU metadata USI column for merging with FASST output table
## Create a function to extract the datasetID and the last segment (filename)

ReDU_USI <- function(USI) {
  USI <- gsub("/", ":", USI)
  USI <- sub("\\.[^\\.]*$", "", USI)
  parts <- unlist(strsplit(USI, ":"))
  combined <- paste(parts[2], parts[length(parts)], sep = ":")
  return(combined)
}

# Apply the function to each row of the fxlename column in the ReDU output table
redu_metadata$USI <- vapply(redu_metadata$USI, ReDU_USI, FUN.VALUE = character(1))

# Merge the ReDU metadata table and the FASST MASST output table
ReDU_MASST <- left_join(molecules_interest_filtered, redu_metadata, by = "USI", relationship = "many-to-many")

# Once both data tables are merged, ones can filter the table which based on the research question
## To note: not all publicly available files have associated metadata and we strongly encourage scientists to make 
### their data available with a very detailed metadata (sample information)
#### As more data are being deposited in repositories more matches will be uncovered and more results will be embedded in heatmaps 

# Standardize the body parts and Health Status
ReDU_MASST_standardize <- ReDU_MASST |> 
  dplyr::mutate(
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'skin of trunk|skin of pes|head or neck skin|axilla skin|skin of manus|arm skin|skin of leg', 'skin'),
    UBERONBodyPartName = str_replace_all(UBERONBodyPartName, 'blood plasma|blood serum', 'blood'),
    HealthStatus = str_replace(HealthStatus, 'Chronic Illness', 'chronic illness'),
    HealthStatus = str_replace(HealthStatus, 'Healthy', 'healthy'))

# Separate humans and rodents from the merged data table
df_humans <- ReDU_MASST_standardize |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

analyze_counts <- function(df, column_interest) {
  
# Create a list of all unique entries in the column of interest and create a df
  df_body_parts <- df |>  distinct(across(all_of(column_interest)))
  
  # Count occurrences of each entry in the column of interest
  df_BodyPartName_counts <- df |> 
    count(across(all_of(column_interest)), name = "Counts_fastMASST")
  
  # Aggregate the number and list of unique Compounds for each entry
  compounds <- df |> 
    group_by(across(all_of(column_interest))) |> 
    summarise(Compounds = n_distinct(Compound),
              CompoundsList = toString(unique(Compound))) |> 
    ungroup()
  
  # Merge all the data into a single data frame
  combined <- df_body_parts |> 
    left_join(df_BodyPartName_counts, by = column_interest) |> 
    left_join(compounds, by = column_interest)
  
  return(combined)
}

# Get a glimpse of the number of counts per organ
body_counts_humans <- analyze_counts(df_humans, "UBERONBodyPartName")
head(body_counts_humans)

# Create a function to pivot the table for data visualization
prepare_pivot_table <- function(df, column_interest, compound) {
  
  grouped_df <- df |> 
    group_by(across(all_of(c(compound, column_interest)))) |> 
    summarise(Count = n(), .groups = 'drop')
  
  pivot_table <- grouped_df |> 
    pivot_wider(names_from = all_of(compound), values_from = Count, values_fill = list(Count = 0))
  
  return(pivot_table)
}

# Define the variables based on your research question 
## Here we are interesting in organ distribution in humans and rodents of the molecule of interest
variable <- 'UBERONBodyPartName'
pivot_table_humans <- prepare_pivot_table(df_humans, variable, 'Compound')
pivot_table_rodents <- prepare_pivot_table(df_rodents, variable, 'Compound')

# Prepare the table to be compatible with pheatmap package
humans_molecules_counts_by_bodypart <- pivot_table_humans |> 
  dplyr::arrange(UBERONBodyPartName) |> 
  tibble::column_to_rownames("UBERONBodyPartName")

clean_molecule_names <- function(nm) {
  nm <- sub("_known.*$", "", nm, ignore.case = TRUE)
  nm <- gsub("(?i)^[0-9]+_(?=5\\s*[-_]?ASA)", "", nm, perl = TRUE)
  nm <- gsub("(?i)^[0-9]+_(?=5\\s*[-_]?aminosalicylic)", "", nm, perl = TRUE)
  nm <- gsub("(?i)[_-]acid\\b", " acid", nm, perl = TRUE)
  nm
}


# Convert all columns to numeric for the humans df
humans_molecules_counts_by_bodypart_mat <- humans_molecules_counts_by_bodypart %>%
  rename_with(clean_molecule_names) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

# Define your chosen colors
colors_version <- c("#FFFFFF", "#C7D6F0", "#EBB0A6")
# Creating the gradient function
color_gradient <- colorRampPalette(colors_version)
# Generate 30 discrete colors from this gradient
gradient_colors <- color_gradient(30)

# The users can log scale or not the data
log_humans_molecules_counts_by_bodypart <- log2(1 + humans_molecules_counts_by_bodypart_mat)
#write.csv(log_humans_molecules_counts_by_bodypart, 
#          file = "log_humans_molecules_counts_by_bodypart.csv", 
#          row.names = TRUE)

# Organ distribution in humans
## Use heatmap for data visualization or organ distribution - humans
### If one MS/MS spectrum is used in reverse metabolomics, the cluster_rows and cluster_cols should be set to FALSE
Organ_humans <- pheatmap(log_humans_molecules_counts_by_bodypart,
                         color = gradient_colors,
                         cluster_rows = FALSE,
                         cluster_cols = TRUE,
                         angle_col = "90",
                         main = "Organ distribution in humans",
                         fontsize = 10,
                         cellwidth = 15,
                         cellheight = 15,
                         treeheight_row = 100,
                         fontsize_row = 12,
                         fontsize_col = 12,
                         #legend_fontsize = 10,
                         border_color = NA)
Organ_humans
#ggsave("Organ_distribution_in_humans_5ASA.svg", plot = Organ_humans, width = 10, height = 10, dpi = 900, bg     = "transparent"  )
#getwd()

# Health phenotype association
## Filter for human information in the ReDU metadata
df_redu_humans <- redu_metadata |>  
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens")

# Humans - filter for lifestage, DOIDCommonName, HealthStatus, BiologicalSex
## For LifeStage
human_ReDU_LifeStage <- df_redu_humans |> 
  dplyr::count(LifeStage) |> 
  dplyr::rename(LifeStage_counts = n, LifeStage = LifeStage)
human_ReDU_LifeStage$LifeStage_counts <- as.numeric(human_ReDU_LifeStage$LifeStage_counts)
# For DOIDCommonName
human_ReDU_DOIDCommonName <- df_redu_humans |> 
  dplyr::count(DOIDCommonName) |> 
  dplyr::rename(DOIDCommonName_counts = n, DOIDCommonName = DOIDCommonName)
human_ReDU_DOIDCommonName$DOIDCommonName_counts <- as.numeric(human_ReDU_DOIDCommonName$DOIDCommonName_counts)
  # For HealthStatus
human_ReDU_HealthStatus <- df_redu_humans |> 
  dplyr::count(HealthStatus) |> 
  dplyr::rename(HealthStatus_counts = n, HealthStatus = HealthStatus)
human_ReDU_HealthStatus$HealthStatus_counts <- as.numeric(human_ReDU_HealthStatus$HealthStatus_counts)
# For BiologicalSex
human_ReDU_BiologicalSex <- df_redu_humans |> 
  dplyr::count(BiologicalSex) |> 
  dplyr::rename(BiologicalSex_counts = n, BiologicalSex = BiologicalSex)
human_ReDU_BiologicalSex$BiologicalSex_counts <- as.numeric(human_ReDU_BiologicalSex$BiologicalSex_counts)

df_humans_filled <- df_humans %>%
  dplyr::mutate(
    DOIDCommonName = as.character(DOIDCommonName),
    HealthStatus   = as.character(HealthStatus),
    DOIDCommonName = if_else(
      str_to_lower(str_trim(DOIDCommonName)) %in% c("missing value", "if indicated"),
      HealthStatus,
      DOIDCommonName))

grouped_df_humans <- df_humans_filled %>%
  group_by(Compound, DOIDCommonName, HealthStatus) %>%
  summarise(Count = n(), .groups = "drop") |> 
  ungroup() |>
  dplyr::select(-HealthStatus) #|> dplyr::filter(Count >= 2)

grouped_df_humans_unique <- grouped_df_humans %>%
  dplyr::group_by(Compound, DOIDCommonName) %>%
  dplyr::summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")


grouped_df_humans_pivot_table <- grouped_df_humans_unique |>
  tidyr::pivot_wider(names_from = Compound, values_from = Count, values_fill = list(Count = 0))

mat_counts <- grouped_df_humans_pivot_table %>%
  filter(DOIDCommonName != "missing value") %>%
  column_to_rownames("DOIDCommonName") %>%
  rename_with(function(nm) {
    nm <- sub("_known.*$", "", nm, ignore.case = TRUE)
    nm <- gsub("(?i)^[0-9]+_(?=5\\s*[-_]?ASA)", "", nm, perl = TRUE)
    nm <- gsub("(?i)^[0-9]+_(?=5\\s*[-_]?aminosalicylic)", "", nm, perl = TRUE)
    nm <- gsub("(?i)[_-]acid\\b", " acid", nm, perl = TRUE)
    nm
  }) %>% as.matrix()

#write_csv(mat_counts, "5_ASA_Humans_Diseases_molecules_counts.csv")

storage.mode(mat_counts) <- "numeric"

# drop all-zero rows/cols (optional but helps visibility)
mat_counts <- mat_counts[rowSums(mat_counts, na.rm = TRUE) > 0, , drop = FALSE]
mat_counts <- mat_counts[, colSums(mat_counts, na.rm = TRUE) > 0, drop = FALSE]

# 2) Per-column frequency in [0,1]: value / column sum
col_den <- colSums(mat_counts, na.rm = TRUE)
# avoid division by zero just in case (already dropped above)
col_den[col_den == 0] <- 1
mat_freq <- sweep(mat_counts, 2, col_den, `/`)
mat_freq[!is.finite(mat_freq)] <- 0  # replace NaN/Inf with 0

# 3) Palette & breaks: 0 = white, (0,1] = dark blue → red (no white in between)
pal <- c("#FFFFFF", colorRampPalette(c("#C7D6F0", "#EBB0A6"))(255))
brks <- seq(0, 1, length.out = length(pal) + 1)   # 256 colors → 257 breaks


Diseases_humans <- pheatmap::pheatmap(
  mat_freq,
  color = pal,
  breaks = brks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  angle_col = "90",
  fontsize = 6,
  cellwidth = 7,
  cellheight = 7,
  treeheight_row = 10,
  treeheight_col = 8,
  border_color = "white",
  na_col = "#F2F2F2")
  #main = "Health phenotype association (column frequency 0–1)")
Diseases_humans

# Export the data 
#ggsave("Diseases_humans_5ASA_derivatives.svg", plot = Diseases_humans, width = 10, height = 10, dpi = 900, bg = "transparent")
#getwd()




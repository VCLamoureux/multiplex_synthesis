setwd("/Users/vincentlamoureux/Library/CloudStorage/OneDrive-UniversityofCalifornia,SanDiegoHealth/Abubaker/")

# load packages
library(tidyverse)
library(data.table) 
library(pheatmap)    
library(duckplyr)

# import the quant table
feature_table <- readr::read_csv("urine_samples/Urine_sample_quant.csv") %>%
  dplyr::rename_with(~ gsub(" Peak area", "", .x, fixed = TRUE))

# import the annotation table
annotations <- data.table::fread("urine_samples/6817a5709ddf4ef2b0ef70c12057a27e-merged_results_with_gnps.tsv") |>
  as.data.frame() |>
  dplyr::rename(Feature = `#Scan#`) |>
  dplyr::mutate(Feature = as.numeric(Feature))

# transpose the table
data_transpose <- feature_table |>
  tibble::column_to_rownames("row ID") |>
  dplyr::select(tidyselect::contains(".mzML")) |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column("filename")
data_transpose$filename <- gsub(".mzML", "", data_transpose$filename)

# change the feature ID with the name of the conjugates
id_to_name <- c(
  "10014" = "carnitine-ibuprofen",
  "8780"  = "Ibuprofen-lysine",
  "7133"  = "Ibuprofen-histidine",
  "7539"  = "Ibuprofen-histidine_2",
  "10567" = "Ibuprofen-glutamine",
  "11425" = "Ibuprofen",
  "10968" = "Carboxyibuprofen")

# keep Ibuprofen derivatives
keep_ids <- c("10014", "8780", "7133", "7539", "10567", "11425", "10968")

data_synthesis <- data_transpose %>%
  dplyr::select(filename, all_of(keep_ids)) %>%
  dplyr::filter(!str_detect(filename, "BLANK")) %>%
  dplyr::rename_with(~ unname(replace(.x, .x %in% names(id_to_name), id_to_name[.x])), .cols = all_of(keep_ids)) %>%
  dplyr::mutate(filename = str_replace(filename, "^AB_(\\d+)$", "Subject\\1"))

row_order <- paste0("Subject", 1:9)
mat_raw <- data_synthesis %>%
  dplyr::filter(filename %in% row_order) %>%
  dplyr::mutate(filename = factor(filename, levels = row_order)) %>%
  arrange(filename) %>%
  column_to_rownames("filename") %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

col_max <- apply(mat_raw, 2, function(x) if(all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
mat_norm <- sweep(mat_raw, 2, col_max, "/")
mat_norm[!is.finite(mat_norm)] <- NA
mat_norm <- pmin(pmax(mat_norm, 0), 1)

pal  <- c("#FFFFFF", colorRampPalette(c("#C7D6F0", "#EBB0A6"))(255))
brks <- seq(0, 1, length.out = length(pal) + 1)

# Generate the scale heatmap with raw relative intensity
HM_synthesis <- pheatmap(
  mat_norm,
  color            = pal,
  breaks           = brks,
  cluster_rows     = FALSE,
  cluster_cols     = TRUE,
  angle_col        = "90",
  fontsize         = 6,
  cellwidth        = 10,
  cellheight       = 10,
  treeheight_row   = 10,
  treeheight_col   = 8,
  border_color     = "white",
  na_col           = "#F2F2F2")

HM_synthesis
#ggsave("Pregnant_women_urine.pdf", plot = HM_synthesis, width = 10, height = 10, dpi = 900, bg = "transparent")
#getwd()


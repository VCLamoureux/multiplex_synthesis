setwd("/Users/vincentlamoureux/Library/CloudStorage/")

# import packages
library(tidyverse)
library(ggrepel)

# Import tables
feature_table <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/MZmine_3_01292025/fbmn_iimn_fbmn_quant.csv")
feature_table <- feature_table |> 
  dplyr::rename_with(~gsub(" Peak area", "", .))
colnames(feature_table)

# import retention time drift and peak area drift table
drift <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/MZmine_3_01292025/fbmn_quant_mzmine.csv")
colnames(drift)[1] <- "Feature"

# import ReDU formatted metadata
metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/ReDU_Diet_Nov19_Template.csv")
metadata <- metadata |> 
  dplyr::select(-MassiveID) 

# import clinical metadata and format them
clinical_metadata <- read_csv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/clinical_data_20patients_updated_may20_final_response_plus_diet_score.csv")
clinical_metadata_renamed_1 <- clinical_metadata |> 
  dplyr::mutate(Patient = paste0("F", Patient))
clinical_metadata_renamed_1$Patient <- gsub("FD010d-0", "FD010d0", clinical_metadata_renamed_1$Patient)
clinical_metadata_renamed_1$Patient <- gsub("FD001D0", "FD001d0", clinical_metadata_renamed_1$Patient)
clinical_metadata_renamed_plasma <- clinical_metadata |> 
  dplyr::mutate(Patient = paste0("P", Patient))
clinical_metadata_renamed_plasma$Patient <- gsub("PD010d-0", "PD010d0", clinical_metadata_renamed_plasma$Patient)
clinical_metadata_renamed_plasma$Patient <- gsub("PD001D0", "PD001d0", clinical_metadata_renamed_plasma$Patient)
clinical_metadata_renamed_2 <- clinical_metadata_renamed_1 |> 
  dplyr::mutate(Patient = str_replace(Patient, "d0", "_T2"),
                Patient = str_replace(Patient, "d14", "_T3"))
clinical_metadata_renamed_plasma_2 <- clinical_metadata_renamed_plasma |> 
  dplyr::mutate(Patient = str_replace(Patient, "d0", "_T2"),
                Patient = str_replace(Patient, "d14", "_T3"))
clinical_metadata_renamed_2 <- clinical_metadata_renamed_2 |> 
  dplyr::rename(filename = Patient)
clinical_metadata_renamed_plasma_2 <- clinical_metadata_renamed_plasma_2 |> 
  dplyr::rename(filename = Patient)

# import annotation table - multiplex library
annotation <- read_tsv("OneDrive-UniversityofCalifornia,SanDiegoHealth/Postdoc_UCSD/Postdoc_projects/Marta_Guma/annotation_multiplex_Oct19.tsv")
annotation$`#Scan#` <- as.character(annotation$`#Scan#`)
annotation <- annotation |> 
  dplyr::rename(Feature = `#Scan#`)
annotation$Feature <- as.character(annotation$Feature)

# get info feature
info_feature <- feature_table[, 1:3]
info_feature <- info_feature |> 
  dplyr::rename(Feature = `row ID`)
info_feature$Feature <- as.character(info_feature$Feature)
colnames(info_feature) <- c("Feature", "mz", "RT")
info_feature_name <- info_feature |> 
 left_join(annotation, by = "Feature")

# transpose data
data_transpose <- feature_table |>
  column_to_rownames("row ID") |>    
  dplyr::select(contains(".mzXML")) |> 
  t() |>  
  as.data.frame() |> 
  rownames_to_column("filename")

# Check retention drift on internal standard Sulfamethizole
feature_row <- drift |>  
  dplyr::filter(Feature == 2394)
rt_columns <- feature_row |> 
  dplyr::select(contains("Feature RT"))
rt_data <- tidyr::gather(rt_columns, key = "filename", value = "RT") 

rt_data_cleaned <- rt_data |> 
  dplyr::filter(!str_detect(filename, "Blank"))
rt_data_cleaned$RT <- as.numeric(rt_data_cleaned$RT)

rt_stats <- rt_data_cleaned |>
  summarise(
    Mean_RT = mean(RT, na.rm = TRUE),
    SD_RT = sd(RT, na.rm = TRUE),
    Variance_RT = var(RT, na.rm = TRUE))

ggplot(rt_data_cleaned, aes(x = filename, y = RT)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "Retention Time (RT)", title = "Retention Time Drift for Feature 2394")

# Check for intensity drift on internal standard Sulfamethizole
feature_row_peakarea <- drift |> 
  dplyr::filter(Feature == 2394)

peakarea_columns <- feature_row_peakarea |> 
  dplyr::select(contains("Peak area"))
peakarea_data <- tidyr::gather(peakarea_columns, key = "filename", value = "Peak")
peakarea_data_cleaned <- peakarea_data |> 
  dplyr::filter(!str_detect(filename, "Blank"))

peakarea_stats <- peakarea_data_cleaned |>
  summarise(
    Mean_peakarea = mean(Peak, na.rm = TRUE),
    SD_peakarea = sd(Peak, na.rm = TRUE),
    Variance_peakarea = var(Peak, na.rm = TRUE))
ggplot(peakarea_data_cleaned, aes(x = filename, y = Peak)) +
  geom_point() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = "Filename", y = "peakarea", title = "peakarea Drift for Feature 2394")


# Filter for QC sixmix
data_sixmix <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "6Mix", filename))
data_sixmix$filename <- gsub(".mzXML", "", data_sixmix$filename)
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix |>  column_to_rownames("filename") |> colMeans(),
                                  SD_sixmix =  data_sixmix |>  column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) |>
  dplyr::filter(Mean_sixmix > 0) |> arrange(desc(Mean_sixmix))
sixmix_feature_info$Feature <- as.character(sixmix_feature_info$Feature)
sixmix_feature_info_name <- left_join(sixmix_feature_info, info_feature_name, by = "Feature" ) |> 
  dplyr::select("Feature", "Mean_sixmix", "SD_sixmix", "CV_sixmix", "mz", "RT", "Compound_Name")

# Filter for blank
data_blank <- data_transpose |>  
  dplyr::filter(str_detect(pattern = "Blank", filename))
blank_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                 Mean_blank = data_blank |> column_to_rownames("filename") |> colMeans(),
                                 SD_blank =  data_blank |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) |>
  dplyr::filter(Mean_blank > 0) |> arrange(desc(Mean_blank))

blank_feature_info_name <- left_join(blank_feature_info, info_feature_name, by = "Feature" )

# Filter for samples
data_sample <- data_transpose |> 
  dplyr::filter(str_detect(filename, "PD|almond|black|bread|chia|FD|flaxseeds|ginger|miso|oats|sesame|tahini|turmeric|vinegar"))
sample_feature_info <- data.frame(Feature = colnames(data_sample)[-1],
                                  Mean_sample = data_sample |> column_to_rownames("filename") |> colMeans(),
                                  SD_sample =  data_sample |> column_to_rownames("filename") |> apply(2, sd)) |>
  dplyr::mutate(CV_sample = SD_sample/Mean_sample) |>
  dplyr::filter(Mean_sample > 0) |> arrange(desc(Mean_sample))
sample_feature_info_name <- left_join(sample_feature_info, info_feature_name, by = "Feature" ) 

specified_features <- c("5552", "2324", "4471","2394", "2872", "7140")

# Features to be removed samples/QCmix < 5
feature_to_remove_qcmix <- sixmix_feature_info |> 
  left_join(sample_feature_info) |> 
  dplyr::filter(Mean_sixmix > 0) |> 
  dplyr::mutate(sample_Mix = Mean_sample/Mean_sixmix) |> 
  dplyr::filter(sample_Mix < 5 | is.na(sample_Mix)) |> 
  dplyr::bind_rows(blank_feature_info |> 
                     dplyr::filter(Feature %in% specified_features))

# Data with QCmix removal, blanks
data_clean2 <- data_transpose |> 
  dplyr::select(-all_of(feature_to_remove_qcmix$Feature)) |> 
  dplyr::filter(!str_detect(filename, "6Mix|Blank"))

data_clean_transpose <- data_clean2 |>
  column_to_rownames("filename") |>
  t() |> 
  as.data.frame() |> 
  rownames_to_column("Feature")

merged <- left_join(data_clean_transpose, annotation, by = "Feature")
merged$Feature <- as.character(merged$Feature)

merged_mzml_mean_all <- merged %>%
  dplyr::select(
    any_of(c("Compound_Name", "Feature", "LibraryName")), 
    matches("_T1|_T2|_T3"))

mzml_columns <- grep("T1|T2|T3", names(merged_mzml_mean_all), value = TRUE)

df_melted <- reshape2::melt(merged_mzml_mean_all, 
                            id.vars = "Feature", 
                            measure.vars = mzml_columns,
                            variable.name = "Sample_Name",  
                            value.name = "Peak_Area")

df_melted_feces <- df_melted %>%
  dplyr::filter(str_detect(Sample_Name, "FD"))
df_melted_feces$Sample_Name <- gsub(".mzXML", "", df_melted_feces$Sample_Name)

df_melted_feces_clinical <- df_melted_feces |>
  left_join(clinical_metadata_renamed_2, by = c("Sample_Name" = "filename"))

# re-label the feature by compound_name
name_map <- c(
  "2428"  = "5-Aminosalicylic acid_propanoic-acid",
  "5455"  = "5-aminosalicylic acid_phenylpropionic acid",
  "10709" = "5-Aminosalicylic acid_palmitic-acid",
  "10832" = "5-Aminosalicylic acid_oleic-acid_iso1",
  "9457"  = "5-Aminosalicylic acid_oleic-acid_iso2",
  "4338"  = "5-Aminosalicylic acid_isovaleric_acid",
  "8419"  = "5-Aminosalicylic acid_linolenic-acid_iso1",
  "10023" = "5-Aminosalicylic acid_linolenic-acid_iso2",
  "1952"  = "5-Aminosalicylic acid_succinic_acid")

# keep only features you map (optional but safer)
feature_ids <- names(name_map)

df_melted_metadata_annotation <- df_melted_feces_clinical |>
  dplyr::filter(Feature %in% feature_ids) |>
  dplyr::filter(!stringr::str_detect(Sample_Name, "T1|T2")) |>
  dplyr::select(Feature, Sample_Name, Peak_Area, Pain_50_improvement, DMARD, dplyr::everything()) |>
  dplyr::mutate(
    sulfasalazine = dplyr::if_else(Sample_Name %in% c("FD003_T3", "FD028_T3"), "user", "non-user"),
    Peak_Area_Log = log10(Peak_Area + 1)) |>
  dplyr::left_join(annotation, by = "Feature") |>
  dplyr::mutate(
  Feature = dplyr::recode(as.character(Feature), !!!name_map, .default = as.character(Feature)),
  Feature_Compound = paste(Feature, Compound_Name, sep = "_")) |>
  dplyr::select(Feature_Compound, dplyr::everything())

# export the formatted table, if needed
#write_csv(df_melted_metadata_annotation, "df_melted_metadata_annotation_RA_COHORT.csv")




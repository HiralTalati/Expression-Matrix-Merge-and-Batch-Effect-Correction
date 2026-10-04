#################################################################
# STEP 8: MASTER MERGE & AUTO-CLINICAL MAPPING
#################################################################

library(openxlsx)
library(dplyr)
library(purrr)
library(tibble)

# --- 1. DEFINE PATHS ---
folders <- c(
  "GSE63060"  = "GSE63060_Final_Results_20260314_171513", 
  "GSE63061"  = "GSE63061_Final_Results_20260314_171611", 
  "GSE140829" = "GSE140829_Final_Results_20260314_180745"
)

cat("[*] Starting Master Merge with Auto-Mapping...\n")

# --- 2. LOAD MATRICES AND EXTRACT GROUPS ---
mat_list <- list()
full_group_map <- data.frame()

for (id in names(folders)) {
  f_path <- folders[[id]]
  
  # A. Load the Expression Matrix
  m_file <- list.files(f_path, pattern = "_06_EXPRESSION_MATRIX.xlsx", full.names = TRUE)
  df_mat <- read.xlsx(m_file)
  mat_list[[id]] <- df_mat
  
  # B. AUTO-RECOVERY OF GROUPS
  # We look for the 'ANNOTATED_ALL' file which contains 'Significance' 
  # and other markers we can use to verify groups.
  a_file <- list.files(f_path, pattern = "_04_ANNOTATED_ALL.xlsx", full.names = TRUE)
  # We read just the first sheet to get the sample list
  df_annot <- read.xlsx(a_file, sheet = 1)
  
  # Every column that starts with 'GSM' is a sample.
  samples_in_this_batch <- colnames(df_mat)[-1]
  
  # IMPORTANT: In your specific studies, the groups are usually identified 
  # in the GEO metadata. Since we already standardized them in the first script, 
  # we will pull them from the column names or the batch structure.
  
  # Let's create a temporary map for this batch
  temp_map <- data.frame(
    SampleID = samples_in_this_batch,
    Batch = id,
    stringsAsFactors = FALSE
  )
  
  # --- LOGIC TO IDENTIFY AD/CTL ---
  # If the GSM IDs follow a known pattern or if you have a pheno file:
  # Since we are running this post-pipeline, let's pull the 'group_clean' 
  # by briefly looking at the mean expression or a marker gene if needed,
  # BUT the best way is to re-access the GEO 'pData'.
  
  library(GEOquery)
  gse_obj <- getGEO(id, GSEMatrix = TRUE, AnnotGPL = FALSE)
  p_data <- pData(gse_obj[[1]])
  
  # Standardizing based on your specific dataset columns
  col_name <- ifelse(id == "GSE140829", "diagnosis:ch1", "status:ch1")
  
  p_data_clean <- p_data %>%
    rownames_to_column("SampleID") %>%
    select(SampleID, group_clean = !!sym(col_name)) %>%
    mutate(group_clean = case_when(
      group_clean %in% c("Control", "CTL", "control") ~ "CTL",
      group_clean %in% c("MCI", "borderline MCI") ~ "MCI",
      group_clean %in% c("AD", "Alzheimer's Disease") ~ "AD",
      TRUE ~ "Other"
    ))
  
  full_group_map <- rbind(full_group_map, p_data_clean)
  cat(" -> Successfully mapped groups for", id, "\n")
}

# --- 3. MERGE EXPRESSION ---
merged_data <- mat_list %>% reduce(inner_join, by = "EntrezID")

# --- 4. CREATE FINAL BATCH METADATA ---
batch_info <- data.frame(SampleID = colnames(merged_data)[-1]) %>%
  left_join(full_group_map, by = "SampleID") %>%
  mutate(Batch = case_when(
    SampleID %in% colnames(mat_list[[1]]) ~ "GSE63060",
    SampleID %in% colnames(mat_list[[2]]) ~ "GSE63061",
    SampleID %in% colnames(mat_list[[3]]) ~ "GSE140829"
  ))

# --- 5. EXPORT ---
write.xlsx(merged_data, "07_MASTER_MERGED_EXPRESSION.xlsx")
write.xlsx(batch_info, "08_BATCH_METADATA.xlsx")

cat("\n✅ SUCCESS: 'group_clean' is no longer pending! Ready for Batch Correction.\n")
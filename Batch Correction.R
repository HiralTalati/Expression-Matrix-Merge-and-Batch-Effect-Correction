#################################################################
# STEP 9: PROTECTED BATCH AUDIT & CORRECTION
#################################################################

library(openxlsx)
library(limma)
library(sva)
library(ggplot2)
library(tibble)

# --- 1. SETUP OUTPUT DIRECTORY ---
output_folder <- "After_Batch_Correction_Results"
if (!dir.exists(output_folder)) dir.create(output_folder)

cat("[*] Initializing Protected Batch Correction Phase...\n")

# --- 2. LOAD THE GENERATED FILES ---
expr_data <- read.xlsx("07_MASTER_MERGED_EXPRESSION.xlsx") %>% column_to_rownames("EntrezID")
batch_meta <- read.xlsx("08_BATCH_METADATA.xlsx")

# Ensure column order in expression matrix matches metadata rows
expr_data <- expr_data[, batch_meta$SampleID]

# --- 3. PRE-CORRECTION AUDIT ---
cat("\n[*] STEP 2: Auditing Pre-Correction Batch Variance...")
pca_pre <- prcomp(t(expr_data), scale. = TRUE)
pc1_pre <- pca_pre$x[,1]
r_squared_pre <- summary(lm(pc1_pre ~ as.factor(batch_meta$Batch)))$r.squared * 100

cat(sprintf("\n    -> Pre-Correction: Batch explains %.2f%% of PC1 variance.", r_squared_pre))

# Save Pre-Correction Plot
pca_df_pre <- data.frame(PC1=pca_pre$x[,1], PC2=pca_pre$x[,2], Batch=batch_meta$Batch)
p1 <- ggplot(pca_df_pre, aes(PC1, PC2, color=Batch)) + 
  geom_point(size=3, alpha=0.7) + theme_minimal() + 
  labs(title="PCA: Pre-Batch Correction", subtitle=paste("Batch Influence:", round(r_squared_pre, 2), "%"))
ggsave(file.path(output_folder, "09_PCA_Pre_Correction.png"), p1, width=7, height=5)

# --- 4. EXECUTE PROTECTED CORRECTION ---
# This Model Matrix protects your biological groups during the noise removal
mod <- model.matrix(~as.factor(group_clean), data = batch_meta)

cat("\n[*] Applying Protected ComBat Correction (preserving AD/MCI/CTL signals)...")
final_matrix <- ComBat(dat = as.matrix(expr_data), 
                       batch = batch_meta$Batch, 
                       mod = mod, 
                       par.prior = TRUE)

# --- 5. POST-CORRECTION AUDIT ---
cat("\n[*] STEP 5: Auditing Post-Correction Results...")
pca_post <- prcomp(t(final_matrix), scale. = TRUE)
pc1_post <- pca_post$x[,1]
r_squared_post <- summary(lm(pc1_post ~ as.factor(batch_meta$Batch)))$r.squared * 100

cat(sprintf("\n    -> Post-Correction: Batch explains %.2f%% of PC1 variance.", r_squared_post))

# Save Post-Correction Plot
pca_df_post <- data.frame(PC1=pca_post$x[,1], PC2=pca_post$x[,2], Batch=batch_meta$Batch, Group=batch_meta$group_clean)
p2 <- ggplot(pca_df_post, aes(PC1, PC2, color=Batch, shape=Group)) + 
  geom_point(size=3, alpha=0.7) + theme_minimal() + 
  labs(title="PCA: Post-Batch Correction", subtitle=sprintf("Remaining Batch Influence: %.2f%%", r_squared_post))
ggsave(file.path(output_folder, "10_PCA_Post_Correction.png"), p2, width=7, height=5)

# --- 6. SAVE FINAL FILES ---
cat("\n[*] Saving final modeling data to folder...")
final_df <- as.data.frame(final_matrix) %>% rownames_to_column("EntrezID")

write.xlsx(final_df, file.path(output_folder, "11_FINAL_BATCH_CORRECTED_MATRIX.xlsx"))

# Save a comparison log for your thesis
audit_log <- data.frame(
  Metric = c("Pre-Correction Variance", "Post-Correction Variance"),
  Value = c(paste0(round(r_squared_pre, 2), "%"), paste0(round(r_squared_post, 2), "%")),
  Status = c("Severe Batch Effect", "Cleaned/Harmonized")
)
write.xlsx(audit_log, file.path(output_folder, "12_Batch_Audit_Comparison.xlsx"))

cat(paste("\n✅ SUCCESS: All results saved in:", output_folder, "\n"))
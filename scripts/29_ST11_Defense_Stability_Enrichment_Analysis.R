#  ST11 Defense and Stability Enrichment Analysis
#  Identifies gene gain (fixation) and gene loss (pruning) in the 
# East Asian ST11 lineage compared to global sporadic isolates.

# --- 1. Environment Setup ----------------------------------------------------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")

library(data.table)
library(dplyr)
library(tidyr)

# --- 2. Configuration and Paths ----------------------------------------------
# Define paths to pangenome and metadata files
PANGENOME_PATH <- "gene_presence_absence.csv"
GLOBAL_META    <- "st_comparison_report_globalST11s_29isolates.csv"
EA_META        <- "st_comparison_report_East_Asian_275_isolates.csv"

# --- 3. Data Loading and Pre-processing --------------------------------------
# Helper: Clean GCF/ID suffixes
clean_ids <- function(x) gsub("\\..*$", "", x)

# Load metadata and extract ST11 IDs
global_ids <- clean_ids(fread(GLOBAL_META)[st_kleborate == 11, id])
ea_ids     <- clean_ids(fread(EA_META)[st_kleborate == 11, id])

# Read pangenome headers to identify relevant isolate columns
headers <- colnames(fread(PANGENOME_PATH, nrows = 0))
header_bases <- clean_ids(headers)
isolate_cols <- headers[header_bases %in% c(global_ids, ea_ids)]

# Load pangenome data for targeted isolates
pan_data <- fread(PANGENOME_PATH, select = c("Gene", "Non-unique Gene name", "Annotation", isolate_cols))

# --- 4. Defense System Mining ------------------------------------------------
# Define broad keywords to capture R-M, TA, CRISPR, and Stability systems
keywords <- paste0("\\b(restriction|modification|methyltransferase|toxin|antitoxin|",
                   "crispr|cas[0-9]|abortive|brex|pifm|par[ab]|vapC|hipA|relE|hok|pemK|pemI)\\b")

is_defense <- grepl(keywords, pan_data$Annotation, ignore.case = TRUE) |
              grepl(keywords, pan_data$`Non-unique Gene name`, ignore.case = TRUE)
defense_genes <- pan_data[is_defense, ]

# --- 5. Statistical Analysis (Fisher's Exact Test) ---------------------------
# Map column indices to groups
iso_names <- names(defense_genes)[-(1:3)]
iso_clean <- clean_ids(iso_names)
ea_idx    <- which(iso_clean %in% ea_ids)
gl_idx    <- which(iso_clean %in% global_ids)

res_list <- lapply(1:nrow(defense_genes), function(i) {
  row_vals <- as.character(defense_genes[i, -(1:3), with = FALSE])
  binary   <- ifelse(is.na(row_vals) | row_vals == "", 0L, 1L)
  
  # Group counts
  ea_pres <- sum(binary[ea_idx]); ea_n <- length(ea_idx)
  gl_pres <- sum(binary[gl_idx]); gl_n <- length(gl_idx)
  
  # 2x2 Contingency Table
  mat <- matrix(c(ea_pres, ea_n - ea_pres, gl_pres, gl_n - gl_pres), nrow = 2, byrow = TRUE)
  
  # Fisher's Exact Test for p-value
  p_val <- fisher.test(mat)$p.value
  
  # Odds Ratio with Haldane-Anscombe Correction (handles zeros)
  mat_h <- mat + 0.5
  or    <- (mat_h[1,1] * mat_h[2,2]) / (mat_h[1,2] * mat_h[2,1])
  se_log_or <- sqrt(sum(1/mat_h))
  
  # 95% Confidence Intervals
  or_lo <- exp(log(or) - 1.96 * se_log_or)
  or_hi <- exp(log(or) + 1.96 * se_log_or)

  data.frame(
    gene_id    = defense_genes$Gene[i],
    gene_name  = defense_genes$`Non-unique Gene name`[i],
    annotation = defense_genes$Annotation[i],
    ea_freq    = round(ea_pres/ea_n, 3),
    gl_freq    = round(gl_pres/gl_n, 3),
    p_value    = p_val,
    odds_ratio = round(or, 3),
    or_lo      = round(or_lo, 3),
    or_hi      = round(or_hi, 3),
    haldane    = any(mat == 0)
  )
})

# Compile results and apply False Discovery Rate (FDR) correction
results <- bind_rows(res_list) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  arrange(padj)

# --- 6. Lineage-Specific Filters ---------------------------------------------
# Positive Selection: Fixed in EA, sporadic in Global
fixed_gatekeepers <- results %>% filter(padj < 0.05, ea_freq > 0.70, odds_ratio > 4)

# Negative Selection: Present in Global, pruned in EA
pandemic_losses <- results %>% filter(padj < 0.05, gl_freq > 0.35, ea_freq < 0.10)

# --- 7. Export Results -------------------------------------------------------
fwrite(results, "ST11_defense_enrichment_full.csv")
fwrite(fixed_gatekeepers, "ST11_fixed_gatekeepers.csv")
fwrite(pandemic_losses, "ST11_pandemic_losses.csv")

message("Analysis Complete. Results saved for fixed gatekeepers and pandemic losses.")

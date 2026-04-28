# ============================================================
# Mendelian Randomization Analysis
# Methods: IVW, Weighted Median, MR-Egger, RAPS, CMM
# ============================================================

library(TwoSampleMR)
library(mr.raps)
library(MendelianRandomization)

# ========================= 1. Load Data =========================
INPUT_PATH  <- "./data/your_gwas_data.csv"   # Path to harmonized instrument variable CSV
OUTPUT_PATH <- "./results/mr/your_results.csv"   # Path to save results CSV

dat <- read.csv(INPUT_PATH, fileEncoding = "UTF-8",
                stringsAsFactors = FALSE, sep = "\t")

# Fallback: retry with comma separator if tab-separated read failed
if (ncol(dat) <= 1) {
  dat <- read.csv(INPUT_PATH, fileEncoding = "UTF-8",
                  stringsAsFactors = FALSE)
}

dat$mr_keep          <- TRUE
dat$mr_keep.exposure <- TRUE
dat$mr_keep.outcome  <- TRUE

cat("Total rows:", nrow(dat), "\n")

# ========================= 2. Initialize Results =========================
all_results <- data.frame()

# ========================= 3. Run MR Per Outcome =========================
outcomes <- unique(dat$id.outcome)

for (oc in outcomes) {
  
  dat_sub       <- dat[dat$id.outcome == oc, ]
  id_exp        <- dat_sub$id.exposure[1]
  id_out        <- dat_sub$id.outcome[1]
  outcome_name  <- dat_sub$outcome[1]
  exposure_name <- dat_sub$exposure[1]
  n_snp         <- nrow(dat_sub)
  
  cat("\nProcessing:", outcome_name, "(n =", n_snp, ")\n")
  
  # --- IVW, Weighted Median, MR-Egger ---
  mr_res <- mr(dat_sub, method_list = c(
    "mr_ivw",
    "mr_weighted_median",
    "mr_egger_regression"
  ))
  
  for (i in 1:nrow(mr_res)) {
    b  <- mr_res$b[i]
    se <- mr_res$se[i]
    all_results <- rbind(all_results, data.frame(
      id.exposure = id_exp,
      id.outcome  = id_out,
      outcome     = outcome_name,
      exposure    = exposure_name,
      method      = mr_res$method[i],
      nsnp        = mr_res$nsnp[i],
      b           = b,
      se          = se,
      pval        = mr_res$pval[i],
      lo_ci       = b - 1.96 * se,
      up_ci       = b + 1.96 * se,
      or          = exp(b),
      or_lci95    = exp(b - 1.96 * se),
      or_uci95    = exp(b + 1.96 * se),
      stringsAsFactors = FALSE
    ))
  }
  
  # --- RAPS ---
  raps_res <- mr.raps(dat_sub$beta.exposure, dat_sub$beta.outcome,
                      dat_sub$se.exposure, dat_sub$se.outcome)
  b_raps  <- raps_res$beta.hat
  se_raps <- raps_res$beta.se
  p_raps  <- pnorm(-abs(b_raps / se_raps)) * 2
  
  all_results <- rbind(all_results, data.frame(
    id.exposure = id_exp,
    id.outcome  = id_out,
    outcome     = outcome_name,
    exposure    = exposure_name,
    method      = "Robust adjusted profile score (RAPS)",
    nsnp        = n_snp,
    b           = b_raps,
    se          = se_raps,
    pval        = p_raps,
    lo_ci       = b_raps - 1.96 * se_raps,
    up_ci       = b_raps + 1.96 * se_raps,
    or          = exp(b_raps),
    or_lci95    = exp(b_raps - 1.96 * se_raps),
    or_uci95    = exp(b_raps + 1.96 * se_raps),
    stringsAsFactors = FALSE
  ))
  
  # --- Contamination Mixture Method ---
  mr_input_obj <- mr_input(
    bx   = dat_sub$beta.exposure,
    bxse = dat_sub$se.exposure,
    by   = dat_sub$beta.outcome,
    byse = dat_sub$se.outcome
  )
  
  cmm_res <- mr_conmix(mr_input_obj)
  b_cmm   <- cmm_res$Estimate
  se_cmm  <- (cmm_res$CIUpper - cmm_res$CILower) / (2 * 1.96)
  
  all_results <- rbind(all_results, data.frame(
    id.exposure = id_exp,
    id.outcome  = id_out,
    outcome     = outcome_name,
    exposure    = exposure_name,
    method      = "Contamination mixture method",
    nsnp        = n_snp,
    b           = b_cmm,
    se          = se_cmm,
    pval        = cmm_res$Pvalue,
    lo_ci       = cmm_res$CILower,
    up_ci       = cmm_res$CIUpper,
    or          = exp(b_cmm),
    or_lci95    = exp(cmm_res$CILower),
    or_uci95    = exp(cmm_res$CIUpper),
    stringsAsFactors = FALSE
  ))
}

# ========================= 4. Export Results =========================
write.csv(all_results, OUTPUT_PATH, row.names = FALSE)

cat("\n========================================\n")
cat("Results saved to:", OUTPUT_PATH, "\n")
cat("Total rows:", nrow(all_results), "\n")
cat("========================================\n")
print(all_results)
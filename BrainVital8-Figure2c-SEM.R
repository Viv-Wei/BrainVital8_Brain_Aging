# ============================================================
# BrainVital8 Second-Order Structural Equation Model
# ============================================================
#
# Structure: 6 first-order latent variables + 2 manifest indicators
#            -> second-order latent variable BrainVital8
#
# Domain structure:
#   Glycaemic   <- Glucose, HbA1c
#   Grip        <- GripR, GripL
#   Physical    <- Stairs, WalkDist, WalkDaily
#   Social      <- Income, CommDiff
#   Emotional   <- Unhappiness, EmoImpact, Interpersonal
#   Respiratory <- FEV1, FVC
#   Sleep_z     <- manifest indicator (directly into 2nd-order)
#   Smoke_rev   <- manifest indicator (directly into 2nd-order)
# ============================================================

library(lavaan)
library(dplyr)
library(tidyr)
library(tibble)

# ============================================================
# 0) Path Configuration
# ============================================================
BASE_DIR  <- "./data/SEM_Analysis"
INPUT_UKB <- file.path(BASE_DIR, "UKB.csv")
INPUT_COV <- "./data/fourbrain_behavior_cov.csv"
OUT_DIR   <- file.path('./results/', "SEM_path_results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1) Read and Merge Two Datasets
# ============================================================
dat_ukb <- read.csv(INPUT_UKB, check.names = FALSE, fileEncoding = "UTF-8")
dat_cov <- read.csv(INPUT_COV, check.names = FALSE, fileEncoding = "UTF-8")

names(dat_ukb) <- make.names(names(dat_ukb), unique = TRUE)
names(dat_cov) <- make.names(names(dat_cov), unique = TRUE)

dat <- left_join(dat_ukb, dat_cov, by = "eid")
cat("Data merged:", nrow(dat), "rows,", ncol(dat), "columns\n")

# ============================================================
# 2) Utility: Z-score Standardization (preserves NA)
# ============================================================
z <- function(x) as.numeric(scale(x))

# ============================================================
# 3) Variable Preprocessing + Direction Alignment
#    (higher = healthier for all indicators)
# ============================================================
dat <- dat %>%
  mutate(
    # ---------- Glycaemic Control ----------
    Glu_rev   = -z(Glucose.levels),           # Higher is worse -> reverse
    HbA1c_rev = -z(HbA1c),                    # Higher is worse -> reverse
    
    # ---------- Grip Strength ----------
    GripR_z = z(Hand.Grip.Strength.R),        # Higher is better
    GripL_z = z(Hand.Grip.Strength.L),        # Higher is better
    
    # ---------- Smoking (manifest, directly into 2nd-order) ----------
    Smoke_rev = -z(as.numeric(Smoke.status)), # More smoking is worse -> reverse
    
    # ---------- Sleep Health (manifest, directly into 2nd-order) ----------
    Sleep_z = z(Sleep.Periods),
    
    # ---------- Physical Function and Activity ----------
    Stairs_z     = z(as.numeric(Liking.for.Stairs.x)),
    WalkDist_rev = -z(as.numeric(              # Higher is worse -> reverse
      Difficulty.with.distance.walking)),
    WalkDaily_rev = -z(as.numeric(             # Higher is worse -> reverse
      Diff..Day.to.Day.Work)),
    
    # ---------- Social and Economic Context ----------
    Income_z     = z(Avg.HH.Income.Pre.Tax),   # Higher is better
    CommDiff_rev = -z(as.numeric(              # Higher is worse -> reverse
      Difficulty.with.joining.in.community.activities)),
    
    # ---------- Emotional Well-being ----------
    Unhappy_rev    = -z(Unhappiness.Score),    # Higher is worse -> reverse
    EmoImpact_rev  = -z(as.numeric(            # Higher is worse -> reverse
      Diff..Emotional.Impact..Health.)),
    Confide_z      = z(as.numeric(             # Higher is better
      Good.interpersonal.relationships)),
    
    # ---------- Respiratory Function ----------
    FEV1_z = z(FEV1),                          # Higher is better
    FVC_z  = z(FVC)                            # Higher is better
  )

# ============================================================
# 4) Missing Rate Check
# ============================================================
key_vars <- c(
  "Glu_rev", "HbA1c_rev",
  "GripR_z", "GripL_z",
  "Smoke_rev", "Sleep_z",
  "Stairs_z", "WalkDist_rev", "WalkDaily_rev",
  "Income_z", "CommDiff_rev",
  "Unhappy_rev", "EmoImpact_rev", "Confide_z",
  "FEV1_z", "FVC_z"
)

miss_check <- sapply(key_vars, function(v) {
  if (!v %in% names(dat)) return("ERROR: column not found")
  sprintf("%.1f%% missing", mean(is.na(dat[[v]])) * 100)
})

cat("\n--- Missing Rate Check for Key Variables ---\n")
print(data.frame(Variable = names(miss_check),
                 Status   = unname(miss_check)), row.names = FALSE)

if (any(grepl("ERROR", miss_check))) {
  stop("Missing columns detected. Please verify column names and re-run.")
}

# ============================================================
# 5) Second-Order Structural Equation Model
#    Sleep_z and Smoke_rev enter as manifest indicators
#    (lavaan fixes their residual variance to 0 by default)
# ============================================================
model_2nd <- '
  # ======== First-Order Measurement Model (6 latent variables) ========

  # 1. Glycaemic Control
  Glycaemic   =~ Glu_rev + HbA1c_rev

  # 2. Grip Strength
  Grip        =~ GripR_z + GripL_z

  # 3. Physical Function and Activity
  Physical    =~ Stairs_z + WalkDist_rev + WalkDaily_rev

  # 4. Social and Economic Context
  Social      =~ Income_z + CommDiff_rev

  # 5. Emotional Well-being
  Emotional   =~ Unhappy_rev + EmoImpact_rev + Confide_z

  # 6. Respiratory Function
  Respiratory =~ FEV1_z + FVC_z

  # ======== Second-Order Structural Model ========
  # 6 first-order latent variables + 2 manifest indicators -> BrainVital8
  BrainVital8 =~ Glycaemic + Grip + Physical + Social +
                 Emotional + Respiratory + Sleep_z + Smoke_rev
'

fit_2nd <- cfa(
  model_2nd,
  data      = dat,
  estimator = "MLR",    # Robust to non-normality
  missing   = "fiml",   # Full information maximum likelihood for missing data
  std.lv    = TRUE      # Standardize latent variable variances to 1
)

# Full output
summary(fit_2nd, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)

# Key fit indices
cat("\n--- Model Fit Indices ---\n")
print(fitMeasures(fit_2nd, c(
  "cfi", "tli",
  "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
  "srmr", "chisq", "df", "pvalue"
)))

# ============================================================
# 6) Extract Parameters
# ============================================================
pe      <- parameterEstimates(fit_2nd, standardized = TRUE)
all_ov  <- lavNames(fit_2nd, type = "ov")

first_order_domains <- c(
  "Glycaemic", "Grip", "Physical",
  "Social", "Emotional", "Respiratory"
)

# ============================================================
# 7) Table A: Second-Order Loadings
#    (Domains / Sleep_z / Smoke_rev -> BrainVital8)
# ============================================================
order2_loadings <- pe %>%
  filter(op == "=~", lhs == "BrainVital8") %>%
  transmute(
    Domain       = rhs,
    Loading_unst = round(est, 3),
    SE           = round(se, 3),
    Z            = round(z, 3),
    Pvalue       = round(pvalue, 4),
    Loading_std  = round(std.all, 3),
    Sig = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

cat("\n--- Second-Order Loadings (Domains / Sleep_z / Smoke_rev -> BrainVital8) ---\n")
print(order2_loadings)

# ============================================================
# 8) Table B: First-Order Loadings (Indicators -> Domains) + R-squared
# ============================================================
order1_loadings <- pe %>%
  filter(op == "=~", lhs %in% first_order_domains) %>%
  transmute(
    Domain       = lhs,
    Indicator    = rhs,
    Loading_unst = round(est, 3),
    SE           = round(se, 3),
    Z            = round(z, 3),
    Pvalue       = round(pvalue, 4),
    Loading_std  = round(std.all, 3),
    Sig = case_when(
      pvalue < 0.001 ~ "***",
      pvalue < 0.01  ~ "**",
      pvalue < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )

r2_df <- lavInspect(fit_2nd, what = "r2") %>%
  as.data.frame() %>%
  rownames_to_column("Indicator") %>%
  rename(R2 = ".")

order1_full <- order1_loadings %>%
  left_join(r2_df, by = "Indicator") %>%
  mutate(R2 = round(R2, 3))

cat("\n--- First-Order Loadings (Indicators -> Domains) + R-squared ---\n")
print(order1_full)

# ============================================================
# 9) Table C: Residual Variances
# ============================================================
resid_tbl <- pe %>%
  filter(op == "~~", lhs == rhs,
         lhs %in% c(all_ov, first_order_domains)) %>%
  transmute(
    Variable     = lhs,
    Type         = if_else(lhs %in% all_ov, "Observed", "Latent(1st-order)"),
    ResidVar     = round(est, 3),
    SE           = round(se, 3),
    ResidVar_std = round(std.all, 3)
  )

cat("\n--- Residual Variances ---\n")
print(resid_tbl)

# ============================================================
# 10) Path Coefficient Summary Table (for path diagrams)
# ============================================================
# Second-order paths
path_2nd <- pe %>%
  filter(op == "=~", lhs == "BrainVital8") %>%
  transmute(
    From    = lhs,
    To      = rhs,
    std.all = round(std.all, 3),
    est     = round(est, 3),
    pvalue  = round(pvalue, 4)
  )

# First-order paths
path_1st <- pe %>%
  filter(op == "=~", lhs %in% first_order_domains) %>%
  transmute(
    From    = lhs,
    To      = rhs,
    std.all = round(std.all, 3),
    est     = round(est, 3),
    pvalue  = round(pvalue, 4)
  )

# Combined
all_paths <- bind_rows(
  path_2nd %>% mutate(Level = "2nd-order"),
  path_1st %>% mutate(Level = "1st-order")
) %>%
  select(Level, From, To, std.all, est, pvalue)

cat("\n--- Full Path Coefficient Table ---\n")
print(all_paths)

# ============================================================
# 11) Export Individual Scores
# ============================================================
lv_scores <- as.data.frame(lavPredict(fit_2nd, type = "lv"))
out_df    <- bind_cols(dat, lv_scores)

score_path <- file.path(OUT_DIR, "UKB_BrainVital8_2ndOrder_scores.csv")
write.csv(out_df, score_path, row.names = FALSE, fileEncoding = "UTF-8")
cat("\nIndividual scores exported:", score_path, "\n")

# ============================================================
# 12) Cross-Cohort Alignment File
# ============================================================
ukb_align <- bind_cols(
  dat %>% select(eid),
  lv_scores
) %>%
  mutate(
    Sleep_manifest = dat$Sleep_z,
    Smoke_manifest = dat$Smoke_rev
  )

align_path <- file.path(OUT_DIR, "UKB_BrainVital8_2ndOrder_alignment.csv")
write.csv(ukb_align, align_path, row.names = FALSE, fileEncoding = "UTF-8")
cat("Cross-cohort alignment file exported\n")

# ============================================================
# 13) Export All Result Tables
# ============================================================
write.csv(order2_loadings,
          file.path(BASE_DIR, "BrainVital8_2ndOrder_loadings.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Second-order loadings table exported\n")

write.csv(order1_full,
          file.path(BASE_DIR, "BrainVital8_1stOrder_loadings_R2.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("First-order loadings + R-squared table exported\n")

write.csv(resid_tbl,
          file.path(BASE_DIR, "BrainVital8_residual_variances.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Residual variances table exported\n")

write.csv(all_paths,
          file.path(BASE_DIR, "BrainVital8_all_path_coefficients.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Path coefficient table exported\n")

fit_idx <- as.data.frame(t(fitMeasures(
  fit_2nd,
  c("cfi", "tli", "rmsea", "rmsea.ci.lower",
    "rmsea.ci.upper", "srmr", "chisq", "df", "pvalue")
)))
write.csv(fit_idx,
          file.path(BASE_DIR, "BrainVital8_model_fit_indices.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Model fit indices exported\n")

# ============================================================
# 14) Export File Summary
# ============================================================
cat("\n", strrep("=", 65), "\n")
cat("All exported files:\n\n")
cat("  [Individual scores]     ", score_path, "\n")
cat("  [Alignment file]        ", align_path, "\n")
cat("  [2nd-order loadings]    ", file.path(BASE_DIR, "BrainVital8_2ndOrder_loadings.csv"), "\n")
cat("  [1st-order loadings+R2] ", file.path(BASE_DIR, "BrainVital8_1stOrder_loadings_R2.csv"), "\n")
cat("  [Residual variances]    ", file.path(BASE_DIR, "BrainVital8_residual_variances.csv"), "\n")
cat("  [Path coefficients]     ", file.path(BASE_DIR, "BrainVital8_all_path_coefficients.csv"), "\n")
cat("  [Model fit indices]     ", file.path(BASE_DIR, "BrainVital8_model_fit_indices.csv"), "\n")
cat(strrep("=", 65), "\n")
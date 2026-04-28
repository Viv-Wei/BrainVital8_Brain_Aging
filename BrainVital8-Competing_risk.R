library(tidyverse)
library(survival)
library(tidycmprsk)
options(scipen = 999)

# ========================= 1. Load Data =========================
df <- read.csv("./data/SEM_Analysis/All_cohort_final_BrainVital8_mortality.csv",
               header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(cohort = factor(cohort, levels = c("ukb", "share", "charls", "elsa", "KLOSA", "HRS")))

cat("Sample size:", nrow(df), "\n")

# ========================= 2. Missing Value Imputation =========================
df <- df %>%
  group_by(cohort) %>%
  mutate(
    across(c(age, sex, bmi, drink, hypertension, BrainVital8),
           ~ifelse(is.na(.), median(., na.rm = TRUE), .))
  ) %>%
  ungroup()

df <- df %>% drop_na(time, status, death_time, death_status,
                     BrainVital8, age, sex, bmi, drink, hypertension)

cat("Final sample size:", nrow(df), "\n")

# ========================= 3. Competing Risk Variables =========================
df <- df %>%
  mutate(
    ftime = pmin(time, death_time, na.rm = TRUE),
    fstatus = case_when(
      status == 1 & time <= death_time ~ "dementia",
      death_status == 1                ~ "death",
      TRUE                             ~ "censored"
    ) %>% factor(levels = c("censored", "dementia", "death"))
  )

# ========================= 4. Within-Cohort Quartiles =========================
df <- df %>%
  group_by(cohort) %>%
  mutate(BrainVital8_Q = ntile(BrainVital8, 4)) %>%
  ungroup() %>%
  mutate(BrainVital8_Q = factor(BrainVital8_Q, 1:4,
                                labels = c("Q1", "Q2", "Q3", "Q4")) %>%
           relevel("Q1"))

# ========================= 5. Fine-Gray Models Per Cohort =========================
cohort_list <- levels(df$cohort)

res_cont_list  <- list()
res_quart_list <- list()

for (coh in cohort_list) {
  
  sub <- df %>% filter(cohort == coh)
  cat("\n==========", coh, "(n =", nrow(sub), ") ==========\n")
  cat("Event distribution:\n")
  print(table(sub$fstatus))
  
  # --- Continuous Fine-Gray ---
  tryCatch({
    fit_c <- tidycmprsk::crr(
      Surv(ftime, fstatus) ~ BrainVital8 + age + sex + bmi + drink + hypertension,
      data = sub
    )
    tmp_c <- tidy(fit_c, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "BrainVital8") %>%
      select(term, HR = estimate, lower = conf.low, upper = conf.high, p = p.value) %>%
      mutate(cohort = coh, across(where(is.numeric), ~round(., 3)))
    res_cont_list[[coh]] <- tmp_c
    print(tmp_c)
  }, error = function(e) {
    cat("  Continuous model failed:", conditionMessage(e), "\n")
  })
  
  # --- Quartile Fine-Gray ---
  tryCatch({
    fit_q <- tidycmprsk::crr(
      Surv(ftime, fstatus) ~ BrainVital8_Q + age + sex + bmi + drink + hypertension,
      data = sub
    )
    tmp_q <- tidy(fit_q, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(str_detect(term, "Q[2-4]")) %>%
      mutate(term = c("Q2 vs Q1", "Q3 vs Q1", "Q4 vs Q1"),
             cohort = coh) %>%
      select(cohort, term, HR = estimate, lower = conf.low, upper = conf.high, p = p.value) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    res_quart_list[[coh]] <- tmp_q
    print(tmp_q)
  }, error = function(e) {
    cat("  Quartile model failed:", conditionMessage(e), "\n")
  })
}

# ========================= 6. Export Results =========================
result_cont_all  <- bind_rows(res_cont_list)
result_quart_all <- bind_rows(res_quart_list)

write.csv(result_cont_all,  "./results/FineGray_continuous_by_cohort.csv",
          row.names = FALSE, fileEncoding = "utf-8")
write.csv(result_quart_all, "./results/FineGray_quartile_by_cohort.csv",
          row.names = FALSE, fileEncoding = "utf-8")

cat("\n\n===== Continuous Variable Summary =====\n")
print(result_cont_all)
cat("\n===== Quartile Summary =====\n")
print(result_quart_all)

# ========================= 7. Descriptive Statistics =========================
desc_table <- df %>%
  group_by(cohort, BrainVital8_Q) %>%
  summarise(
    N = n(),
    case = sum(fstatus == "dementia"),
    death = sum(fstatus == "death"),
    .groups = "drop"
  )

desc_total <- df %>%
  group_by(cohort) %>%
  summarise(
    BrainVital8_Q = "Total",
    N = n(),
    case = sum(fstatus == "dementia"),
    death = sum(fstatus == "death"),
    .groups = "drop"
  )

desc_all <- bind_rows(
  desc_table %>% mutate(BrainVital8_Q = as.character(BrainVital8_Q)),
  desc_total
) %>%
  arrange(cohort, BrainVital8_Q)

cat("\n===== Quartile Descriptive Statistics by Cohort =====\n")
print(desc_all, n = Inf)

write.csv(desc_all, "./results/Quartile_descriptive_by_cohort.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nDone.\n")
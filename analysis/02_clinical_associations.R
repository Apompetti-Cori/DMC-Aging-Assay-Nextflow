suppressPackageStartupMessages({
  library(duckdb)
  library(jsonlite)
  library(tidyverse)
  library(emmeans)
  library(dbplyr)
  library(DBI)
  library(here)
})

conn <- dbConnect(
  duckdb::duckdb(),
  dbdir = here("db/0.1.1/db.duckdb"),
  read_only = TRUE
)
on.exit(dbDisconnect(conn, shutdown = TRUE))

# 1. Fetch and clean data
binary_vars <- c(
  "stroke.tia",
  "chf",
  "alcohol",
  "peripheral_vascular_disease",
  "smoking",
  "renal_disease",
  "copd_or_asthma",
  "afib",
  "autoimmune",
  "thyroid",
  "cancer",
  "anticoagulant_or_antiplatelet",
  "diabetes",
  "hld",
  "cad",
  "hbp",
  "connective_tissue_disease",
  "dementia"
)

data <- left_join(
  tbl(conn, "pred_data"),
  tbl(conn, "metadata") %>% distinct(sample, meta, race, sex),
  by = "sample"
) %>%
  collect() %>%
  filter(cohort == "CUH") %>%
  mutate(meta = map(meta, ~ as.data.frame(fromJSON(.x)))) %>%
  unnest(meta) %>%
  mutate(
    sex = factor(sex, levels = c("F", "M"), labels = c("Female", "Male")),
    race = factor(
      race,
      levels = c("Caucasian", "African American or Black", "Asian", "Other")
    ),
    smoking_history = factor(smoking_history, levels = 1:3),
    weight_cat = factor(
      if_else(bmi >= 30, "Obese", "Non-obese", missing = "Unknown"),
      levels = c("Non-obese", "Obese", "Unknown")
    ),
    across(all_of(binary_vars), ~ factor(.x, levels = 0:1))
  ) %>%
  select(
    sample,
    age,
    pred_age.jsd,
    age_error.jsd,
    Sex = sex,
    Race = race,
    `Smoking History` = smoking_history,
    Weight = weight_cat,
    `Cerebrovascular Disease` = stroke.tia,
    `Congestive Heart Failure` = chf,
    `Alcohol abuse` = alcohol,
    `Peripheral Vascular Disease` = peripheral_vascular_disease,
    Smoking = smoking,
    `Renal Disease` = renal_disease,
    `COPD/Asthma` = copd_or_asthma,
    `Atrial fibrillation` = afib,
    `Autoimmune Disease` = autoimmune,
    `Thryroid Disease` = thyroid,
    `Cancer History` = cancer,
    `Anticoagulant/Antiplatelet Usage` = anticoagulant_or_antiplatelet,
    `Diabetes` = diabetes,
    Hyperlipidemia = hld,
    `Coronary Artery Disease` = cad,
    Hypertension = hbp,
    `Connective Tissue Disease` = connective_tissue_disease,
    Dementia = dementia
  )

# 2. Variable filter map
var_names <- names(data)[5:ncol(data)]
filter_rules <- if_else(
  var_names %in% c("Sex", "Race", "Smoking History", "Weight"),
  "all",
  "over50"
)

# 3. Model execution
res <- map2_dfr(var_names, filter_rules, function(var, flt) {
  df_sub <- data %>%
    {
      if (flt == "over50") filter(., age > 50) else .
    } %>%
    filter(!is.na(.data[[var]]))

  counts <- count(df_sub, group = as.character(.data[[var]]))

  fit <- lm(
    reformulate(paste0("`", var, "`"), response = "age_error.jsd"),
    data = df_sub
  )
  rg <- emmeans(fit, specs = var)

  df_res <- as.data.frame(rg) %>%
    rename(group = !!sym(var)) %>%
    mutate(variable_name = var) %>%
    left_join(counts, by = "group")

  if (nrow(df_res) > 2 && var != "Weight") {
    contrast_df <- joint_tests(fit) %>%
      filter(`model term` == var) %>%
      transmute(
        variable_name = `model term`,
        group = counts$group[1],
        estimate = NA_real_,
        p.value
      )
  } else {
    contrast_df <- contrast(rg, method = "trt.vs.ctrl", ref = 1) %>%
      as.data.frame() %>%
      mutate(
        group = gsub(paste0(var, "| - .*|\\(|\\)"), "", contrast),
        variable_name = var
      ) %>%
      select(
        contrast,
        estimate,
        contrast_SE = SE,
        p.value,
        variable_name,
        group
      )
  }

  df_res %>%
    left_join(contrast_df, by = c("group", "variable_name")) %>%
    mutate(estimate = if_else(row_number() == 1, emmean, estimate))
})

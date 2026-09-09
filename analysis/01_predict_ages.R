suppressPackageStartupMessages({
  library(tidyverse)
  library(duckdb)
  library(here)
  library(dbplyr)
  library(DBI)
})

# 1. Setup target directory and duplicate 0.1.0 DB to 0.1.1
src_db_path <- here("db/0.1.0/db.duckdb")
target_dir <- here("db/0.1.1")
dest_db_path <- here("db/0.1.1/db.duckdb")

if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

# Duplicate the existing 0.1.0 DB so 0.1.1 retains all original tables
file.copy(from = src_db_path, to = dest_db_path, overwrite = TRUE)

# Connect to the destination database in read/write mode
conn <- dbConnect(
  duckdb::duckdb(),
  dbdir = dest_db_path,
  read_only = FALSE
)

# 2. Reference tables in DuckDB
collapsed_tbl <- tbl(conn, "collapsed_data") %>% complete(sample, target)
metadata_tbl <- tbl(conn, "metadata")

target_filter <- c(
  "R5434",
  "R8436",
  "B6_151",
  "R23",
  "Ks07",
  "C13_194",
  "T1_200"
)

# Build and filter base dataset lazily in DuckDB
preds_db <- left_join(
  collapsed_tbl %>%
    select(sample, target, jsd.mean, meth.mean, entropy.relative.mean),
  metadata_tbl %>% distinct(sample, age),
  by = "sample"
) %>%
  mutate(
    cohort = case_when(
      starts_with(sample, "ND") ~ "NINDS",
      starts_with(sample, "CUH") ~ "CUH",
      TRUE ~ "Other"
    )
  ) %>%
  filter(target %in% target_filter) %>%
  group_by(sample) %>%
  filter(n() > (!!length(target_filter) - 2)) %>%
  ungroup()

# 3. Compute GLMs in R memory (fit on remote data)
fit_data <- preds_db %>%
  filter(sample != "ND37325") %>%
  collect()

fit_jsd <- glm(
  jsd.mean ~ age + target + cohort,
  family = "gaussian",
  data = fit_data
)
fit_meth <- glm(
  meth.mean ~ age + target + cohort,
  family = "gaussian",
  data = fit_data
)
fit_entropy <- glm(
  entropy.relative.mean ~ age + target + cohort,
  family = "gaussian",
  data = fit_data
)

# Extract predictions to impute missing values lazily via temporary SQL table
impute_df <- fit_data %>%
  mutate(
    pred_jsd = predict(fit_jsd, newdata = .),
    pred_meth = predict(fit_meth, newdata = .),
    pred_entropy = predict(fit_entropy, newdata = .)
  ) %>%
  select(sample, target, pred_jsd, pred_meth, pred_entropy)

dbWriteTable(
  conn,
  "temp_imputed_preds",
  impute_df,
  temporary = TRUE,
  overwrite = TRUE
)
impute_tbl <- tbl(conn, "temp_imputed_preds")

# Join imputed estimates and aggregate at sample level inside DuckDB
preds_imputed <- preds_db %>%
  left_join(impute_tbl, by = c("sample", "target")) %>%
  mutate(
    jsd.mean = coalesce(jsd.mean, pred_jsd),
    meth.mean = coalesce(meth.mean, pred_meth),
    entropy.relative.mean = coalesce(entropy.relative.mean, pred_entropy)
  )

preds_summary_db <- preds_imputed %>%
  group_by(sample, age, cohort) %>%
  summarise(
    jsd.mean = mean(jsd.mean, na.rm = TRUE),
    meth.mean = mean(meth.mean, na.rm = TRUE),
    entropy.relative.mean = mean(entropy.relative.mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sample)

# 4. Extract sample-level models for age calculations
summary_df <- preds_summary_db %>% collect()

fit_jsd_ninds <- lm(
  jsd.mean ~ age,
  data = filter(summary_df, startsWith(sample, "ND") & sample != "ND37325")
)
fit_jsd_cuh <- lm(
  jsd.mean ~ age,
  data = filter(summary_df, startsWith(sample, "CUH"))
)
fit_meth_ninds <- lm(
  meth.mean ~ age,
  data = filter(summary_df, startsWith(sample, "ND") & sample != "ND37325")
)
fit_meth_cuh <- lm(
  meth.mean ~ age,
  data = filter(summary_df, startsWith(sample, "CUH"))
)
fit_entropy_ninds <- lm(
  entropy.relative.mean ~ age,
  data = filter(summary_df, startsWith(sample, "ND") & sample != "ND37325")
)
fit_entropy_cuh <- lm(
  entropy.relative.mean ~ age,
  data = filter(summary_df, startsWith(sample, "CUH"))
)

# Extract model coefficients as scalar constants
b0_jsd_nd <- coef(fit_jsd_ninds)["(Intercept)"] |> unname()
b1_jsd_nd <- coef(fit_jsd_ninds)["age"] |> unname()
b0_jsd_cu <- coef(fit_jsd_cuh)["(Intercept)"] |> unname()
b1_jsd_cu <- coef(fit_jsd_cuh)["age"] |> unname()

b0_meth_nd <- coef(fit_meth_ninds)["(Intercept)"] |> unname()
b1_meth_nd <- coef(fit_meth_ninds)["age"] |> unname()
b0_meth_cu <- coef(fit_meth_cuh)["(Intercept)"] |> unname()
b1_meth_cu <- coef(fit_meth_cuh)["age"] |> unname()

b0_ent_nd <- coef(fit_entropy_ninds)["(Intercept)"] |> unname()
b1_ent_nd <- coef(fit_entropy_ninds)["age"] |> unname()
b0_ent_cu <- coef(fit_entropy_cuh)["(Intercept)"] |> unname()
b1_ent_cu <- coef(fit_entropy_cuh)["age"] |> unname()

# 5. Compute final predicted ages inside DuckDB and write directly to pred_data
final_preds_db <- preds_summary_db %>%
  collect() %>%
  mutate(
    pred_age.jsd = case_when(
      startsWith(sample, "ND") ~ (jsd.mean - !!b0_jsd_nd) / !!b1_jsd_nd,
      startsWith(sample, "CUH") ~ (jsd.mean - !!b0_jsd_cu) / !!b1_jsd_cu
    ),
    pred_age.meth = case_when(
      startsWith(sample, "ND") ~ (meth.mean - !!b0_meth_nd) / !!b1_meth_nd,
      startsWith(sample, "CUH") ~ (meth.mean - !!b0_meth_cu) / !!b1_meth_cu
    ),
    pred_age.entropy = case_when(
      startsWith(sample, "ND") ~ (entropy.relative.mean - !!b0_ent_nd) /
        !!b1_ent_nd,
      startsWith(sample, "CUH") ~ (entropy.relative.mean - !!b0_ent_cu) /
        !!b1_ent_cu
    ),
    age_error.jsd = pred_age.jsd - age,
    age_error.meth = pred_age.meth - age,
    age_error.entropy = pred_age.entropy - age
  )

# Materialize query results straight into target database table
dbWriteTable(
  conn,
  name = "pred_data",
  value = final_preds_db,
  temporary = FALSE,
  overwrite = TRUE
)

# 6. Index and close
dbExecute(conn, "CREATE INDEX idx_pred_data ON pred_data (sample)")
dbDisconnect(conn, shutdown = TRUE)

library(fastverse)
library(dplyr)
library(tidyr)
library(tibble)
library(duckdb)

db <- dbConnect(
  duckdb::duckdb(),
  dbdir = "/home/apompetti/data/2025-09-10_Aging_Assay/data/db/aging_assay/0.1.1/db.duckdb",
  read_only = TRUE
)

dbListTables(db)

# Create table to be used for predicting ages from jsd values
preds <- left_join(
  tbl(db, "collapsed_data") %>% complete(sample, target),
  tbl(db, "metadata") %>% distinct(sample, age),
  by = "sample"
) %>%
  collect() %>%
  filter(
    target %in%
      c(
        "R5434",
        "R8436",
        "B6_151",
        "R23",
        "Ks07",
        "C13_194",
        "T1_200"
      )
  ) %>%
  mutate(
    cohort = case_when(
      startsWith(sample, "ND") ~ "NINDS",
      startsWith(sample, "CUH") ~ "CUH",
      TRUE ~ "Other"
    )
  ) %>%
  select(
    sample,
    target,
    age,
    cohort,
    jsd.mean,
    meth.mean,
    entropy.relative.mean
  )

# Filter samples with more than one target missing
preds <- preds %>%
  group_by(
    sample
  ) %>%
  filter(
    n() > (n_distinct(preds$target) - 2)
  ) %>%
  ungroup()

preds_wide <- preds %>%
  select(sample, age, cohort, target, jsd.mean) %>%
  pivot_wider(
    id_cols = c(sample, age, cohort),
    names_from = target,
    values_from = c(jsd.mean)
  )

# Impute missing values with mean jsd at that target's age and cohort
fit.jsd <- glm(
  as.formula("jsd.mean ~ age + target + cohort"),
  family = "gaussian",
  data = preds %>% filter(sample != "ND37325")
)

fit.meth <- glm(
  as.formula("meth.mean ~ age + target + cohort"),
  family = "gaussian",
  data = preds %>% filter(sample != "ND37325")
)

(sum(is.na(preds$jsd.mean)) / nrow(preds)) * 100

# Fill missing values with imputated values
preds <- preds %>%
  mutate(
    jsd.mean = ifelse(
      is.na(jsd.mean),
      predict(fit.jsd, newdata = .),
      jsd.mean
    ),
    meth.mean = ifelse(
      is.na(meth.mean),
      predict(fit.meth, newdata = .),
      meth.mean
    )
  )

# Average values for each sample
preds <- preds %>%
  group_by(sample, age, cohort) %>%
  summarise(
    jsd.mean = mean(jsd.mean),
    meth.mean = mean(meth.mean)
  )

# Create fit for reverse regression (jsd ~ age and then predict age)
fit.jsd.ninds <- lm(
  jsd.mean ~ age,
  data = preds %>%
    filter(startsWith(sample, "ND")) %>%
    filter(sample != "ND37325")
)
fit.jsd.cuh <- lm(
  jsd.mean ~ age,
  data = preds %>% filter(startsWith(sample, "CUH"))
)
fit.meth.ninds <- lm(
  meth.mean ~ age,
  data = preds %>%
    filter(startsWith(sample, "ND")) %>%
    filter(sample != "ND37325")
)
fit.meth.cuh <- lm(
  meth.mean ~ age,
  data = preds %>% filter(startsWith(sample, "CUH"))
)

preds <- preds %>%
  mutate(
    pred_age.jsd = case_when(
      startsWith(sample, "ND") ~ (jsd.mean -
        fit.jsd.ninds$coefficients[["(Intercept)"]]) /
        fit.jsd.ninds$coefficients[["age"]],
      startsWith(sample, "CUH") ~ (jsd.mean -
        fit.jsd.cuh$coefficients[["(Intercept)"]]) /
        fit.jsd.cuh$coefficients[["age"]]
    ),
    pred_age.meth = case_when(
      startsWith(sample, "ND") ~ (meth.mean -
        fit.meth.ninds$coefficients[["(Intercept)"]]) /
        fit.meth.ninds$coefficients[["age"]],
      startsWith(sample, "CUH") ~ (meth.mean -
        fit.meth.cuh$coefficients[["(Intercept)"]]) /
        fit.meth.cuh$coefficients[["age"]]
    ),
    age_error.jsd = pred_age.jsd - age,
    age_error.meth = pred_age.meth - age
  )

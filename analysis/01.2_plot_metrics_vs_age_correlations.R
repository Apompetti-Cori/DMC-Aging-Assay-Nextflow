suppressPackageStartupMessages({
  library(fastverse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(duckdb)
})

conn <- dbConnect(
  duckdb::duckdb(),
  dbdir = here("db/0.1.1/db.duckdb"),
  read_only = TRUE
)

data <- left_join(
  tbl(conn, "collapsed_data"),
  tbl(conn, "metadata") %>% distinct(sample, age),
  by = "sample"
) %>%
  mutate(
    cohort = case_when(
      starts_with(sample, "ND") ~ "NINDS",
      starts_with(sample, "CUH") ~ "CUH",
      TRUE ~ "Other"
    )
  ) %>%
  filter(cohort == "NINDS") %>%
  collect() %>%
  group_by(
    target
  ) %>%
  mutate(
    meth.mean = ifelse(target == "R3988", 1 - meth.mean, meth.mean),
    target = paste0(
      "bold('",
      target,
      "')",
      " *' (n = ",
      n(),
      ")'"
    ),
    target = as.factor(target)
  )

xlim <- c(18, max(data$age))
ylim <- c(0, 1)

p1 <- data %>%
  ggplot(
    .,
    aes(
      x = age,
      y = jsd.mean
    )
  ) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  ggpubr::stat_cor(
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.05,
    vjust = 1.5,
    size = 5,
    geom = "label",
    output.type = "text",
    fill = "white",
    linewidth = 0
  ) +
  facet_wrap(
    ~target,
    labeller = label_parsed
  ) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 16, color = "white"),
    strip.background = element_rect(colour = "black", fill = "black")
  ) +
  labs(
    x = "Chronological Age",
    y = "Jensen-Shannon Distance"
  )

p2 <- data %>%
  ggplot(
    .,
    aes(
      x = age,
      y = meth.mean
    )
  ) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    labels = scales::unit_format(scale = 100, suffix = "")
  ) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  ggpubr::stat_cor(
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.05,
    vjust = 1.5,
    size = 5,
    geom = "label",
    output.type = "text",
    fill = "white",
    linewidth = 0
  ) +
  facet_wrap(
    ~target,
    labeller = label_parsed
  ) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 16, color = "white"),
    strip.background = element_rect(colour = "black", fill = "black")
  ) +
  labs(
    x = "Chronological Age",
    y = "DNA Methylation (%)"
  )

data <- left_join(
  tbl(conn, "collapsed_data"),
  tbl(conn, "metadata") %>% distinct(sample, age, cohort),
  by = "sample"
) %>%
  filter(cohort == "CUH") %>%
  collect() %>%
  group_by(
    target
  ) %>%
  mutate(
    meth.mean = ifelse(target == "R3988", 1 - meth.mean, meth.mean),
    target = paste0(
      "bold('",
      target,
      "')",
      " *' (n = ",
      n(),
      ")'"
    ),
    target = as.factor(target)
  )

xlim <- c(18, max(data$age))
ylim <- c(0, 1)

p1 <- data %>%
  ggplot(
    .,
    aes(
      x = age,
      y = jsd.mean
    )
  ) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), limit = c(0, 1)) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  ggpubr::stat_cor(
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.05,
    vjust = 1.5,
    size = 5,
    geom = "label",
    output.type = "text",
    fill = "white",
    linewidth = 0
  ) +
  facet_wrap(
    ~target,
    labeller = label_parsed
  ) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 16, color = "white"),
    strip.background = element_rect(colour = "black", fill = "black")
  ) +
  labs(
    x = "Chronological Age",
    y = "Jensen-Shannon Distance"
  )

p2 <- data %>%
  ggplot(
    .,
    aes(
      x = age,
      y = meth.mean
    )
  ) +
  geom_point(size = 1) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1),
    labels = scales::unit_format(scale = 100, suffix = "")
  ) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  ggpubr::stat_cor(
    label.x = -Inf,
    label.y = Inf,
    hjust = -0.05,
    vjust = 1.5,
    size = 5,
    geom = "label",
    output.type = "text",
    fill = "white",
    linewidth = 0
  ) +
  facet_wrap(
    ~target,
    labeller = label_parsed
  ) +
  #geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 16, color = "white"),
    strip.background = element_rect(colour = "black", fill = "black")
  ) +
  labs(
    x = "Chronological Age",
    y = "DNA Methylation (%)"
  )

suppressPackageStartupMessages({
  library(fastverse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(duckdb)
  library(patchwork)

  plotting_function <- function(
    df,
    y,
    x,
    plot_label,
    y_label,
    x_label,
    outlier_gap = FALSE
  ) {
    p <- ggplot(
      df,
      aes(
        x = !!sym(x),
        y = !!sym(y)
      )
    ) +
      geom_point(size = 1.5) +
      annotate(
        "label",
        label = plot_label,
        x = -Inf,
        y = Inf,
        hjust = -0.05,
        vjust = 1.2,
        size = 5,
        fill = "white",
        linewidth = 0
      ) +
      scale_x_continuous(breaks = seq(0, 100, 25)) +
      scale_y_continuous(breaks = seq(0, max(df[[y]]), 25)) +
      coord_cartesian(xlim = c(0, 100)) +
      theme_bw() +
      theme(
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 14, color = "black"),
        strip.text = element_text(size = 16, color = "white"),
        strip.background = element_rect(colour = "black", fill = "black")
      ) +
      labs(
        x = x_label,
        y = y_label
      )
    if (outlier_gap) {
      dimy1 <- c(
        floor((max(p$data[[y]], na.rm = TRUE)) / 10) * 10,
        max(p$data[[y]]) + 5
      )
      dimy2 <- p$data[[y]]
      dimy2 <- dimy2[!(dimy2 %in% max(p$data[[y]]))]

      p1 <- p +
        coord_cartesian(ylim = dimy1, xlim = c(0, 100)) +
        scale_y_continuous(breaks = seq(0, max(dimy1), 10))
      p1@layers$annotate <- NULL

      p2 <- p + coord_cartesian(ylim = c(0, max(dimy2)), xlim = c(0, 100))

      p <- p1 /
        p2 +
        plot_layout(
          axis_title = "collect",
          axes = "collect",
          heights = c(.1, .9)
        )
    }

    return(p)
  }
})

conn <- dbConnect(
  duckdb::duckdb(),
  dbdir = here("db/0.1.1/db.duckdb"),
  read_only = TRUE
)

data <- tbl(conn, "pred_data") %>%
  collect() %>%
  filter(
    cohort == "NINDS"
  )

cor_label <- data %>%
  filter(sample != "ND37325") %$%
  cor.test(age, pred_age.jsd) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  filter(sample != "ND37325") %>%
  pull(age_error.jsd) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p1 <- plotting_function(
  df = data,
  y = "pred_age.jsd",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (JSD-Derived)",
  x_label = "Chronological Age",
  outlier_gap = TRUE
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 1
  )


cor_label <- data %>%
  filter(sample != "ND37325") %$%
  cor.test(age, pred_age.meth) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  filter(sample != "ND37325") %>%
  pull(age_error.meth) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p2 <- plotting_function(
  df = data,
  y = "pred_age.meth",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (Methylation-Derived)",
  x_label = "Chronological Age",
  outlier_gap = TRUE
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 1
  )

cor_label <- data %$%
  cor.test(age, pred_age.jsd) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  #filter(sample != "ND37325") %>%
  pull(age_error.jsd) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p1 <- plotting_function(
  df = data,
  y = "pred_age.jsd",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (JSD-Derived)",
  x_label = "Chronological Age",
  outlier_gap = TRUE
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 1
  )


cor_label <- data %$%
  cor.test(age, pred_age.meth) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  #filter(sample != "ND37325") %>%
  pull(age_error.meth) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p2 <- plotting_function(
  df = data,
  y = "pred_age.meth",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (Methylation-Derived)",
  x_label = "Chronological Age",
  outlier_gap = TRUE
)

data <- tbl(conn, "pred_data") %>%
  collect() %>%
  filter(
    cohort == "CUH"
  )

cor_label <- data %$%
  cor.test(age, pred_age.jsd) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  pull(age_error.jsd) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p1 <- plotting_function(
  df = data,
  y = "pred_age.jsd",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (JSD-Derived)",
  x_label = "Chronological Age",
  outlier_gap = F
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 1
  )

cor_label <- data %$%
  cor.test(age, pred_age.meth) %$%
  {
    p <- ifelse(
      p.value < 2.2e-16,
      "p < 2.2e-16",
      paste0("p = ", round(p.value, 2))
    )
    paste0("R = ", round(estimate, 2), ", ", p)
  }

mae_label <- data %>%
  pull(age_error.meth) %>%
  abs() %>%
  median() %>%
  {
    paste0("Median Absolute Error = ", round(., 2))
  }

plot_label <- paste0(
  cor_label,
  "\n",
  mae_label
)

p2 <- plotting_function(
  df = data,
  y = "pred_age.meth",
  x = "age",
  plot_label = plot_label,
  y_label = "Predicted Age (Methylation-Derived)",
  x_label = "Chronological Age",
  outlier_gap = F
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray40",
    linewidth = 1
  )

dbDisconnect(conn)
rm(conn)

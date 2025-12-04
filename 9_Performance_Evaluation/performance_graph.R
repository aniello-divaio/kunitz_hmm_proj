sessionInfo()

setwd("C:\Users\nello\Desktop\Project_Lab_Redo") # Change to your wd 
getwd()

# Install the required package
install.packages("ggplot2")
install.packages("tidyverse")

# Load required libraries
library(tidyverse)
library(readr)

# Step 1: Read data
df1 <- read_tsv("set_1_performance.tsv")
df2 <- read_tsv("set_2_performance.tsv")

# Step 2: Preprocess both dataframes
df1 <- df1 %>%
  mutate(
    mode = ifelse(fullseq, "fullseq", "domain"),
    threshold_numeric = as.numeric(threshold),
    minus_log10_threshold = -log10(threshold_numeric),
    dataset = "set_1"
  )

df2 <- df2 %>%
  mutate(
    mode = ifelse(fullseq, "fullseq", "domain"),
    threshold_numeric = as.numeric(threshold),
    minus_log10_threshold = -log10(threshold_numeric),
    dataset = "set_2"
  )

# Step 3: Combine datasets
df_all <- bind_rows(df1, df2)

# Step 4: Plot MCC vs -log10(E-value), colored by dataset, faceted by mode
ggplot(df_all, aes(x = minus_log10_threshold, y = mcc, color = dataset)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_color_manual(values = c("set_1" = "purple", "set_2" = "darkgreen")) +
  scale_x_continuous(
    breaks = 1:12,
    labels = paste0(1:12),
    name = "-log10(e-value)"
  ) +
  scale_y_continuous(limits = c(0.92, 1)) +
  ylab("MCC") +
  facet_wrap(~mode, ncol = 1) +
  labs(
    title = "MCC vs E-value Threshold",
    color = "Dataset:"
  ) +
  theme_minimal(base_size = 14)

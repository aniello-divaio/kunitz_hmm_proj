# Load required libraries
library(ggplot2)
library(dplyr)

# Set working directory (update this path!)
setwd("Users/nello/Desktop/Università/Bioinformatica/8_Lab_Bio/Project")

# Load the results
data <- read.csv("hmm_evaluation_results.csv")

# Convert threshold to numeric for plotting on a log scale
data$Threshold_num <- as.numeric(gsub("1e-", "1e-", data$Threshold))

# Plot: All in one, log x-axis
ggplot(data, aes(x = Threshold_num, y = MCC, color = EvalType, linetype = Set, group = interaction(EvalType, Set))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_log10(
    trans = "log10",
    breaks = c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12),
    labels = scales::label_scientific()
  ) +
  scale_y_continuous(limits = c(0.9, 1.0)) +
  scale_color_manual(values = c("fullseq" = "darkgreen", "bestdomain" = "purple")) +
  labs(
    title = "MCC across E-value thresholds",
    x = "E-value threshold",
    y = "Matthews Correlation Coefficient (MCC)",
    color = "Evaluation type",
    linetype = "Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


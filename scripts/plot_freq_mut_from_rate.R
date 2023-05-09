# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)
library(ggrepel)
library(ggpubr)

# Read the CSV file
data <- read.csv("./SIM_RESULTS_STEPHEN/12_determine_experimental_mutrate_from_freq/12_experimental_mutrate_search_freq_mut_from_rate.csv", sep=',', na.strings=c('NaN','Inf'))

# Filter out rows with NA/NaN/Inf values in 'freq_mut_outcome' column
data <- data[!is.na(data$freq_mut_outcome) & is.finite(data$freq_mut_outcome), ]

filtered_data <- data %>% filter(mut_rate >= 0 & mut_rate <= 0.45)

# Calculate mean and confidence intervals
mean_ci <- data %>% group_by(mut_rate) %>% summarize(mean = mean(freq_mut_outcome), ci = 1.96 * sd(freq_mut_outcome) / sqrt(n()))
mean_ci2 <- data.frame(mean_ci)

# Perform logarithmic regression on filtered data
log_model <- lm(freq_mut_outcome ~ log(mut_rate), data = filtered_data)
log_eq <- expression(y == exp(a) * x^b)

# Calculate Nagelkerke's pseudo R-squared
y_hat <- predict(log_model, newdata = filtered_data)
y <- filtered_data$freq_mut_outcome
null_model <- lm(freq_mut_outcome ~ 1, data = filtered_data)
p_r_squared <- 1 - exp((logLik(log_model) - logLik(null_model)) / length(y))

# Format the equation of the line of best fit
eq_label <- paste("y =", format(coef(log_model)[1], digits = 4), "*", "x^", format(coef(log_model)[2], digits = 4))

# Create scatterplot with mean values and logarithmic regression line
plot <- ggplot(data, aes(x = mut_rate, y = freq_mut_outcome)) +
  geom_point(size = 0.5) +
  geom_line(data = mean_ci, aes(x = mut_rate, y = mean), color = "red", linewidth=0.5) +
  geom_smooth(data = filtered_data, method = "lm", formula = y ~ log(x), se = FALSE, color = "blue", linewidth=0.5) +
  geom_text(x = 0.4, y = 0.1, label = eq_label, color = "blue", size = 4.5) +
  labs(x = "Mutation rate", y = "Frequency of mutation in the outcome") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  theme_classic() +
  theme(text = element_text(size = 15, color='black'),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.line=element_line(linewidth=0.5))

# Save the plot as a PNG file
ggsave("SIM_RESULTS_STEPHEN/12_determine_experimental_mutrate_from_freq/freq_mut_mutrate_plot.png", plot, width = 8, height = 6, dpi = 300)

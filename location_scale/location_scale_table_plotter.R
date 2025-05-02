# Load necessary libraries
library(readxl)
library(dplyr)
library(flextable)
library(xtable)

# Load the multiple aggregated data into DataFrames
df1 <- read_excel("location_scale/dp_pf/summary/dp_pf_location_scale_sim_summary.xlsx")
df2 <- read_excel("location_scale/rejection_abc/summary/rejection_abc_location_scale_sim_summary.xlsx")

# Prepare the data for the first dataset
table_data1 <- df1 %>%
  group_by(epsilon, sample_size) %>%
  summarize(
    mu_mean_post = round(mean(post_mean_mu), 3),  # Mean posterior
    mu_sd_post = round(sd(post_mean_mu), 3),       # Standard deviation
    sigma_mean_post = round(mean(post_mean_sigma), 3),  # Mean posterior
    sigma_sd_post = round(sd(post_mean_sigma), 3),       # Standard deviation
    time_per_ess = round(mean(time_per_ess), 5)
  ) %>%
  mutate(
    mu_mean_sd = paste(mu_mean_post, "(", mu_sd_post, ")", sep = " "),
    sigma_mean_sd = paste(sigma_mean_post, "(", sigma_sd_post, ")", sep = " ")
  ) %>%
  dplyr::select(epsilon, sample_size, mu_mean_sd, sigma_mean_sd, time_per_ess) %>%
  arrange(epsilon, sample_size)

# Prepare the data for the second dataset
table_data2 <- df2 %>%
  group_by(epsilon, sample_size) %>%
  summarize(
    mu_mean_post = round(mean(post_mean_mu), 3),  # Mean posterior
    mu_sd_post = round(sd(post_mean_mu), 3),       # Standard deviation
    sigma_mean_post = round(mean(post_mean_sigma), 3),  # Mean posterior
    sigma_sd_post = round(sd(post_mean_sigma), 3),       # Standard deviation
    time_per_ess = round(mean(time_per_ess), 5)
  ) %>%
  mutate(
    mu_mean_sd = paste(mu_mean_post, "(", mu_sd_post, ")", sep = " "),
    sigma_mean_sd = paste(sigma_mean_post, "(", sigma_sd_post, ")", sep = " ")
  ) %>%
  dplyr::select(epsilon, sample_size, mu_mean_sd, sigma_mean_sd, time_per_ess) %>%
  arrange(epsilon, sample_size)

# Combine the data from both datasets by joining on epsilon and sample_size
combined_data <- full_join(
  table_data1 %>% rename(mu_mean_sd_1 = mu_mean_sd, sigma_mean_sd_1 = sigma_mean_sd, time_per_ess_1 = time_per_ess),
  table_data2 %>% rename(mu_mean_sd_2 = mu_mean_sd, sigma_mean_sd_2 = sigma_mean_sd, time_per_ess_2 = time_per_ess),
  by = c("epsilon", "sample_size")
) %>%
  arrange(epsilon, sample_size)

xtab <- xtable(combined_data)

# Customize the xtable and save as a LaTeX file
tex_filename <- "location_scale/tables/location_scale_comparison_table_test.tex"
print(xtab, type = "latex", file = tex_filename,
      include.rownames = FALSE, 
      sanitize.text.function = function(x) x, 
      add.to.row = list(pos = list(-1), command = c(paste0("\\hline\n\\multicolumn{2}{|c|}{} & \\multicolumn{3}{c|}{ DP-PF } & \\multicolumn{3}{c|}{ DP-Reject-ABC } \\\\\n\\hline\n"))))
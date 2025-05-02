# Load necessary libraries
library(readxl)
library(dplyr)
library(xtable)

# Load the multiple aggregated data into DataFrames
df1 <- read_excel("linear_regression/dp_pf/conjugate_prior/summary/dp_pf_conj_rand_sdp_sim_summary.xlsx")
df2 <- read_excel("linear_regression/dp_data_augmentation_mcmc/conjugate_prior/summary/mcmc_conj_sim_summary.xlsx")
# df3 <- read_excel("linear_regression/dp_rejection_abc/summary/rejection_abc_sim_summary.xlsx")

# Prepare the data for the first dataset
table_data1 <- df1 %>%
  group_by(epsilon, sample_size) %>%
  summarize(
    mu_mean_post = round(mean(post_mu_beta_1), 3),  # Mean posterior
    mu_sd_post = round(sd(post_mu_beta_1), 3),       # Standard deviation for error bars
    time_per_ess = mean(seconds_per_ess)
  ) %>%
  mutate(
    mu_mean_sd = paste(mu_mean_post, "(", mu_sd_post, ")", sep = " ")
  ) %>%
  dplyr::select(epsilon, sample_size, , mu_mean_sd, time_per_ess) %>%
  arrange(epsilon, sample_size)

# Prepare the data for the first dataset
table_data2 <- df2 %>%
  group_by(epsilon, sample_size) %>%
  summarize(
    mu_mean_post = round(mean(mean_estimate_1), 3),  # Mean posterior
    mu_sd_post = round(sd(mean_estimate_1), 3),       # Standard deviation for error bars
    time_per_ess = mean(time_per_ess)
  ) %>%
  mutate(
    mu_mean_sd = paste(mu_mean_post, "(", mu_sd_post, ")", sep = " ")
  ) %>%
  dplyr::select(epsilon, sample_size, mu_mean_sd, time_per_ess) %>%
  arrange(epsilon, sample_size)

# Prepare the data for the third dataset
# table_data3 <- df3 %>%
#   group_by(epsilon, sample_size) %>%
#   summarize(
#     mu_mean_post = round(mean(post_beta_1), 3),  # Mean posterior
#     mu_sd_post = round(sd(post_beta_1), 3),       # Standard deviation for error bars
#     time_per_ess = mean(time / ESS)
#   ) %>%
#   mutate(
#     mu_mean_sd = paste(mu_mean_post, "(", mu_sd_post, ")", sep = " ")
#   ) %>%
#   dplyr::select(epsilon, sample_size, mu_mean_sd, time_per_ess) %>%
#   arrange(epsilon, sample_size)

# # Combine the data from both datasets by joining on epsilon and sample_size
combined_data <- full_join(
  table_data1 %>% rename(mu_mean_sd_1 = mu_mean_sd, time_per_ess_1 = time_per_ess),
  table_data2 %>% rename(mu_mean_sd_2 = mu_mean_sd, time_per_ess_2 = time_per_ess),
  by = c("epsilon", "sample_size")
) %>%
#   # full_join(
#   #   table_data3 %>% rename(mu_mean_sd_3 = mu_mean_sd, time_per_ess_3 = time_per_ess),
#   #   by = c("epsilon", "sample_size")
#   # ) %>%
  arrange(epsilon, sample_size)

xtab <- xtable(combined_data)

# Customize the xtable and save as a LaTeX file
tex_filename <- "linear_regression/visuals/rand_sdp_conjugate_comparison_table.tex"
print(xtab, type = "latex", file = tex_filename,
      include.rownames = FALSE,
      sanitize.text.function = function(x) x,
      add.to.row = list(pos = list(-1), command = c(paste0("\\hline\n\\multicolumn{2}{|c|}{} & \\multicolumn{3}{c|}{ DP-PF } & \\multicolumn{3}{c|}{ DP-DA-MCMC } \\\\\n\\hline\n"))))
      # add.to.row = list(pos = list(-1), command = c(paste0("\\hline\n\\multicolumn{2}{|c|}{} & \\multicolumn{3}{c|}{ DP-PF } & \\multicolumn{3}{c|}{ DP-DA-MCMC } & \\multicolumn{3}{c|}{ DP-Reject-ABC } \\\\\n\\hline\n"))))

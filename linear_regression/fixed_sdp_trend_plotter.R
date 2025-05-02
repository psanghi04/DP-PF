# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(patchwork)

max_coverage = 0.95
df1 <- read_excel("linear_regression/dp_pf/conjugate_prior/summary/dp_pf_conj_fixed_sdp_trend_sim_summary.xlsx")
df2 <- read_excel("linear_regression/dp_pf/non_conjugate_prior/summary/dp_pf_nonconj_fixed_sdp_trend_sim_summary.xlsx")

plot_data_1 <- df1 %>%
  group_by(epsilon, particles) %>%
  summarize(
    mean_coverage = mean(beta_coverage_1),
    sd_post = (mean_coverage * (1 - mean_coverage) / 100) ^ 0.5
  )

plot_data_2 <- df2 %>%
  group_by(epsilon, particles) %>%
  summarize(
    mean_coverage = mean(beta_coverage_1),
    sd_post = (mean_coverage * (1 - mean_coverage) / 100) ^ 0.5
  )

plot1 <- ggplot(plot_data_1, aes(x = as.factor(particles), y = mean_coverage, group = epsilon)) +
  geom_hline(yintercept = max_coverage, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(aes(color = as.factor(epsilon), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(color = as.factor(epsilon), linewidth = ifelse(epsilon == 1, 1, 0.5), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5)) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  labs(
    title = "Conjugate Prior",
    x = expression("Particles:" ~ N),
    y = expression("Coverage of" ~ beta[1]),
    color = expression("Epsilon:" ~ epsilon)
  ) +
  theme_light(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ylim(0.8, 1)

plot2 <- ggplot(plot_data_2, aes(x = as.factor(particles), y = mean_coverage, group = epsilon)) +
  geom_hline(yintercept = max_coverage, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(aes(color = as.factor(epsilon), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5), size = 3) +
  geom_line(aes(color = as.factor(epsilon), linewidth = ifelse(epsilon == 1, 1, 0.5), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5)) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  labs(
    title = "Non-Conjugate Prior",
    x = expression("Particles:" ~ N),
    y = NULL,
    color = expression("Epsilon:" ~ epsilon)
  ) +
  theme_light(base_size = 14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  ylim(0.8, 1)

combined_plot <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")

ggsave("linear_regression/visuals/dp_pf_fixed_sdp_coverage_trend", plot = combined_plot, width = 12, height = 6, units = "in", dpi = 300)

# Uncomment the following lines if you want to plot the confidence interval width instead of coverage

# plot_data_1 <- df1 %>%
#   group_by(epsilon, particles) %>%
#   summarize(
#     mean_coverage = round(mean(beta_CI_1), 3)
#   )
#
# plot_data_2 <- df2 %>%
#   group_by(epsilon, particles) %>%
#   summarize(
#     mean_coverage = round(mean(beta_CI_1), 3)
#   )
#
# plot1 <- ggplot(plot_data_1, aes(x = as.factor(particles), y = mean_coverage, group = epsilon)) +
#   geom_point(aes(color = as.factor(epsilon), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5), size = 3) +
#   geom_line(aes(color = as.factor(epsilon), linewidth = ifelse(epsilon == 1, 1, 0.5), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5)) +
#   scale_alpha_identity() +
#   scale_linewidth_identity() +
#   labs(
#     title = "Conjugate Prior",
#     x = expression("Particles:" ~ N),
#     y = expression("Confidence Interval Width for" ~ beta[1]),
#     color = expression("Epsilon:" ~ epsilon)
#   ) +
#   theme_light(base_size = 14) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()) +
#   ylim(0, 0.6)
#
# plot2 <- ggplot(plot_data_2, aes(x = as.factor(particles), y = mean_coverage, group = epsilon)) +
#   geom_point(aes(color = as.factor(epsilon), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5), size = 3) +
#   geom_line(aes(color = as.factor(epsilon), linewidth = ifelse(epsilon == 1, 1, 0.5), alpha = ifelse(epsilon == 1, 1, 0.5)), position = position_dodge(width = 0.5)) +
#   scale_alpha_identity() +
#   scale_linewidth_identity() +
#   labs(
#     title = "Non-Conjugate Prior",
#     x = expression("Particles:" ~ N),
#     y = NULL,
#     color = expression("Epsilon:" ~ epsilon)
#   ) +
#   theme_light(base_size = 14) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank()) +
#   ylim(0, 0.6)
#
# combined_plot <- plot1 + plot2 + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
#
# ggsave("linear_regression/visuals/dp_pf_fixed_sdp_ci_trend", plot = combined_plot, width = 12, height = 6, units = "in", dpi = 300)

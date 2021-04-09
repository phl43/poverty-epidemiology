library(tidyverse)
library(lubridate)
library(zoo)
library(rstan)
library(matrixStats)
library(cowplot)
library(EpiEstim)
library(modelsummary)
library(gt)
library(fixest)

# NOTE: all the results presented in the post are based to data up to April 4th, but the code below will use the latest
# data available, so if you want to replicate the results in the post you have to exclude data after that date

# plot charts that summarize a model fit and compare it to the actual data
plot_model_fit_summary <- function(actual_data, stan_data, stan_output, file_path, file_suffix) {
  predicted_R_t <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "Total",
    R = colMeans(stan_output$Rt)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$Rt, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$Rt, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_R_h <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "Other variants",
    R = colMeans(stan_output$Rt_h)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$Rt_h, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$Rt_h, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_R_v <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "B.1.1.7",
    R = colMeans(stan_output$Rt_v)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$Rt_v, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$Rt_v, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_R_model_only <- predicted_R_t %>%
    bind_rows(predicted_R_h) %>%
    bind_rows(predicted_R_v)
  
  predicted_R_model_actual <- predicted_R_t %>%
    mutate(type = "Model") %>%
    bind_rows(actual_R)
  
  predicted_incidence_h <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "Other variants",
    incidence = colMeans(stan_output$expected_C_h)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$expected_C_h, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$expected_C_h, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_incidence_v <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "B.1.1.7",
    incidence = colMeans(stan_output$expected_C_v)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$expected_C_v, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$expected_C_v, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_incidence_t <- tibble(
    date = seq(ymd("2021-01-01"), max(actual_data$date), by = "1 day"),
    type = "Total",
    incidence = colMeans(stan_output$expected_C_h + stan_output$expected_C_v)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    lower = colQuantiles(stan_output$expected_C_h + stan_output$expected_C_v, probs = 0.025)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)],
    upper = colQuantiles(stan_output$expected_C_h + stan_output$expected_C_v, probs = 0.975)[(stan_data$seed_length + 1):(stan_data$seed_length + stan_data$prediction_length)]
  )
  
  predicted_incidence <- actual_data %>%
    mutate(type = "Actual") %>%
    rename(incidence = smoothed_cases) %>%
    bind_rows(predicted_incidence_t %>% mutate(type = "Model"))
  
  ggplot(predicted_R_model_only, aes(x = date, y = R, group = type, color = type)) +
    geom_line(size = 0.75) +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.1, show.legend = FALSE) +
    theme_bw() +
    ggtitle("Effective reproduction number in France predicted by the model (B.1.1.7, other variants and total)") +
    xlab("Date") +
    ylab("R") +
    scale_x_date(
      labels = scales::date_format("%m/%d"),
      date_breaks = "7 day"
    ) +
    ylim(0, max(predicted_R_model_only$upper)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
    ggsave(
      paste0(
        file_path,
        "Effective reproduction number in France predicted by the model (B.1.1.7, other variants and total) - ",
        file_suffix,
        ".png"
        ),
      width = 12,
      height = 6
      )
  
  R_plot_model_actual <- ggplot(predicted_R_model_actual, aes(x = date, y = R, group = type, color = type)) +
    geom_line(size = 0.75) +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.1, show.legend = FALSE) +
    theme_bw() +
    ggtitle("Effective reproduction number in France (actual and predicted by the model)") +
    xlab("Date") +
    ylab("R") +
    scale_x_date(
      labels = scales::date_format("%m/%d"),
      date_breaks = "7 day"
    ) +
    ylim(0, max(predicted_R_model_actual$upper)) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")
  
  incidence_plot <- ggplot(predicted_incidence, aes(x = date, y = incidence, group = type, color = type)) +
    geom_line(size = 0.75) +
    geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.1, show.legend = FALSE) +
    theme_bw() +
    ggtitle("Daily number of COVID-19 cases in France (actual and predicted by the model)") +
    xlab("Date") +
    ylab("Daily number of COVID-19 cases") +
    scale_x_date(
      labels = scales::date_format("%m/%d"),
      date_breaks = "7 day"
    ) +
    scale_y_continuous(
      labels = scales::comma
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank()
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)")
  
  plot_grid(
    R_plot_model_actual,
    incidence_plot,
    labels = c("", ""),
    ncol = 1
  ) +
    ggsave(
      paste0(
        file_path,
        "Daily number of COVID-19 cases by variant in France (actual and predicted by the model) - ",
        file_suffix,
        ".png"
        ),
      width = 12,
      height = 12
      )
}

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#########################################################################################################
#                                       DATA PREPARATION                                                #
#########################################################################################################

# I only include data from metropolitan France and exclude overseas territories
dom <- c(
  "971",
  "972",
  "973",
  "974",
  "975",
  "976",
  "977",
  "978"
)

url_department_test_data <- "https://www.data.gouv.fr/fr/datasets/r/406c6a23-e283-4300-9484-54e78c8ae675"

national_incidence_data <- read_delim(url_department_test_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% dom)) %>%
  mutate(
    date = ymd(jour)
  ) %>%
  group_by(date) %>%
  summarize(
    P = sum(P)
  ) %>%
  mutate(
    smoothed_cases = rollmean(P, 7, fill = c(0, 0, 0), align = "right")
  ) %>%
  select(
    date,
    smoothed_cases
  )

# see https://www.medrxiv.org/content/10.1101/2020.11.17.20231548v2 for the serial interval's mean and standard deviation
actual_R <- estimate_R(
  national_incidence_data$smoothed_cases,
  method = "parametric_si",
  config = make_config(list(mean_si = 5.6, std_si = 4.2))
)$R %>%
  filter (t_end >= which(national_incidence_data$date == ymd("2021-01-01"))) %>%
  add_column(
    date = seq(ymd("2021-01-01"), last(national_incidence_data$date), by = "1 day"),
    type = "Actual"
    ) %>%
  rename(R = `Mean(R)`) %>%
  select(
    date,
    R,
    type
  )

# plot the effective reproduction number in France since the beginning of the year
ggplot(actual_R, aes(x = date, y = R)) +
  geom_line(size = 0.75, color = "steelblue") +
  theme_bw() +
  ggtitle("Effective reproduction number in France") +
  xlab("Date") +
  ylab("R") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  ylim(0, 1.5) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/Effective reproduction number in France.png", width = 12, height = 6)

actual_data <- national_incidence_data %>%
  inner_join(actual_R, by = "date") %>%
  filter(date >= ymd("2021-01-01"))

# those are the only data points used in Gaymard et al.'s paper
gaymard_variant_data <- tibble(
  date = c(ymd("2021-01-08"), ymd("2021-01-27"), ymd("2021-02-12"), ymd("2021-02-18")),
  prevalence = c(0.033, 0.13, 0.37, 0.49),
  type = c("Training", "Training", "Test", "Test")
)

url_voc_data <- "https://www.data.gouv.fr/fr/datasets/r/16f4fd03-797f-4616-bca9-78ff212d06e8"

# create a data frame with all the publicly available data on variants of concern in France
national_voc_data <- read_delim(url_voc_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% dom)) %>%
  mutate(
    # the prevalence of B.1.1.7 in the dataset is the average for d to d + 7, so
    # I assign the estimate for a given week to d + 4
    date = ymd(str_extract(semaine, "^[:digit:]{4}-[:digit:]{2}-[:digit:]{2}")) + 3,
  ) %>%
  group_by(date) %>%
  summarize(
    prevalence_b117 = sum(Nb_susp_501Y_V1) / sum(Nb_tests_PCR_TA_crible)
  )

incidence_by_variant <- national_voc_data %>%
  inner_join(national_incidence_data, by = "date") %>%
  mutate(
    incidence_h = round((1 - prevalence_b117) * smoothed_cases),
    incidence_v = round(prevalence_b117 * smoothed_cases)
  ) %>%
  bind_rows(
    tibble(
      date = c(ymd("2021-01-08"), ymd("2021-01-27")),
      incidence_h = round((1 - gaymard_variant_data$prevalence[gaymard_variant_data$date %in% c(ymd("2021-01-08"), ymd("2021-01-27"))]) * actual_data$smoothed_cases[actual_data$date %in% c(ymd("2021-01-08"), ymd("2021-01-27"))]),
      incidence_v = round(gaymard_variant_data$prevalence[gaymard_variant_data$date %in% c(ymd("2021-01-08"), ymd("2021-01-27"))] * actual_data$smoothed_cases[actual_data$date %in% c(ymd("2021-01-08"), ymd("2021-01-27"))])
    )
  ) %>%
  mutate(
    prevalence_v = incidence_v / (incidence_v + incidence_h),
    type = ifelse(date %in% c(ymd("2021-01-08"), ymd("2021-01-27")), "Training", "Test")
  ) %>%
  arrange(date) %>%
  select(
    date,
    incidence_h,
    incidence_v,
    prevalence_v,
    type
  )

national_incidence_data <- national_incidence_data %>%
  filter(date >= ymd("2021-01-01"))

seed_length <- 10
prediction_length <- as.integer(max(national_incidence_data$date) - ymd("2021-01-01")) + 1

# LIST OF ESTIMATES OF THE GENERATION TIME DISTRIBUTION I FOUND IN THE LITERATURE WITH THE VALUE OF THE PARAMETERS
# mean = 5.20, sd = 1.72 (https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257)
# mean = 3.95, sd = 1.51 (https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257)
# mean = 5.5, sd = 1.8 but 5.5 and 2.1 when Gamma distribution is assumed (https://www.medrxiv.org/content/10.1101/2020.09.04.20188516v2)
# median = 5, sd = 1.9 (https://science.sciencemag.org/content/368/6491/eabb6936)
# mean = 3.99, sd = 2.96 (https://www.sciencedirect.com/science/article/pii/S2468042720300634)
# median = 3.71, sd = 2.79 (https://elifesciences.org/articles/57149, I obtained the standard deviation from the authors through personal communication)
# median = 2.82, sd = 2.49 (https://elifesciences.org/articles/57149, I obtained the standard deviation from the authors through personal communication)

mean_gt_range <- seq(3, 7, by = 0.1)
sd_gt_range <- seq(1.5, 3, by = 0.1)
R_h_range <- c(0.9, 1, 1.1)

dir.create("Figures")
dir.create("Figures/Econometric Analysis")
dir.create("Figures/Epidemiological Modeling")
dir.create("Figures/Epidemiological Modeling/Model 1")
dir.create("Figures/Epidemiological Modeling/Model 2")
dir.create("Figures/Epidemiological Modeling/Model 3")

#########################################################################################################
#                               FIRST MODEL (REPLICATION OF GAYMARD ET AL.)                             #
#########################################################################################################

model1 <- stan_model("Stan/model1.stan")

advantage_estimates1 <- tibble()

for (mu in mean_gt_range) {
  for (sigma in sd_gt_range) {
    for (R in R_h_range) {
      stan_data1 <- list(
        N = 2,
        t = c(8, 27),
        mu = mu,
        sigma = sigma,
        R_h = R,
        N_t = c(round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-08")]), round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-27")])),
        C_v = c(
          round(gaymard_variant_data$prevalence[gaymard_variant_data$date == ymd("2021-01-08")] * actual_data$smoothed_cases[actual_data$date == ymd("2021-01-08")]),
          round(gaymard_variant_data$prevalence[gaymard_variant_data$date == ymd("2021-01-27")] * actual_data$smoothed_cases[actual_data$date == ymd("2021-01-27")])
        ),
        prediction_length = prediction_length
      )

      fit1 <- sampling(
        model1,
        data = stan_data1,
        iter = 4000,
        warmup = 2000,
        chains = 4,
        thin = 1,
        control = list(adapt_delta = 0.95, max_treedepth = 15)
      )

      output1 <- rstan::extract(fit1)

      advantage_estimates1 <- advantage_estimates1 %>%
        bind_rows(
          tibble(
            mean_gt = mu,
            sd_gt = sigma,
            R_non_b117 = R,
            estimate = mean(output1$R_v / R) - 1,
            lower = quantile(output1$R_v / R, probs = 0.025) - 1,
            upper = quantile(output1$R_v / R, probs = 0.975) - 1
          )
        )
    }
  }
}

ggplot(advantage_estimates1, aes(x = estimate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  theme_bw() +
  ggtitle(
    "Distribution of transmissibility advantage estimates with a model similar to that used by Gaymard et al.",
    subtitle = "(R of historical lineage = {0.9, 1, 1.1}, mean of generation time spanning [3,7] by increment of 0.1, standard deviation spanning [1.5,3] by increment of 0.1)"
    ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  ylab("Density") +
  xlab("Estimate") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  ggsave(
    "Figures/Epidemiological Modeling/Model 1/Distribution of transmissibility advantage estimates - Model 1.png",
    width = 12,
    height = 6
  )

# plot charts showing the result with the same value of R_h, mu and sigma as in Gaymard et al.'s central scenario

stan_data1 <- list(
  N = 2,
  t = c(
    as.integer(ymd("2021-01-08") - ymd("2021-01-01") + 1),
    as.integer(ymd("2021-01-27") - ymd("2021-01-01") + 1)
    ),
  mu = 6.5,
  sigma = 6.5 * 0.62,
  R_h = 1,
  N_t = c(round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-08")]), round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-27")])),
  C_v = c(
    incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-08")],
    incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-27")]
  ),
  prediction_length = prediction_length
)

fit1 <- sampling(
  model1,
  data = stan_data1,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

output1 <- rstan::extract(fit1)

b117_transmissibility_advantage_replication_gaymard1 <- mean(output1$alpha)

data_for_plot_prevalence1 <- tibble(
  date = seq(ymd("2021-01-01"), max(national_incidence_data$date), by = "1 day"),
  type = "Model",
  prevalence = colMeans(output1$expected_p_t),
  lower = colQuantiles(output1$expected_p_t, probs = 0.025),
  upper = colQuantiles(output1$expected_p_t, probs = 0.975)
) %>%
  bind_rows(
    incidence_by_variant %>%
      select(date, prevalence_v) %>%
      rename(prevalence = prevalence_v) %>%
      mutate(type = "Actual")
    ) %>%
  filter(date <= max(incidence_by_variant$date))

ggplot(data_for_plot_prevalence1, aes(x = date, y = prevalence, group = type, color = type)) +
  geom_line(data = data_for_plot_prevalence1 %>% filter(type == "Model"), size = 0.75) +
  geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  geom_point(data = data_for_plot_prevalence1 %>% filter(type == "Actual"), size = 3) +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "Prevalence of B.1.1.7 in France showing up to date data - Actual and predicted by a model similar to that used in Gaymard et al.'s paper",
    subtitle = "(mean of generation time = 6.5, standard deviation of generation time = 4, R of non-B.1.1.7 variants = 1)"
  ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("Date") +
  scale_color_discrete(
    name = NULL,
    guide = guide_legend(
      override.aes = list(linetype = c(0, 1), shape = c(16, NA))
      )
    ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave(
    "Figures/Epidemiological Modeling/Model 1/Prevalence of B.1.1.7 in France showing up to date data - Replication of Gaymard et al. (mean of generation time = 6.5, standard deviation of generation time = 4, R of non-B.1.1.7 variants = 1).png",
    width = 12,
    height = 6
  )

predicted_incidence_t1 <- tibble(
  date = seq(ymd("2021-01-01"), max(national_incidence_data$date), by = "1 day"),
  type = "Model",
  incidence = colQuantiles(output1$expected_C_h + output1$expected_C_v, probs = 0.5),
  lower = colQuantiles(output1$expected_C_h + output1$expected_C_v, probs = 0.025),
  upper = colQuantiles(output1$expected_C_h + output1$expected_C_v, probs = 0.975)
)

predicted_incidence1 <- national_incidence_data %>%
  mutate(type = "Actual") %>%
  rename(incidence = smoothed_cases) %>%
  bind_rows(predicted_incidence_t1)

ggplot(predicted_incidence1, aes(x = date, y = incidence, group = type, color = type)) +
  geom_line(size = 0.75) +
  geom_ribbon(aes(ymin = lower, ymax = upper), linetype = 0, alpha = 0.1, show.legend = FALSE) +
  theme_bw() +
  ggtitle(
    "Daily number of COVID-19 cases in France - Actual and predicted by a model similar to that used in Gaymard et al.'s paper",
    subtitle = "(mean of generation time = 6.5, standard deviation of generation time = 4, R of non-B.1.1.7 variants = 1)"
  ) +
  xlab("Date") +
  ylab("Daily number of COVID-19 cases") +
  scale_x_date(
    labels = scales::date_format("%m/%d"),
    date_breaks = "7 day"
  ) +
  scale_y_continuous(
    labels = scales::comma
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_blank()
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave(
    "Figures/Epidemiological Modeling/Model 1/Daily number of COVID-19 cases in France - Replication of Gaymard et al. (mean of generation time = 6.5, standard deviation of generation time = 4, R of non-B.1.1.7 variants = 1).png",
    width = 12,
    height = 6
  )

# estimate B.1.1.7's transmissibility advantage with the same model as Gaymard et al. but using the estimate
# of the generation time distribution from Challen et al.'s paper

stan_data1 <- list(
  N = 2,
  t = c(
    as.integer(ymd("2021-01-08") - ymd("2021-01-01") + 1),
    as.integer(ymd("2021-01-27") - ymd("2021-01-01") + 1)
  ),
  mu = 4.8,
  sigma = 1.7,
  R_h = 1,
  N_t = c(round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-08")]), round(actual_data$smoothed_cases[actual_data$date == ymd("2021-01-27")])),
  C_v = c(
    incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-08")],
    incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-27")]
  ),
  prediction_length = prediction_length
)

fit1 <- sampling(
  model1,
  data = stan_data1,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

output1 <- rstan::extract(fit1)

b117_transmissibility_advantage_replication_gaymard2 <- mean(output1$alpha)

#########################################################################################################
#                                       SECOND MODEL                                                    #
#########################################################################################################

model2 <- stan_model("Stan/model2.stan")

advantage_estimates2 <- tibble()

for (mu in mean_gt_range) {
  for (sigma in sd_gt_range) {
    stan_data2 <- list(
      N = nrow(incidence_by_variant),
      t = as.integer(incidence_by_variant$date - ymd("2021-01-01") + 1),
      C_h = incidence_by_variant$incidence_h,
      C_v = incidence_by_variant$incidence_v,
      mu = mu,
      sigma = sigma,
      prediction_length = prediction_length
    )

    fit2 <- sampling(
      model2,
      data = stan_data2,
      iter = 4000,
      warmup = 2000,
      chains = 4,
      thin = 1,
      control = list(adapt_delta = 0.95, max_treedepth = 15)
    )

    output2 <- rstan::extract(fit2)

    advantage_estimates2 <- advantage_estimates2 %>%
      bind_rows(
        tibble(
          mean_gt = mu,
          sd_gt = sigma,
          estimate = mean(output2$R_v / output2$R_h) - 1,
          lower = quantile(output2$R_v / output2$R_h, probs = 0.025) - 1,
          upper = quantile(output2$R_v / output2$R_h, probs = 0.975) - 1
        )
      )
  }
}

ggplot(advantage_estimates2, aes(x = estimate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  theme_bw() +
  ggtitle(
    "Distribution of transmissibility advantage estimates",
    subtitle = "(mean of generation time spanning [3,7] by increment of 0.1, standard deviation spanning [1.5,3] by increment of 0.1)"
    ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  ylab("Density") +
  xlab("Estimate") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  ggsave(
    "Figures/Epidemiological Modeling/Model 2/Distribution of transmissibility advantage estimates - Model 2.png",
    width = 12,
    height = 6
  )

# plot charts showing the result when using the estimate of the generation time distribution from Challen et al.'s paper

mu <- 4.8
sd <- 1.7

stan_data2 <- list(
  N = nrow(incidence_by_variant),
  t = as.integer(incidence_by_variant$date - ymd("2021-01-01") + 1),
  C_h = incidence_by_variant$incidence_h,
  C_v = incidence_by_variant$incidence_v,
  mu = mu,
  sigma = sd,
  prediction_length = prediction_length
)

fit2 <- sampling(
  model2,
  data = stan_data2,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

output2 <- rstan::extract(fit2)

b117_transmissibility_advantage_replication_gaymard3 <- mean(output2$R_v / output2$R_h) - 1

#########################################################################################################
#                                       THIRD MODEL                                                     #
#########################################################################################################

interventions_start <- c(
  ymd("2021-01-16"),
  ymd("2021-02-06"),
  ymd("2021-02-13"),
  ymd("2021-02-20")
)

interventions_end <- c(
  ymd("2021-05-01"),
  ymd("2021-02-21"),
  ymd("2021-02-28"),
  ymd("2021-03-07")
)

X <- matrix(nrow = seed_length + prediction_length, ncol = 4)

# create the features matrix
for (i in 1:(seed_length + prediction_length)) {
  for (j in 1:4) {
    if (i >= as.integer(interventions_start[j] - ymd("2021-01-01")) + seed_length + 1 & i <= as.integer(interventions_end[j] - ymd("2021-01-01")) + seed_length + 1) {
      X[i, j] = 1
    } else {
      X[i, j] = 0
    }
  }
}

# I use Challen et al.'s estimate of the generation time distribution
stan_data3 <- list(
  N = nrow(incidence_by_variant),
  t = as.integer(incidence_by_variant$date - ymd("2021-01-01") + 1),
  P = 4,
  seed_length = seed_length,
  prediction_length = prediction_length,
  X = X,
  C_h = incidence_by_variant$incidence_h,
  C_v = incidence_by_variant$incidence_v,
  mu = 4.8,
  sigma = 1.7
)

model3 <- stan_model("Stan/model3.stan")

fit3 <- sampling(
  model3,
  data = stan_data3,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

output3 <- rstan::extract(fit3)

b117_transmissibility_advantage_replication_gaymard4 <- mean(output3$R0_v / output3$R0_h) - 1

plot_model_fit_summary(actual_data, stan_data3, output3, "Figures/Epidemiological Modeling/Model 3/", "Model 3")

#########################################################################################################
#                                       ECONOMETRIC ANALYSIS                                            #
#########################################################################################################

# this is the assumption about the mean of the generation time distribution used to estimate R from incidence growth
# rates below, which I varied to see how it affected the results of the analysis, as reported in the post
mean_gt <- 4.8

department_voc_data <- read_delim(url_voc_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% dom)) %>%
  mutate(
    date = ymd(str_extract(semaine, "^[:digit:]{4}-[:digit:]{2}-[:digit:]{2}")),
    dep = factor(dep),
    prevalence_b117 = Prc_susp_501Y_V1 / 100
  )

national_growth_prevalence_data <- read_delim(url_department_test_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% dom)) %>%
  mutate(
    date = ymd(jour)
  ) %>%
  group_by(date) %>%
  summarize(
    P = sum(P)
  ) %>%
  mutate(
    weekly_cases_total = rollsumr(P, 7, fill = rep(0, 6), align = "left")
  ) %>%
  inner_join(national_voc_data %>% mutate(date = date - 3), by = "date") %>%
  mutate(
    weekly_cases_b117 = weekly_cases_total * prevalence_b117,
    weekly_cases_others = weekly_cases_total * (1 - prevalence_b117),
    weekly_growth_factor_b117 = weekly_cases_b117 / lag(weekly_cases_b117, 7),
    weekly_growth_factor_others = weekly_cases_others / lag(weekly_cases_others, 7),
    R_b117 = weekly_growth_factor_b117 ^ (mean_gt / 7),
    R_others = weekly_growth_factor_others ^ (mean_gt / 7),
    advantage = R_b117 / R_others - 1,
    prevalence_b117 = lag(prevalence_b117, 7)
  ) %>%
  filter(
    date %in% seq(ymd("2021-01-04"), ymd("2021-03-29"), by = "7 day")
  ) %>%
  mutate(
    period = paste0("Week ", as.integer(date - 7 - ymd("2021-01-04")) / 7 + 1, " to week ", as.integer(date - ymd("2021-01-04")) / 7 + 1)
  ) %>%
  na.omit() %>%
  bind_rows(
    tibble(
      period = c("01/08 to 01/27", "01/27 to 02/15"),
      prevalence_b117 = c(
        gaymard_variant_data$prevalence[gaymard_variant_data$date == ymd("2021-01-08")],
        gaymard_variant_data$prevalence[gaymard_variant_data$date == ymd("2021-01-27")]
        ),
      R_b117 = c(
        (incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-27")] / incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-08")]) ^ (mean_gt / as.integer(ymd("2021-01-27") - ymd("2021-01-08"))),
        (incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-02-15")] / incidence_by_variant$incidence_v[incidence_by_variant$date == ymd("2021-01-27")]) ^ (mean_gt / as.integer(ymd("2021-02-15") - ymd("2021-01-27")))
      ),
      R_others = c(
        (incidence_by_variant$incidence_h[incidence_by_variant$date == ymd("2021-01-27")] / incidence_by_variant$incidence_h[incidence_by_variant$date == ymd("2021-01-08")]) ^ (mean_gt / as.integer(ymd("2021-01-27") - ymd("2021-01-08"))),
        (incidence_by_variant$incidence_h[incidence_by_variant$date == ymd("2021-02-15")] / incidence_by_variant$incidence_h[incidence_by_variant$date == ymd("2021-01-27")]) ^ (mean_gt / as.integer(ymd("2021-02-15") - ymd("2021-01-27")))
      ),
      advantage = R_b117 / R_others - 1
    )
  ) %>%
  arrange(period) %>%
  select(
    period,
    prevalence_b117,
    weekly_growth_factor_b117,
    weekly_growth_factor_others,
    R_b117,
    R_others,
    advantage
  )

national_growth_prevalence_data$period <- fct_relevel(
  national_growth_prevalence_data$period,
  c(
    "01/08 to 01/27",
    "01/27 to 02/15",
    "Week 7 to week 8",
    "Week 8 to week 9",
    "Week 9 to week 10",
    "Week 10 to week 11",
    "Week 11 to week 12",
    "Week 12 to week 13"
  )
)

ggplot(national_growth_prevalence_data, aes(x = prevalence_b117, y = advantage)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "B.1.1.7's transmissibility advantage vs. prevalence of B.1.1.7 at the national level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("B.1.1.7's transmissibility advantage") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/B.1.1.7's transmissibility advantage vs. prevalence of B.1.1.7 at the national level.png", width = 12, height = 6)

ggplot(national_growth_prevalence_data, aes(x = prevalence_b117, y = R_b117)) +
  geom_point(aes(group = period, color = period)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "B.1.1.7's effective reproduction number vs. prevalence of B.1.1.7 at the national level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("B.1.1.7's effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/B.1.1.7's effective reproduction number vs. prevalence of B.1.1.7 at the national level.png", width = 12, height = 6)

department_growth_prevalence_data <- read_delim(url_department_test_data, delim = ";") %>%
  filter(cl_age90 == "0" & !(dep %in% dom)) %>%
  mutate(
    date = ymd(jour),
    weekly_cases_total = rollsumr(P, 7, fill = rep(0, 6), align = "left")
  ) %>%
  inner_join(department_voc_data, by = c("date", "dep")) %>%
  group_by(dep) %>%
  mutate(
    weekly_cases_b117 = weekly_cases_total * prevalence_b117,
    weekly_cases_others = weekly_cases_total * (1 - prevalence_b117),
    weekly_growth_factor_total = weekly_cases_total / lag(weekly_cases_total, 7),
    weekly_growth_factor_b117 = weekly_cases_b117 / lag(weekly_cases_b117, 7),
    weekly_growth_factor_others = weekly_cases_others / lag(weekly_cases_others, 7),
    R_total = weekly_growth_factor_total ^ (mean_gt / 7),
    R_b117 = weekly_growth_factor_b117 ^ (mean_gt / 7),
    R_others = weekly_growth_factor_others ^ (mean_gt / 7),
    advantage = R_b117 / R_others - 1,
    prevalence_b117 = lag(prevalence_b117, 7)
  ) %>%
  ungroup() %>%
  filter(
    date %in% seq(ymd("2021-01-04"), ymd("2021-03-29"), by = "7 day")
  ) %>%
  mutate(
    week = paste0("Week ", as.integer(date - 7 - ymd("2021-01-04")) / 7 + 1, " to week ", as.integer(date - ymd("2021-01-04")) / 7 + 1)
  ) %>%
  select(
    week,
    dep,
    prevalence_b117,
    weekly_growth_factor_total,
    weekly_growth_factor_b117,
    weekly_growth_factor_others,
    R_total,
    R_b117,
    R_others,
    advantage
  ) %>%
  na.omit()

department_growth_prevalence_data$week <- fct_relevel(
  department_growth_prevalence_data$week,
  c(
    "Week 7 to week 8",
    "Week 8 to week 9",
    "Week 9 to week 10",
    "Week 10 to week 11",
    "Week 11 to week 12",
    "Week 12 to week 13"
  )
)

econometric_model1 <- femlm(
  log(R_total) ~ prevalence_b117,
  data = department_growth_prevalence_data,
  cluster = c("dep", "week"),
  family = "gaussian"
)

econometric_model2 <- femlm(
  log(R_total) ~ prevalence_b117 | week,
  data = department_growth_prevalence_data,
  cluster = c("dep", "week"),
  family = "gaussian"
)

econometric_model3 <- femlm(
  log(R_total) ~ prevalence_b117 | dep,
  data = department_growth_prevalence_data,
  cluster = c("dep", "week"),
  family = "gaussian"
)

econometric_model4 <- femlm(
  log(R_total) ~ prevalence_b117 | week + dep,
  data = department_growth_prevalence_data,
  cluster = c("dep", "week"),
  family = "gaussian"
)

econometric_models <- list()
econometric_models[["Model with no fixed effect"]] <- econometric_model1
econometric_models[["Model with period fixed effect"]] <- econometric_model2
econometric_models[["Model with department fixed effect"]] <- econometric_model3
econometric_models[["Model with both period and department fixed effects"]] <- econometric_model4

cm <- c(
  "prevalence_b117" = "% of B.1.1.7"
)

msummary(
  econometric_models,
  output = "gt",
  title = "Econometric analysis of the relationship between the effective reproduction number and the prevalence of B.1.1.7 in France",
  coef_map = cm,
  gof_omit = "FE|R2|Std",
  statistic = "std.error",
  stars = TRUE,
  exponentiate = TRUE
  ) %>%
  tab_footnote(
    footnote = "Standard errors are clustered by department and period.",
    locations = cells_body(rows = 2, columns = 2:5)
  ) %>%
  gtsave("Figures/Econometric Analysis/Econometric analysis of the relationship between the effective reproduction number and the prevalence of B.1.1.7 in France.png")

coefficient_plot_data <- tribble(
  ~model, ~estimate, ~lower, ~upper,
  "Model with no fixed effect", exp(econometric_model1$coefficients[2]) - 1, exp(confint(econometric_model1)[2,1]) - 1, exp(confint(econometric_model1)[2,2]) - 1,
  "Model with period fixed effect", exp(econometric_model2$coefficients[1]) - 1, exp(confint(econometric_model2)[1,1]) - 1, exp(confint(econometric_model2)[1,2]) - 1,
  "Model with department fixed effect", exp(econometric_model3$coefficients[1]) - 1, exp(confint(econometric_model3)[1,1]) - 1, exp(confint(econometric_model3)[1,2]) - 1,
  "Model with both period and department fixed effects", exp(econometric_model4$coefficients[1]) - 1, exp(confint(econometric_model4)[1,1]) - 1, exp(confint(econometric_model4)[1,2]) - 1
)

coefficient_plot_data$model <- fct_relevel(
  coefficient_plot_data$model,
  c(
    "Model with no fixed effect",
    "Model with period fixed effect",
    "Model with department fixed effect",
    "Model with both period and department fixed effects"
  )
)

ggplot(coefficient_plot_data, aes(x = model, y = estimate)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = .1)) +
  geom_point(size = 2, color = "steelblue") +
  theme_bw() +
  ggtitle("Estimate of B.1.1.7's transmissibility advantage with 95% confidence interval according to various econometric models") +
  xlab("") +
  ylab("B.1.1.7's transmissibility advantage") +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
  ) +
  scale_x_discrete() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggsave("Figures/Econometric Analysis/Estimate of B.1.1.7's transmissibility advantage according to various econometric models.png", width = 12, height = 6)

# I remove the estimate for the Creuse (department number 23) for week 11 because due to the fact that only
# a tiny number of cases were recorded that week and the previous one in that department (it has a very small
# population), the estimate of the advantage is ridiculously large for that week and, while this doesn't affect
# the fit, it makes the graph more difficult to read by increasing the range of the y-axis a lot
ggplot(department_growth_prevalence_data %>% filter(week != "Week 10 to week 11" | dep != "23"), aes(x = prevalence_b117, y = advantage)) +
  geom_point(aes(group = week, color = week)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "B.1.1.7's transmissibility advantage vs. prevalence of B.1.1.7 at the department level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("B.1.1.7's transmissibility advantage") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/B.1.1.7's transmissibility advantage vs. prevalence of B.1.1.7 at the department level.png", width = 12, height = 6)

ggplot(department_growth_prevalence_data, aes(x = prevalence_b117, y = R_b117)) +
  geom_point(aes(group = week, color = week)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "B.1.1.7's effective reproduction number vs. prevalence of B.1.1.7 at the department level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("B.1.1.7's effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/B.1.1.7's effective reproduction number vs. prevalence of B.1.1.7 at the department level.png", width = 12, height = 6)

ggplot(department_growth_prevalence_data, aes(x = prevalence_b117, y = R_total)) +
  geom_point(aes(group = week, color = week)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of B.1.1.7 at the department level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("Effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/Econometric Analysis/Effective reproduction number vs. prevalence of B.1.1.7 at the department level.png", width = 12, height = 6)

ggplot(department_growth_prevalence_data, aes(x = prevalence_b117, y = R_total)) +
  geom_point(color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x, color = "red") +
  facet_wrap(~ week, ncol = 2) +
  scale_x_continuous(labels = scales::percent) +
  theme_bw() +
  ggtitle(
    "Effective reproduction number vs. prevalence of B.1.1.7 at the department level",
    subtitle = "(growth rates are converted growth into effective reproduction numbers by assuming the generation time distribution has a mean of 4.8 days)"
    ) +
  xlab("Prevalence of B.1.1.7") +
  ylab("Effective reproduction number") +
  scale_color_discrete(name = "Period") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(caption = "Source: Santé publique France - Computation and chart by Philippe Lemoine (@phl43)") +
  ggsave("Figures/Econometric Analysis/Effective reproduction number vs. prevalence of B.1.1.7 at the department level by week.png", width = 12, height = 6)

#########################################################################################################
#                                        MATHEMATICAL NOTES                                             #
#########################################################################################################

# For the econometric analysis, I use the formula R ~= (1 + r)^g to convert incidence growth rates into effective
# reproduction numbers. This is a cruder approximation than what I use to replicate Gaymard et al.'s results, which unlike
# the latter doesn't involve the standard deviation of the generation time distribution, but it's good enough and often
# used in the literature. Here is the sketch of a proof that it holds approximately. As equation (2) in this
# paper (https://www.medrxiv.org/content/medrxiv/early/2020/03/30/2020.03.27.20045575.full.pdf) states, in a model
# where all infectors are assumed to transmit the virus exactly g units of time after getting infected, r = log(R) / g.
# It follows that R_0 = exp(r * g). Since log(1 + r) ~= r and exp(log(x)) = x, 1 + r ~= exp(r), so (1 + r)^g. So
# R ~= (1 + r)^g in a model where all infectors are assumed to transmit the virus exactly g units of time after getting
# infected.

# While we're at it, it's easy to show that, in such a model, a longer generation time implies a higher transmissibility
# advantage for B.1.1.7. Here is a quick proof. Let r_v be the growth rate of B.1.1.7, r_h the growth rate of the
# historical lineage and g’ = g + a. By the previous result, R_v ~= (1 + r)^g and R_h = (1 + r)^g if g is the mean
# of the generation time distribution and R’_v ~= (1 + r)^(g + a) and R’_h = (1 + r)^(g + a) if g’ is the mean of
# the generation time distribution, where R_v is the R of B.1.1.7 and R_h is the R of the historical lineage. Therefore
# (R’_v / R’_h) / (R_v / R_h) = ((1 + r_v) / (1 + r_h))^a, which provided that r_v > r_h, so that (1 + r_v) / (1 + r_h) > 1,
# is greater than 1 if and only if a > 0. Conversely, it seems that other things being equal a larger standard deviation
# of the generation time distribution implies a lower transmissibility advantage for B.1.1.7, but proving that would
# require assuming a more complicated model and I haven't tried.
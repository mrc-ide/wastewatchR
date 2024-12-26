# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2)

# Run branching process
population <- 1000000
model_output <- simulate_branching_process(initial_mn_offspring = 0.15,
                                           disp_offspring = 1.0,
                                           generation_time_dist = function(n) { rgamma(n, shape = 12, rate = 2) },
                                           prob_symptomatic = 0.8,
                                           infection_to_onset_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           prob_severe = 0.35,
                                           prob_seek_healthcare_non_severe = 0.9,
                                           prob_seek_healthcare_severe = 0.9,
                                           onset_to_healthcare_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           prob_seroconvert_asymptomatic = 0.8,
                                           prob_seroconvert_severe = 0.8,
                                           prob_seroconvert_non_severe = 0.8,
                                           infection_to_seroconversion_dist = function(n) { rgamma(n, shape = 56, rate = 2) },
                                           seroconversion_to_seroreversion_dist = function(n) { rgamma(n, shape = 730, rate = 2) },
                                           prob_beneficial_mutation = 0,
                                           beneficial_mutation_effect_dist = function(n) {rexp(n, rate = 10) },
                                           max_mn_offspring = 1.5,
                                           annual_spillover_rate = 2,
                                           spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
                                           initial_immune = 0,
                                           t0 = 0,
                                           tf = Inf,
                                           population = population, # always make sure this is a pretty large number unless deliberately want susceptible depletion
                                           check_final_size = 10000,
                                           max_num_outbreaks = 200,
                                           seed = 2025)

infections <- generate_infections_time_series(branching_process_output = model_output)
shedding_dist <- EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16) / max(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16))
# note that shedding dist is normalised to highest values to match the approach used in Hewitt et al.
# can think of this quantity as "the effective number of people shedding at their peak shedding time"
shedding <- generate_number_shedding_time_series(branching_process_output = model_output,
                                                 shedding_dist = shedding_dist,
                                                 shedding_relative_SC2 = 1)
plot(shedding$day, shedding$shedding_value)
plot(shedding$day, shedding$new_infections)

plot(infections$day, infections$new_infections)

ggplot(infections, aes(x = day, y = new_infections)) +
  geom_line(col = "#8075FF") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 10)) +
  labs(x = "Time (Days)", y = "Clinical Visits Per Day", title = "Clinical Surveillance") +
  theme_bw()

ggplot(shedding, aes(x = day, y = shedding_value)) +
  geom_line(col = "#1DBD83") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 10)) +
  labs(x = "Time (Days)", y = "Effective Number of People Shedding", title = "Wastewater Surveillance") +
  theme_bw()


calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 14,
                         sampling_method = "autosampler",
                         detection_approach = "threshold",
                         detection_params = list(population = 10000,
                                                 threshold_limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 80, 150, 250)))

calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 14,
                         sampling_method = "moore_swab",
                         detection_approach = "threshold",
                         detection_params = list(population = 10000,
                                                 threshold_limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 80, 150, 250)))

calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 1,
                         sampling_method = "autosampler",
                         detection_approach = "probit_curve",
                         detection_params = list(population = 100000,
                                                 logistic_beta_0 = -1.229996,
                                                 logistic_beta_1 = 0.258775,
                                                 limit_of_detection = 3,
                                                 seed = 1000))

calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 1,
                         sampling_method = "moore_swab",
                         detection_approach = "probit_curve",
                         detection_params = list(population = 100000,
                                                 probit_beta_0 = -3,
                                                 probit_beta_1 = 1,
                                                 seed = 100))


x <- seq(0, 15, 0.01)
y <- plogis(-5 + x * 1)
plot(x, y, type = "l", xlab = "Effective Number Shedding", ylab = "Prob. Detection", las = 1, lwd = 4)


cowplot::plot_grid(a, b, c, nrow = 3)

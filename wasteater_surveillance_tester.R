# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2)

# Run branching process
population <- 1000000
model_output <- simulate_branching_process(initial_mn_offspring = 0.85,
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

infections <- generate_infections_time_series(branching_process_output = model_output, population = population)
shedding_dist <- EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16)
shedding <- generate_number_shedding_time_series(branching_process_output = model_output, shedding_dist = shedding_dist, population = population)

## something important for me to work out here is whether I want shedding_dist to sum to 1 (so proportion of total shed over time)
## or whether I want to normalise it to peak shedding (i.e. make peak shedding the max). This is unclear to me currently.
## I think in practice it doesn't matter too much as long as I make sure I'm using the same quantity to pin down the sensitivity.

ggplot(infections, aes(x = day, y = new_infections)) +
  geom_line(col = "#8075FF") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 10)) +
  labs(x = "Time (Days)", y = "Clinical Visits Per Day", title = "Clinical Surveillance") +
  theme_bw()

ggplot(shedding, aes(x = day, y = shedding_value)) +
  geom_line(col = "#1DBD83") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 10)) +
  labs(x = "Time (Days)", y = "Log(Gene Copies per ml)", title = "Wastewater Surveillance") +
  theme_bw()


calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 14,
                         sampling_method = "autosampler",
                         detection_approach = "threshold",
                         detection_params = list(population = 100000,
                                                 threshold_limits = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 40, 80, 150, 250)))

calculate_wastewater_ttd(wastewater_number_shedding_time_series = shedding,
                         sampling_frequency = 1,
                         sampling_method = "autosampler",
                         detection_approach = "probit_curve",
                         detection_params = list(population = 100000,
                                                 probit_beta_0 = -3,
                                                 probit_beta_1 = 1))


calculate_wastewater_ttd <- function(wastewater_number_shedding_time_series,
                                     sampling_frequency,
                                     sampling_method,
                                     detection_approach,
                                     detection_params) {

  ## Checking that sampling_frequency is an integer
  if ((sampling_frequency / floor(sampling_frequency)) != 1) {
    stop("sampling_frequency must be an integer")
  }
  if (sampling_frequency < 1) {
    stop("sampling_frequency must be greater than or equal to 1 (corresponding to daily sampling)")
  }

  ## Checking that the user has specified a suitable smpling method
  if (!(sampling_method %in% c("autosampler", "grab", "moore_swab"))) {
    stop("Error - sampling_method must be one of autosampler, grab or moore_swab")
  }
  if (sampling_method %in% c("autosampler", "grab")) {
    warning("Note that as we don't do sub-daily time-resolution atm, there is no difference in our approach to representing autosampling and grab")
  }

  ## Checking that the user has specified a suitable detection type
  if (!(detection_approach %in% c("threshold", "probit_curve", "per_person_probability"))) {
    stop("detection_approach must be one of threshold, probit_curve or per_person_probability")
  }

  ## Checking that detection_params is a list
  if (!is.list(detection_params)) {
    stop("detection_params must be a list containing detection_approach-specific parameters")
  }

  ## For detection_approach == "threshold", ttd is the first time at which the effective number
  ## of shedding individuals eclipses said threshold
  if (detection_approach == "threshold") {

    if (!("threshold_limits" %in% names(detection_params))) {
      stop("detection_params must contain a named element called threshold_limits which contains reals that specify sensitivity per 100,000 population")
    }
    if (!("population" %in% names(detection_params))) {
      stop("detection_params must contain a named element called population which contains the size of the population")
    }
    if (sum(!is.numeric(detection_params$threshold_limits)) > 0) {
      stop("threshold_limits must only contain numerics")
    }
    wastewater_shedding_ttd <- tibble(threshold = detection_params$threshold_limits) %>%
      rowwise() %>%
      mutate(wastewater_first_day = {
        filtered_data <- wastewater_number_shedding_time_series %>%
          filter(day %% sampling_frequency == 0) %>%
          filter((100000 * shedding_value / detection_params$population) >= threshold)
          if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
        }) %>%
      ungroup()
  }

  ## For detection_approach == "probit_curve", the effective number of shedding at each sampling timepoint is converted
  ## to a probability and a draw done from a bernoulli. ttd is the first time at which the bernoulli is successful.
  if (detection_approach == "probit_curve") {

    # Check that the necessary probit parameters exist
    # (Rename these as needed, e.g. "probit_intercept", "probit_slope", etc.)
    if (!all(c("probit_beta_0", "probit_beta_1") %in% names(detection_params))) {
      stop("For detection_approach == 'probit_curve', detection_params must contain probit_beta_0 and probit_beta_1.")
    }

    # Filter to the sampling days
    sampled_data <- wastewater_number_shedding_time_series %>%
      filter(day %% sampling_frequency == 0) %>%
      mutate(
        # Convert 'shedding_value' to a probability of detection via probit
        prob_detect = pnorm(
          detection_params$probit_beta_0 +
            detection_params$probit_beta_1 *
            (100000 * shedding_value / detection_params$population)),
        # Draw once from a Bernoulli with this probability
        detect_draw = rbinom(n = n(), size = 1, prob = prob_detect)
      )

    # The time-to-detection is the first sampled day at which detect_draw == 1
    detection_day <- sampled_data %>%
      ungroup() %>%
      filter(detect_draw == 1) %>%
      summarize(first_day = ifelse(n() == 0, NA_real_, min(day))) %>%
      pull(first_day)

    # Return as a tibble with a single row
    wastewater_shedding_ttd <- tibble(wastewater_first_day = detection_day)

  }

  return(wastewater_shedding_ttd)

}

cowplot::plot_grid(a, b, c, nrow = 3)

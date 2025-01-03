#' Utility function to calculate shedding given infection time-series and shedding distribution
#'
#' @param day The day being considered
#' @param infection_counts A time series of the number of individuals infected on each day
#' @param shedding_dist The distribution of amount shed over time following infection, normalised so the amount of shedding on the day
#' with the maximum amount is 1.
#' @param shedding_relative_SC2 Amount of shedding that occurs relative to SARS-CoV-2 (required given our pinning to the SC2 data from Hewitt et al)
#' @family utils
#' @export
calculate_shedding <- function(day, infection_counts, shedding_dist, shedding_relative_SC2) {
  # Shift infection counts based on the day difference
  shedding_contributions <- sapply(1:length(shedding_dist), function(i) {
    lagged_day <- day - (i - 1)
    if (lagged_day >= 0) {
      return(infection_counts$new_infections[lagged_day + 1] * shedding_dist[i] * shedding_relative_SC2)
    } else {
      return(0)
    }
  })
  # Sum up all the contributions
  return(sum(shedding_contributions))
}

#' Convert stochastic realisation of branching process into number shedding time-series
#'
#' This function converts a branching process output into time-series of number shedding,
#' weighted by the shedding distribution.
#'
#' @param branching_process_output Output from simulate_branching_process
#' @param shedding_dist The distribution of amount shed over time following infection, normalised so the amount of shedding on the day
#' with the maximum amount is 1.
#' @param shedding_relative_SC2 Amount of shedding that occurs relative to SARS-CoV-2 (required given our pinning to the SC2 data from Hewitt et al)
#'
#' @family post-processing
#' @export
generate_number_shedding_time_series <- function(branching_process_output, shedding_dist, shedding_relative_SC2) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)
  infection_counts <- generate_infections_time_series(branching_process_output)

  # Apply the function for each day
  shedding_results <- tibble(day = days) |>
    rowwise() |>
    mutate(shedding_value = calculate_shedding(day, infection_counts, shedding_dist, shedding_relative_SC2))
  shedding_results <- shedding_results %>%
    left_join(infection_counts, by = "day")

  return(shedding_results)

}

#' @export
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
  if (sampling_method == "moore_swab") {
    wastewater_number_shedding_time_series <- wastewater_number_shedding_time_series %>%
      ungroup() %>%
      mutate(shedding_value = rollapply(data = shedding_value,
                                        width = 7,
                                        FUN = mean,
                                        align = "right",
                                        partial = TRUE))
  }
  # if (sampling_method %in% c("autosampler", "grab")) {
  #   warning("Note that as we don't do sub-daily time-resolution atm, there is no difference in our approach to representing autosampling and grab")
  # }

  ## Checking that the user has specified a suitable detection type
  if (!(detection_approach %in% c("threshold", "logistic_curve", "per_person_probability"))) {
    stop("detection_approach must be one of threshold, logistic_curve or per_person_probability")
  }

  ## Checking that detection_params is a list
  if (!is.list(detection_params)) {
    stop("detection_params must be a list containing detection_approach-specific parameters")
  }

  ## For detection_approach == "threshold", ttd is the first time at which the effective number
  ## of shedding individuals eclipses said threshold
  if (detection_approach == "threshold") {

    # Check that the necessary logistic parameters exist
    if (!all(c("threshold_limits", "population") %in% names(detection_params))) {
      stop("For detection_approach == 'threshold', detection_params must contain threshold_limits and population")
    }
    if (sum(!is.numeric(detection_params$threshold_limits)) > 0) {
      stop("threshold_limits must only contain numerics")
    }
    calculated_wastewater_shedding_ttd <- tibble(threshold = detection_params$threshold_limits) %>%
      rowwise() %>%
      mutate(wastewater_first_day = {
        filtered_data <- wastewater_number_shedding_time_series %>%
          filter(day %% sampling_frequency == 0) %>%
          filter((100000 * shedding_value / detection_params$population) >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
      ungroup()
  }

  ## For detection_approach == "logistic_curve", the effective number of shedding at each sampling timepoint is converted
  ## to a probability and a draw done from a bernoulli. ttd is the first time at which the bernoulli is successful.
  if (detection_approach == "logistic_curve") {

    # Check that the necessary logistic parameters exist
    if (!all(c("logistic_beta_0", "logistic_beta_1", "seed", "limit_of_detection", "population") %in% names(detection_params))) {
      stop("For detection_approach == 'logistic_curve', detection_params must contain logistic_beta_0, logistic_beta_1, limit_of_detection and a seed")
    }

    # Filter to the sampling days
    set.seed(detection_params$seed)
    sampled_data <- wastewater_number_shedding_time_series %>%
      filter(day %% sampling_frequency == 0) %>%
      mutate(
        # Convert 'shedding_value' to a probability of detection via probit
        prob_detect = ifelse(shedding_value < detection_params$limit_of_detection, 0,
                             plogis(
                               detection_params$logistic_beta_0 +     # -1.229996 from Hewitt et al Fig 5B
                               detection_params$logistic_beta_1 *     # 0.258775 from Hewitt et al Fig 5B
                                  (100000 * shedding_value / detection_params$population))),
        # Draw once from a Bernoulli with this probability
        detect_draw = rbinom(n = n(), size = 1, prob = prob_detect))

    # The time-to-detection is the first sampled day at which detect_draw == 1
    detection_day <- sampled_data %>%
      ungroup() %>%
      dplyr::filter(detect_draw == 1) %>%
      dplyr::summarize(first_day = ifelse(n() == 0, NA_real_, min(day))) %>%
      dplyr::pull(first_day)

    # Return as a tibble with a single row
    calculated_wastewater_shedding_ttd <- tibble(wastewater_first_day = detection_day)

  }

  ## For detection_approach == "per_person_probability", we look at new_infections on each sampling day, do a
  ## binomial draw with size = new_infections and prob = detection_params$per_infection_probability_detection.
  ## The time to detection is the first day that draw > 0.
  ## Note: We do this based on the timing of the infection and so the timing isn't quite right - we would have do
  ## something weird with the shedding dist to fully capture this approach.
  if (detection_approach == "per_person_probability") {
    # Check that the necessary parameter is present
    if (!all(c("per_infection_probability_detection") %in% names(detection_params))) {
      stop("For detection_approach == 'per_person_probability', detection_params must contain per_infection_probability_detection")
    }

    # Filter by sampling frequency
    sampled_data <- wastewater_number_shedding_time_series %>%
      filter(day %% sampling_frequency == 0) %>%
      mutate(
        detect_draw = rbinom(
          n = n(),
          size = new_infections,
          prob = detection_params$per_infection_probability_detection
        )
      )

    # TTD is the first day the binomial draw is > 0
    detection_day <- sampled_data %>%
      ungroup() %>%
      dplyr::filter(detect_draw > 0) %>%
      dplyr::summarize(first_day = ifelse(n() == 0, NA_real_, min(day))) %>%
      dplyr::pull(first_day)

    calculated_wastewater_shedding_ttd <- tibble(wastewater_first_day = detection_day)
  }

  return(calculated_wastewater_shedding_ttd)

}

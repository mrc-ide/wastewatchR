#' Utility function to calculate shedding given infection time-series and shedding distribution
#'
#' @param day The day being considered
#' @param infection_counts A time series of the number of individuals infected on each day
#' @param method The method to use - either "effective_n_shedders" to use the effective number of shedders based on the shedding profile "shedding_dist" (see below) (equivalent to Hewitt Model 3) or "n_shedders" to use the number of shedders (equivalent to Hewitt Model 2) on any one day
#' @param shedding_dist The distribution of amount shed over time following infection, normalised so the amount of shedding on the day
#' with the maximum amount is 1. Used if method = "effective_n_shedders"
#' @param duration duration of shedding, if method = "n_shedders". Should be NA if using shedding_dist
#' @param shedding_relative_SC2 Amount of shedding that occurs relative to SARS-CoV-2 (required given our pinning to the SC2 data from Hewitt et al)
#' @family utils
#' @export
calculate_shedding <- function(day, infection_counts, method, shedding_dist, duration, shedding_relative_SC2) {
  # Shift infection counts based on the day difference

  if(method == "effective_n_shedders"){

    shedding_contributions <- sapply(1:length(shedding_dist), function(i) {
      lagged_day <- day - (i - 1)
      if (lagged_day >= 0) {
        return(infection_counts$new_infections[lagged_day + 1] * shedding_dist[i] * shedding_relative_SC2)
      } else {
        return(0)
      }
    })
  }

  if(method == "n_shedders"){

    shedding_contributions <- sapply(1:duration, function(i) {
      lagged_day <- day - (i - 1)
      if (lagged_day >= 0) {
        return(infection_counts$new_infections[lagged_day + 1] * shedding_relative_SC2)
      } else {
        return(0)
      }
    })
  }

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
#' @param method The method to use - either "effective_n_shedders" to use the effective number of shedders based on the shedding profile "shedding_dist" (equivalent to Hewitt Model 3) or "n_shedders" to use the number of shedders (equivalent to Hewitt Model 2) on any one day
#' @param duration duration of shedding, if method = "n_shedders". Should be NA if using shedding_dist
#' @family post-processing
#' @export
generate_number_shedding_time_series <- function(branching_process_output, shedding_dist, shedding_relative_SC2, method, duration) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)
  infection_counts <- generate_infections_time_series(branching_process_output)

  # Apply the function for each day
  shedding_results <- tibble(day = days) |>
    rowwise() |>
    mutate(shedding_value = calculate_shedding(day, infection_counts, method, shedding_dist, duration, shedding_relative_SC2))
  shedding_results <- shedding_results %>%
    left_join(infection_counts, by = "day")

  return(shedding_results)

}

#' @param wastewater_number_shedding_time_series Output from generate_number_shedding_time_series
#' @param sampling_frequency The frequency of wastewate sampling, with 1 being daily, 7 being weekly, 14 being fortnightly
#' @param sampling_method The method used, either "autosampler", "grab" or "moore_swab"
#' @param detection_approach Approach to translating shedding into probability of detection, either "threshold", "logistic_curve" or "per_person_probability"
#' @param detection_params A named list. If detection_approach == "threshold", must contain "threshold_limits" and "population". If detection_approach == "logistic_curve", must contain "logistic_beta_0", "logistic_beta_1", "seed", "limit_of_detection", "population". If detection_approach == "per_person_probability", must contain "per_infection_probability_detection". Additionally, if "sampling_method" == "moore_swab", "detection_params must also contain "duration" in days (>=1).
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

  ## Checking that the user has specified a suitable sampling method
  if (!(sampling_method %in% c("autosampler", "grab", "moore_swab"))) {
    stop("Error - sampling_method must be one of autosampler, grab or moore_swab")
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

  ## Checking that duration is speciifed if the user has selected moore_swab as sampling_method
  if (sampling_method == "moore_swab" ){
    stop("if sampling method is moore_swab,'detection_params$duration' must be a number of days >=1")
  }

  ## If sampling method = moore_swab, calculate mean shedding over moore swab window
  if (sampling_method == "moore_swab") {
    wastewater_number_shedding_time_series <- wastewater_number_shedding_time_series %>%
      ungroup() %>%
      mutate(shedding_value = rollapply(data = shedding_value,
                                        width = detection_params$duration,
                                        FUN = mean,
                                        align = "right",
                                        partial = TRUE))
  }

  ## To ensure first date of sampling is random with respect to the beginning of an outbreak
  ## set a random nudge >=0 and < sampling_frequency
  set.seed(detection_params$seed)

  x <- if (sampling_frequency == 1) { # ensuring daily sampling (sampling_frequency == 1) starts on day 0
    0
  } else {
    sample(0:(sampling_frequency - 1), 1)
  }

  ## For detection_approach == "threshold", ttd is the first time at which the effective number
  ## of shedding individuals eclipses said threshold
  ## NOTE - we do not currently output probability of detection at the population level as this approach does not lend itself to this as prob is either 1 or 0
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
          filter((day-x) %% sampling_frequency == 0) %>%
          filter((100000 * shedding_value / detection_params$population) >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
      ungroup()

    sampled_data <- NA # for consistency across detection_approaches
  }

  ## For detection_approach == "logistic_curve", the effective number of shedding at each sampling timepoint is converted
  ## to a probability and a draw done from a bernoulli. ttd is the first time at which the bernoulli is successful.
  if (detection_approach == "logistic_curve") {

    # Check that the necessary logistic parameters exist
    if (!all(c("logistic_beta_0", "logistic_beta_1", "seed", "limit_of_detection", "population") %in% names(detection_params))) {
      stop("For detection_approach == 'logistic_curve', detection_params must contain logistic_beta_0, logistic_beta_1, limit_of_detection, population, and a seed")
    }

    # Filter to the sampling days
    sampled_data <- wastewater_number_shedding_time_series %>%
      mutate(
        # Convert 'shedding_value' to a probability of detection via logistic curve
        prob_detect = ifelse(shedding_value < detection_params$limit_of_detection, 0,
                             plogis(
                               detection_params$logistic_beta_0 +     # -1.229996 from Hewitt et al Fig 5B for Model 3
                               detection_params$logistic_beta_1 *     # 0.258775 from Hewitt et al Fig 5B for Model 3
                                  (log10(100000 * (shedding_value+1e-3) / detection_params$population)))), # 100000 = population used by Hewitt et. al.
        sampled = ifelse((day-x) %% sampling_frequency == 0, "Yes", "No"),
        detect_draw = ifelse((day-x) %% sampling_frequency == 0, rbinom(n = n(), size = 1, prob = prob_detect), 0))# Draw once from a Bernoulli with this probability

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
  ## NOTE - we do not currently output probability of detection at the population level as this approach does not lend itself to this
  if (detection_approach == "per_person_probability") {
    # Check that the necessary parameter is present
    if (!all(c("per_infection_probability_detection") %in% names(detection_params))) {
      stop("For detection_approach == 'per_person_probability', detection_params must contain per_infection_probability_detection")
    }

    # Filter by sampling frequency
    sampled_data <- wastewater_number_shedding_time_series %>%
      filter(day-x %% sampling_frequency == 0) %>%
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

  return(list(ttd = calculated_wastewater_shedding_ttd,
              sampled_data = sampled_data))

}

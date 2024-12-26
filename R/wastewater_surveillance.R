#' Utility function to calculate shedding given infection time-series and shedding distribution
#'
#' @param day The day being considered
#' @param infection_counts A time series of the number of individuals infected on each day
#' @param shedding_dist The distribution of amount shed over time following infection.
#'
#' @family utils
#' @export
calculate_shedding <- function(day, infection_counts, shedding_dist) {
  # Shift infection counts based on the day difference
  shedding_contributions <- sapply(1:length(shedding_dist), function(i) {
    lagged_day <- day - (i - 1)
    if (lagged_day >= 0) {
      return(infection_counts$new_infections[lagged_day + 1] * shedding_dist[i])
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
#' @param shedding_dist The distribution of amount shed over time following infection.
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @export
generate_number_shedding_time_series <- function(branching_process_output, shedding_dist, population) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)
  infection_counts <- generate_infections_time_series(branching_process_output, population)

  # Apply the function for each day
  shedding_results <- tibble(day = days) |>
    rowwise() |>
    mutate(shedding_value = calculate_shedding(day, infection_counts, shedding_dist)) %>%
    mutate(shedding_value_per_thousand = 1000 * shedding_value / population)
  return(shedding_results)

}


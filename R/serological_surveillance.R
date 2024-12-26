#' Convert stochastic realisation of branching process into hospitalisation time-series
#'
#' This function converts a branching process output into time-series of hospitalisations
#' and hospital occupancy.
#'
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @export
generate_seropositivity_timeseries <- function(branching_process_output, population) {

  ## Calculating daily incidence of hospitalisations
  seroconversion_df <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
    filter(seroconvert == 1) %>%
    mutate(time_seroconversion_floor = floor(time_seroconversion)) |>
    group_by(time_seroconversion_floor) %>%
    summarise(n_seroconverted = n()) %>%
    complete(time_seroconversion_floor = seq(0, max(time_seroconversion_floor, na.rm = TRUE), by = 1),
             fill = list(n_seroconverted = 0)) |>
    mutate(cumulative_seroconversions = cumsum(n_seroconverted))

  seroreversion_df <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
    filter(seroconvert == 1) %>%
    mutate(time_seroreversion_floor = floor(time_seroreversion)) |>
    group_by(time_seroreversion_floor) %>%
    summarise(n_seroreverted = n()) %>%
    complete(time_seroreversion_floor = seq(0, max(time_seroreversion_floor, na.rm = TRUE), by = 1),
             fill = list(n_seroreverted = 0)) |>
    mutate(cumulative_seroreversions = cumsum(n_seroreverted))

  max_seroconversion_value <- max(seroconversion_df$cumulative_seroconversions)
  sero_df <- seroconversion_df %>%
    full_join(seroreversion_df, by = c("time_seroconversion_floor" = "time_seroreversion_floor")) %>%
    rename(time = time_seroconversion_floor) %>%
    complete(time = seq(0, max(time, na.rm = TRUE), by = 1),
             fill = list(n_seroconverted = 0)) %>%
    complete(time = seq(0, max(time, na.rm = TRUE), by = 1),
             fill = list(cumulative_seroconversions = max_seroconversion_value)) %>%
    mutate(seropositive_abs = cumulative_seroconversions - cumulative_seroreversions) %>% ## NOTE - CHECK THIS IS CORRECT!!
    mutate(seropositive_prop_pop = seropositive_abs / population)

  return(sero_df)
}

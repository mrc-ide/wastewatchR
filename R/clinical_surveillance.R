#' Convert stochastic realisation of branching process into symptom onset time-series
#' This function converts a branching process output into time-series of symptom onsets
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @import dplyr
#' @import tidyr
#' @export
generate_symptom_onset_time_series <- function(branching_process_output) {

  ## Generating base time-series to then modify
  max_day <- max(c(branching_process_output$time_infection,
                   branching_process_output$time_symptom_onset,
                   branching_process_output$time_seek_healthcare),
                 na.rm = TRUE)
  symptom_onset_incidence_base <- data.frame(day = 1:floor(max_day), incidence_symptom_onset = 0)

  ## Calculating daily incidence of symptom onsets
  symptom_onset_incidence <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
    filter(symptomatic == 1)

  ## If no incidence, just return dataframe of 0s
  if (nrow(symptom_onset_incidence) == 0) {
    symptom_onset_incidence <- symptom_onset_incidence_base

  ## Else, modify the base df as appropriate
  } else {
    symptom_onset_incidence <- symptom_onset_incidence |>
      mutate(time_symptom_onset_floor = floor(time_symptom_onset)) |>
      group_by(time_symptom_onset_floor) |>
      summarise(incidence_symptom_onset = n()) |>
      complete(time_symptom_onset_floor = seq(0, max(time_symptom_onset_floor, na.rm = TRUE), by = 1),
               fill = list(incidence_symptom_onset = 0)) |>
      rename(day = time_symptom_onset_floor)
    index_to_replace <- symptom_onset_incidence_base$day %in% symptom_onset_incidence$day
    symptom_onset_incidence_base$incidence_symptom_onset[index_to_replace] <- symptom_onset_incidence$incidence_symptom_onset
    symptom_onset_incidence <- symptom_onset_incidence_base
  }
  return(symptom_onset_incidence)
}

#' Convert stochastic realisation of branching process into hospitalisation time-series
#'
#' This function converts a branching process output into time-series of hospitalisations
#' and hospital occupancy.
#'
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#'
#' @family post-processing
#' @export
generate_healthcare_seeking_time_series <- function(branching_process_output) {

  ## Generating base time-series to then modify
  max_day <- max(c(branching_process_output$time_infection,
                   branching_process_output$time_symptom_onset,
                   branching_process_output$time_seek_healthcare),
                 na.rm = TRUE)
  healthcare_seeking_incidence_base <- data.frame(day = 1:floor(max_day), incidence_seek_healthcare = 0)

  ## Calculating daily incidence of healthcare seeking infections
  healthcare_seeking_incidence <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
    filter(seek_healthcare == 1)

  ## If no incidence, just return dataframe of 0s
  if (nrow(healthcare_seeking_incidence) == 0) {
    healthcare_seeking_incidence <- healthcare_seeking_incidence_base

    ## Else, modify the base df as appropriate
  } else {
    healthcare_seeking_incidence <- healthcare_seeking_incidence |>
      mutate(time_seek_healthcare_floor = floor(time_seek_healthcare)) |>
      group_by(time_seek_healthcare_floor) |>
      summarise(incidence_seek_healthcare = n()) |>
      rename(day = time_seek_healthcare_floor)
    index_to_replace <- healthcare_seeking_incidence_base$day %in% healthcare_seeking_incidence$day
    healthcare_seeking_incidence_base$incidence_seek_healthcare[index_to_replace] <- healthcare_seeking_incidence$incidence_seek_healthcare
    healthcare_seeking_incidence <- healthcare_seeking_incidence_base
  }
  return(healthcare_seeking_incidence)
}


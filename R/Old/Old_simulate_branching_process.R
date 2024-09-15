#' Generate a stochastic realisation from a stochastic branching process
#'
#' This function creates a named list containing the output of simulating a stochastic branching
#' process.
#'
#' @family simulation
#' @export
simulate_branching_process <- function(

    ## Offspring distribution parameters
    mn_offspring = 3,
    disp_offspring = 1.1,

    ## Natural history parameters
    generation_time_dist = function(n) { rgamma(n, shape = 12, rate = 2) }, ## poss need to reverse rate and shape - TBD

    ## Disease severity parameters
    prob_symptomatic = 0.8,
    infection_to_onset_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
    prob_severe = 0.35,
    prob_seek_healthcare_non_severe = 0.9,
    prob_seek_healthcare_severe = 0.9,
    onset_to_healthcare_dist = function(n) { rgamma(n, shape = 6, rate = 2) },

    ## Serology related parameters
    prob_seroconvert_asymptomatic = 0.8,
    prob_seroconvert_severe = 0.8,
    prob_seroconvert_non_severe = 0.8,
    infection_to_seroconversion_dist = function(n) { rgamma(n, shape = 56, rate = 2) },
    seroconversion_to_seroreversion_dist = function(n) { rgamma(n, shape = 730, rate = 2) },

    ## Simulation parameters
    initial_immune = 0,
    seeding_cases = 10,
    t0 = 0,
    tf = Inf,
    population = 100000,
    check_final_size = 10000

) {

  ## Setting up the Offspring Distribution
  if (disp_offspring <= 1.0) {
    offspring_fun <- function(n, susc) {
      new_mn <- mn_offspring * susc/population
      rpois(n = n, lambda = new_mean)
    }
  } else {
    offspring_fun <- function(n, susc) {
      new_mn <- mn_offspring * susc/population
      size <- new_mn/(disp_offspring - 1)
      rnbinom(n = n, size = size, mu = new_mn)
    }
  }

  # Pre-allocate a dataframe with the maximum size
  tdf <- data.frame(
    infection_id = integer(check_final_size),
    infection_ancestor = integer(check_final_size),
    infection_generation = integer(check_final_size),
    time_infection = NA_real_,
    symptomatic = integer(check_final_size),
    time_symptom_onset = numeric(check_final_size),
    severe = integer(check_final_size),
    seek_healthcare = integer(check_final_size),
    time_seek_healthcare = numeric(check_final_size),
    seroconvert = integer(check_final_size),
    time_seroconversion = numeric(check_final_size),
    time_seroreversion = numeric(check_final_size),
    n_offspring = integer(check_final_size),
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )

  # Initialising the model with the seeding cases
  seeding_cases_time_infection <- t0 + seq(from = 0, to = 0.01, length.out = seeding_cases) ## Time of infection for seeding cases

  ## Generating symptom status of seeding cases
  seeding_cases_symptomatic <- sample(c(0, 1), seeding_cases, replace = TRUE,               ## Determining whether seeding cases are symptomatic
                                      prob = c(1 - prob_symptomatic, prob_symptomatic))
  number_seeding_cases_symptomatic <- sum(seeding_cases_symptomatic)                        ## Calculating number of symptomatic index cases
  index_seeding_cases_symptomatic <- which(seeding_cases_symptomatic == 1)                  ## Which seeding cases are symptomatic
  index_seeding_cases_asymptomatic <- which(seeding_cases_symptomatic == 0)                  ## Which seeding cases are symptomatic
  seeding_cases_time_symptom_onset <- rep(NA, seeding_cases)
  seeding_cases_time_symptom_onset[index_seeding_cases_symptomatic] <- seeding_cases_time_infection[index_seeding_cases_symptomatic] + infection_to_onset_dist(number_seeding_cases_symptomatic)

  ## Generating disease severity status of seeding cases
  prob_severe_if_symptomatic <- prob_severe / prob_symptomatic  ## Probability of having severe disease conditional on being symptomatic
  if (prob_severe_if_symptomatic >= 1) {
    stop("prob_severe / prob_symptomatic must be <= 1")
  }
  seeding_cases_severe <- rep(0, seeding_cases)
  seeding_cases_severe[index_seeding_cases_symptomatic] <- sample(c(0, 1), number_seeding_cases_symptomatic, replace = TRUE,                  ## Determining whether symptomatic seeding cases have severe disease
                                                                        prob = c(1 - prob_severe_if_symptomatic, prob_severe_if_symptomatic))
  number_seeding_cases_severe <- sum(seeding_cases_severe)                        ## Calculating number of severe disease seeding cases
  index_seeding_cases_severe <- which(seeding_cases_symptomatic == 1 & seeding_cases_severe == 1)                  ## Which seeding cases are severe disease
  index_seeding_cases_not_severe <- which(seeding_cases_symptomatic == 1  & seeding_cases_severe != 1)

  ## Generating time of seeking healthcare for each seeding case (stratified by disease severity)
  seeding_case_seek_healthcare <- rep(0, seeding_cases)
  seeding_case_seek_healthcare[index_seeding_cases_not_severe] <- sample(c(0, 1), number_seeding_cases_symptomatic - number_seeding_cases_severe,
                                                                         replace = TRUE, prob = c(1 - prob_seek_healthcare_non_severe, prob_seek_healthcare_non_severe))
  seeding_case_seek_healthcare[index_seeding_cases_severe] <- sample(c(0, 1), number_seeding_cases_severe,
                                                                     replace = TRUE, prob = c(1 - prob_seek_healthcare_severe, prob_seek_healthcare_severe))
  number_seeding_cases_seek_healthcare <- sum(seeding_case_seek_healthcare)
  index_seeding_case_seek_healthcare <- which(seeding_case_seek_healthcare == 1)
  seeding_case_time_seek_healthcare <- rep(NA, seeding_cases)
  seeding_case_time_seek_healthcare[index_seeding_case_seek_healthcare] <- seeding_cases_time_symptom_onset[index_seeding_case_seek_healthcare] + onset_to_healthcare_dist(number_seeding_cases_seek_healthcare)

  ## Generating time of seroconversion/seroreversion for each seeding case (stratified by disease severity)
  seeding_case_seroconvert <- rep(0, seeding_cases)
  seeding_case_seroconvert[index_seeding_cases_asymptomatic] <- sample(c(0, 1), seeding_cases - number_seeding_cases_symptomatic,
                                                                       replace = TRUE, prob = c(1 - prob_seroconvert_asymptomatic, prob_seroconvert_asymptomatic))
  seeding_case_seroconvert[index_seeding_cases_not_severe] <- sample(c(0, 1), number_seeding_cases_symptomatic - number_seeding_cases_severe,
                                                                     replace = TRUE, prob = c(1 - prob_seroconvert_non_severe, prob_seroconvert_non_severe))
  seeding_case_seroconvert[index_seeding_cases_severe] <- sample(c(0, 1), number_seeding_cases_severe,
                                                                 replace = TRUE, prob = c(1 - prob_seroconvert_severe, prob_seroconvert_severe))
  number_seeding_cases_seroconvert <- sum(seeding_case_seroconvert)
  index_seeding_case_seroconvert <- which(seeding_case_seroconvert == 1)
  seeding_case_time_seroconvert <- rep(NA, seeding_cases)
  seeding_case_time_seroconvert[index_seeding_case_seroconvert] <- seeding_cases_time_infection[index_seeding_case_seroconvert] + infection_to_seroconversion_dist(number_seeding_cases_seroconvert)
  seeding_case_time_serorevert <- rep(NA, seeding_cases)
  seeding_case_time_serorevert[index_seeding_case_seroconvert] <- seeding_case_time_seroconvert[index_seeding_case_seroconvert] + seroconversion_to_seroreversion_dist(number_seeding_cases_seroconvert)

  ## Initialize the dataframe with the seeding cases
  tdf[1:seeding_cases, ] <- data.frame(
    infection_id = seq_len(seeding_cases),
    infection_ancestor = NA_integer_,
    infection_generation = 1L,
    time_infection = seeding_cases_time_infection,
    symptomatic = seeding_cases_symptomatic,
    time_symptom_onset = seeding_cases_time_symptom_onset,
    severe = seeding_cases_severe,
    seek_healthcare = seeding_case_seek_healthcare,
    time_seek_healthcare = seeding_case_time_seek_healthcare,
    seroconvert = seeding_case_seroconvert,
    time_seroconversion = seeding_case_time_seroconvert,
    time_seroreversion = seeding_case_time_serorevert,
    n_offspring = NA_integer_,
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )

  ## This loops through each infection individually and generates offspring, and the status (symptoms, hospitalisation etc) of all of these offspring
  ## Note that we have check_final_size to cap how long we simulate for, but that when you reach this threshold you might have infections for which
  ## you haven't generated offspring for.
  susc <- population - initial_immune - 1L
  while (susc > 0 & nrow(tdf) <= check_final_size) {

    # Selecting the earliest infection for which we haven't generated offspring
    time_infection_index <- min(tdf$time_infection[tdf$offspring_generated == 0 & !is.na(tdf$time_infection)])

    # Extracting information on the "parent" infection whose offspring are being generated
    parent_idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1] # get the id of the earliest unsimulated infection
    parent_time_infection <- tdf$time_infection[parent_idx]                                       # infection time of the earliest unsimulated infection
    parent_gen <- tdf$infection_generation[parent_idx]                                            # generation of the earliest unsimulated infection
    parent_symptomatic <- tdf$symptomatic[parent_idx]                                             # whether the earliest unsimulated infection is symptomatic
    parent_time_symptom_onset <- tdf$time_symptom_onset[parent_idx]                               # if symptomatic, the time when the earliest unsimulated infection develops symptoms

    # Calculate total number of infections in the dataframe currently (so we can figure out how to label the new infections)
    current_max_id <- max(tdf$infection_id)

    # Calculating number of offspring generated
    n_offspring <- offspring_fun(1, susc)
    tdf$n_offspring[parent_idx] <- n_offspring
    tdf$offspring_generated[parent_idx] <- TRUE
    offspring_time_infection <- parent_time_infection + generation_time_dist(n_offspring)

    ## Generating symptom and hospitalisation status/timing of offspring
    if (n_offspring > 0) {

      ## Symptom status
      offspring_symptomatic <- sample(c(0, 1), n_offspring, replace = TRUE, prob = c(1 - prob_symptomatic, prob_symptomatic))
      number_offspring_symptomatic <- sum(offspring_symptomatic)
      index_offspring_symptomatic <- which(offspring_symptomatic == 1)
      index_offspring_asymptomatic <- which(offspring_symptomatic == 0)
      offspring_time_symptom_onset <- rep(NA, n_offspring)
      offspring_time_symptom_onset[index_offspring_symptomatic] <- offspring_time_infection[index_offspring_symptomatic] + infection_to_onset_dist(number_offspring_symptomatic)

      ## Severe disease
      offspring_severe <- rep(0, n_offspring)
      offspring_severe[index_offspring_symptomatic] <- sample(c(0, 1), number_offspring_symptomatic, replace = TRUE,
                                                              prob = c(1 - prob_severe_if_symptomatic, prob_severe_if_symptomatic))
      number_offspring_severe <- sum(offspring_severe)
      index_offspring_severe <- which(offspring_symptomatic == 1 & offspring_severe == 1)                  ## Which offspring are severe disease
      index_offspring_not_severe <- which(offspring_symptomatic == 1  & offspring_severe != 1)

      ## Seeking healthcare
      offspring_seek_healthcare <- rep(0, n_offspring)
      offspring_seek_healthcare[index_offspring_not_severe] <- sample(c(0, 1), number_offspring_symptomatic - number_offspring_severe,
                                                                      replace = TRUE, prob = c(1 - prob_seek_healthcare_non_severe, prob_seek_healthcare_non_severe))
      offspring_seek_healthcare[index_offspring_severe] <- sample(c(0, 1), number_offspring_severe,
                                                                  replace = TRUE, prob = c(1 - prob_seek_healthcare_severe, prob_seek_healthcare_severe))
      number_offspring_seek_healthcare <- sum(offspring_seek_healthcare)
      index_offspring_seek_healthcare <- which(offspring_seek_healthcare == 1)
      offspring_time_seek_healthcare <- rep(NA, n_offspring)
      offspring_time_seek_healthcare[index_offspring_seek_healthcare] <- offspring_time_symptom_onset[index_offspring_seek_healthcare] + onset_to_healthcare_dist(number_offspring_seek_healthcare)

      ## Time of seroconversion/seroreversion for offspring (stratified by disease severity)
      offspring_seroconvert <- rep(0, n_offspring)
      offspring_seroconvert[index_offspring_asymptomatic] <- sample(c(0, 1), n_offspring - number_offspring_symptomatic,
                                                                    replace = TRUE, prob = c(1 - prob_seroconvert_asymptomatic, prob_seroconvert_asymptomatic))
      offspring_seroconvert[index_offspring_not_severe] <- sample(c(0, 1), number_offspring_symptomatic - number_offspring_severe,
                                                                  replace = TRUE, prob = c(1 - prob_seroconvert_non_severe, prob_seroconvert_non_severe))
      offspring_seroconvert[index_offspring_severe] <- sample(c(0, 1), number_offspring_severe,
                                                              replace = TRUE, prob = c(1 - prob_seroconvert_severe, prob_seroconvert_severe))
      number_offspring_seroconvert <- sum(offspring_seroconvert)
      index_offspring_seroconvert <- which(offspring_seroconvert == 1)
      offspring_time_seroconvert <- rep(NA, n_offspring)
      offspring_time_seroconvert[index_offspring_seroconvert] <- offspring_time_infection[index_offspring_seroconvert] + infection_to_seroconversion_dist(number_offspring_seroconvert)
      offspring_time_serorevert <- rep(NA, n_offspring)
      offspring_time_serorevert[index_offspring_seroconvert] <- offspring_time_seroconvert[index_offspring_seroconvert] + seroconversion_to_seroreversion_dist(number_offspring_seroconvert)

      ## Appending this information on the offspring to the table
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_id"] <- c(current_max_id + seq_len(n_offspring))
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_ancestor"] <- parent_idx
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_generation"] <- parent_gen + 1L
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_infection"] <- offspring_time_infection
      tdf[(current_max_id+1):(current_max_id+n_offspring), "symptomatic"] <- offspring_symptomatic
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_symptom_onset"] <- offspring_time_symptom_onset
      tdf[(current_max_id+1):(current_max_id+n_offspring), "severe"] <- offspring_severe
      tdf[(current_max_id+1):(current_max_id+n_offspring), "seek_healthcare"] <- offspring_seek_healthcare
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_seek_healthcare"] <- offspring_time_seek_healthcare
      tdf[(current_max_id+1):(current_max_id+n_offspring), "seroconvert"] <- offspring_seroconvert
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_seroconversion"] <- offspring_time_seroconvert
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_seroreversion"] <- offspring_time_serorevert
      tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
      tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
    }
    susc <- susc - n_offspring
  }
  tdf <- tdf[order(tdf$time_infection, tdf$infection_id), ]
  return(tdf)
}

#' Utility function to calculate shedding given infection time-series and shedding distribution

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


#' Convert stochastic realisation of branching process into number of new infections time-series
#'
#' @family post-processing
#' @export
generate_infections_time_series <- function(branching_process_output, population) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)

  # Create a dataframe with counts of infections per day
  infection_counts <- branching_process_output |>
    mutate(day = floor(time_infection)) |>
    group_by(day) |>
    summarise(new_infections = n(), .groups = 'drop') |>
    complete(day = seq(0, max_day, by = 1),
             fill = list(new_infections = 0)) |>
    mutate(new_infections_per_thousand = 1000 * new_infections / population)
  return(infection_counts)

}

#' Convert stochastic realisation of branching process into number shedding time-series
#'
#' This function converts a branching process output into time-series of number shedding,
#' weighted by the shedding distribution.
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

#' Convert stochastic realisation of branching process into symptom onset time-series
#'
#' This function converts a branching process output into time-series of symptom onsets
#'
#' @family post-processing
#' @import dplyr
#' @import tidyr
#' @export
generate_symptom_onset_time_series <- function(branching_process_output, population) {

  ## Calculating daily incidence of symptom onsets
  symptom_onset_incidence <- branching_process_output |>
    filter(symptomatic == 1) |>
    mutate(time_symptom_onset_floor = floor(time_symptom_onset)) |>
    group_by(time_symptom_onset_floor) |>
    summarise(incidence_symptom_onset = n()) |>
    complete(time_symptom_onset_floor = seq(0, max(time_symptom_onset_floor, na.rm = TRUE), by = 1),
             fill = list(incidence_symptom_onset = 0)) |>
    rename(day = time_symptom_onset_floor) %>%
    mutate(incidence_symptom_onset_per_thousand = 1000 * incidence_symptom_onset / population)
  return(symptom_onset_incidence)

}

#' Convert stochastic realisation of branching process into hospitalisation time-series
#'
#' This function converts a branching process output into time-series of hospitalisations
#' and hospital occupancy.
#'
#' @family post-processing
#' @export
generate_healthcare_seeking_time_series <- function(branching_process_output, population) {

  ## Calculating daily incidence of hospitalisations
  healthcare_seeking_incidence <- branching_process_output |>
    filter(seek_healthcare == 1) |>
    mutate(time_seek_healthcare_floor = floor(time_seek_healthcare)) |>
    group_by(time_seek_healthcare_floor) |>
    summarise(incidence_seek_healthcare = n()) |>
    complete(time_seek_healthcare_floor = seq(0, max(time_seek_healthcare_floor, na.rm = TRUE), by = 1),
             fill = list(incidence_seek_healthcare = 0)) |>
    rename(day = time_seek_healthcare_floor) %>%
    mutate(incidence_seek_healthcare_per_thousand = 1000 * incidence_seek_healthcare / population)

  ## OLD Calculating number of individuals in hospital on any given day
  # max_day <- ceiling(max(branching_process_output$time_leave_hospital, na.rm = TRUE))
  # days <- seq(0, max_day, by = 1)
  # hospital_occupancy <- sapply(days, function(t) {
  #   sum(branching_process_output$time_hospitalised < t & branching_process_output$time_leave_hospital > t, na.rm = TRUE)
  # })
  # hospital_occupancy <- data.frame(day = days, hospitalised = hospital_occupancy)
  ## Merging dataframes together
  # overall_hospital_df <- hospitalisation_incidence |>
  #   left_join(hospital_occupancy, by = "day") |>
  #   mutate(incidence_hospitalisation_per_thousand = 1000 * incidence_hospitalisation / population,
  #          hospitalised_per_thousand = 1000 * hospitalised / population)

  return(healthcare_seeking_incidence)
}

#' Convert stochastic realisation of branching process into hospitalisation time-series
#'
#' This function converts a branching process output into time-series of hospitalisations
#' and hospital occupancy.
#'
#' @family post-processing
#' @export
generate_seropositivity_timeseries <- function(branching_process_output, population) {

  ## Calculating daily incidence of hospitalisations
  seroconversion_df <- branching_process_output |>
    filter(seroconvert == 1) %>%
    mutate(time_seroconversion_floor = floor(time_seroconversion)) |>
    group_by(time_seroconversion_floor) %>%
    summarise(n_seroconverted = n()) %>%
    complete(time_seroconversion_floor = seq(0, max(time_seroconversion_floor, na.rm = TRUE), by = 1),
             fill = list(n_seroconverted = 0)) |>
    mutate(cumulative_seroconversions = cumsum(n_seroconverted))

  seroreversion_df <- branching_process_output |>
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

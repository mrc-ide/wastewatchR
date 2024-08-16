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
    prob_hosp = 0.35,
    prob_seek_healthcare = 0.9,
    hospitalisation_delay_dist = function(n) { rgamma(n, shape = 24, rate = 2) },
    hospitalisation_duration_dist = function(n) { rgamma(n, shape = 12, rate = 2) },

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
    hospitalised = integer(check_final_size),
    time_hospitalised = numeric(check_final_size),
    time_leave_hospital = numeric(check_final_size),
    n_offspring = integer(check_final_size),
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )

  # Initialising the model with the seeding cases
  seeding_cases_time_infection <- t0 + seq(from = 0, to = 0.01, length.out = seeding_cases) ## Time of infection for seeding cases

  ## Generating symptom status of seeding cases
  seeding_cases_symptomatic <- sample(c(0, 1), seeding_cases, replace = TRUE,               ## Determining whether seeding cases are symptomatic
                                      prob = c(1-prob_symptomatic, prob_symptomatic))
  number_seeding_cases_symptomatic <- sum(seeding_cases_symptomatic)                        ## Calculating number of symptomatic index cases
  index_seeding_cases_symptomatic <- which(seeding_cases_symptomatic == 1)                  ## Which seeding cases are symptomatic
  seeding_cases_time_symptom_onset <- rep(NA, seeding_cases)
  seeding_cases_time_symptom_onset[index_seeding_cases_symptomatic] <- seeding_cases_time_infection[index_seeding_cases_symptomatic] + infection_to_onset_dist(number_seeding_cases_symptomatic)

  ## Generating hospitalised status of seeding cases
  prob_hosp_if_symptomatic <- prob_hosp / prob_symptomatic  ## Probability of being hospitalised conditional on being symptomatic
  if (prob_hosp_if_symptomatic >= 1) {
    stop("prob_hosp / prob_symptomatic must be <= 1")
  }
  seeding_cases_hospitalised <- rep(0, seeding_cases)
  seeding_cases_hospitalised[index_seeding_cases_symptomatic] <- sample(c(0, 1), number_seeding_cases_symptomatic, replace = TRUE,               ## Determining whether seeding cases are symptomatic
                                                                        prob = c(1-(prob_hosp_if_symptomatic * prob_seek_healthcare),
                                                                                 (prob_hosp_if_symptomatic * prob_seek_healthcare)))
  number_seeding_cases_hospitalised <- sum(seeding_cases_hospitalised)                        ## Calculating number of hospitalised index cases
  index_seeding_cases_hospitalised <- which(seeding_cases_hospitalised == 1)                  ## Which seeding cases are hospitalised
  seeding_cases_time_hospitalised <- rep(NA, seeding_cases)
  seeding_cases_time_hospitalised[index_seeding_cases_hospitalised] <- seeding_cases_time_symptom_onset[index_seeding_cases_hospitalised] + hospitalisation_delay_dist(number_seeding_cases_hospitalised)
  seeding_cases_time_leave_hospital <- rep(NA, seeding_cases)
  seeding_cases_time_leave_hospital[index_seeding_cases_hospitalised] <- seeding_cases_time_hospitalised[index_seeding_cases_hospitalised] + hospitalisation_duration_dist(number_seeding_cases_hospitalised)

  ## Initialize the dataframe with the seeding cases
  tdf[1:seeding_cases, ] <- data.frame(
    infection_id = seq_len(seeding_cases),
    infection_ancestor = NA_integer_,
    infection_generation = 1L,
    time_infection = seeding_cases_time_infection,
    symptomatic = seeding_cases_symptomatic,
    time_symptom_onset = seeding_cases_time_symptom_onset,
    hospitalised = seeding_cases_hospitalised,
    time_hospitalised = seeding_cases_time_hospitalised,
    time_leave_hospital = seeding_cases_time_leave_hospital,
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
      offspring_symptomatic <- sample(c(0, 1), n_offspring, replace = TRUE, prob = c(1-prob_symptomatic, prob_symptomatic))
      number_offspring_symptomatic <- sum(offspring_symptomatic)
      index_offspring_symptomatic <- which(offspring_symptomatic == 1)
      offspring_time_symptom_onset <- rep(NA, n_offspring)
      offspring_time_symptom_onset[index_offspring_symptomatic] <- offspring_time_infection[index_offspring_symptomatic] + infection_to_onset_dist(number_offspring_symptomatic)

      ## Hospitalisation status/timing
      offspring_hospitalised <- rep(0, n_offspring)
      offspring_hospitalised[index_offspring_symptomatic] <- sample(c(0, 1), number_offspring_symptomatic, replace = TRUE,
                                                                    prob = c(1-(prob_hosp_if_symptomatic * prob_seek_healthcare),
                                                                             (prob_hosp_if_symptomatic * prob_seek_healthcare)))
      number_offspring_hospitalised <- sum(offspring_hospitalised)
      index_offspring_hospitalised <- which(offspring_hospitalised == 1)
      offspring_time_hospitalised <- rep(NA, n_offspring)
      offspring_time_hospitalised[index_offspring_hospitalised] <- offspring_time_symptom_onset[index_offspring_hospitalised] + hospitalisation_delay_dist(number_offspring_hospitalised)
      offspring_time_leave_hospital <- rep(NA, n_offspring)
      offspring_time_leave_hospital[index_offspring_hospitalised] <- offspring_time_hospitalised[index_offspring_hospitalised] + hospitalisation_duration_dist(number_offspring_hospitalised)

      ## Appending this information on the offspring to the table
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_id"] <- c(current_max_id + seq_len(n_offspring))
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_ancestor"] <- parent_idx
      tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_generation"] <- parent_gen + 1L
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_infection"] <- offspring_time_infection
      tdf[(current_max_id+1):(current_max_id+n_offspring), "symptomatic"] <- offspring_symptomatic
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_symptom_onset"] <- offspring_time_symptom_onset
      tdf[(current_max_id+1):(current_max_id+n_offspring), "hospitalised"] <- offspring_hospitalised
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_hospitalised"] <- offspring_time_hospitalised
      tdf[(current_max_id+1):(current_max_id+n_offspring), "time_leave_hospital"] <- offspring_time_leave_hospital
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
generate_hospitalisation_time_series <- function(branching_process_output, population) {

  ## Calculating daily incidence of hospitalisations
  hospitalisation_incidence <- branching_process_output |>
    filter(hospitalised == 1) |>
    mutate(time_hospitalised_floor = floor(time_hospitalised)) |>
    group_by(time_hospitalised_floor) |>
    summarise(incidence_hospitalisation = n()) |>
    complete(time_hospitalised_floor = seq(0, max(time_hospitalised_floor, na.rm = TRUE), by = 1),
             fill = list(incidence_hospitalisation = 0)) |>
    rename(day = time_hospitalised_floor)

  ## Calculating number of individuals in hospital on any given day
  max_day <- ceiling(max(branching_process_output$time_leave_hospital, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)
  hospital_occupancy <- sapply(days, function(t) {
    sum(branching_process_output$time_hospitalised < t & branching_process_output$time_leave_hospital > t, na.rm = TRUE)
  })
  hospital_occupancy <- data.frame(day = days, hospitalised = hospital_occupancy)

  ## Merging dataframes together
  overall_hospital_df <- hospitalisation_incidence |>
    left_join(hospital_occupancy, by = "day") |>
    mutate(incidence_hospitalisation_per_thousand = 1000 * incidence_hospitalisation / population,
           hospitalised_per_thousand = 1000 * hospitalised / population)
  return(overall_hospital_df)
}

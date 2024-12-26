#' Generate a stochastic realisation from a stochastic branching process
#'
#' This function creates a named list containing the output of simulating a stochastic branching
#' process.
#' @param initial_mn_offspring The initial value for the mean of the offspring distribution (inital as it can be modified by evolutionary forces, if these are switched on)
#' @param disp_offspring The overdisperion of the offspring distribution. Must be >= 1. When set to 1, equivalent to a Poisson distribution; >1 is a Negative Binomial distribution.
#' @param generation_time_dist The generation time distribution
#' @param prob_symptomatic The probability that an infection develops symptoms
#' @param infection_to_onset_dist The infection to symptom onset distribution
#' @param prob_severe The probability that an infection develops severe disease. Note this is NOT conditional on being symptomatic (this conditional probability is calculated internally within the model)
#' @param prob_seek_healthcare_non_severe The probability that a non-severe infection seeks healthcare and is tested and diagnosed successfully.
#' @param prob_seek_healthcare_severe The probability that a severe infection seeks healthcare and is tested and diagnosed successfully.
#' @param onset_to_healthcare_dist The sympmtom onset to seeking healthcare delay distribution.
#' @param prob_seroconvert_asymptomatic The probability that an asymptomatic infection seroconverts.
#' @param prob_seroconvert_severe The probability that an infection with severe disease seroconverts.
#' @param prob_seroconvert_non_severe The probability that an infection with non-severe disease seroconverts.
#' @param infection_to_seroconversion_dist The infection to seroconversion delay distribution.
#' @param seroconversion_to_seroreversion_dist The seroconversion to seroreversion delay distribution.
#' @param prob_beneficial_mutation The probability that any given infection acquires a mutation that modifies the mean of their offspring distribution (i.e. increases transmissibility).
#' @param beneficial_mutation_effect_dist The distribution of effect sizes for each beneficial mutation (i.e. how much they increase transmissibility by).
#' @param max_mn_offspring An upper cap on how high the mean of the offspring distribution can be.
#' @param annual_spillover_rate Rate of zoonotic spillover per year
#' @param spillover_seeding_cases_dist The distribution of initial seeding cases associated with each zoonotic spillover.
#' @param t0 The time to start the simulation from.
#' @param tf The time to finish the simulation at. Not functionally used much as usually we hit check_final_size or max_num_outbreaks first.
#' @param population Size of the overall population.
#' @param initial_immune The number of individuals initially immune.
#' @param check_final_size The maximum number of infected individuals that will be simulated.
#' @param max_num_outbreaks The maximum number of distinct outbreaks that will be simulated.
#' @param seed The stochastic seed used in the simulation.
#'
#'
#' @family simulation
#' @export
simulate_branching_process <- function(

    ## Offspring distribution parameters
    initial_mn_offspring = 0.9,
    disp_offspring = 0.9,

    ## Natural history parameters
    generation_time_dist = function(n) { rgamma(n, shape = 12, rate = 2) }, ## poss need to reverse rate and shape - TBD

    ## Disease severity parameters
    prob_symptomatic = 0.8,
    infection_to_onset_dist = function(n) { rgamma(n, shape = 6, rate = 2) }, ## poss need to reverse rate and shape - TBD
    prob_severe = 0.35,
    prob_seek_healthcare_non_severe = 0.9,
    prob_seek_healthcare_severe = 0.9,
    onset_to_healthcare_dist = function(n) { rgamma(n, shape = 6, rate = 2) }, ## poss need to reverse rate and shape - TBD

    ## Serology related parameters
    prob_seroconvert_asymptomatic = 0.8,
    prob_seroconvert_severe = 0.8,
    prob_seroconvert_non_severe = 0.8,
    infection_to_seroconversion_dist = function(n) { rgamma(n, shape = 56, rate = 2) }, ## poss need to reverse rate and shape - TBD
    seroconversion_to_seroreversion_dist = function(n) { rgamma(n, shape = 730, rate = 2) }, ## poss need to reverse rate and shape - TBD

    ## Evolution related parameters
    prob_beneficial_mutation = 0.75,
    beneficial_mutation_effect_dist = function(n) {rexp(n, rate = 10) },
    max_mn_offspring = 4,

    ## Simulation parameters
    annual_spillover_rate = 2,
    spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
    t0 = 0,
    tf = Inf,
    population = 100000,
    initial_immune = 0,
    check_final_size = 10000,
    max_num_outbreaks = 10,
    seed = 98) {

  ## Setting the seed
  set.seed(seed)

  ## Setting up the Offspring Distribution
  if (disp_offspring <= 1.0) {
    offspring_fun <- function(n, susc, mn_offspring) {
      new_mn <- mn_offspring * susc/population
      rpois(n = n, lambda = new_mn)
    }
  } else {
    offspring_fun <- function(n, susc, mn_offspring) {
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
    outbreak = integer(check_final_size),
    mutated = integer(check_final_size),
    effect_size = numeric(check_final_size),
    infection_mn_offspring = numeric(check_final_size),
    seeding = NA_character_,
    offspring_generated = FALSE,
    stringsAsFactors = FALSE
  )
  susc <- population - initial_immune - 1L ## do we need to be account for seroreversion/waning immunity for this - possibly...
  outbreak_index <- 1

  while (outbreak_index <= max_num_outbreaks & nrow(tdf) <= check_final_size & sum(is.na(tdf$time_infection)) != 0) {

    # Getting the first unfilled spot in the overall dataframe
    tdf_start_index <- min(which(is.na(tdf$time_infection)))

    # Generating seeding cases for the initial outbreak
    seeding_cases <- spillover_seeding_cases_dist()
    if (seeding_cases == 0) {
      seeding_cases <- 1
    }
    time_spillover <- ifelse(outbreak_index == 1, 0, min(tdf$time_infection[tdf$outbreak == (outbreak_index - 1)]) + rexp(1, rate = annual_spillover_rate / 365))

    # Initialising the model with the seeding cases
    seeding_cases_time_infection <- time_spillover + seq(from = 0, to = 0.01, length.out = seeding_cases) ## Time of infection for seeding cases

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
    current_max_id <- max(tdf$infection_id)
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "infection_id"] <- (current_max_id + 1):(current_max_id + seeding_cases)
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "infection_ancestor"] <- NA_integer_
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "infection_generation"] <- 1L
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "time_infection"] <- seeding_cases_time_infection
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "symptomatic"] <- seeding_cases_symptomatic
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "time_symptom_onset"] <- seeding_cases_time_symptom_onset
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "severe"] <- seeding_cases_severe
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "seek_healthcare"] <- seeding_case_seek_healthcare
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "time_seek_healthcare"] <- seeding_case_time_seek_healthcare
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "seroconvert"] <- seeding_case_seroconvert
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "time_seroconversion"] <- seeding_case_time_seroconvert
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "time_seroreversion"] <- seeding_case_time_serorevert
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "n_offspring"] <- NA_integer_
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "outbreak"] <- outbreak_index
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "mutated"] <- 0
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "effect_size"] <- 0
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "infection_mn_offspring"] <- initial_mn_offspring
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "seeding"] <- "seeding"
    tdf[tdf_start_index:(tdf_start_index + seeding_cases - 1), "offspring_generated"] <- FALSE

    ## This loops through each infection individually and generates offspring, and the status (symptoms, hospitalisation etc) of all of these offspring
    ## Note that we have check_final_size to cap how long we simulate for, but that when you reach this threshold you might have infections for which
    ## you haven't generated offspring for.
    while (any(tdf$offspring_generated[tdf$outbreak == outbreak_index] == FALSE) & nrow(tdf) <= check_final_size) {

      # Selecting the earliest infection for which we haven't generated offspring
      time_infection_index <- min(tdf$time_infection[tdf$offspring_generated == FALSE & !is.na(tdf$time_infection) & tdf$outbreak == outbreak_index])

      # Extracting information on the "parent" infection whose offspring are being generated
      parent_idx <- which(tdf$time_infection == time_infection_index & !tdf$offspring_generated)[1] # get the id of the earliest unsimulated infection
      parent_outbreak <- tdf$outbreak[parent_idx]                                                   # which outbreak does the parent belong to
      parent_time_infection <- tdf$time_infection[parent_idx]                                       # infection time of the earliest unsimulated infection
      parent_gen <- tdf$infection_generation[parent_idx]                                            # generation of the earliest unsimulated infection
      parent_symptomatic <- tdf$symptomatic[parent_idx]                                             # whether the earliest unsimulated infection is symptomatic
      parent_time_symptom_onset <- tdf$time_symptom_onset[parent_idx]                               # if symptomatic, the time when the earliest unsimulated infection develops symptoms
      parent_mn_offspring <- tdf$infection_mn_offspring[parent_idx]

      # Calculate total number of infections in the dataframe currently (so we can figure out how to label the new infections)
      current_max_id <- max(tdf$infection_id)

      # Calculating number of offspring generated
      n_offspring <- offspring_fun(1, susc, parent_mn_offspring)
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

        ## Evolution of transmissibility
        mn_offspring_increase <- sample(c(0, 1), n_offspring, replace = TRUE, prob = c(1 - prob_beneficial_mutation, prob_beneficial_mutation))
        index_mn_offspring_increase <- which(mn_offspring_increase == 1)
        number_mn_offspring_increase <- sum(mn_offspring_increase)
        mn_offspring_increase_effect_size <- rep(0, n_offspring)
        mn_offspring_increase_effect_size[index_mn_offspring_increase] <- beneficial_mutation_effect_dist(number_mn_offspring_increase)
        mn_offspring_increase_effect_size[mn_offspring_increase >= max_mn_offspring] <- max_mn_offspring

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
        tdf[(current_max_id+1):(current_max_id+n_offspring), "outbreak"] <- parent_outbreak
        tdf[(current_max_id+1):(current_max_id+n_offspring), "n_offspring"] <- NA
        tdf[(current_max_id+1):(current_max_id+n_offspring), "outbreak"] <- outbreak_index
        tdf[(current_max_id+1):(current_max_id+n_offspring), "mutated"] <- mn_offspring_increase
        tdf[(current_max_id+1):(current_max_id+n_offspring), "effect_size"] <- mn_offspring_increase_effect_size
        tdf[(current_max_id+1):(current_max_id+n_offspring), "infection_mn_offspring"] <- parent_mn_offspring + mn_offspring_increase_effect_size
        tdf[(current_max_id+1):(current_max_id+n_offspring), "seeding"] <- "onwards_infections"
        tdf[(current_max_id+1):(current_max_id+n_offspring), "offspring_generated"] <- FALSE
      }
      susc <- susc - n_offspring
    }
    outbreak_index <- outbreak_index + 1

  }
  tdf <- tdf[order(tdf$time_infection, tdf$infection_id), ]
  return(tdf)
}

#' Convert stochastic realisation of branching process into number of new infections time-series
#'
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @export
generate_infections_time_series <- function(branching_process_output) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)

  # Create a dataframe with counts of infections per day
  infection_counts <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
    mutate(day = floor(time_infection)) |>
    group_by(day) |>
    summarise(new_infections = n(), .groups = 'drop') |>
    complete(day = seq(0, max_day, by = 1),
             fill = list(new_infections = 0))
  return(infection_counts)
}

#' Convert stochastic realisation of branching process into number of new infections time-series
#'
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @export
generate_outbreak_size <- function(branching_process_output, population) {

  max_day <- ceiling(max(branching_process_output$time_infection, na.rm = TRUE))
  days <- seq(0, max_day, by = 1)

  ## Calculate outbreak time
  outbreak_times <- branching_process_output |>
    filter(outbreak != 0) %>%
    group_by(outbreak) %>%
    summarise(outbreak_start = min(time_infection, na.rm = TRUE),
              outbreak_end = max(c(time_infection, time_symptom_onset, time_seek_healthcare), na.rm = TRUE))

  ## For infections
  infection_counts_per_outbreak_including_seeding <- branching_process_output |>
    filter(outbreak != 0) %>%
    group_by(outbreak) %>%
    summarise(total_infection_size_with_seeding = n())

  infection_counts_per_outbreak_excluding_seeding <- branching_process_output |>
    filter(seeding != "seeding") %>%
    filter(outbreak != 0) %>%
    group_by(outbreak) %>%
    summarise(total_infection_size_without_seeding = n())

  infection_overall_outbreak_size_df <- infection_counts_per_outbreak_including_seeding %>%
    left_join(infection_counts_per_outbreak_excluding_seeding, by = "outbreak")

  ## For symptom onsets
  symptom_onset_counts_per_outbreak_including_seeding <- branching_process_output |>
    filter(outbreak != 0) %>%
    filter(symptomatic == 1) %>%
    group_by(outbreak) %>%
    summarise(total_onset_size_with_seeding = n())

  symptom_onset_counts_per_outbreak_excluding_seeding <- branching_process_output |>
    filter(symptomatic == 1) %>%
    filter(outbreak != 0) %>%
    filter(seeding != "seeding") %>%
    group_by(outbreak) %>%
    summarise(total_onset_size_without_seeding = n())

  onsets_overall_outbreak_size_df <- symptom_onset_counts_per_outbreak_including_seeding %>%
    left_join(symptom_onset_counts_per_outbreak_excluding_seeding, by = "outbreak")

  ## For healthcare seeking
  healthcare_counts_per_outbreak_including_seeding <- branching_process_output |>
    filter(seek_healthcare == 1) |>
    filter(outbreak != 0) %>%
    group_by(outbreak) %>%
    summarise(total_healthcare_seek_size_with_seeding = n())

  healthcare_counts_per_outbreak_excluding_seeding <- branching_process_output |>
    filter(seek_healthcare == 1) |>
    filter(seeding != "seeding") %>%
    filter(outbreak != 0) %>%
    group_by(outbreak) %>%
    summarise(total_healthcare_seek_size_without_seeding = n())

  healthcare_seek_overall_outbreak_size_df <- healthcare_counts_per_outbreak_including_seeding %>%
    left_join(healthcare_counts_per_outbreak_excluding_seeding, by = "outbreak")

  overall_df <- infection_overall_outbreak_size_df %>%
    left_join(onsets_overall_outbreak_size_df, by = "outbreak") %>%
    left_join(healthcare_seek_overall_outbreak_size_df, by = "outbreak") %>%
    left_join(outbreak_times, by = "outbreak")

  return(overall_df)
}

#' Convert stochastic realisation of branching process into symptom onset time-series
#' This function converts a branching process output into time-series of symptom onsets
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#' @family post-processing
#' @import dplyr
#' @import tidyr
#' @export
generate_symptom_onset_time_series <- function(branching_process_output, population) {

  ## NOTE - NEED TO CHANGE THIS TO MATCH THE NEW STRUCTURE OF generate_healthcare_seeking_time_series
  ## Calculating daily incidence of symptom onsets
  symptom_onset_incidence <- branching_process_output |>
    filter(!is.na(time_infection)) %>%
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
#' @param branching_process_output Output from simulate_branching_process
#' @param population The population size used in the simulate_branching_process output
#'
#'
#' @family post-processing
#' @export
generate_healthcare_seeking_time_series <- function(branching_process_output, population) {

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

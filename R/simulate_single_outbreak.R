#' Generate a realisation of a single outbreak from a stochastic
#' branching process
#'
#' These functions specify a minimal model for simulating a single outbreak
#' of an infectious disease following a spillover event.
#'
#' This version of the outbreak simulater is designed to run quickly and simply
#' and cannot accommodate dependence of infectees characteristics on infectors
#' characteristics (e.g. mutations --> changes in R0, or decreases in R0 based on
#' decreases in the proportion of the population susceptible over time)
#'
#' This function simulates a branching process
#'
#' @param mn_offspring The mean of the offspring distribution (R0)
#' @param disp_offspring The overdispersion of the offspring distribution. Must be >= 1. When set to 1, equivalent to a Poisson distribution; >1 is a Negative Binomial distribution.
#' @param max_gen The maximum number of generations of transmission to simulate - default is Inf but with mn_offspring <1 transmission dies out eventually.
#' @param index_cases The number of initial seeding cases associated with the zoonotic spillover.
#' @param initial_immune The proportion of the population initially immune.
#'
#' @family simulation
#' @export
sim_minimal <- function(mn_offspring = 0.90,
                        disp_offspring = 1,
                        max_gen = Inf,
                        index_cases = 1,
                        initial_immune = 0){

  Z <- list()
  Z[[1]] <- index_cases
  i <- 1

  if (disp_offspring <= 1.0) {

    while(sum(Z[[i]]) > 0 && i <= max_gen) {
      Z[[i+1]] <- rpois(n = sum(Z[[i]]),
                        lambda = mn_offspring * (1 - initial_immune))
      i <- i+1
    }

  } else {

    while(sum(Z[[i]]) > 0 && i <= max_gen) {

      Z[[i+1]] <- rnbinom(n = sum(Z[[i]]),
                          size =  (1 - initial_immune) *
                            mn_offspring/(disp_offspring - 1),
                          mu = mn_offspring)
      i <- i+1

    }
  }

  return(Z)

}

#' This function takes the branching process output of sim_minimal and creates
#' a linelist consisting of a dataframe with 1 row per infected individual with
#' assigned clinical characteristics based on specified probability distributions.
#'
#' @param spillover_day time that spillover occurred (taken from output of spillover function)
#' @param index_case_ID ID of the index case (the spillover)
#' @param mn_offspring The mean of the offspring distribution (R0)
#' @param disp_offspring The overdisperion of the offspring distribution. Must be >= 1. When set to 1, equivalent to a Poisson distribution; >1 is a Negative Binomial distribution.
#' @param max_gen The maximum number of generations of transmission to simulate - default is Inf but with mn_offspring <1 transmission dies out eventually.
#' @param index_cases The number of initial seeding cases associated with the zoonotic spillover.
#' @param generation_time_dist The generation time distribution
#' @param prob_symptomatic The probability that an infected individual develops symptoms
#' @param infection_to_onset_dist The infection to symptom onset distribution
#' @param prob_severe The probability that an infection develops severe disease. Note this IS CONDITIONAL on being symptomatic.
#' @param prob_seek_healthcare_non_severe The probability that a non-severe infection seeks healthcare.
#' @param prob_seek_healthcare_severe The probability that a severe infection seeks healthcare.
#' @param onset_to_healthcare_dist The symptom onset to seeking healthcare delay distribution.
#' @param prob_diagnosis The probability of diagnosis given healthcare sought. Note this is not modeled as being dependent on severity, beyond the difference in probability of seeking care.
#' @param healthcare_to_diagnosis_dist The distribution of time from healthcare seeking to diagnosis
#' @param initial_immune The proportion of the population initially immune.
#'
#' @family simulation
#' @export

  sim_single_outbreak <- function(spillover_day,
                                  index_case_ID,
                                  mn_offspring = 0.90,
                     disp_offspring = 1,
                     max_gen = Inf,
                     index_cases = 1,
                     generation_time_dist = function(n)
                       { rgamma(n, shape = 12, rate = 2) },
                     prob_symptomatic = 0.80,
                     infection_to_onset_dist = function(n)
                       { rgamma(n, shape = 6, rate = 2) },
                     prob_severe = 0.30,
                     prob_seek_healthcare_non_severe = 0.50,
                     prob_seek_healthcare_severe = 0.95,
                     onset_to_healthcare_dist = function(n)
                       { rgamma(n, shape = 6, rate = 2) },
                     prob_diagnosis = 0.8,
                     healthcare_to_diagnosis_dist = function(n)
                       { rgamma(n, shape = 6, rate = 2) },
                     initial_immune = 0,
                      ...){

  #-----------------------------------------------------------------------------
  # simulate branching process -------------------------------------------------

    bp <- sim_minimal(mn_offspring = mn_offspring,
                      disp_offspring = disp_offspring,
                      max_gen = max_gen,
                      index_cases = index_cases,
                      initial_immune = initial_immune)


  #-----------------------------------------------------------------------------
  # format output of sim_minimal into dataframe w/ 1 row per infected individual

    tmp <- melt(bp)%>%
    mutate(number = 1, # number for counting generation size
           infection_generation = L1-1)%>% # index case gen 1
    filter(infection_generation>=1)%>%
    rename(n_offspring = value)%>%
    group_by(infection_generation)%>%
    mutate(node = cumsum(number))%>% # count generation size
    ungroup()%>%
    select(-number, -L1)%>%
    mutate(id = paste0(infection_generation, "-", node))

  ## look up infectors
  tmp$infector <- "animal"

  if(dim(tmp)[1]>=index_cases+1){
    tmp$infector[(index_cases + 1):dim(tmp)[1]] <- tmp%>% uncount(n_offspring) %>%
      pull(id)}

  ## order more intuitively
  tmp <- tmp %>%
    relocate(id, infection_generation, n_offspring, infector) %>%
    select(-node)

  #-----------------------------------------------------------------------------
  # assign characteristics (symptoms & healthcare seeking) ---------------------

  tmp <- tmp %>%
    mutate(time_inf_rel = if_else(infector == "animal", 0,
                                  generation_time_dist(nrow(.))),
           symptomatic = sample(c(0, 1), nrow(.),
                                replace = TRUE,
                                prob = c(1 - prob_symptomatic,
                                         prob_symptomatic))) %>%
    mutate(infection_to_onset = if_else(symptomatic == 0,
                                        NA, infection_to_onset_dist(nrow(.))),
           severe = if_else(symptomatic == 1,
                            sample(c(0, 1), nrow(.),
                                   replace = TRUE,
                                   prob = c(1-prob_severe, prob_severe)), 0),
           seek_healthcare = case_when((symptomatic == 1 & severe == 1) ~
                                         sample(c(0, 1), nrow(.),replace = TRUE,
                                                prob = c(1 - prob_seek_healthcare_severe,
                                                         prob_seek_healthcare_severe)),
                                       (symptomatic == 1 & severe == 0) ~
                                         sample(c(0, 1), nrow(.),
                                                replace = TRUE,
                                                prob = c(1 - prob_seek_healthcare_non_severe,
                                                         prob_seek_healthcare_non_severe)),
                                       symptomatic == 0 ~ 0),
           onset_to_healthcare = if_else(seek_healthcare == 1,
                                         onset_to_healthcare_dist(nrow(.)), NA),
           diagnosis = if_else(seek_healthcare == 1,
                               sample(c(0,1), nrow(.),
                                      replace = TRUE,
                                      prob = c(1-prob_diagnosis,
                                               prob_diagnosis)), 0),
           healthcare_to_diagnosis = if_else(diagnosis == 1,
                                             healthcare_to_diagnosis_dist(nrow(.)), NA))


  #-----------------------------------------------------------------------------
  # anchor to time. Note the day of the spillover that initiated the outbreak
  # replaces time = 0. If using this as a stand-alone function outside of the
  # wrapper, spillover_day can be set as zero or loaded from column "t" of the
  # dataframe output of the "spillover" function.

  tmp$time_infection <- spillover_day

  ## time_infection = time_infection of infector + generation time

  if(dim(tmp)[1]>index_cases+1){
    infectors <- unique(tmp$infector)
    tmp2 <- vector(mode = "list", length = length(infectors))
    tmp2[[1]] <- tmp %>% filter(infector == "animal")
    tmp3 <- bind_rows(tmp2)

    for(i in 2:(dim(tmp)[1])){
      tmp2[[i]] <- tmp %>% filter(infector == infectors[i]) %>%
        mutate(time_infection = tmp3 %>%
                 filter(id == infectors[i]) %>%
                 pull(time_infection) + time_inf_rel)
      tmp3 <- bind_rows(tmp2)
    }
    tmp <- tmp3
  }


  linelist <- tmp %>%
    mutate(time_symptom_onset = time_infection + infection_to_onset) %>%
    mutate(time_seek_healthcare = time_symptom_onset + onset_to_healthcare)%>%
    mutate(time_diagnosis = time_seek_healthcare + healthcare_to_diagnosis)%>%
    mutate(spillover_ID = index_case_ID)

  return(linelist) }



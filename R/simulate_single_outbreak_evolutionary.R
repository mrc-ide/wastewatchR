#' Simulate one outbreak following a spillover
#'
#' Simulates onward transmission from one spillover event using an individual-level
#' branching process. Offspring means are adjusted for the susceptible population,
#' and beneficial mutations can increase the transmissibility inherited by descendants.
#' Clinical and serological characteristics are assigned to every infection.
#'
#' @param spillover_day Day on which spillover occurred.
#' @param index_case_ID Identifier for the spillover event, stored as `spillover_ID`.
#' @param initial_mn_offspring Initial mean offspring number (R0 before susceptible adjustment).
#' @param disp_offspring Offspring overdispersion. Values <= 1 use Poisson; values > 1 use negative binomial.
#' @param generation_time_dist Function taking `n` and returning `n` generation times.
#' @param prob_symptomatic Probability an infection is symptomatic.
#' @param infection_to_onset_dist Function returning infection-to-onset delays.
#' @param prob_severe Probability of severe disease conditional on being symptomatic.
#' @param prob_seek_healthcare_non_severe Probability a symptomatic non-severe case seeks healthcare.
#' @param prob_seek_healthcare_severe Probability a severe case seeks healthcare.
#' @param onset_to_healthcare_dist Function returning onset-to-healthcare delays.
#' @param prob_seroconvert_asymptomatic Seroconversion probability for asymptomatic infections.
#' @param prob_seroconvert_severe Seroconversion probability for severe infections.
#' @param prob_seroconvert_non_severe Seroconversion probability for symptomatic non-severe infections.
#' @param infection_to_seroconversion_dist Function returning infection-to-seroconversion delays.
#' @param seroconversion_to_seroreversion_dist Function returning seroconversion-to-seroreversion durations.
#' @param prob_beneficial_mutation Probability an infection acquires a beneficial mutation.
#' @param beneficial_mutation_effect_dist Function returning additive offspring-mean effects for mutated infections.
#' @param max_mn_offspring Maximum inherited offspring mean after mutation.
#' @param population Total population size used for susceptible adjustment.
#' @param initial_immune Number of people immune before this outbreak.
#' @param check_final_size Maximum number of infections to simulate, including seed cases.
#' @param seed Optional random seed. If `NULL`, the current R random stream is used.
#' @param prob_diagnosis Probability of diagnosis conditional on seeking healthcare.
#' @param healthcare_to_diagnosis_dist Function returning delays from healthcare seeking to diagnosis.
#'
#' @return A data frame with one row per simulated infection. It contains infection and
#' ancestor IDs, infection generation and times, clinical and serological outcomes,
#' offspring counts, mutation effects and inherited offspring means, plus `spillover_ID`.
#' Healthcare-seeking infections include diagnosis status and diagnosis time when diagnosed.
#' If `check_final_size` is reached, infections already created are returned even if
#' their offspring have not yet been generated.
#'
#' @family simulation
#' @export
simulate_single_outbreak_evolutionary <- function(
  spillover_day,
  index_case_ID,
  initial_mn_offspring = 0.9,
  disp_offspring = 1,
  generation_time_dist = function(n) rgamma(n, shape = 12, rate = 2),
  prob_symptomatic = 0.8,
  infection_to_onset_dist = function(n) rgamma(n, shape = 6, rate = 2),
  prob_severe = 0.35,
  prob_seek_healthcare_non_severe = 0.9,
  prob_seek_healthcare_severe = 0.9,
  onset_to_healthcare_dist = function(n) rgamma(n, shape = 6, rate = 2),
  prob_seroconvert_asymptomatic = 0.8,
  prob_seroconvert_severe = 0.8,
  prob_seroconvert_non_severe = 0.8,
  infection_to_seroconversion_dist = function(n) rgamma(n, shape = 56, rate = 2),
  seroconversion_to_seroreversion_dist = function(n) rgamma(n, shape = 730, rate = 2),
  prob_beneficial_mutation = 0.75,
  beneficial_mutation_effect_dist = function(n) rexp(n, rate = 10),
  max_mn_offspring = 4,
  population = 100000,
  initial_immune = 0,
  check_final_size = 10000,
  seed = NULL,
  prob_diagnosis = 0.8,
  healthcare_to_diagnosis_dist = function(n) rgamma(n, shape = 6, rate = 2)) {

  probabilities <- c(prob_symptomatic, prob_severe, prob_seek_healthcare_non_severe,
    prob_seek_healthcare_severe, prob_diagnosis, prob_seroconvert_asymptomatic,
    prob_seroconvert_severe, prob_seroconvert_non_severe, prob_beneficial_mutation)
  if (any(!is.finite(probabilities) | probabilities < 0 | probabilities > 1)) {
    stop("All probability arguments must be between 0 and 1.")
  }
  if (length(check_final_size) != 1L || is.na(check_final_size) ||
      check_final_size < 1) stop("check_final_size must be at least 1.")
  if (!is.finite(population) || population <= 0 || initial_immune < 0 || initial_immune >= population) {
    stop("population must be positive and initial_immune must be in [0, population).")
  }
  if (!is.finite(initial_mn_offspring) || initial_mn_offspring < 0 ||
      !is.finite(max_mn_offspring) || max_mn_offspring < 0) stop("Offspring means must be non-negative.")
  if (!is.finite(disp_offspring) || disp_offspring < 1) stop("disp_offspring must be at least 1.")
  if (!is.null(seed)) set.seed(seed)

  draw_bernoulli <- function(n, p) if (n == 0L) integer() else rbinom(n, 1, p)
  draw_offspring <- function(mean_offspring, susceptible) {
    mu <- max(0, mean_offspring * susceptible / population)
    if (mu == 0) return(0L)
    if (disp_offspring <= 1) rpois(1, mu) else rnbinom(1, size = mu / (disp_offspring - 1), mu = mu)
  }
  draw_delays <- function(fun, n) {
    if (!n) return(numeric())
    values <- fun(n)
    if (length(values) != n || any(!is.finite(values)) || any(values < 0)) {
      stop("Delay and effect distribution functions must return n finite, non-negative values.")
    }
    values
  }
  make_cases <- function(ids, ancestors, generations, times, inherited_mean, seeding) {
    n <- length(ids)
    symptomatic <- draw_bernoulli(n, prob_symptomatic)
    onset <- rep(NA_real_, n)
    si <- which(symptomatic == 1L)
    onset[si] <- times[si] + draw_delays(infection_to_onset_dist, length(si))
    severe <- integer(n)
    if (length(si)) severe[si] <- draw_bernoulli(length(si), prob_severe)
    seek <- integer(n)
    sev <- which(severe == 1L)
    nonsev <- which(symptomatic == 1L & severe == 0L)
    seek[sev] <- draw_bernoulli(length(sev), prob_seek_healthcare_severe)
    seek[nonsev] <- draw_bernoulli(length(nonsev), prob_seek_healthcare_non_severe)
    time_seek <- rep(NA_real_, n)
    hi <- which(seek == 1L)
    time_seek[hi] <- onset[hi] + draw_delays(onset_to_healthcare_dist, length(hi))
    diagnosis <- integer(n)
    diagnosis[hi] <- draw_bernoulli(length(hi), prob_diagnosis)
    time_diagnosis <- rep(NA_real_, n)
    di <- which(diagnosis == 1L)
    time_diagnosis[di] <- time_seek[di] +
      draw_delays(healthcare_to_diagnosis_dist, length(di))

    sero_p <- ifelse(symptomatic == 0L, prob_seroconvert_asymptomatic,
      ifelse(severe == 1L, prob_seroconvert_severe, prob_seroconvert_non_severe))
    seroconvert <- if (n) rbinom(n, 1, sero_p) else integer()
    time_seroconversion <- time_seroreversion <- rep(NA_real_, n)
    sci <- which(seroconvert == 1L)
    time_seroconversion[sci] <- times[sci] + draw_delays(infection_to_seroconversion_dist, length(sci))
    time_seroreversion[sci] <- time_seroconversion[sci] + draw_delays(seroconversion_to_seroreversion_dist, length(sci))

    data.frame(infection_id = ids, infection_ancestor = ancestors,
      infection_generation = generations, time_infection = times,
      symptomatic = symptomatic, time_symptom_onset = onset,
      severe = severe, seek_healthcare = seek, time_seek_healthcare = time_seek,
      diagnosis = diagnosis, time_diagnosis = time_diagnosis,
      seroconvert = seroconvert, time_seroconversion = time_seroconversion,
      time_seroreversion = time_seroreversion, n_offspring = NA_integer_,
      outbreak = 1L, mutated = 0L, effect_size = 0,
      infection_mn_offspring = inherited_mean, seeding = seeding,
      offspring_generated = FALSE, spillover_ID = index_case_ID,
      stringsAsFactors = FALSE)
  }

  cases <- make_cases(1L, NA_integer_, 1L, spillover_day,
    initial_mn_offspring, "seeding")
  susceptible <- max(0, population - initial_immune - 1L)
  next_id <- 2L
  parent_idx <- 1L

  while (parent_idx <= nrow(cases) && nrow(cases) < check_final_size) {
    parent_id <- cases$infection_id[parent_idx]
    parent_mean <- cases$infection_mn_offspring[parent_idx]
    n_offspring <- draw_offspring(parent_mean, susceptible)
    cases$n_offspring[parent_idx] <- n_offspring
    cases$offspring_generated[parent_idx] <- TRUE
    if (n_offspring > 0L) {
      n_add <- min(as.integer(n_offspring), as.integer(check_final_size - nrow(cases)))
      offspring_times <- cases$time_infection[parent_idx] + draw_delays(generation_time_dist, n_add)
      mutated <- draw_bernoulli(n_add, prob_beneficial_mutation)
      effect <- numeric(n_add)
      mi <- which(mutated == 1L)
      effect[mi] <- draw_delays(beneficial_mutation_effect_dist, length(mi))
      offspring_mean <- pmin(max_mn_offspring, parent_mean + effect)
      children <- make_cases(next_id:(next_id + n_add - 1L), rep(parent_id, n_add),
        rep(cases$infection_generation[parent_idx] + 1L, n_add), offspring_times,
        offspring_mean, rep("onwards_infections", n_add))
      children$mutated <- mutated
      children$effect_size <- effect
      cases <- rbind(cases, children)
      next_id <- next_id + n_add
      susceptible <- max(0, susceptible - n_add)
    }
    parent_idx <- parent_idx + 1L
  }

  cases <- cases[order(cases$time_infection, cases$infection_id), , drop = FALSE]
  rownames(cases) <- NULL
  cases
}

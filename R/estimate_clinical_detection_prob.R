#' Map infection-to-diagnosis delays to clinical detection probabilities
#'
#' Calculates the probability that an individual infected at a given lag
#' will be clinically diagnosed, conditional on the delay from infection to
#' diagnosis and the probabilities of being symptomatic, seeking healthcare,
#' and receiving a diagnosis.
#'
#' This is an internal helper function used by
#' [calculate_clinical_detection_prob()].
#'
#' @param lag Numeric vector of lags, in days, from infection to the
#'   potential day of diagnosis.
#' @param max_lag Maximum lag, in days, represented by `inf_to_diag`.
#' @param prob_symptomatic Probability that an infected individual
#'   develops symptoms.
#' @param prob_seek_healthcare Probability that a symptomatic
#'   individual seeks healthcare.
#' @param prob_diagnosis Probability that an individual seeking
#'   healthcare receives a diagnosis.
#' @param inf_to_diag Numeric vector giving the probability distribution of
#'   the delay from infection to diagnosis. The element at index `k + 1`
#'   corresponds to a delay of `k` days.
#'
#' @return A numeric vector containing the probability of clinical diagnosis
#'   for each value of `lag`. Values outside the range `0` to `max_lag` are
#'   assigned a probability of zero.
#'
#' @keywords internal
#'
get_q <- function(lag, max_lag, prob_symptomatic, prob_seek_healthcare, prob_diagnosis, inf_to_diag) {
  out <- numeric(length(lag))
  valid <- lag >= 0 & lag <= max_lag
  out[valid] <- prob_symptomatic * prob_seek_healthcare * prob_diagnosis * inf_to_diag[lag[valid] + 1]
  out
}
#' Calculate probability of clinical detection
#'
#' Calculates the probability of clinical detection over time based on
#' wastewater-derived infections and distributions describing delays from
#' infection to symptom onset, healthcare seeking, and diagnosis.
#'
#' The delay distributions are normalised before being convolved to obtain
#' the distribution of the total delay from infection to diagnosis. The
#' resulting probabilities are then used to calculate the probability of
#' observing at least one clinical detection for each day.
#'
#' @param nts Output of [generate_number_shedding_time_series()]
#' @param det Output of [calculate_wastewater_ttd()]
#' @param delay_onset Numeric vector giving the probability distribution of
#'   the delay from infection to symptom onset.
#' @param delay_seek Numeric vector giving the probability distribution of
#'   the delay from symptom onset to seeking healthcare.
#' @param delay_diag Numeric vector giving the probability distribution of
#'   the delay from seeking healthcare to diagnosis.
#' @param prob_symptomatic Numeric probability that an infected individual
#'   develops symptoms.
#' @param prob_seek_healthcare Probability that a symptomatic
#'   individual seeks healthcare.
#' @param prob_diagnosis Probability that an individual seeking
#'   healthcare is diagnosed.
#'
#' @return The object supplied to `det$sampled_data`,
#'   with an additional `prob_detection_clinical` element containing the
#'   probability of clinical detection for each day.
#'
#' @export
#'
calculate_clinical_detection_prob <- function(nts,
                                              det,
                                              delay_onset,
                                              delay_seek,
                                              delay_diag,
                                              prob_symptomatic,
                                              prob_seek_healthcare,
                                              prob_diagnosis){

  # normalise delay distributions since we truncate - tail should be captured but just in case
  delay_onset  <- delay_onset / sum(delay_onset)
  delay_seek  <- delay_seek / sum(delay_seek)
  delay_diag  <- delay_diag / sum(delay_diag)

  # convolution to get total delay inf_to_diag[k+1] = P(diagnosis k days after onset)
  ## (note "rev" is needed as stats::convolve actually does correlation)
  conv12 <- convolve(delay_onset, rev(delay_seek), type = "open")  # infection --> seek
  inf_to_diag <- convolve(conv12, rev(delay_diag), type = "open") # infection--> diagnosis
  inf_to_diag <- inf_to_diag / sum(inf_to_diag)  # normalise just incase
  max_lag <- length(inf_to_diag) - 1   # max lag (days) from infection --> diagnosis

  # map probability that those infected on days s (before t) are diagnosed on day t
  days <- nts$day
  n_days <- length(days)
  lag_mat <- outer(days, days, FUN = function(t, s) t - s)
  ## matrix of prob that a person infected on day s (cols) is diagnosed on day t (rows)
  q_mat <- matrix(get_q(as.vector(lag_mat), max_lag = max_lag,
                        prob_symptomatic, prob_seek_healthcare, prob_diagnosis, inf_to_diag),
                  nrow = n_days, ncol = n_days)

  # for each detection day t, compute prob of >=1 detection
  ## using log scale for stability
  N_s <- det$sampled_data$new_infections

  log_prob_none <- sapply(1:n_days, function(i) {
    q_vec <- q_mat[i, ]           # contribs for detection day i from each infection day s
    # if q_vec are zero for many entries this is fine
    sum(N_s * log1p(-q_vec))     # log1p(-q) is numerically stable for small q
  })

  ### attach to input
  det$sampled_data$prob_detection_clinical <- 1 - exp(log_prob_none)


  return(det)
}

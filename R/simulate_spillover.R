#' Generate spillover cases
#'
#' These functions define the daily probability of spillover, which can be time varying.
#' They then simulate the occurrence of spillovers over a time period of interest
#'
#' @param t The day (base = 1, min = 1, max = Inf)
#' @param tmax The day of the year that the variable peaks (min = 1, max = 365)
#' @param P The period of seasonality in days (default = 365)
#' @param b The baseline of the variable outside of the peak
#' @param d The seasonal forcing (amplitude increase rel to baseline, default = 0)
#' @param sigma The spread of the peak (sd of underlying normal distribution)

circular_dist <- function(t, tmax, P = 365) {

  # support function needed to make the Gaussian periodic
  # (to wrap back around each year)
  # Returns the shortest signed distance in range -P/2 to P/2
  dist <- (t - tmax + P/2) %% P - P/2
  return(dist)
}

#' Gaussian forcing function for spillover
#' @export
gaussian_forcing <- function(t, tmax, P = 365, b = 0, d = 0, sigma = 30) {

  # Where t is the day and tmax is the peak (between 0-P)
  dist <- circular_dist(t, tmax, P)

  # Gaussian seasonal forcing form
  f <- b * (1 + d * exp(- (dist)^2 / (2 * sigma^2)))
  return(f)
}

gaussian_mean <- function(tmax, P = 365, sigma = 30){

  t <- 1:P
  dist <- circular_dist(t, tmax, P)

  mean(exp(-(dist^2) / (2 * sigma^2)))
}

#' Calculate the baseline rate
#'
#' Calculates the baseline rate `b` required to achieve a specified
#' mean spillover rate given the seasonal forcing parameters.
#'
#' @param x Desired mean spillover rate.
#' @param d Seasonal forcing amplitude.
#' @param tmax Day of the year at which seasonal forcing peaks.
#' @param P Period of seasonality in days.
#' @param sigma Spread of the seasonal peak.
#'
#' @return The baseline rate `b`.
#'
#' @export
solve_b <- function(x, d, tmax, P = 365, sigma = 30){

  # this function gives you the value of b that corresponds to a desired rate
  # of spillover, for a given d and tmax

  gbar <- gaussian_mean(tmax, P, sigma)

  b <- x / (1 + d * gbar)
  b
}


#' Simulating spillovers over time
#' @param time The time period you would like to run the model over (days)
#' @param specify The parameters you would like to specify: either "spillover_rate" to specify spillover rate itself  OR "swiss_cheese" to specify factors contributing to spillover rate (prevalence in animals, contact rate between animals and humans, and probability of infection upon contact)
#' Used if specify == "spillover_rate"
#' @param spillover_rate_pars A list of: "tmax" - the day of the year (1-365) that spillover rate peaks, "b" - the baseline spillover rate, outside of a peak, "d" - the annual seasonal forcing (amplitude increase relative to baseline, default = 0), and "sigma" - the spread of the peak.
#' Used if specify == "swiss_cheese"
#' @param prevalence_pars A list of: "tmax" the day of the year (1-365) that prevalence peaks, "b" - the baseline prevalence of infection in the reservoir, outside of a peak, "d" - the annual seasonal forcing (amplitude increase relative to baseline, default = 0), and "sigma" - the spread of the peak.
#' @param contact_rate_pars A list of: "tmax" - the day of the year (1-365) that human-reservoir contact rate peaks, "b" - the baseline daily contact rate between people and the reservoir, "d" - the annual seasonal forcing (amplitude increase relative to baseline, default = 0), and "sigma" - the spread of the peak.
#' @param p_infection The probability of a person getting infected given contact with an infected animal. This is assumed to be constant.
#'
#' @family simulation
#' @export
#'

spillover <- function(time, specify, seasonal_period,
                      spillover_rate_pars, prevalence_pars, contact_rate_pars, p_infection) {

  # Check specify has been used correctly
  if(specify != "spillover_rate" & specify != "swiss_cheese"){
    stop("Error: 'specify' must be set to either 'spillover_rate' or 'swiss_cheese'")
  }

  t <- 1:time
  if(specify == "spillover_rate"){

      # Check that the necessary parameters exist
      if (!all(c("tmax", "b", "d", "sigma") %in% names(spillover_rate_pars))) {
        stop("For specify == 'spillover_rate', spillover_rate_pars must contain tmax, b, d and sigma")
      }

    s <- gaussian_forcing(t = t, tmax = spillover_rate_pars$tmax, P = seasonal_period,
                          b = spillover_rate_pars$b, d = spillover_rate_pars$d, sigma = spillover_rate_pars$sigma) # vectorised
    prevalence <- NA
    contact_rate <- NA
    p_infection <- NA

  } else if(specify == "swiss_cheese"){

    # Check that the necessary parameters exist
    if (!all(c("tmax", "b", "d", "sigma") %in% names(prevalence_pars)) &
        !all(c("tmax", "b", "d", "sigma") %in% names(contact_rate_pars))) {
      stop("For specify == 'swiss_cheese', prevalence_pars and contact_rate_pars must both contain tmax, b, d and sigma")
    }

    # Compute prevalence each day
    prevalence <- gaussian_forcing(t = t, tmax = prevalence_pars$tmax, P = seasonal_period,
                                   b = prevalence_pars$b, d = prevalence_pars$d, sigma = prevalence_pars$sigma) # vectorised
    # Compute contact rate each day
    contact_rate <- gaussian_forcing(t = t, tmax = contact_rate_pars$tmax, P = seasonal_period,
                                     b = contact_rate_pars$b, d = contact_rate_pars$d, sigma = contact_rate_pars$sigma) # vectorised
    # Compute subsequent daily spillover rate, s
    s <- prevalence * contact_rate * p_infection
  }

  # Simulate spillovers based on daily spillover rate
  spillovers <- rpois(time, s) # vectorised

  # Include parameters to be returned for sanity checks
  dat <- data.frame(t = t,
                    spillovers = spillovers, spillover_rate = s,
                    prevalence = prevalence,
                    contact_rate = contact_rate,
                    p_infection = p_infection)
  return(dat)
}

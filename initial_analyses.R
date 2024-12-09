# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2); library(parallel)

# Parameters for both archetypes
R0_scan <- c(0.75, 0.95, 1.25)
annual_spillover_rate_scan <- 5
prob_symptomatic_scan <- seq(0.05, 1, 0.05)
prob_seek_healthcare_scan <- c(0.2, 0.8)
iterations_scan <- 1:10
params <- expand_grid(R0 = R0_scan, annual_spillover_rate = annual_spillover_rate_scan,
                      prob_symptomatic = prob_symptomatic_scan, prob_seek_healthcare = prob_seek_healthcare_scan,
                      iteration = iterations_scan)

# SC2 specific params
params_SC2 <- params
params_SC2$Tg <- list(function(n) { rgamma(n, shape = 12, rate = 2) })
params_SC2$infection_to_onset_dist <- list(function(n) { rgamma(n, shape = 6, rate = 2) })
params_SC2$onset_to_healthcare_dist <- list(function(n) { rgamma(n, shape = 24, rate = 2) })
params_SC2$pathogen <- "SC2"
params_SC2$shedding_dist <- list(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16))

# EBV specific params
params_EBV <- params
params_EBV$Tg <- list(function(n) { rgamma(n, shape = 24, rate = 2) })
params_EBV$infection_to_onset_dist <- list(function(n) { rgamma(n, shape = 24, rate = 2) })
params_EBV$onset_to_healthcare_dist <- list(function(n) { rgamma(n, shape = 30, rate = 2) })
params_EBV$pathogen <- "EBV"
params_EBV$shedding_dist <- list(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16))

# Overall
overall_params <- rbind(params_SC2, params_EBV)
overall_params$index <- 1:nrow(overall_params)

# Create cluster
cl <- makeCluster(parallel::detectCores() - 2)
clusterEvalQ(cl, {
  library(wastewatchR)
  library(tidyr)
  library(dplyr)
})
clusterExport(cl, varlist = c("overall_params"))

# Run branching process
results <- parLapply(cl, 1:nrow(overall_params), function(i) {

  ## Running the branching process
  population <- 100000
  model_output <- simulate_branching_process(## Parameters shaping the outputs we'll be using
                                             initial_mn_offspring = overall_params$R0[i],
                                             disp_offspring = 1.0,
                                             generation_time_dist = overall_params$Tg[[i]],
                                             prob_symptomatic = overall_params$prob_symptomatic[i],
                                             infection_to_onset_dist = overall_params$infection_to_onset_dist[[i]],
                                             prob_severe =  overall_params$prob_symptomatic[i] / 2,
                                             prob_seek_healthcare_non_severe = overall_params$prob_seek_healthcare[i],
                                             prob_seek_healthcare_severe = overall_params$prob_seek_healthcare[i],
                                             onset_to_healthcare_dist = overall_params$onset_to_healthcare_dist[[i]],
                                             annual_spillover_rate = overall_params$annual_spillover_rate[i],
                                             spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
                                             t0 = 0,
                                             tf = Inf,
                                             population = population,
                                             check_final_size = 10000,
                                             max_num_outbreaks = 200,
                                             initial_immune = 0,
                                             seed = 2025,

                                             ## Not concerned with serology or evolution so these parameters / their outputs aren't used
                                             prob_seroconvert_asymptomatic = 1,
                                             prob_seroconvert_severe = 1,
                                             prob_seroconvert_non_severe = 1,
                                             infection_to_seroconversion_dist = function(n) { rep(100, n) },
                                             seroconversion_to_seroreversion_dist = function(n) { rep(100, n) },
                                             prob_beneficial_mutation = 0,
                                             beneficial_mutation_effect_dist = function(n) { rep(0, n) },
                                             max_mn_offspring = 10)

  ## Calculating number of healthcare seeking individuals over time
  healthcare_seeking_thresholds <- seq(1, 10, 1)
  healthcare_seeking <- generate_healthcare_seeking_time_series(branching_process_output = model_output, population = population)
  healthcare_seeking_ttd <- tibble(threshold = healthcare_seeking_thresholds) %>%
    rowwise() %>%
    mutate(
      healthcare_seeking_first_day = {
        filtered_data <- healthcare_seeking %>% filter(incidence_seek_healthcare >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
    ungroup()
  healthcare_seeking_ttd$index <- i

  ## Calculating effective number of individuals shedding over time
  wastewate_shedding_thresholds <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)
  wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = model_output,
                                                              shedding_dist = overall_params$shedding_dist[[i]],
                                                              population = population)
  wastewater_shedding_ttd <- tibble(threshold = wastewate_shedding_thresholds) %>%
    rowwise() %>%
    mutate(
      wastewater_first_day = {
        filtered_data <- wastewater_shedding %>% filter((100000 * shedding_value / population)  >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
    ungroup()
  wastewater_shedding_ttd$index <- i

  ## Returning overall outputs
  overall_outputs <- healthcare_seeking %>%
    left_join(wastewater_shedding, by = "day") %>%
    mutate(index = i)

  print(i)

  return(list(output = overall_outputs,
              time_detection_wastewater = wastewater_shedding_ttd,
              time_detection_healthcare = healthcare_seeking_ttd))

})

stopCluster(cl)



output <- vector(mode = "list", length = nrow(overall_params))
time_detection_wastewater <- vector(mode = "list", length = nrow(overall_params))
time_detection_healthcare <- vector(mode = "list", length = nrow(overall_params))

for (i in 1:nrow(overall_params)) {

  ## Running the branching process
  population <- 1000000
  model_output <- simulate_branching_process(## Parameters shaping the outputs we'll be using
    initial_mn_offspring = overall_params$R0[i],
    disp_offspring = 1.0,
    generation_time_dist = overall_params$Tg[[i]],
    prob_symptomatic = overall_params$prob_symptomatic[i],
    infection_to_onset_dist = overall_params$infection_to_onset_dist[[i]],
    prob_severe =  overall_params$prob_symptomatic[i] / 2,
    prob_seek_healthcare_non_severe = overall_params$prob_seek_healthcare[i],
    prob_seek_healthcare_severe = overall_params$prob_seek_healthcare[i],
    onset_to_healthcare_dist = overall_params$onset_to_healthcare_dist[[i]],
    annual_spillover_rate = overall_params$annual_spillover_rate[i],
    spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
    t0 = 0,
    tf = Inf,
    population = population,
    check_final_size = 10000,
    max_num_outbreaks = 200,
    initial_immune = 0,
    seed = 2025,

    ## Not concerned with serology or evolution so these parameters / their outputs aren't used
    prob_seroconvert_asymptomatic = 1,
    prob_seroconvert_severe = 1,
    prob_seroconvert_non_severe = 1,
    infection_to_seroconversion_dist = function(n) { rep(100, n) },
    seroconversion_to_seroreversion_dist = function(n) { rep(100, n) },
    prob_beneficial_mutation = 0,
    beneficial_mutation_effect_dist = function(n) { rep(0, n) },
    max_mn_offspring = 10)

  ## Calculating number of healthcare seeking individuals over time
  healthcare_seeking_thresholds <- seq(1, 10, 1)
  healthcare_seeking <- generate_healthcare_seeking_time_series(branching_process_output = model_output, population = population)
  healthcare_seeking_ttd <- tibble(threshold = healthcare_seeking_thresholds) %>%
    rowwise() %>%
    mutate(
      healthcare_seeking_first_day = {
        filtered_data <- healthcare_seeking %>% filter(incidence_seek_healthcare >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
    ungroup()
  healthcare_seeking_ttd$index <- i

  ## Calculating effective number of individuals shedding over time
  wastewate_shedding_thresholds <- seq(1, 10, 1)
  wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = model_output,
                                                              shedding_dist = overall_params$shedding_dist[[i]],
                                                              population = population)
  wastewater_shedding_ttd <- tibble(threshold = wastewate_shedding_thresholds) %>%
    rowwise() %>%
    mutate(
      wastewater_first_day = {
        filtered_data <- wastewater_shedding %>% filter(shedding_value >= threshold)
        if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
      }) %>%
    ungroup()
  wastewater_shedding_ttd$index <- i

  ## Returning overall outputs
  overall_outputs <- healthcare_seeking %>%
    left_join(wastewater_shedding, by = "day") %>%
    mutate(index = i)


  output[[i]] <- overall_outputs
  time_detection_wastewater[[i]] <- wastewater_shedding_ttd
  time_detection_healthcare[[i]] <- healthcare_seeking_ttd

  print(i)
}

## Scrap testing
x <- results[[1]]
y <- generate_outbreak_size(x, 1000000)
z <- generate_infections_time_series(x, 1000000)
plot(z$day, z$new_infections, type = "l")
a <- generate_healthcare_seeking_time_series(x, 1000000)
plot(a$day, a$incidence_seek_healthcare)
wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = x,
                                                            shedding_dist = EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16),
                                                            population = 1000000)
plot(shedding$day, shedding$shedding_value)

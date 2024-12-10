# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2); library(parallel); library(pbapply)

# Parameters for both archetypes
R0_scan <- c(0.5, 0.75, 0.95, 1.25)
annual_spillover_rate_scan <- 5
prob_symptomatic_scan <- seq(0.05, 1, 0.05)
prob_seek_healthcare_scan <- c(0.1, 0.9)
iterations_scan <- 1:25
iterations_seed <- data.frame(iteration = iterations_scan,
                              seed = ceiling(runif(n = length(iterations_scan), min = 1, max = 10^7)))
params <- expand_grid(R0 = R0_scan, annual_spillover_rate = annual_spillover_rate_scan,
                      prob_symptomatic = prob_symptomatic_scan, prob_seek_healthcare = prob_seek_healthcare_scan)

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
overall_params <- expand_grid(overall_params, iteration = iterations_scan) %>%
  left_join(iterations_seed, by = "iteration")
overall_params$index2 <- 1:nrow(overall_params)

fresh_run <- FALSE
if (fresh_run) {

  # Create cluster
  cl <- makeCluster(parallel::detectCores() - 2)
  clusterEvalQ(cl, {
    library(wastewatchR)
    library(tidyr)
    library(dplyr)
  })
  clusterExport(cl, varlist = c("overall_params"))

  # Run branching process
  results <- pblapply(X = 1:nrow(overall_params), FUN = function(i) {

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
      spillover_seeding_cases_dist = function() { rpois(1, lambda = 2) },
      t0 = 0,
      tf = Inf,
      population = population,
      check_final_size = 10000,
      max_num_outbreaks = 200,
      initial_immune = 0,
      seed = overall_params$seed[i],

      ## Not concerned with serology or evolution so these parameters / their outputs aren't used
      prob_seroconvert_asymptomatic = 1,
      prob_seroconvert_severe = 1,
      prob_seroconvert_non_severe = 1,
      infection_to_seroconversion_dist = function(n) { rep(100, n) },
      seroconversion_to_seroreversion_dist = function(n) { rep(100, n) },
      prob_beneficial_mutation = 0,
      beneficial_mutation_effect_dist = function(n) { rep(0, n) },
      max_mn_offspring = 10)

    outbreak_info <- generate_outbreak_size(branching_process_output = model_output, population = population)

    ## Calculating number of healthcare seeking individuals over time
    healthcare_seeking_thresholds <- seq(1, 10, 1)

    ### Time to detection based on clinical surveillance
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

    ### Number of outbreaks before clinical surveillance detects
    outbreak_detection_healthcare <- tibble(threshold = healthcare_seeking_thresholds) %>%
      rowwise() %>%
      mutate(
        healthcare_seeking_outbreak_number = {
          filtered_data <- outbreak_info %>% filter(total_healthcare_seek_size_with_seeding >= threshold)
          if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$outbreak)
        }) %>%
      ungroup()
    outbreak_detection_healthcare$index <- i
    outbreak_detection_healthcare <- outbreak_detection_healthcare %>%
      left_join(healthcare_seeking_ttd, by = c("threshold", "index")) %>%
      select(index, threshold, healthcare_seeking_outbreak_number, healthcare_seeking_first_day)

    ## Calculating effective number of individuals shedding over time
    wastewate_shedding_thresholds <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)

    ### Time to detection based on wastewater surveillance
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

    ### Number of outbreaks before wastewater surveillance detects
    outbreak_detection_wastewater <- wastewater_shedding_ttd %>%
      rowwise() %>%
      mutate(wastewater_outbreak_number = {
        candidate <- outbreak_info %>%
          filter(!is.na(wastewater_first_day),
                 outbreak_start <= wastewater_first_day,
                 wastewater_first_day < outbreak_end) %>%
          arrange(outbreak_start)
        if (nrow(candidate) == 0) NA_integer_ else candidate$outbreak[1]
      }) %>%
      ungroup()
    outbreak_detection_wastewater$index <- i
    outbreak_detection_wastewater <- outbreak_detection_wastewater %>%
      select(index, threshold, wastewater_outbreak_number, wastewater_first_day)


    ## Returning overall outputs
    overall_outputs <- healthcare_seeking %>%
      left_join(wastewater_shedding, by = "day") %>%
      mutate(index = i)


    return(list(output = overall_outputs,
                outbreak_info = outbreak_info,
                time_detection_wastewater = outbreak_detection_wastewater,
                time_detection_healthcare = outbreak_detection_healthcare))

  }, cl = cl)
  stopCluster(cl)
  saveRDS(list(results = results,
               overall_params = overall_params), "inst/temp_results_10thDec2024.rds")

} else {

  overall <- readRDS("inst/temp_results_10thDec2024.rds")
  overall_params <- overall[[2]]
  results <- overall[[1]]

}


## Plotting time to detection, healthcare
combined_df_time_detection_healthcare <- results %>%
  lapply(function(x) x$time_detection_healthcare[1, ]) %>%
  bind_rows() %>%
  select(index, healthcare_seeking_first_day)

combined_df_time_detection_healthcare_summary <- overall_params %>%
  left_join(combined_df_time_detection_healthcare, by = c("index2" = "index")) %>%
  group_by(R0, annual_spillover_rate, prob_symptomatic, prob_seek_healthcare, pathogen) %>%
  mutate(healthcare_seeking_first_day = ifelse(!is.na(healthcare_seeking_first_day), healthcare_seeking_first_day, max(healthcare_seeking_first_day, na.rm = TRUE))) %>%
  summarise(healthcare_seeking_iteration_success = sum(!is.na(healthcare_seeking_first_day)) / length(iterations_scan),
            healthcare_seeking_first_day_median = median(healthcare_seeking_first_day, na.rm = TRUE),
            healthcare_seeking_first_day_lower = quantile(healthcare_seeking_first_day, 0.25, na.rm = TRUE),
            healthcare_seeking_first_day_upper = quantile(healthcare_seeking_first_day, 0.75, na.rm = TRUE))

ggplot(combined_df_time_detection_healthcare_summary,
       aes(x = prob_symptomatic, y = healthcare_seeking_first_day_median,
           col = factor(prob_seek_healthcare), fill = factor(prob_seek_healthcare))) +
  geom_line() +
  geom_ribbon(aes(ymin = healthcare_seeking_first_day_lower,
                  ymax = healthcare_seeking_first_day_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Probability of Being Symptomatic", y = "Time to Detection",
       color = "Prob.\nSeek Healthcare",
       fill = "Prob.\nSeek Healthcare")

## Plotting outbreak detection number, healthcare
combined_df_time_detection_healthcare_outbreak <- results %>%
  lapply(function(x) x$time_detection_healthcare[1, ]) %>%
  bind_rows() %>%
  select(index, healthcare_seeking_outbreak_number )

combined_df_time_detection_healthcare_outbreak_summary <- overall_params %>%
  left_join(combined_df_time_detection_healthcare_outbreak, by = c("index2" = "index")) %>%
  group_by(R0, annual_spillover_rate, prob_symptomatic, prob_seek_healthcare, pathogen) %>%
  mutate(healthcare_seeking_outbreak_number  = ifelse(!is.na(healthcare_seeking_outbreak_number ), healthcare_seeking_outbreak_number, max(healthcare_seeking_outbreak_number))) %>%
  summarise(healthcare_seeking_iteration_success = sum(!is.na(healthcare_seeking_outbreak_number)) / length(iterations_scan),
            healthcare_seeking_outbreak_median = median(healthcare_seeking_outbreak_number, na.rm = TRUE),
            healthcare_seeking_outbreak_lower = quantile(healthcare_seeking_outbreak_number, 0.25, na.rm = TRUE),
            healthcare_seeking_outbreak_upper = quantile(healthcare_seeking_outbreak_number, 0.75, na.rm = TRUE))

ggplot(combined_df_time_detection_healthcare_outbreak_summary,
       aes(x = prob_symptomatic, y = healthcare_seeking_outbreak_median,
           col = factor(prob_seek_healthcare), fill = factor(prob_seek_healthcare))) +
  geom_line() +
  geom_ribbon(aes(ymin = healthcare_seeking_outbreak_lower,
                  ymax = healthcare_seeking_outbreak_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Probability of Being Symptomatic", y = "Outbreak Detection Occurs At",
       color = "Prob.\nSeek Healthcare",
       fill = "Prob.\nSeek Healthcare")

## Plotting time to detection, healthcare
combined_df_time_detection_wastewater <- results %>%
  lapply(function(x) x$time_detection_wastewater) %>%
  bind_rows() %>%
  select(index, threshold, wastewater_first_day)

combined_df_time_detection_wastewater_summary <- combined_df_time_detection_wastewater %>%
  left_join(overall_params, by = c("index" = "index2")) %>%
  mutate(wastewater_first_day = ifelse(!is.na(wastewater_first_day), wastewater_first_day, max(wastewater_first_day, na.rm = TRUE))) %>%
  group_by(R0, threshold, annual_spillover_rate, pathogen) %>%
  summarise(wastewater_first_day_median = median(wastewater_first_day, na.rm = TRUE),
            wastewater_first_day_lower = quantile(wastewater_first_day, 0.25, na.rm = TRUE),
            wastewater_first_day_upper = quantile(wastewater_first_day, 0.75, na.rm = TRUE))

ggplot(combined_df_time_detection_wastewater_summary,
       aes(x = threshold, y = wastewater_first_day_median)) +
  geom_line() +
  geom_ribbon(aes(ymin = wastewater_first_day_lower,
                  ymax = wastewater_first_day_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Sensitivity Relative to SC2 Surveillance", y = "Time to Detection") +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
                labels = paste0(rev(c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)), "x")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 1.2))

## Plotting time to detection, healthcare
combined_df_outbreak_wastewater <- results %>%
  lapply(function(x) x$time_detection_wastewater) %>%
  bind_rows() %>%
  select(index, threshold, wastewater_first_day, wastewater_outbreak_number )

combined_df_outbreak_wastewater_summary <- combined_df_outbreak_wastewater %>%
  mutate(wastewater_outbreak_number = ifelse(is.na(wastewater_outbreak_number) & !is.na(wastewater_first_day), 1, wastewater_outbreak_number)) %>%
  left_join(overall_params, by = c("index" = "index2")) %>%
  mutate(wastewater_outbreak_number = ifelse(!is.na(wastewater_outbreak_number), wastewater_outbreak_number, max(wastewater_outbreak_number, na.rm = TRUE))) %>%
  group_by(R0, threshold, annual_spillover_rate, pathogen) %>%
  summarise(wastewater_outbreak_median = median(wastewater_outbreak_number, na.rm = TRUE),
            wastewater_outbreak_lower = quantile(wastewater_outbreak_number, 0.25, na.rm = TRUE),
            wastewater_outbreak_upper = quantile(wastewater_outbreak_number, 0.75, na.rm = TRUE))

ggplot(combined_df_outbreak_wastewater_summary,
       aes(x = threshold, y = wastewater_outbreak_median)) +
  geom_line() +
  geom_ribbon(aes(ymin = wastewater_outbreak_lower,
                  ymax = wastewater_outbreak_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Sensitivity Relative to SC2 Surveillance", y = "Outbreak Detection Occurs At") +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
                labels = paste0(rev(c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)), "x")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 1.2))


## Final set of plots for now -
colnames(combined_df_time_detection_healthcare_summary)
colnames(combined_df_time_detection_wastewater_summary)

combined_overall <- combined_df_time_detection_healthcare_summary %>%
  left_join(combined_df_time_detection_wastewater_summary,
            by = c("R0", "annual_spillover_rate", "pathogen")) %>%
  mutate(detection = ifelse(healthcare_seeking_first_day_median > wastewater_first_day_median,
                            "wastewater", "clinical")) %>%
  select(R0, annual_spillover_rate, pathogen, prob_symptomatic, prob_seek_healthcare, threshold,
         healthcare_seeking_first_day_median, wastewater_first_day_median, detection) %>%
  filter(prob_seek_healthcare == 0.9)

ggplot(combined_overall,
       aes(x = threshold, y = prob_symptomatic, fill = detection)) +
  geom_tile() +
  facet_grid(pathogen ~ R0) +
  scale_fill_manual(values = c("wastewater" = "#80ADA0", "clinical" = "#EED2CC")) +
  theme_bw() +
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000),
                labels = paste0(rev(c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)), "x"))

colnames(combined_overall)

# output <- vector(mode = "list", length = nrow(overall_params))
# time_detection_wastewater <- vector(mode = "list", length = nrow(overall_params))
# time_detection_healthcare <- vector(mode = "list", length = nrow(overall_params))
#
# for (i in 441:nrow(overall_params)) {
#
#   ## Running the branching process
#   population <- 1000000
#   model_output <- simulate_branching_process(## Parameters shaping the outputs we'll be using
#     initial_mn_offspring = overall_params$R0[i],
#     disp_offspring = 1.0,
#     generation_time_dist = overall_params$Tg[[i]],
#     prob_symptomatic = overall_params$prob_symptomatic[i],
#     infection_to_onset_dist = overall_params$infection_to_onset_dist[[i]],
#     prob_severe =  overall_params$prob_symptomatic[i] / 2,
#     prob_seek_healthcare_non_severe = overall_params$prob_seek_healthcare[i],
#     prob_seek_healthcare_severe = overall_params$prob_seek_healthcare[i],
#     onset_to_healthcare_dist = overall_params$onset_to_healthcare_dist[[i]],
#     annual_spillover_rate = overall_params$annual_spillover_rate[i],
#     spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
#     t0 = 0,
#     tf = Inf,
#     population = population,
#     check_final_size = 10000,
#     max_num_outbreaks = 200,
#     initial_immune = 0,
#     seed = 2025,
#
#     ## Not concerned with serology or evolution so these parameters / their outputs aren't used
#     prob_seroconvert_asymptomatic = 1,
#     prob_seroconvert_severe = 1,
#     prob_seroconvert_non_severe = 1,
#     infection_to_seroconversion_dist = function(n) { rep(100, n) },
#     seroconversion_to_seroreversion_dist = function(n) { rep(100, n) },
#     prob_beneficial_mutation = 0,
#     beneficial_mutation_effect_dist = function(n) { rep(0, n) },
#     max_mn_offspring = 10)
#
#   ## Calculating number of healthcare seeking individuals over time
#   healthcare_seeking_thresholds <- seq(1, 10, 1)
#   healthcare_seeking <- generate_healthcare_seeking_time_series(branching_process_output = model_output, population = population)
#   healthcare_seeking_ttd <- tibble(threshold = healthcare_seeking_thresholds) %>%
#     rowwise() %>%
#     mutate(
#       healthcare_seeking_first_day = {
#         filtered_data <- healthcare_seeking %>% filter(incidence_seek_healthcare >= threshold)
#         if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
#       }) %>%
#     ungroup()
#   healthcare_seeking_ttd$index <- i
#
#   ## Calculating effective number of individuals shedding over time
#   wastewate_shedding_thresholds <- seq(1, 10, 1)
#   wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = model_output,
#                                                               shedding_dist = overall_params$shedding_dist[[i]],
#                                                               population = population)
#   wastewater_shedding_ttd <- tibble(threshold = wastewate_shedding_thresholds) %>%
#     rowwise() %>%
#     mutate(
#       wastewater_first_day = {
#         filtered_data <- wastewater_shedding %>% filter(shedding_value >= threshold)
#         if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
#       }) %>%
#     ungroup()
#   wastewater_shedding_ttd$index <- i
#
#   ## Returning overall outputs
#   overall_outputs <- healthcare_seeking %>%
#     left_join(wastewater_shedding, by = "day") %>%
#     mutate(index = i)
#
#   output[[i]] <- overall_outputs
#   time_detection_wastewater[[i]] <- wastewater_shedding_ttd
#   time_detection_healthcare[[i]] <- healthcare_seeking_ttd
#
#   print(i)
# }
## Scrap testing
# x <- results[[1]]
# y <- generate_outbreak_size(x, 1000000)
# z <- generate_infections_time_series(x, 1000000)
# plot(z$day, z$new_infections, type = "l")
# a <- generate_healthcare_seeking_time_series(x, 1000000)
# plot(a$day, a$incidence_seek_healthcare)
# wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = x,
#                                                             shedding_dist = EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16),
#                                                             population = 1000000)
# plot(shedding$day, shedding$shedding_value)

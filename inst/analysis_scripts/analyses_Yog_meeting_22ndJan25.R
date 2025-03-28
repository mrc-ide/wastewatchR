# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2); library(parallel); library(pbapply)

# Parameters for both archetypes
wastewater_shedding_relative_SC2 <- seq(log10(0.1), log10(15), 0.025)
R0_scan <- c(0.50, 0.75, 0.95)
annual_spillover_rate_scan <- 5
prob_symptomatic_scan <- seq(0.05, 1, 0.025)
prob_seek_healthcare_scan <- c(0.1, 0.5, 0.9)
iterations_scan <- 1:30
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
params_SC2$shedding_dist <- list(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16) / max(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16)))

# EBV specific params
params_EBV <- params
params_EBV$Tg <- list(function(n) { rgamma(n, shape = 24, rate = 2) })
params_EBV$infection_to_onset_dist <- list(function(n) { rgamma(n, shape = 24, rate = 2) })
params_EBV$onset_to_healthcare_dist <- list(function(n) { rgamma(n, shape = 30, rate = 2) })
params_EBV$pathogen <- "EBV"
params_EBV$shedding_dist <- list(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16) / max(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16)))

# Overall
overall_params <- rbind(params_SC2, params_EBV)
overall_params$index <- 1:nrow(overall_params)
overall_params <- expand_grid(overall_params, iteration = iterations_scan) %>%
  left_join(iterations_seed, by = "iteration")
overall_params$index2 <- 1:nrow(overall_params)

fresh_run <- FALSE
if (fresh_run) {

  # Create cluster
  num_cores <- parallel::detectCores() - 2 # 4
  cl <- makeCluster(num_cores)
  clusterEvalQ(cl, {
    library(wastewatchR)
    library(tidyr)
    library(dplyr)
  })
  clusterExport(cl, varlist = c("overall_params", "wastewater_shedding_relative_SC2"))

  # Run branching process
  results <- pblapply(X = 1:nrow(overall_params), FUN = function(i) {

    ## Running the branching process
    population <- 100000
    model_output <- simulate_branching_process(

      ## Parameters relevant to the model outputs in this particular set of analyses
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

    ## Generating information on each of the outbreaks
    outbreak_info <- generate_outbreak_size(branching_process_output = model_output, population = population)

    ############################
    # Clinical Surveillance
    ############################

    # Extracting time to detection based on clinical surveillance and the properties of the outbreak that led to detection
    ## Note - poss alternative approach to this to consider next time is just find the first person who seeks healthcare
    ##        in the raw results dataframe, and then get outbreak id from their row etc.
    healthcare_seeking_time_series <- generate_healthcare_seeking_time_series(branching_process_output = model_output)

    ## Time-to-detection
    healthcare_seeking_ttd <- tibble(threshold = 1) %>%
      mutate(
        healthcare_seeking_first_day = {
          filtered_data <- healthcare_seeking_time_series %>% filter(incidence_seek_healthcare >= threshold)
          if (nrow(filtered_data) == 0) NA_real_ else min(filtered_data$day)
        }) %>%
      select(healthcare_seeking_first_day)

    ## Outbreak info
    healthcare_seeking_outbreak_info <- extract_outbreak_characteristics_time(outbreak_info = outbreak_info,
                                                                              outbreak_id = healthcare_seeking_ttd$healthcare_seeking_first_day,
                                                                              id_type = "time") %>%
      select(outbreak, total_infection_size_with_seeding)

    healthcare_seeking_results <- cbind(healthcare_seeking_ttd, healthcare_seeking_outbreak_info) %>%
      mutate(index = i)

    ############################
    # Wastewater Surveillance
    ############################

    ## Creating dataframe to store results for detection based on wastewater surveillance and different shedding amounts
    wastewater_results <- tibble(wastewater_shedding_relative_SC2 = 10^wastewater_shedding_relative_SC2,
                                 wastewater_first_day = NA_real_,
                                 outbreak = NA_real_,
                                 total_infection_size_with_seeding = NA_real_,
                                 index = i)

    ## Looping over the different shedding amounts and calculating time-to-detection, extracting outbreak characteristics etc
    for (j in 1:length(wastewater_shedding_relative_SC2)) {

      ## Generate time-series of effective number of individuals shedding over time
      wastewater_shedding <- generate_number_shedding_time_series(branching_process_output = model_output,
                                                                  shedding_dist = overall_params$shedding_dist[[i]],
                                                                  shedding_relative_SC2 = 10^(wastewater_shedding_relative_SC2[j]))

      ## Calculating time-to-detection for wastewater surveillance
      temp_wastewater_ttd <- calculate_wastewater_ttd(wastewater_number_shedding_time_series = wastewater_shedding,
                                                      sampling_frequency = 14,
                                                      sampling_method = "autosampler",
                                                      detection_approach = "logistic_curve",
                                                      detection_params = list(population = 100000,
                                                                              logistic_beta_0 = -1.229996,
                                                                              logistic_beta_1 = 0.258775,
                                                                              limit_of_detection = 0.5,
                                                                              seed = 1000))
      wastewater_results$wastewater_first_day[j] <- temp_wastewater_ttd$wastewater_first_day

      ## Outbreak info (only if detection successfully occurs)
      if (!is.na(temp_wastewater_ttd$wastewater_first_day)) {
        wastewater_outbreak_info <- extract_outbreak_characteristics_time(outbreak_info = outbreak_info,
                                                                          outbreak_id = temp_wastewater_ttd$wastewater_first_day,
                                                                          id_type = "time") %>%
          select(outbreak, total_infection_size_with_seeding)
        wastewater_results$outbreak[j] <- wastewater_outbreak_info$outbreak
        wastewater_results$total_infection_size_with_seeding[j] <- wastewater_outbreak_info$total_infection_size_with_seeding
      }
    }

    return(list(wastewater_results = wastewater_results,
                healthcare_seeking_results = healthcare_seeking_results))

  }, cl = cl)
  stopCluster(cl)
  saveRDS(list(results = results,
               overall_params = overall_params), "inst/temp_results_30thDec2024.rds")

} else {

  overall <- readRDS("inst/temp_results_30thDec2024.rds")
  overall_params <- overall[[2]]
  results <- overall[[1]]

}
overall_params$index2 <- 1:nrow(overall_params)

## Plotting time to detection, healthcare
combined_df_time_detection_healthcare <- results %>%
  lapply(function(x) x$healthcare_seeking_results) %>%
  bind_rows() %>%
  rename(healthcare_seeking_outbreak = outbreak,
         healthcare_seeking_outbreak_size = total_infection_size_with_seeding)

combined_df_time_detection_healthcare_summary <- overall_params %>%
  left_join(combined_df_time_detection_healthcare, by = c("index2" = "index")) %>%
  group_by(R0, annual_spillover_rate, prob_symptomatic, prob_seek_healthcare, pathogen) %>%
  mutate(healthcare_seeking_first_day = ifelse(!is.na(healthcare_seeking_first_day), healthcare_seeking_first_day, max(healthcare_seeking_first_day, na.rm = TRUE))) %>%
  mutate(healthcare_seeking_outbreak  = ifelse(!is.na(healthcare_seeking_outbreak), healthcare_seeking_outbreak, max(healthcare_seeking_outbreak, na.rm = TRUE))) %>%
  mutate(healthcare_seeking_outbreak_size  = ifelse(!is.na(healthcare_seeking_outbreak_size), healthcare_seeking_outbreak_size, max(healthcare_seeking_outbreak_size, na.rm = TRUE))) %>%
  summarise(healthcare_seeking_iteration_success = sum(!is.na(healthcare_seeking_first_day)) / length(iterations_scan),
            healthcare_seeking_first_day_median = median(healthcare_seeking_first_day, na.rm = TRUE),
            healthcare_seeking_first_day_lower = quantile(healthcare_seeking_first_day, 0.25, na.rm = TRUE),
            healthcare_seeking_first_day_upper = quantile(healthcare_seeking_first_day, 0.75, na.rm = TRUE),
            healthcare_seeking_outbreak_median = median(healthcare_seeking_outbreak, na.rm = TRUE),
            healthcare_seeking_outbreak_lower = quantile(healthcare_seeking_outbreak, 0.25, na.rm = TRUE),
            healthcare_seeking_outbreak_upper = quantile(healthcare_seeking_outbreak, 0.75, na.rm = TRUE),
            healthcare_seeking_outbreak_size_median = median(healthcare_seeking_outbreak_size, na.rm = TRUE),
            healthcare_seeking_outbreak_size_lower = quantile(healthcare_seeking_outbreak_size, 0.25, na.rm = TRUE),
            healthcare_seeking_outbreak_size_upper = quantile(healthcare_seeking_outbreak_size, 0.75, na.rm = TRUE))

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

ggplot(combined_df_time_detection_healthcare_summary,
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

ggplot(combined_df_time_detection_healthcare_summary,
       aes(x = prob_symptomatic, y = healthcare_seeking_outbreak_size_median,
           col = factor(prob_seek_healthcare), fill = factor(prob_seek_healthcare))) +
  geom_line() +
  geom_ribbon(aes(ymin = healthcare_seeking_outbreak_size_lower,
                  ymax = healthcare_seeking_outbreak_size_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Probability of Being Symptomatic", y = "Outbreak Size Detection Occurs At",
       color = "Prob.\nSeek Healthcare",
       fill = "Prob.\nSeek Healthcare")

## Plotting time to detection, healthcare
combined_df_time_detection_wastewater <- results %>%
  lapply(function(x) x$wastewater_results) %>%
  bind_rows() %>%
  rename(wastewater_outbreak = outbreak,
         wastewater_outbreak_size = total_infection_size_with_seeding)

combined_df_time_detection_wastewater_summary <- combined_df_time_detection_wastewater %>%
  left_join(overall_params, by = c("index" = "index2")) %>%
  mutate(wastewater_first_day = ifelse(!is.na(wastewater_first_day), wastewater_first_day, max(wastewater_first_day, na.rm = TRUE))) %>%
  mutate(wastewater_outbreak  = ifelse(!is.na(wastewater_outbreak ), wastewater_outbreak , max(wastewater_outbreak, na.rm = TRUE))) %>%
  mutate(wastewater_outbreak_size  = ifelse(!is.na(wastewater_outbreak_size), wastewater_outbreak_size, max(wastewater_outbreak_size, na.rm = TRUE))) %>%
  group_by(R0, wastewater_shedding_relative_SC2, annual_spillover_rate, pathogen) %>%
  summarise(wastewater_first_day_median = median(wastewater_first_day, na.rm = TRUE),
            wastewater_first_day_lower = quantile(wastewater_first_day, 0.25, na.rm = TRUE),
            wastewater_first_day_upper = quantile(wastewater_first_day, 0.75, na.rm = TRUE),
            wastewater_outbreak_median = median(wastewater_outbreak, na.rm = TRUE),
            wastewater_outbreak_lower = quantile(wastewater_outbreak, 0.25, na.rm = TRUE),
            wastewater_outbreak_upper = quantile(wastewater_outbreak, 0.75, na.rm = TRUE),
            wastewater_outbreak_size_median = median(wastewater_outbreak_size, na.rm = TRUE),
            wastewater_outbreak_size_lower = quantile(wastewater_outbreak_size, 0.25, na.rm = TRUE),
            wastewater_outbreak_size_upper = quantile(wastewater_outbreak_size, 0.75, na.rm = TRUE))

ggplot(combined_df_time_detection_wastewater_summary,
       aes(x = wastewater_shedding_relative_SC2, y = wastewater_first_day_median)) +
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
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100),
                labels = paste0(rev(c(1e-3, 1e-2, 1e-1, 1, 10, 100)), "x")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 1.2))

ggplot(combined_df_time_detection_wastewater_summary,
       aes(x = wastewater_shedding_relative_SC2, y = wastewater_outbreak_median)) +
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
  scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10, 100),
                labels = paste0(c(1e-3, 1e-2, 1e-1, 1, 10, 100), "x")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 1.2))


ggplot(combined_df_time_detection_wastewater_summary,
       aes(x = wastewater_shedding_relative_SC2, y = wastewater_outbreak_size_median)) +
  geom_line() +
  geom_ribbon(aes(ymin = wastewater_outbreak_size_lower,
                  ymax = wastewater_outbreak_size_upper,
                  col = NULL),
              alpha = 0.2) +
  facet_wrap(~ pathogen + R0, , scales = "free_y",
             nrow = 2, ncol = 4,
             labeller = labeller(R0 = function(x) paste0("R0 = ", x))) +
  theme_bw() +
  labs(x = "Sensitivity Relative to SC2 Surveillance", y = "Outbreak Size Detection Occurs At") +
  scale_x_log10(breaks = c(1e-2, 1e-1, 1, 10, 100),
                labels = paste0(c(1e-2, 1e-1, 1, 10, 100), "x")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 1.2))

## Final set of plots for now -
colnames(combined_df_time_detection_healthcare_summary)
colnames(combined_df_time_detection_wastewater_summary)

combined_overall <- combined_df_time_detection_healthcare_summary %>%
  left_join(combined_df_time_detection_wastewater_summary,
            by = c("R0", "annual_spillover_rate", "pathogen")) %>%
  mutate(detection = ifelse(healthcare_seeking_first_day_median > wastewater_first_day_median,
                            "Wastewater\nSurveillance", "Clinical\nSurveillance")) %>%
  select(R0, annual_spillover_rate, pathogen, prob_symptomatic, prob_seek_healthcare, wastewater_shedding_relative_SC2,
         healthcare_seeking_first_day_median, wastewater_first_day_median, detection) %>%
  filter(R0 == 0.95)


my_pathogen_labels <- c("EBV" = "Ebola Virus-Like",
                        "SC2" = "SARS-CoV-2-Like")
my_prob_seek_labels <- c("0.1" = "10% Seek Healthcare",
                         "0.5" = "50% Seek Healthcare",
                         "0.9" = "90% Seek Healthcare")
ggplot(combined_overall,
       aes(x = wastewater_shedding_relative_SC2, y = prob_symptomatic, fill = detection)) +
  geom_tile() +
  facet_grid(prob_seek_healthcare ~ pathogen,
             labeller = labeller(
               pathogen = my_pathogen_labels,
               prob_seek_healthcare = my_prob_seek_labels
             )) +
  labs(y = "Probability of Symptoms", x = "Shedding Relative to SC2",
       fill = "Detects First") +
  scale_fill_manual(values = c("Wastewater\nSurveillance" = "#EB4B98",
                               "Clinical\nSurveillance" = "#28C6D5")) +
  theme_bw() +
  scale_x_log10(breaks = c(1e-1, 1, 10, 100),
                labels = paste0(c(1e-1, 1, 10, 100), "x"))

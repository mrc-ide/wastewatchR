# Compare the simple, evolutionary single-outbreak, and multi-outbreak simulators.
# Run from the package root with Rscript inst/analysis_scripts/compare_outbreak_simulators.R

suppressPackageStartupMessages({
  library(devtools)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(cowplot)
})
devtools::load_all(".", quiet = TRUE)

horizon_days <- 3 * 365
n_replicates <- 100L
catchment_size <- 100000
spillover_rate_per_day <- 0.1
diagnosis_probability <- 0.1
max_simulated_infections <- 50000L
parameter_sets <- data.frame(
  parameter_set = factor(paste0("R0 = ", c(0.20, 0.40, 0.60, 0.80, 0.95)),
                         levels = paste0("R0 = ", c(0.20, 0.40, 0.60, 0.80, 0.95))),
  R0 = c(0.20, 0.40, 0.60, 0.80, 0.95)
)
methods <- c("Simple single outbreak", "Evolutionary single outbreak",
             "Multi-outbreak simulator")
shedding_profile <- readRDS("misc/ww_sensitivity_data/fecal_shedding.rds")$shedding
detection_fit <- readRDS("misc/ww_sensitivity_data/ww_detection_params.rds")
beta0 <- unname(detection_fit["alpha", "mean"])
beta1 <- unname(detection_fit["beta", "mean"])
lod <- 0.01

# Severe disease probability is conditional on symptoms in both single-outbreak
# functions and absolute among infections in simulate_branching_process.
prob_symptomatic <- 0.20
prob_severe_absolute <- 0.10
prob_severe_given_symptoms <- prob_severe_absolute / prob_symptomatic
p_seek_non_severe <- 0.10
p_seek_severe <- 0.90

draw_spillovers <- function(seed) {
  set.seed(seed)
  daily <- spillover(
    time = horizon_days, specify = "spillover_rate", seasonal_period = 365,
    spillover_rate_pars = list(tmax = 183, b = spillover_rate_per_day,
                               d = 0, sigma = 63),
    prevalence_pars = NULL, contact_rate_pars = NULL, p_infection = NA_real_
  )
  event_days <- rep(daily$t, daily$spillovers)
  data.frame(index_case_ID = seq_along(event_days), spillover_day = event_days)
}

simulate_single_events <- function(method, R0, events, seed) {
  set.seed(seed)
  output <- vector("list", nrow(events))
  for (i in seq_len(nrow(events))) {
    common <- list(
      spillover_day = events$spillover_day[i],
      index_case_ID = events$index_case_ID[i],
      disp_offspring = 2,
      generation_time_dist = function(n) rgamma(n, shape = 12, rate = 2),
      prob_symptomatic = prob_symptomatic,
      infection_to_onset_dist = function(n) rgamma(n, shape = 6, rate = 2),
      prob_seek_healthcare_non_severe = p_seek_non_severe,
      prob_seek_healthcare_severe = p_seek_severe,
      onset_to_healthcare_dist = function(n) rgamma(n, shape = 9, rate = 1.5),
      initial_immune = 0
    )
    if (method == "Simple single outbreak") {
      output[[i]] <- do.call(sim_single_outbreak, c(common, list(
        mn_offspring = R0, max_gen = 100,
        prob_severe = prob_severe_given_symptoms,
        prob_diagnosis = 1,
        healthcare_to_diagnosis_dist = function(n) rep(0, n)
      )))
    } else {
      output[[i]] <- do.call(simulate_single_outbreak_evolutionary, c(common, list(
        initial_mn_offspring = R0, prob_severe = prob_severe_given_symptoms,
        prob_diagnosis = 1,
        healthcare_to_diagnosis_dist = function(n) rep(0, n),
        prob_seroconvert_asymptomatic = 0.8,
        prob_seroconvert_severe = 0.8,
        prob_seroconvert_non_severe = 0.8,
        infection_to_seroconversion_dist = function(n) rgamma(n, shape = 56, rate = 2),
        seroconversion_to_seroreversion_dist = function(n) rgamma(n, shape = 730, rate = 2),
        prob_beneficial_mutation = 0, max_mn_offspring = 4,
        population = catchment_size, check_final_size = max_simulated_infections
      )))
    }
  }
  dplyr::bind_rows(output)
}

simulate_multi_outbreak <- function(R0, seed) {
  set.seed(seed)
  simulate_branching_process(
    initial_mn_offspring = R0, disp_offspring = 2,
    generation_time_dist = function(n) rgamma(n, shape = 12, rate = 2),
    prob_symptomatic = prob_symptomatic,
    infection_to_onset_dist = function(n) rgamma(n, shape = 6, rate = 2),
    prob_severe = prob_severe_absolute,
    prob_seek_healthcare_non_severe = p_seek_non_severe,
    prob_seek_healthcare_severe = p_seek_severe,
    onset_to_healthcare_dist = function(n) rgamma(n, shape = 9, rate = 1.5),
    prob_seroconvert_asymptomatic = 0.8,
    prob_seroconvert_severe = 0.8,
    prob_seroconvert_non_severe = 0.8,
    infection_to_seroconversion_dist = function(n) rgamma(n, shape = 56, rate = 2),
    seroconversion_to_seroreversion_dist = function(n) rgamma(n, shape = 730, rate = 2),
    prob_beneficial_mutation = 0, max_mn_offspring = 4,
    annual_spillover_rate = spillover_rate_per_day * 365,
    spillover_seeding_cases_dist = function() 1L,
    t0 = 0, tf = horizon_days, population = catchment_size,
    initial_immune = 0, check_final_size = max_simulated_infections,
    max_num_outbreaks = 250, seed = seed
  )
}

detect_wastewater <- function(cases, seed) {
  set.seed(seed)
  infection_days <- floor(cases$time_infection)
  infection_days <- infection_days[infection_days >= 0 & infection_days <= horizon_days]
  incidence <- tabulate(infection_days + 1L, nbins = horizon_days + 1L)
  shedding <- numeric(horizon_days + 1L)
  for (lag in seq_along(shedding_profile)) {
    offset <- lag - 1L
    if (offset <= horizon_days) {
      idx <- (offset + 1L):(horizon_days + 1L)
      shedding[idx] <- shedding[idx] +
        incidence[seq_len(horizon_days + 1L - offset)] * shedding_profile[lag]
    }
  }
  phase <- sample.int(7L, 1L) - 1L
  sample_days <- seq.int(phase, horizon_days, by = 7L)
  sampled <- shedding[sample_days + 1L]
  p <- ifelse(sampled < lod, 0,
              plogis(beta0 + beta1 * log10(100000 * (sampled + 1e-3) / catchment_size)))
  any(rbinom(length(p), 1, p) == 1L)
}

run_trial <- function(param_i, replicate, R0) {
  seed <- 202611170 + param_i * 10000 + replicate * 10
  events <- draw_spillovers(seed + 1L)
  methods_cases <- list(
    "Simple single outbreak" = simulate_single_events(
      "Simple single outbreak", R0, events, seed + 2L),
    "Evolutionary single outbreak" = simulate_single_events(
      "Evolutionary single outbreak", R0, events, seed + 3L),
    "Multi-outbreak simulator" = simulate_multi_outbreak(R0, seed + 4L)
  )
  dplyr::bind_rows(lapply(seq_along(methods_cases), function(m) {
    method <- names(methods_cases)[m]
    all_cases <- methods_cases[[m]]
    cases <- all_cases[!is.na(all_cases$time_infection) &
                         all_cases$time_infection >= 0 &
                         all_cases$time_infection <= horizon_days, , drop = FALSE]
    seeking <- which(cases$seek_healthcare == 1 &
                       !is.na(cases$time_seek_healthcare) &
                       cases$time_seek_healthcare <= horizon_days)
    set.seed(seed + 100L + m)
    clinical <- length(seeking) > 0 &&
      any(rbinom(length(seeking), 1, diagnosis_probability) == 1L)
    ww <- detect_wastewater(cases, seed + 200L + m)
    spill_count <- if (method == "Multi-outbreak simulator") {
      sum(!is.na(all_cases$time_infection) &
            all_cases$infection_generation == 1L &
            all_cases$time_infection <= horizon_days)
    } else nrow(events)
    data.frame(
      parameter_set = parameter_sets$parameter_set[param_i],
      R0 = R0, replicate = replicate, method = method,
      spillovers = spill_count, infections_in_window = nrow(cases),
      wastewater_detected = ww, clinical_detected = clinical
    )
  }))
}

# 5 parameter sets x 100 replicates x 3 methods.
results <- vector("list", nrow(parameter_sets) * n_replicates)
k <- 1L
for (p in seq_len(nrow(parameter_sets))) {
  message("Running ", parameter_sets$parameter_set[p], "...")
  for (replicate in seq_len(n_replicates)) {
    results[[k]] <- run_trial(p, replicate, parameter_sets$R0[p])
    k <- k + 1L
  }
}
replicate_results <- dplyr::bind_rows(results)
replicate_results$method <- factor(replicate_results$method, levels = methods)

summarise_detection <- function(x) {
  rate <- mean(x)
  n <- length(x)
  se <- sqrt(rate * (1 - rate) / n)
  c(rate = rate, lower = max(0, rate - 1.96 * se),
    upper = min(1, rate + 1.96 * se), n = n)
}
summary_rows <- replicate_results |>
  dplyr::group_by(parameter_set, R0, method) |>
  dplyr::summarise(
    wastewater = list(summarise_detection(wastewater_detected)),
    clinical = list(summarise_detection(clinical_detected)),
    .groups = "drop"
  )
detection_summary <- dplyr::bind_rows(lapply(seq_len(nrow(summary_rows)), function(i) {
  ww <- summary_rows$wastewater[[i]]
  cl <- summary_rows$clinical[[i]]
  data.frame(parameter_set = summary_rows$parameter_set[i],
    R0 = summary_rows$R0[i], method = summary_rows$method[i],
    wastewater_rate = ww["rate"], wastewater_lower = ww["lower"],
    wastewater_upper = ww["upper"], clinical_rate = cl["rate"],
    clinical_lower = cl["lower"], clinical_upper = cl["upper"], n = ww["n"])
}))

theme_clean <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.title.x = element_blank())
p_size <- ggplot(replicate_results,
  aes(x = parameter_set, y = infections_in_window, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.82), outlier.alpha = 0.15,
               width = 0.72) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Infections during the three-year window",
       y = "Number of infections (log scale)", fill = NULL) +
  theme_clean + theme(legend.position = "none")

plot_detection <- function(rate, lower, upper, title) {
  ggplot(detection_summary, aes(x = parameter_set, y = .data[[rate]],
                                color = method, group = method)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = .data[[lower]], ymax = .data[[upper]]), width = 0.1) +
    scale_color_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = title, y = "Replicates detected", color = NULL) + theme_clean
}
p_ww <- plot_detection("wastewater_rate", "wastewater_lower",
                       "wastewater_upper", "Wastewater detection frequency")
p_clinical <- plot_detection("clinical_rate", "clinical_lower",
                             "clinical_upper", "Clinical detection frequency")
legend <- cowplot::get_legend(p_ww)
figure <- cowplot::plot_grid(
  p_size, p_ww + theme(legend.position = "none"),
  p_clinical + theme(legend.position = "none"),
  ncol = 1, align = "v", rel_heights = c(1.2, 1, 1)
)
figure <- cowplot::plot_grid(figure, legend, ncol = 1, rel_heights = c(1, 0.08))

output_dir <- "inst/analysis_scripts/outbreak_simulator_comparison"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(replicate_results,
          file.path(output_dir, "replicate_detection_results.csv"), row.names = FALSE)
write.csv(detection_summary,
          file.path(output_dir, "detection_frequency_summary.csv"), row.names = FALSE)
ggsave(file.path(output_dir, "outbreak_simulator_comparison.png"),
       figure, width = 11, height = 12, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "outbreak_simulator_comparison.pdf"),
       figure, width = 11, height = 12, bg = "white")
writeLines(c(
  "Three-year spillover and observation window; 100 replicates per R0 set and method.",
  "Flat spillover rate: 0.1/day (36.5/year), as agreed for this comparison.",
  "Both single-outbreak functions share the same daily Poisson spillover schedule. The multi-outbreak simulator draws exponential spillover intervals internally; its rate is set to the equivalent annual rate, so event times are not paired.",
  paste0("Absolute P(severe) = ", prob_severe_absolute,
         "; both single-outbreak functions receive P(severe | symptoms) = ",
         prob_severe_given_symptoms, "."),
  "Both single-outbreak functions take prob_severe conditional on symptoms; simulate_branching_process takes an absolute probability among all infections. The script passes 0.5 to both single-outbreak functions and 0.1 to the multi-outbreak function, with P(symptomatic)=0.2.",
  "Both single-outbreak functions have prob_diagnosis and healthcare_to_diagnosis_dist. The multi-outbreak function only records healthcare seeking and its time, without diagnosis probability or diagnosis time.",
  paste0("For a common clinical endpoint, the script sets the single-outbreak functions' internal diagnosis probability to 1 and their delays to 0, then applies an external P(diagnosis | healthcare sought) = ",
         diagnosis_probability, " draw to every method. Clinical detection time is not saved in this run."),
  "The evolutionary single-outbreak and multi-outbreak functions model seroconversion and seroreversion; sim_single_outbreak does not. Serology was assigned common settings but was not used in the detection comparison.",
  "Beneficial mutation probability is zero in all functions. The same negative-binomial dispersion (2) and natural-history delay distributions were used.",
  "Wastewater: weekly grab samples, example shedding profile, and logistic coefficient posterior means.",
  "The simple simulator keeps its susceptible fraction fixed (initial_immune is a proportion). The evolutionary and multi-outbreak functions decrement an absolute susceptible count as infections accrue. The evolutionary function resets that count on each spillover call; simulate_branching_process shares it across outbreaks in one run.",
  "The simple simulator is capped at 100 generations; the evolutionary and multi-outbreak simulations are capped at 50,000 infections.",
  "This run emitted three recycling warnings from sim_single_outbreak's infector assignment. These occurred with the finite generation cap; treat those simple-simulator chains as potentially truncated or mislinked.",
  "simulate_branching_process currently does not use tf to stop generation. This script allows up to 250 outbreaks and filters infections and detection to the three-year window.",
  "Outbreak size is infections occurring within the three-year window."
), file.path(output_dir, "comparison_notes.txt"))
print(detection_summary)
message("Saved results in ", output_dir)

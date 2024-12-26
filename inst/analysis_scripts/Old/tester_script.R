# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2)

# Run branching process
population <- 1000000
model_output <- simulate_branching_process(initial_mn_offspring = 0.75,
                                           disp_offspring = 1.0,
                                           generation_time_dist = function(n) { rgamma(n, shape = 12, rate = 2) },
                                           prob_symptomatic = 0.8,
                                           infection_to_onset_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           prob_severe = 0.35,
                                           prob_seek_healthcare_non_severe = 0.9,
                                           prob_seek_healthcare_severe = 0.9,
                                           onset_to_healthcare_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           prob_seroconvert_asymptomatic = 0.8,
                                           prob_seroconvert_severe = 0.8,
                                           prob_seroconvert_non_severe = 0.8,
                                           infection_to_seroconversion_dist = function(n) { rgamma(n, shape = 56, rate = 2) },
                                           seroconversion_to_seroreversion_dist = function(n) { rgamma(n, shape = 730, rate = 2) },
                                           prob_beneficial_mutation = 0.05,
                                           beneficial_mutation_effect_dist = function(n) {rexp(n, rate = 10) },
                                           max_mn_offspring = 1.5,
                                           annual_spillover_rate = 2,
                                           spillover_seeding_cases_dist = function() { rpois(1, lambda = 5) },
                                           initial_immune = 0,
                                           t0 = 0,
                                           tf = Inf,
                                           population = population, # always make sure this is a pretty large number unless deliberately want susceptible depletion
                                           check_final_size = 10000,
                                           max_num_outbreaks = 200,
                                           seed = 2025)
x <- generate_outbreak_size(model_output, 1000000)

# Post-processing
infections <- generate_infections_time_series(branching_process_output = model_output)
plot(infections$day, infections$new_infections, type = "l", ylim = c(0, 50),
     xlab = "Time", ylab = "Infections")

symptom_onsets <- generate_symptom_onset_time_series(branching_process_output = model_output, population = population)
plot(symptom_onsets$day, symptom_onsets$incidence_symptom_onset, type = "l")

healthcare_seeking <- generate_healthcare_seeking_time_series(branching_process_output = model_output)
plot(healthcare_seeking$day, healthcare_seeking$incidence_seek_healthcare, type = "l")

healthcare_seeking <- generate_healthcare_seeking_time_series(branching_process_output = model_output, population = population)
plot(healthcare_seeking$day, healthcare_seeking$incidence_seek_healthcare, type = "l")

seropositivity <- generate_seropositivity_timeseries(branching_process_output = model_output, population = population)
plot(seropositivity$time, seropositivity$seropositive_abs)
plot(seropositivity$time, log10(seropositivity$cumulative_seroconversions + 1), type = "l")

shedding_dist <- EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16) / max(EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16))
shedding <- generate_number_shedding_time_series(branching_process_output = model_output, shedding_dist = shedding_dist, shedding_relative_SC2 = 5)
plot(shedding$day, log(shedding$shedding_value + 1))


d <- ggplot(infections, aes(x = day, y = new_infections)) +
  geom_line(col = "#D49137") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 10)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  labs(x = "Time (Days)", y = "Infection Incidence", title = "Daily Incidence of Infections") +
  theme_bw()
a <- ggplot(healthcare_seeking, aes(x = day, y = incidence_seek_healthcare)) +
  geom_line(col = "#8075FF") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 5)) +
  labs(x = "Time (Days)", y = "Clinical Visits Per Day", title = "Clinical Surveillance") +
  theme_bw()
b <- ggplot(shedding, aes(x = day, y = log(shedding_value + 1))) +
  geom_line(col = "#1DBD83") +
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 4)) +
  labs(x = "Time (Days)", y = "Effective Number of People Shedding", title = "Wastewater Surveillance") +
  theme_bw()
c <- ggplot(seropositivity, aes(x = time, y = 100 * 2500 * seropositive_prop_pop)) +
  geom_line(col = "#D566C5") + # "#B279A7"
  coord_cartesian(xlim = c(0, 365 * 3), ylim = c(0, 100 * 0.25)) +
  labs(x = "Time (Days)", y = "Seropositivity (% of Population)", title = "Serological Surveillance") +
  theme_bw()

cowplot::plot_grid(a, b, c, nrow = 3)


nrow(model_output)
sum(model_output$symptomatic, na.rm = TRUE)
sum(model_output$seek_healthcare, na.rm = TRUE)
sum(model_output$time_symptom_onset > model_output$time_infection, na.rm = TRUE)
sum(model_output$time_seek_healthcare > model_output$time_symptom_onset, na.rm = TRUE)

sum(model_output$symptomatic)
sum(model_output$severe)
sum(model_output$seek_healthcare)
hist(model_output$n_offspring[!is.na(model_output$n_offspring)], na.rm = TRUE)
mean(model_output$n_offspring[!is.na(model_output$n_offspring)], na.rm = TRUE)

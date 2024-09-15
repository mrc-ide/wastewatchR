# Load required libraries
library(wastewatchR); library(dplyr); library(tidyr); library(ggplot2)

# Run branching process
population <- 10^6
model_output <- simulate_branching_process(mn_offspring = 1.75,
                                           check_final_size = 10000,
                                           prob_symptomatic = 0.5,
                                           prob_hosp = 0.1,
                                           prob_seek_healthcare = 0.8,
                                           generation_time_dist = function(n) { rgamma(n, shape = 12, rate = 2) }, ## poss need to reverse rate and shape - TBD
                                           hospitalisation_delay_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           hospitalisation_duration_dist = function(n) { rgamma(n, shape = 6, rate = 2) },
                                           population = population)

# Post-processing
hospitalisations <- generate_hospitalisation_time_series(branching_process_output = model_output, population = population)
symptom_onsets <- generate_symptom_onset_time_series(branching_process_output = model_output, population = population)
shedding_dist <- EpiSewer::get_discrete_gamma(gamma_mean = 4, gamma_sd = 2.4, maxX = 16)
shedding <- generate_number_shedding_time_series(branching_process_output = model_output, shedding_dist = shedding_dist, population = population)
infections <- generate_infections_time_series(branching_process_output = model_output, population = population)

# Plotting
overall <- hospitalisations %>%
  left_join(symptom_onsets, by = "day") %>%
  left_join(shedding, by = "day") %>%
  left_join(infections, by = "day") %>%
  select(day, new_infections, incidence_symptom_onset, incidence_hospitalisation, hospitalised, shedding_value) %>%
  pivot_longer(-day, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  filter(day <= day[which.max(value)])

ggplot(overall, aes(x = day, y = value, col = metric)) +
  geom_line() +
  labs(x = "Day", y = "Number of Individuals") +
  scale_colour_discrete(labels = c("In Hospital", "Hosp. Incidence", "Symptom. Incidence",
                                   "Infection Incidence", "Number Shedding"), name = "") +
  theme_bw()






# plot(infections$day, log(infections$new_infections + 1), type = "l",
#      xlab = "Time (Days)", ylab = "Incidence of Infections")
# lines(symptom_onsets$day, log(symptom_onsets$incidence_symptom_onset + 1), type = "l",
#       xlab = "Time (Days)", ylab = "Symptom Onset Incidence", col = "blue")
# lines(shedding$day, log(shedding$shedding_value + 1), type = "l",
#       xlab = "Time (Days)", ylab = "Symptom Onset Incidence", col = "orange")
# lines(hospitalisations$day, log(hospitalisations$hospitalised + 1), type = "l",
#       xlab = "Time (Days)", ylab = "Incidence", col = "purple", ylim = c(0, max(hospitalisations$hospitalised)))
# lines(hospitalisations$day, log(hospitalisations$incidence_hospitalisation + 1), type = "l",
#       xlab = "Time (Days)", ylab = "Incidence", col = "red", ylim = c(0, max(infections$new_infections)))

# Example

# This script provides an example of how to use wastewatchR to simulate
# spillover, onward transmission, and surveillance for a virus with high
# spillover rate and low R0
# disease

# load libraries
library(tidyverse)
library(reshape2)
library(igraph)
library(ggraph)
library(cowplot)

devtools::install_github("mrc-ide/wastewatchR@dev", force = TRUE)

library(wastewatchR)

# 1. Simulate spillovers over a 1-yr period ------------------------------------
# ------------------------------------------------------------------------------
# Here assume a seasonal spillover rate that varies with an annual periodicity
# peaking on day 183 of the year, with an average daily spillover rate of 0.7
# in our chosen catchment area of 100,000 people

# Note that if you wish to simulate spillover mechanistically, there is also the
# option to specify "swiss_cheese" instead of "spillover_rate".
# In this case you can specify the prevalence of infection in the reservoir,
# contact rate with people, and probability of infection through contact, the
# product of which defines the spillover rate.

catchment_size <- 100000

sp <- spillover(time = 365, specify = "spillover_rate",
                seasonal_period = 365,
                spillover_rate_pars = list("tmax" = 183,
                                           "b" = solve_b(x=0.7,tmax=183,
                                                         sigma=60,d=3),
                                           "d" = 3,
                                           "sigma" = 63))%>%
  filter(spillovers >0)

# account for possibility of multiple spillovers on the same day
sp <- sp %>%
  uncount(spillovers) %>%
  rowid_to_column("ID")

# 2. For each spillover, simulate onward transmission based aon assumed R0 -----
# ------------------------------------------------------------------------------

# List of spillover events that could initiate onward transmission
events <- vector("list", length = dim(sp)[1])

# Specify clinical surveillance parameters
prob_symptomatic <- 0.2
prob_seek_healthcare <- 0.9
prob_diagnosis <- 0.1

## delay infection-->onset
shape_onset <- 6
rate_onset <- 2

## delay onset--> healthcare
shape_hc <- 9
rate_hc <- 1.5

## delay healthcare-->diagnosis
shape_diag <- 6
rate_diag <- 2

# Simulate outbreaks following each spillover event
## Specify offspring distribution and delays:
for(s in 1:(length(events))){
  events[[s]] <- sim_single_outbreak(spillover_day = sp$t[s],
                                     index_case_ID = sp$ID[s],
                                     mn_offspring = 0.2, # low R0 of 0.2
                                     disp_offspring = 2, # >1 --> overdispersion
                                     max_gen = Inf,
                                     index_cases = 1,
                                     generation_time_dist = function(n)
                                     { rgamma(n, shape = 12, rate = 2) },
                                     prob_symptomatic = prob_symptomatic,
                                     infection_to_onset_dist = function(n)
                                     { rgamma(n, shape = shape_onset, rate = rate_onset) },
                                     prob_severe = 1, # given symptomatic
                                     prob_seek_healthcare_non_severe = 0,
                                     prob_seek_healthcare_severe = prob_seek_healthcare,
                                     onset_to_healthcare_dist = function(n)
                                     { rgamma(n, shape = shape_hc, rate = rate_hc) },
                                     prob_diagnosis = prob_diagnosis,
                                     healthcare_to_diagnosis_dist = function(n)
                                     { rgamma(n, shape = shape_diag, rate = rate_diag) },
                                     initial_immune = 0, # no prior immunity
                                     seed = 102)
}

## Merge linelists from all outbreaks and spillovers into 1 mega linelist
linelist <- bind_rows(events)%>%
  mutate(type = if_else(infector == "animal", "Animal-to-human", "Human-to-human"))


# 3. Wastewater detection ------------------------------------------------------
# ------------------------------------------------------------------------------

# Load your assumed shedding distribution here
## as a placeholder we used one based on fecal shedding profile of SARS-CoV-2
shedding_dist <- readRDS("misc/ww_sensitivity_data/fecal_shedding.rds")

# Generate time series of number of individuals shedding each day
## Choose whether the model is based on the total number of people shedding
## (i.e. anyone shedding counts as one shedder), or on the effective number
## of shedders, which accounts for changes in shedding amount over time since
## infection. Under the latter approach, a person on their peak day of shedding
## counts as 1 shedder, while someone shedding at 10% of their peak amount
## counts as 0.1 shedders.
nts <- generate_number_shedding_time_series(linelist,
                                            method = "effective_n_shedders",
                                            shedding_dist = shedding_dist$shedding,
                                            shedding_relative_SC2 = 1)

# Specify your chosen model describing the relationship between your measure of
# shedding (number of shedders or effective shedders) and probabilty of detection
## Here we use a model based on data from SARS-CoV-2 at a quarantine facility in
## New Zealand.
logistic_params <- readRDS("misc/ww_sensitivity_data/ww_detection_params.rds")
detection_params <- list(
  population = catchment_size,
  logistic_beta_0 = logistic_params[1,1],
  logistic_beta_1 = logistic_params[2,1],
  limit_of_detection = 0.01, # Specify LOD in terms of shedders/effective shedders
  seed = 123
)

# Simulate wastewater detection
det <- calculate_wastewater_ttd(wastewater_number_shedding_time_series = nts,
                                sampling_frequency = 7,
                                sampling_method = "grab",
                                detection_approach = "logistic_curve",
                                detection_params)

# 4. Clinical detection probability --------------------------------------------
# ------------------------------------------------------------------------------

# Compute probability of clinical detection each day based on delay distributions
# and assumed probabilities of symptoms/healthcare seeking/diagnosis.
det_clinical <- calculate_clinical_detection_prob(nts = nts,
                                                  det = det,
                                                  delay_onset = dgamma(0:21,
                                                                       shape = shape_onset,
                                                                       rate = rate_onset),
                                                  delay_seek = dgamma(0:21,
                                                                      shape = shape_hc,
                                                                      rate = rate_hc),
                                                  delay_diag = dgamma(0:21,
                                                                      shape = shape_diag,
                                                                      rate = rate_diag),
                                                  prob_symptomatic = prob_symptomatic,
                                                  prob_seek_healthcare = prob_seek_healthcare,
                                                  prob_diagnosis = prob_diagnosis)

# 5. Take a look at results ----------------------------------------------------
# ------------------------------------------------------------------------------

##### A. infections (by type) and viral load shed to wastewater

# Approximate start and mid points for each month (non-leap year)
month_start <- c(1, 32, 60, 91, 121, 152,
                 182, 213, 244, 274, 305, 335, 365)
p1 <- ggplot(linelist) +
  theme_classic() +
  theme(legend.position = "top",
        text = element_text(size = 15),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(hjust = -0.66),) +
  geom_ribbon(data = det$sampled_data, aes(x = day, ymax = shedding_value*1.5, ymin = 0), fill = "#D3D3D3") +
  geom_histogram(aes(
    x = time_infection,
    fill = factor(type, levels = c("Human-to-human", "Animal-to-human"))
  ),
  binwidth = 1
  ) +
  scale_fill_manual(values = c("#826699", "#9CAF88"), name = "Source") +
  xlab(NULL) +
  ylab("Incident infections") +
  scale_x_continuous(
    breaks = month_start, # tick marks at start of each month
    labels = c(month.abb, " "),
    limits = c(1, 365),
    expand = c(0.01,0)
  ) +
  scale_y_continuous(
    # secondary y-axis - rescale to account for shedding rel to sc2
    sec.axis = sec_axis(~ .*(1/1)/1.5, name = "Effective n individuals\nshedding into WW")
  )

##### B. probability of detection in wastewater & stochastic detection events
p2 <- ggplot(det$sampled_data)+
  geom_vline(data = det$sampled_data %>% filter(sampled == "Yes"),
             aes(xintercept = day), linewidth = 0.5, colour = "grey90")+
  geom_ribbon(data = det$sampled_data,
              aes(x = day, ymax = prob_detect, ymin = 0), fill = "salmon4", alpha = 0.6)+
  geom_line(aes(y = prob_detect, x = day), colour = "salmon4", alpha = 0.3)+
  geom_vline(data = det$sampled_data %>% filter(sampled == "Yes" & detect_draw>0),
             aes(xintercept = day), linewidth = 0.5, colour = "coral3")+

  scale_x_continuous(
    breaks = month_start, # tick marks at start of each month
    labels = c(NULL), # labels for months
    limits = c(1, 365),
    expand = c(0.01,0)
  ) +
  xlab(NULL)+
  ylab("Probability detected\nin wastewater")+
  theme_classic()+
  theme(text = element_text(size = 15))

##### C. probability of clinical detection & stochastic detection events

p3 <- ggplot(det_clinical)+
  geom_ribbon(aes(x = day, ymax = prob_detection_clinical, ymin = 0), fill = "#5B7C99", alpha = 0.5)+
  geom_line(aes(y = prob_detection_clinical, x = day), colour = "#5B7C99", alpha = 0.7)+
  geom_vline(data = linelist, aes(xintercept = time_diagnosis), linewidth = 0.5, colour = "coral3")+
  xlab(NULL)+
  scale_x_continuous(
    breaks = month_start,        # tick marks at start of each month
    labels = c(month.abb, " "),           # labels for months
    limits = c(1, 365),
    expand = c(0.01,0)
  ) +
  scale_x_continuous(
    breaks = month_start,                 # tick marks at start of each month
    labels = c(month.abb, " "),
    limits = c(1, 365),
    expand = c(0.01,0)
  ) +
  ylab("Probability detected\nclinically ")+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = -0.66),
        text = element_text(size = 15))

##### D. plot transmission trees

# Create edges and vertices
linelistplot <- linelist %>%
  filter(time_infection<366)%>%
  mutate(
    id_outbreak = paste0(id, "-", spillover_ID),
    id_infector = paste0(infector, "-", spillover_ID),
  ) %>%
  select(id_outbreak, id_infector, symptomatic, diagnosis, type,
         time_infection, infection_generation) %>%
  rename(to = id_outbreak, from = id_infector)
edges <- linelistplot %>%
  select(from, to)

# Make sure vertices$"name" exists and covers ALL from/to
vertices <- linelistplot %>%
  transmute(name = to,
            symptomatic, diagnosis, type, time_infection,
            infection_generation)
# Add missing infector names (animal reservoir)
animal_names <- unique(linelistplot$from[grepl(linelistplot$from,
                                               pattern = "animal")])
animal_res <- linelistplot %>% filter(from %in% animal_names)%>%
  select(from, time_infection)%>%
  mutate(name = from,
         symptomatic = NA,
         diagnosis = NA,
         type = "index",
         infection_generation = 0)
vertices <- bind_rows(vertices, animal_res) %>%
  distinct(name, .keep_all = TRUE)  # remove duplicates if any

# Build edges with attributes of the 'to' node
edges_attr <- edges %>%
  left_join(vertices %>% select(name, type), by = c("to" = "name"))

# Build graph with edge attributes included
g <- graph_from_data_frame(d = edges_attr, vertices = vertices, directed = TRUE)

# Layout to plot by time and generation and add jitter to humans but not animals
vertices <- vertices %>%
  mutate(infection_generation = if_else(grepl(name, pattern = "animal"),
                                        infection_generation,
                                        infection_generation + rnorm(n(), 0, 0.1)))
# Flip the layout so gneration 0 is at the top of the y axis
flipped_layout <- as.matrix(vertices[, c("time_infection", "infection_generation")])
flipped_layout[, 2] <- -flipped_layout[, 2]  # Multiply the y-axis by -1

# Plot
p0 <- ggraph(g, layout = flipped_layout) +
  geom_edge_diagonal(arrow = arrow(length = unit(2, 'mm')),
                     end_cap = circle(1.2, 'mm'),
                     width = 0.5, colour = "grey40") +
  geom_node_point(aes(color = type), size = 4) +
  theme_classic() +
  ylab("Infection generation")+
  xlab(NULL)+
  scale_x_continuous(
    breaks = month_start,        # tick marks at start of each month
    labels = NULL,           # labels for months
    limits = c(1, 365),
    expand = c(0.01,0)
  ) +
  scale_y_continuous(labels = c(5,4,3,2,1,0), breaks = c(-5,-4,-3,-2,-1,0))+
  scale_color_manual(values = c("#9CAF88","#826699","grey40"),
                     name = "Transmission source")+
  guides(
    color = guide_legend(nrow = 1, byrow = TRUE),
    edge_color = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(1.2, "cm"),
        axis.text.x = element_text(hjust = -0.75),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 15))

legend <- get_legend(p1)

p0123<- plot_grid(p0+theme(legend.position = "none"), p1+theme(legend.position = "none"),p2,p3, nrow = 4, align = "v")

ptot <- plot_grid(legend,p0123, nrow = 2, rel_heights = c(1,10))

print(ptot)

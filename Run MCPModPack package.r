setwd("c:/temp")

library(MCPModPack)

####################################################

# MCPMod-based analysis

# Example 1 (normally distributed endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 0.5

# Perform an MCPMod-based analysis of the trial's data
results = MCPModAnalysis(endpoint_type = "Normal", 
		                 models = models, 
		                 dose = normal$dose, 
		                 resp = normal$resp, 
		                 alpha = alpha, 
		                 direction = direction, 
		                 model_selection = model_selection, 
		                 Delta = Delta)

# Simple summary of the MCPMod analysis results
results

# Detailed summary of the MCPMod analysis results
AnalysisReport(results, "MCPMod analysis summary (Normally distributed endpoint)", "MCPMod analysis summary (Normally distributed endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based analysis of the trial's data
AnalysisApp()

####################################################

# MCPMod-based analysis

# Example 2 (binary endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 0.3

# Perform an MCPMod-based analysis of the trial's data
results = MCPModAnalysis(endpoint_type = "Binary", 
                     models = models, 
                     dose = binary$dose, 
                     resp = binary$resp, 
                     alpha = alpha, 
                     direction = direction, 
                     model_selection = model_selection, 
                     Delta = Delta)

# Simple summary of the MCPMod analysis results
results

# Detailed summary of the MCPMod analysis results
AnalysisReport(results, "MCPMod analysis summary (Binary endpoint)", "MCPMod analysis summary (Binary endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based analysis of the trial's data
AnalysisApp()

####################################################

# MCPMod-based analysis

# Example 3 (count endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 5

# Vector of overdisperstion parameters
theta = c(2, 2, 2, 2, 2)

# Perform an MCPMod-based analysis of the trial's data
results = MCPModAnalysis(endpoint_type = "Count", 
                     models = models, 
                     dose = count$dose, 
                     resp = count$resp, 
                     alpha = alpha, 
                     direction = direction, 
                     model_selection = model_selection, 
                     Delta = Delta,
                     theta = theta)

# Simple summary of the MCPMod analysis results
results

# Detailed summary of the MCPMod analysis results
AnalysisReport(results, "MCPMod analysis summary (Count endpoint)", "MCPMod analysis summary (Count endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based analysis of the trial's data
AnalysisApp()

####################################################

# MCPMod-based simulations

# Example 4 (normally distributed endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 0.2

# Select the assumed dose-response model, values of the non-linear model parameters and the standard deviation of the outcome variable in each trial arm (required for normally distributed endpoint)
sim_models = list(emax = 1, 
                  placebo_effect = 0.5, 
                  max_effect = seq(from = 0, to = 0.5, by = 0.1), 
                  sd = c(0.5, 0.5, 0.75, 0.75, 0.95))

# Simulation parameters (number of patients in each trial arm, dose levels, dropout rate, number of simulations, threshold for computing go probabilities)
sim_parameters = list(n = c(40, 40, 40, 40, 40),
                      doses = c(0, 0.05, 0.2, 0.6, 1),
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = 0.6)

# Perform an MCPMod-based simulation
results = MCPModSimulation(endpoint_type = "Normal", 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results
SimulationReport(results, "MCPMod simulation summary (Normally distributed endpoint)", "MCPMod simulation summary (Normally distributed endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based simulation
SimulationApp()

####################################################

# MCPMod-based simulations

# Example 5 (binary endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 0.3

# Select the assumed dose-response model and values of the non-linear model parameters
sim_models = list(emax = 1, 
                  placebo_effect = 0.2, 
                  max_effect = seq(from = 0, to = 0.5, by = 0.1))

# Simulation parameters (number of patients in each trial arm, dose levels, dropout rate, number of simulations, threshold for computing go probabilities)
sim_parameters = list(n = c(40, 40, 40, 40, 40),
                      doses = c(0, 0.05, 0.2, 0.6, 1),
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = 0.4)

# Perform an MCPMod-based simulation
results = MCPModSimulation(endpoint_type = "Binary", 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results
SimulationReport(results, "MCPMod simulation summary (Binary endpoint)", "MCPMod simulation summary (Binary endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based simulation
SimulationApp()

####################################################

# MCPMod-based simulations

# Example 6 (count endpoint)

# Select the candidate dose-response models and initial values of the non-linear model parameters (linear, quadratic, exponential, emax, Logistic and sigemax)
models = list(linear = NA, quadratic = -1, exponential = 2, emax = 0.2, logistic = c(0.1, 1), sigemax = c(0.1, 1))

# One-sided Type I error rate
alpha = 0.025

# Direction of the dose-response relationship (a larger value of the mean treatment difference corresponds to a beneficial treatment effect)
direction = "increasing"

# Model selection criterion
model_selection = "AIC"

# The treatment effect for identifying the target dose 
Delta = 3

# Vector of overdisperstion parameters
theta = c(2, 2, 2, 2, 2)

# Select the assumed dose-response model and values of the non-linear model parameters
sim_models = list(exponential = 1, 
                  placebo_effect = 2, 
                  max_effect = seq(from = 0, to = 5, by = 1))

# Simulation parameters (number of patients in each trial arm, dose levels, dropout rate, number of simulations, threshold for computing go probabilities)
sim_parameters = list(n = c(40, 40, 40, 40, 40),
                      doses = c(0, 0.05, 0.2, 0.6, 1),
                      dropout_rate = 0.05,
                      nsims = 1000,
                      go_threshold = 4)

# Perform an MCPMod-based simulation
results = MCPModSimulation(endpoint_type = "Count", 
                           models = models, 
                           alpha = alpha, 
                           direction = direction, 
                           model_selection = model_selection, 
                           Delta = Delta,
                           theta = theta,
                           sim_models = sim_models,
                           sim_parameters = sim_parameters)

# Simple summary of the MCPMod simulation results
results

# Detailed summary of the MCPMod simulation results
SimulationReport(results, "MCPMod simulation summary (Count endpoint)", "MCPMod simulation summary (Count endpoint).docx") 

# Launch a Shiny application to perform an MCPMod-based simulation
SimulationApp()

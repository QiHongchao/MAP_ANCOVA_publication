##Packages
library(mvtnorm)
library(rjags)
library(extraDistr)

##Simulation scenarios
num_historical <- c(3, 5)
with_trt_eff <- c(0, 1)
heterogeneity_level <- 1:5 ##1: No, 2: Small, 3: Moderate, 4: Substantial, 5: Large.
simulation_scenarios <- expand.grid(num_historical = num_historical, with_trt_eff = with_trt_eff, 
                                    heterogeneity_level = heterogeneity_level)
simulation_scenarios <- simulation_scenarios[order(simulation_scenarios$num_historical, simulation_scenarios$with_trt_eff),]
nrow(simulation_scenarios)
num_scenarios <- nrow(simulation_scenarios)

##Random seeds
set.seed(1234)
seed_scenarios <- sample.int(.Machine$integer.max, size = num_scenarios)

##Simulated parameters
beta <- c(1.06, 1.16)
trt_eff <- -3 ##clinically meaningful difference for ADAS-cog in one year
between_rho_beta <- -0.9
between_sd_beta <- list(No = c(0, 0), Small = c(20/8, 0.72/8), Moderate = c(20/4, 0.72/4), 
                        Substantial = c(20/2, 0.72/2), Large = c(20, 0.72))
num_per_arm <- 60
baseline_mean <- 24.53
baseline_sd <- 9.47
error_sd <- 6.79
num_simulations <- 1000
num_chains <- 4
num_iter <- 5000
perc_burnin <- 0.2

##JAGS data
##Number of fixed effects, i.e. the intercept and the baseline effect
num_fixef <- 2
##Scale parameters of HN priors for between-study standard deviations of the intercept and the baseline effect
sd_sigma_beta0 <- 20/2
sd_sigma_beta1 <- 0.72/2

##Function definition
##ANCOVA model data generation
data_gen <- function(seed, beta, trt_eff, between_sd_beta, between_rho_beta, num_per_arm, num_historical,
                     baseline_mean, baseline_sd, error_sd) {
  set.seed(seed)
  ##Between-study covariance matrix for the regression coefficients, the current results are not stable with high between-study sd if we use cov_beta
  ##to model beta_current
  cov_beta <- matrix(c(between_sd_beta[1]^2,
                       rep(between_rho_beta*between_sd_beta[1]*between_sd_beta[2], 2),
                       between_sd_beta[2]^2), 2)
  
  ##Sample the current parameters
  beta_current <- beta
  data_current <- data.frame(study = "c", 
                             baseline = rnorm(num_per_arm*2, baseline_mean, baseline_sd),
                             group = rep(0:1, each = num_per_arm))
  current_design_matrix <- model.matrix(~ baseline + group, data = data_current)
  data_current$followup <- as.numeric(current_design_matrix %*% c(beta_current, trt_eff) + 
                                        rnorm(num_per_arm*2, 0, error_sd))
  
  ##Historical data generation
  data_historical <- data.frame()
  for (i in 1:num_historical) {
    ##Sample historical parameters
    beta_historical <- rmvnorm(1, beta, cov_beta)
    
    ##Generate separate historical data sets
    data_historical_tmp <- data.frame(study = paste0("h", i), 
                                      baseline = rnorm(num_per_arm, baseline_mean, baseline_sd),
                                      group = rep(0, num_per_arm))
    historical_design_matrix <- model.matrix(~ baseline, data = data_historical_tmp)
    data_historical_tmp$followup <- as.numeric(historical_design_matrix %*% t(beta_historical) + 
                                                 rnorm(num_per_arm, 0, error_sd))
    data_historical <- rbind(data_historical, data_historical_tmp)
  }
  
  ##Return the data sets
  return(list(data_current = data_current, data_historical = data_historical))
}
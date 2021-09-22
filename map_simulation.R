rm(list = ls())

##Simulation setting
source("simulation_settings.R")

##Data sets to save the results
res_map <- list()
res_mac <- list()

start <- Sys.time()
for (i in 1:num_scenarios) {
  ##Generate seeds for data generation process
  set.seed(seed_scenarios[i])
  seed_datagen <- sample.int(.Machine$integer.max, size = num_simulations)
  ##Generate seeds for JAGS model
  set.seed(seed_scenarios[i])
  seed_simulations <- sample.int(.Machine$integer.max, size = num_simulations * num_chains)
  
  ##Save the results per scenario
  res_map[[i]] <- data.frame(simulation = 1:num_simulations, beta0 = NA, beta1 = NA,
                             sigma_beta0 = NA, sigma_beta1 = NA, sd_beta0 = NA, sd_beta1 = NA, rho = NA)
  res_mac[[i]] <- data.frame(simulation = 1:num_simulations, bias = NA, sd = NA, mse = NA, power = NA)
  
  for (j in 1:num_simulations) {
    ##Generate data sets
    data_all <- data_gen(seed = seed_datagen[j], beta = beta, trt_eff = trt_eff * simulation_scenarios$with_trt_eff[i], 
                         between_sd_beta = between_sd_beta[[simulation_scenarios$heterogeneity_level[i]]], 
                         between_rho_beta = between_rho_beta, num_per_arm = num_per_arm, 
                         num_historical = simulation_scenarios$num_historical[i],
                         baseline_mean = baseline_mean, baseline_sd = baseline_sd, error_sd = error_sd)
    
    ##Prepare JAGS data for MAP and MAC
    studyindex <- sort(c((0:(simulation_scenarios$num_historical[i]-1))*num_per_arm + 1,
                         (1:simulation_scenarios$num_historical[i])*num_per_arm))
    jags_data_map <- list(fixef_historical = model.matrix(~ baseline, data = data_all$data_historical), 
                          response_historical = data_all$data_historical$followup,
                          H = simulation_scenarios$num_historical[i], studyindex = studyindex,
                          N = nrow(data_all$data_current), num_fixef = num_fixef,
                          fixef_new = model.matrix(~ baseline + group, data = data_all$data_current), 
                          response_new = data_all$data_current$followup,
                          mean_beta = rep(0, num_fixef), prec_beta = diag(c(1e-4, 1e-2)))
    ##Precision of HN priors for between-study standard deviations of the intercept and the baseline effect
    jags_data_map$prec_sigma_beta0 <- 1/(sd_sigma_beta0^2)
    jags_data_map$prec_sigma_beta1 <- 1/(sd_sigma_beta1^2)
                    
    ##MAP
    map_jags_model <- jags.model(file = "./jags/map_jags_model.txt", data = jags_data_map, 
                                     n.chains = num_chains, n.adapt = 1000, quiet = T,
                                     inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains]))
    )
    
    ##Burnin stage
    update(map_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta_new", "sigma_beta0", "sigma_beta1", "rho")
    
    map_jags <- coda.samples(map_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                                 progress.bar = "none")
    
    map_sample <- as.data.frame(do.call(rbind, map_jags))
    
    ##Operating characteristics
    res_map[[i]]$beta0[j] <- mean(map_sample$'beta_new[1]')
    res_map[[i]]$beta1[j] <- mean(map_sample$'beta_new[2]')
    res_map[[i]]$sigma_beta0[j] <- mean(map_sample$sigma_beta0)
    res_map[[i]]$sigma_beta1[j] <- mean(map_sample$sigma_beta1)
    res_map[[i]]$sd_beta0[j] <- sd(map_sample$'beta_new[1]')
    res_map[[i]]$sd_beta1[j] <- sd(map_sample$'beta_new[2]')
    res_map[[i]]$rho[j] <- mean(map_sample$rho)
    
    ##MAC
    mac_jags_model <- jags.model(file = "./jags/mac_jags_model.txt", data = jags_data_map, 
                                     n.chains = num_chains, n.adapt = 1000, quiet = T,
                                     inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains])))
    ##Burnin stage
    update(mac_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta_new", "lambda")
    
    mac_jags <- coda.samples(mac_jags_model, params, n.iter = num_iter * (1 - perc_burnin), progress.bar = "none")
    
    mac_sample <- as.data.frame(do.call(rbind, mac_jags))
    
    ##Operating characteristics
    res_mac[[i]]$bias[j] <- mean(mac_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i]
    res_mac[[i]]$sd[j] <- sd(mac_sample$lambda)
    res_mac[[i]]$mse[j] <- (mean(mac_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i])^2
    res_mac[[i]]$power[j] <- (quantile(mac_sample$lambda, 0.025) > 0 | 
                                    quantile(mac_sample$lambda, 0.975) < 0)
    
    ##Progress
    if (j%%50 == 0) {
      print(paste("scenario:", i, ", simulation:", j)) 
    }
  }
}
Sys.time() - start

rm(list = ls())

##Simulation setting
source("simulation_settings.R")

##Data sets to save the results
res_current <- list()
res_pooling_hist <- list()
res_pooling <- list()

start <- Sys.time()
for (i in 1:num_scenarios) {
  ##Generate seeds for data generation process
  set.seed(seed_scenarios[i])
  seed_datagen <- sample.int(.Machine$integer.max, size = num_simulations)
  ##Generate seeds for JAGS model
  set.seed(seed_scenarios[i])
  seed_simulations <- sample.int(.Machine$integer.max, size = num_simulations * num_chains)
  
  ##Save the results per scenario
  res_current[[i]] <- data.frame(simulation = 1:num_simulations, bias = NA, sd = NA, mse = NA, power = NA)
  res_pooling_hist[[i]] <- data.frame(simulation = 1:num_simulations, beta0 = NA, beta1 = NA, sd_beta0 = NA, sd_beta1 = NA)
  res_pooling[[i]] <- data.frame(simulation = 1:num_simulations, bias = NA, sd = NA, mse = NA, power = NA)
  
  for (j in 1:num_simulations) {
    ##Generate data sets
    data_all <- data_gen(seed = seed_datagen[j], beta = beta, trt_eff = trt_eff * simulation_scenarios$with_trt_eff[i], 
                         between_sd_beta = between_sd_beta[[simulation_scenarios$heterogeneity_level[i]]], 
                         between_rho_beta = between_rho_beta, num_per_arm = num_per_arm, 
                         num_historical = simulation_scenarios$num_historical[i],
                         baseline_mean = baseline_mean, baseline_sd = baseline_sd, error_sd = error_sd)
    
    ##Prepare JAGS data
    ##No borrowing
    jags_data_current <- list(N = nrow(data_all$data_current), 
                              fixef = model.matrix(~ baseline + group, data = data_all$data_current), 
                              response = data_all$data_current$followup,
                              mean_beta = rep(0, 2), prec_beta = diag(c(1e-4, 1e-2)))
    ##Pooling historical
    jags_data_pooling_hist <- list(N = nrow(data_all$data_historical), 
                                   fixef = model.matrix(~ baseline, data = data_all$data_historical), 
                                   response = data_all$data_historical$followup,
                                   mean_beta = rep(0, 2), prec_beta = diag(c(1e-4, 1e-2)))
    
    ##Pooling
    pooling <- rbind(data_all$data_current, data_all$data_historical)
    jags_data_pooling <- list(N = nrow(pooling), 
                              fixef = model.matrix(~ baseline + group, data = pooling), 
                              response = pooling$followup,
                              mean_beta = rep(0, 2), prec_beta = diag(c(1e-4, 1e-2)))
    
    ##No borrowing
    current_jags_model <- jags.model(file = "./jags/current_pooling_jags_model.txt", data = jags_data_current, 
                                     n.chains = num_chains, n.adapt = 1000, quiet = T,
                                     inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains]))
                                     )
    ##Burnin stage
    update(current_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta", "lambda")
    
    current_jags <- coda.samples(current_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                                 progress.bar = "none")
    
    current_sample <- as.data.frame(do.call(rbind, current_jags))
    
    ##Operating characteristics
    res_current[[i]]$bias[j] <- mean(current_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i]
    res_current[[i]]$sd[j] <- sd(current_sample$lambda)
    res_current[[i]]$mse[j] <- (mean(current_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i])^2
    res_current[[i]]$power[j] <- (quantile(current_sample$lambda, 0.025) > 0 | 
                                    quantile(current_sample$lambda, 0.975) < 0)
    
    ##Pooling historical data sets
    pooling_hist_jags_model <- jags.model(file = "./jags/pooling_hist_jags_model.txt", data = jags_data_pooling_hist, 
                                     n.chains = num_chains, n.adapt = 1000, quiet = T,
                                     inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains]))
    )
    ##Burnin stage
    update(pooling_hist_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta")
    
    pooling_hist_jags <- coda.samples(pooling_hist_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                                      progress.bar = "none")
    
    pooling_hist_sample <- as.data.frame(do.call(rbind, pooling_hist_jags))
    
    ##Operating characteristics
    res_pooling_hist[[i]]$beta0[j] <- mean(pooling_hist_sample$'beta[1]')
    res_pooling_hist[[i]]$beta1[j] <- mean(pooling_hist_sample$'beta[2]')
    res_pooling_hist[[i]]$sd_beta0[j] <- sd(pooling_hist_sample$'beta[1]')
    res_pooling_hist[[i]]$sd_beta1[j] <- sd(pooling_hist_sample$'beta[2]')
    
    ##Pooling
    pooling_jags_model <- jags.model(file = "./jags/current_pooling_jags_model.txt", data = jags_data_pooling, 
                                     n.chains = num_chains, n.adapt = 1000, quiet = T,
                                     inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                  list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains])))
    ##Burnin stage
    update(pooling_jags_model, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta", "lambda")
    
    pooling_jags <- coda.samples(pooling_jags_model, params, n.iter = num_iter * (1 - perc_burnin), 
                                 progress.bar = "none")
    
    pooling_sample <- as.data.frame(do.call(rbind, pooling_jags))
    
    ##Operating characteristics
    res_pooling[[i]]$bias[j] <- mean(pooling_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i]
    res_pooling[[i]]$sd[j] <- sd(pooling_sample$lambda)
    res_pooling[[i]]$mse[j] <- (mean(pooling_sample$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i])^2
    res_pooling[[i]]$power[j] <- (quantile(pooling_sample$lambda, 0.025) > 0 | 
                                    quantile(pooling_sample$lambda, 0.975) < 0)
    
    ##Progress
    if (j%%50 == 0) {
      print(paste("scenario:", i, ", simulation:", j)) 
    }
  }
}
Sys.time() - start

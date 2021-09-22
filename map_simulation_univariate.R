rm(list = ls())

##Simulation setting
source("simulation_settings.R")

##Data sets to save the results
res_map_univariate_common <- list()
res_mac_univariate_common <- list()
res_map_univariate_separate <- list()
res_mac_univariate_separate <- list()

start <- Sys.time()
for (i in 1:num_scenarios) {
  ##Generate seeds for data generation process
  set.seed(seed_scenarios[i])
  seed_datagen <- sample.int(.Machine$integer.max, size = num_simulations)
  ##Generate seeds for JAGS model
  set.seed(seed_scenarios[i])
  seed_simulations <- sample.int(.Machine$integer.max, size = num_simulations * num_chains)
  
  ##Save the results per scenario
  res_map_univariate_common[[i]] <- data.frame(simulation = 1:num_simulations, beta0 = NA,
                                   sigma_beta0 = NA, sd_beta0 = NA)
  res_mac_univariate_common[[i]] <- data.frame(simulation = 1:num_simulations, bias = NA, sd = NA, mse = NA, power = NA)
  res_map_univariate_separate[[i]] <- data.frame(simulation = 1:num_simulations, beta0 = NA,
                                               sigma_beta0 = NA, sd_beta0 = NA)
  res_mac_univariate_separate[[i]] <- data.frame(simulation = 1:num_simulations, bias = NA, sd = NA, mse = NA, power = NA)
  
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
                          response_new = data_all$data_current$followup)
    ##Precision of HN priors for between-study standard deviations of the intercept
    jags_data_map$prec_sigma_beta0 <- 1/(sd_sigma_beta0^2)
    
    ##Univariate, common baseline effect
    ##MAP
    map_jags_model_univariate_common <- jags.model(file = "./jags/map_jags_model_univariate_common.txt", data = jags_data_map, 
                                       n.chains = num_chains, n.adapt = 1000, quiet = T,
                                       inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains]))
    )
    
    ##Burnin stage
    update(map_jags_model_univariate_common, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta0_new", "sigma_beta0")
    
    map_jags_univariate_common <- coda.samples(map_jags_model_univariate_common, params, n.iter = num_iter * (1 - perc_burnin), 
                                   progress.bar = "none")
    
    map_sample_univariate_common <- as.data.frame(do.call(rbind, map_jags_univariate_common))
    
    ##Operating characteristics
    res_map_univariate_common[[i]]$beta0[j] <- mean(map_sample_univariate_common$'beta0_new')
    res_map_univariate_common[[i]]$sigma_beta0[j] <- mean(map_sample_univariate_common$sigma_beta0)
    res_map_univariate_common[[i]]$sd_beta0[j] <- sd(map_sample_univariate_common$'beta0_new')
    
    ##MAC
    mac_jags_model_univariate_common <- jags.model(file = "./jags/mac_jags_model_univariate_common.txt", data = jags_data_map, 
                                       n.chains = num_chains, n.adapt = 1000, quiet = T,
                                       inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                    list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains])))
    ##Burnin stage
    update(mac_jags_model_univariate_common, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta0_new", "lambda")
    
    mac_jags_univariate_common <- coda.samples(mac_jags_model_univariate_common, params, n.iter = num_iter * (1 - perc_burnin), progress.bar = "none")
    
    mac_sample_univariate_common <- as.data.frame(do.call(rbind, mac_jags_univariate_common))
    
    ##Operating characteristics
    res_mac_univariate_common[[i]]$bias[j] <- mean(mac_sample_univariate_common$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i]
    res_mac_univariate_common[[i]]$sd[j] <- sd(mac_sample_univariate_common$lambda)
    res_mac_univariate_common[[i]]$mse[j] <- (mean(mac_sample_univariate_common$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i])^2
    res_mac_univariate_common[[i]]$power[j] <- (quantile(mac_sample_univariate_common$lambda, 0.025) > 0 | 
                                      quantile(mac_sample_univariate_common$lambda, 0.975) < 0)
    
    ##Univariate, separate baseline effects
    ##MAP
    map_jags_model_univariate_separate <- jags.model(file = "./jags/map_jags_model_univariate_separate.txt", data = jags_data_map, 
                                                   n.chains = num_chains, n.adapt = 1000, quiet = T,
                                                   inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains]))
    )
    
    ##Burnin stage
    update(map_jags_model_univariate_separate, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta0_new", "sigma_beta0")
    
    map_jags_univariate_separate <- coda.samples(map_jags_model_univariate_separate, params, n.iter = num_iter * (1 - perc_burnin), 
                                               progress.bar = "none")
    
    map_sample_univariate_separate <- as.data.frame(do.call(rbind, map_jags_univariate_separate))
    
    ##Operating characteristics
    res_map_univariate_separate[[i]]$beta0[j] <- mean(map_sample_univariate_separate$'beta0_new')
    res_map_univariate_separate[[i]]$sigma_beta0[j] <- mean(map_sample_univariate_separate$sigma_beta0)
    res_map_univariate_separate[[i]]$sd_beta0[j] <- sd(map_sample_univariate_separate$'beta0_new')
    
    ##MAC
    mac_jags_model_univariate_separate <- jags.model(file = "./jags/mac_jags_model_univariate_separate.txt", data = jags_data_map, 
                                                   n.chains = num_chains, n.adapt = 1000, quiet = T,
                                                   inits = list(list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 1]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 2]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[(j - 1) * num_chains + 3]),
                                                                list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = seed_simulations[j * num_chains])))
    ##Burnin stage
    update(mac_jags_model_univariate_separate, n.iter = num_iter * perc_burnin, progress.bar = "none")
    
    ##Sampling after burnin
    params <- c("beta0_new", "lambda")
    
    mac_jags_univariate_separate <- coda.samples(mac_jags_model_univariate_separate, params, n.iter = num_iter * (1 - perc_burnin), progress.bar = "none")
    
    mac_sample_univariate_separate <- as.data.frame(do.call(rbind, mac_jags_univariate_separate))
    
    ##Operating characteristics
    res_mac_univariate_separate[[i]]$bias[j] <- mean(mac_sample_univariate_separate$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i]
    res_mac_univariate_separate[[i]]$sd[j] <- sd(mac_sample_univariate_separate$lambda)
    res_mac_univariate_separate[[i]]$mse[j] <- (mean(mac_sample_univariate_separate$lambda) - trt_eff * simulation_scenarios$with_trt_eff[i])^2
    res_mac_univariate_separate[[i]]$power[j] <- (quantile(mac_sample_univariate_separate$lambda, 0.025) > 0 | 
                                                  quantile(mac_sample_univariate_separate$lambda, 0.975) < 0)
    
    ##Progress
    if (j%%50 == 0) {
      print(paste("scenario:", i, ", simulation:", j)) 
    }
  }
}
Sys.time() - start


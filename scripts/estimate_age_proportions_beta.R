#code by Patrick Thompson

library(tidyverse)
library(rstan)
library(tidybayes)

theme_set(theme_bw())

#simulate fake data####
Y <- 30 # years
phi <- 0.1
mu_theta <- 2
sigma_theta <- 0.5
mu_M <- -1.5
sigma_M <- 0.2

N1_mu <- 1e4
N1_phi <- 10
init_N2 <- 500

N1 <- rnbinom(n = Y, mu = N1_mu, size = N1_phi)
N2 <- rep(NA, Y)
N2[1] <- init_N2

theta <- plogis(rnorm(Y, mu_theta, sigma_theta))
hist(theta)

M <- plogis(rnorm(Y, mu_M, sigma_M))

N_lake <- O1 <- O2 <- rep(NA, Y)
for(y in 1:Y){
  if(y > 1){
  N2[y] = N1[y-1] * (1 - theta[y-1]) * M[y-1]
  }
  N_lake[y] = N1[y] + N2[y]
  
  # Outmigrants
  O1[y] = theta[y] * N1[y]
  O2[y] = N2[y]
}

obs_error <- 0.001
N_obs <- rnbinom(Y, mu = N_lake, 1/obs_error)


A_total <- rbinom(n = Y, size = N_obs, prob = 0.02)
A1_obs <- rbinom(n = Y, size = A_total, prob = O1 / (O1 + O2))

plot(N_obs ~ N_lake)

dat_list <- list(Y = Y, 
                 N_obs = N_obs, 
                 A1_obs = A1_obs, 
                 A_total = A_total, 
                 obs_error_prior = 1/0.001,
                 mu_theta_prior = 2,
                 sigma_theta_prior = 0.5,
                 mu_M_prior = -1.5,
                 sigma_M_prior = 0.5,
                 N2_init_prior = 300,
                 N2_init_sigma_prior = 0.2
                 )
  

#run stan model####
model <- stan_model("./stan/age_prop_beta.stan")
fit <- sampling(model, data = dat_list, chains = 4, cores = 4, iter = 1000)

#assess convergence###
worst_Rhat <- summary(fit)$summary %>% 
  as.data.frame() %>% 
  mutate(Rhat = round(Rhat, 3)) %>% 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

head(worst_Rhat, n = 20)

traceplot(fit, pars = rownames(worst_Rhat)[1:20])
pairs(fit, pars = c("mu_theta", "mu_M"))

#explore posterior####
post <- extract(fit)

plot(density(post$mu_M))
abline(v = mu_M)

plot(density(post$mu_theta))
abline(v = mu_theta)

plot(density(post$init_N2))
abline(v = init_N2)

spread_draws(fit, O1[year]) %>% 
  ggplot(aes(x = year, y = O1)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, O1 = O1), color = "red")

spread_draws(fit, O2[year]) %>% 
  ggplot(aes(x = year, y = O2)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, O2 = O2), color = "red")

spread_draws(fit, theta[year]) %>% 
  ggplot(aes(x = year, y = theta)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, theta = theta), color = "red")

spread_draws(fit, M[year]) %>% 
  ggplot(aes(x = year, y = M)) +
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = data.frame(year = 1:Y, M = M), color = "red")
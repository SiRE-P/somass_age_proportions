#code by Patrick Thompson
#patrick.thompson@dfo-mpo.gc.ca

#load packages####
library(rstan)
library(tidybayes)
library(tidyverse)

theme_set(theme_bw())

#load data####
age_data <- readxl::read_excel("./data/Annual Smolt Age Size Comps From Nick 250212 1.xlsx")

age_data$n <- as.numeric(age_data$n)
age_data$age1_pct[is.na(age_data$n)] <- NA

ggplot(age_data, aes(x = smolt_year, y = age1_pct, color = is.na(n)))+
  geom_point()+
  facet_wrap(~lake)

ggplot(age_data, aes(y = age1_pct, x = n))+
  geom_point()+
  facet_wrap(~lake, scales = "free_x")+
  scale_color_viridis_c()

age_data %>% 
  mutate(lake_f = as.numeric(factor(lake))) %>% 
  mutate(age1_pct_adg = ifelse(age1_pct == 1, age1_pct-0.01, age1_pct)) %>% 
  mutate(logit_prop = log(age1_pct_adg/(1-age1_pct_adg))) %>% 
  group_by(lake_f) %>% 
  summarise(age1_mean = mean(logit_prop, na.rm = TRUE), age1_sd = sd(logit_prop, na.rm = TRUE))

dat_list <- list(N = nrow(age_data),
  L = length(unique(age_data$lake)),
  lake = as.numeric(factor(age_data$lake)), 
  observed = as.numeric(!is.na(age_data$n)),
  age_1s = ifelse(is.na(age_data$n), -999, round(age_data$age1_pct * age_data$n)),
  sample_size = ifelse(is.na(age_data$n), -999, age_data$n),
  mu_prior = c(2,5,4)) # did this because the model was having trouble converging with a single mu prior - base this on knowledge of the lake - prior is on logit scale
  
#run stan model####
model <- stan_model("./stan/age_prop.stan")
fit <- sampling(model, data = dat_list, chains = 4, cores = 4)

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

pairs(fit, pars = c("mu[1]", "p_z[1]", "p_z[2]", "p_z[3]", "p_z[4]", "p_z[5]"))
pairs(fit, pars = c("mu[3]", "p_z[70]", "p_z[71]", "p_z[72]"))
pairs(fit, pars = c("mu[2]", "p_z[89]", "p_z[90]", "p_z[91]"))

#plot results
post <- rstan::extract(fit)

spread_draws(fit, p_sigma[lake_n]) %>% 
  left_join(data.frame(lake_n = 1:3, lake = levels(factor(age_data$lake)))) %>% 
  ggplot(aes(x = lake, y = p_sigma))+
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))

spread_draws(fit, mu[lake_n]) %>% 
  left_join(data.frame(lake_n = 1:3, lake = levels(factor(age_data$lake)))) %>% 
  ggplot(aes(x = lake, y = plogis(mu)))+
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))

spread_draws(fit, p[obs]) %>% 
  left_join(age_data %>% mutate(obs = 1:n())) %>% 
  ggplot(aes(x = n, y = p))+
  stat_pointinterval(.width = c(0.5, 0.8, 0.95))+
  geom_point(data = age_data, aes(y = age1_pct), color = "red")+
  scale_x_log10()+
  facet_wrap(~lake)+
  geom_hline(data = data.frame(lake = levels(factor(age_data$lake)), mu = apply(plogis(post$mu), 2, median)), aes(yintercept = mu), lty = 2)+
  ylab("proportion age 1")+
  xlab("sample size")

spread_draws(fit, p[obs]) %>% 
  left_join(age_data %>% mutate(obs = 1:n())) %>% 
  ggplot(aes(x = smolt_year, y = p))+
  stat_pointinterval(.width = c(0.5, 0.8, 0.95), aes(color = !is.na(n)))+
  scale_color_manual(values = c("grey30", "grey60"), "samples")+
  geom_point(data = age_data, aes(y = age1_pct), color = "red")+
  facet_wrap(~lake)+
  geom_hline(data = data.frame(lake = levels(factor(age_data$lake)), mu = apply(plogis(post$mu), 2, median)), aes(yintercept = mu), lty = 2)+
  ylab("proportion age 1")

spread_draws(fit, p_z[obs]) %>% 
  left_join(age_data %>% mutate(obs = 1:n())) %>% 
  ggplot(aes(x = smolt_year, y = p_z))+
  stat_pointinterval(.width = c(0.5, 0.8, 0.95), aes(color = !is.na(n)))+
  scale_color_manual(values = c("grey30", "grey60"), "samples")+
  facet_wrap(~lake)+
  ylab("proportion age 1")


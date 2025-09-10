## meta-regressions on parcel-wise gPPI results
# testing the hypothesis that DMN > Ctrl=LIM and ITI alignment is superior


pacman::p_load(parallel, tidyverse, readxl, glue, data.table, brms, tidybayes, extrafont)


load("/Users/alexdombrovski/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/collected_letters/papers/sceptic_fmri/vmPFC/gppi/2025-09-09-gPPI-l2-betas.Rdata")
d <- gPPI_l2_betas %>% filter(grepl("ppi", l1_cope_name)) %>% 
  group_by(event) %>% mutate(beta = beta*100,
                             se_z = abs(beta/zstat),
                             sd = sqrt(betavar)*100,
                             se = sd/sqrt(8)) %>% select(id, session, l1_cope_name, atlas_value, beta, sd, se, se_z, zstat, event, network) %>%
  ungroup()
str(d)

sdf <- d %>% group_by(network, event, atlas_value) %>% summarize(z = mean(zstat))
setwd("/Users/alexdombrovski/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Documents/collected_letters/papers/sceptic_fmri/vmPFC/gppi/")
pdf("betas_by_network_event_across_parcels.pdf", height = 3, width = 6)
ggplot(sdf, aes(event, z, color = network)) + geom_violin(draw_quantiles = 0.5) + 
  geom_point(alpha = 0.6, position = position_dodge2(width = 0.8)) + theme_minimal() +
  scale_color_brewer(palette = "Dark2") + geom_hline(yintercept = 0)
dev.off()


m1 <- brm(beta|se(se) ~ 1 + network * event + (1|atlas_value) + (1|id), 
          data = d,
          # prior = priors,
          iter = 4000,
          silent = F,
          control = list(
            adapt_delta = 0.8,  # Higher adapt_delta for complex models
            max_treedepth = 10   # Allow deeper trees if needed
          ),
          chains = 4, cores = 4, backend = "cmdstanr", threads = threading(2))
summary(m1)
str(m1$data)

newdata <- expand.grid(network = levels(as.factor(m1$data$network)), 
                       event = levels(as.factor(m1$data$event))) 
newdata$se <- 0  # Use a placeholder value (not used in prediction)
newdata$atlas_value <- 55


# Show the newdata used for prediction
# print(newdata)

# Sample from the posterior for each TMS x network combination
posterior_samples <- posterior_linpred(m1, newdata = newdata, re_formula = NA, transform = TRUE)


# 2. Prepare the data for plotting
plot_data <- newdata %>%
  add_epred_draws(m1, re_formula = NA) # This is the modern tidybayes equivalent


# 3. Create the plot
pdf("gppi_network_event.pdf", height = 4, width = 6)
ggplot(plot_data, aes(x = event, y = .epred, fill = network)) +
  scale_fill_brewer(palette = "Dark2") + 
  stat_halfeye(p_limits = c(0.05, 0.95), orientation = "vertical", 
               position = position_dodge(width = 0.75), alpha = 0.8) +
  labs(x = 'Event', , 
       y = glue::glue("PPI: entropy * hippocampal connectivity\nposterior estimate"),
       fill = 'Network') + 
  geom_hline(yintercept = 0, lty = "dashed") + #facet_grid(period ~ .) 
  theme_minimal() 
dev.off()
  # theme(legend.position = "bottom")
# Reshape to tidy format for plotting

# 
# 
# m1 <- lme4::lmer(value ~ network * event + (1 | id) + (1 | atlas_value), d)
# summary(m1)
# em <- emmeans::emmeans(m1, ~ network * event) %>% as_tibble()
# ggplot(em, aes(event, emmean, color = network)) + geom_point(position = position_dodge2(width = 0.9)) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width = 0.2)) 
# 
# m1_sc <- lme4::lmer(value_sc ~ network * event + (1 | id) + (1 | atlas_value), d)
# summary(m1_sc)
# em <- emmeans::emmeans(m1_sc, ~ network * event) %>% as_tibble()
# ggplot(em, aes(event, emmean, color = network)) + geom_point(position = position_dodge2(width = 0.9)) + 
#   geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width = 0.2)) 
# 
# 
# 
# m2 <- lme4::lmer(value ~ l1_cope_name * network + (1 | id) + (1 | atlas_value), d)
# summary(m2)
# 

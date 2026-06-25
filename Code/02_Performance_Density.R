library(tidyverse)
library(magrittr)
library(readxl)
library(broom)
library(rsample)
library(ggh4x)
library(here)

##########################################################################
##########################################################################
### Analysis -- Model fitting
##########################################################################
##########################################################################

##############
### 0. Load data
##############

load(file = "Data\\CleanData\\performance.rdata")


##############
### 1. Model selection
##############

################################
### Flexible intercepts
### focal performance v.s. competitor density
### performance = log(tot_focal) - log(tot_focal_initial)
################################

cur.data <- dat %>%
  #filter(Remove_density == 0) %>% # Not to remove samples, follow the original experimental setting
  mutate(performance = log(tot_focal) - log(tot_focal_initial),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  select(focal, competitor, density, performance, decay, treatment, inoculum, competition, `focal x comp`)


cur.model.1.sum <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ tidy(lm(performance ~ density, .x)))

cur.model.2.sum <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ tidy(lm(performance ~ density + I(density^2), .x)))

cur.model.compare.F.21 <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ anova(lm(performance ~ density, .x), 
                       lm(performance ~ density + I(density^2), .x)))
# Scenarios with quadratic terms suggested by F-test
## ENGERO (focal) + MACHZU (comp) + dead + sterilized (P-value = 0.04432038)
## MACHZU (focal) + MACHZU (comp) + living + unsterilized (P-value = 0.03623963)

cur.model.compare.F.10 <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ anova(lm(performance ~ 1, .x), 
                       lm(performance ~ density, .x)))
# Significant density dependent effect suggested by F-test
## ENGERO (focal) + ENGERO (comp) + dead + unsterilized (P-value = 0.001222282)

cur.model.compare.F.20 <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ anova(lm(performance ~ 1, .x), 
                       lm(performance ~ density, .x),
                       lm(performance ~ density + I(density^2), .x)))

cur.model.compare.BIC <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ .x %>% 
                 mutate(BIC.1 = BIC(lm(performance ~ density, .x))) %>%
                 mutate(BIC.2 = BIC(lm(performance ~ density + I(density^2), .x))) %>%
                 select(BIC.1, BIC.2) %>%
                 unique %>%
                 mutate(order = case_when(
                   BIC.1 < BIC.2 ~ 1,
                   BIC.1 > BIC.2 ~ 2,
                   TRUE ~ 0
                 )))
# Scenarios with quadratic terms suggested by lower BIC
## ENGERO (focal) + MACHZU (comp) + dead + sterilized
## MACHZU (focal) + MACHZU (comp) + living + unsterilized

cur.model.compare.AIC <- cur.data %>% 
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ .x %>% 
                 mutate(AIC.1 = AIC(lm(performance ~ density, .x))) %>%
                 mutate(AIC.2 = AIC(lm(performance ~ density + I(density^2), .x))) %>%
                 select(AIC.1, AIC.2) %>%
                 unique %>%
                 mutate(order = case_when(
                   AIC.1 < AIC.2 ~ 1,
                   AIC.1 > AIC.2 ~ 2,
                   TRUE ~ 0
                 )))
# Scenarios with quadratic terms suggested by lower AIC
## ENGERO (focal) + MACHZU (comp) + dead + sterilized
## MACHZU (focal) + MACHZU (comp) + living + unsterilized

# Estimated quadratic term
## ENGERO (focal) + MACHZU (comp) + dead + sterilized (quadratic term = -0.224979201) -> negative
## MACHZU (focal) + MACHZU (comp) + living + unsterilized (quadratic term = 0.105957800) -> positive (not consider positive quadratic term)

# setwd("Data\\Model Selection")
# write_csv(cur.model.compare.F.10, file = "density-competition_F-test_1_0.csv")
# write_csv(cur.model.compare.F.21, file = "density-competition_F-test_2_1.csv")
# write_csv(cur.model.compare.F.20, file = "density-competition_F-test_2_0.csv")
# write_csv(cur.model.compare.AIC, file = "density-competition_AIC.csv")
# write_csv(cur.model.compare.BIC, file = "density-competition_BIC.csv")
# setwd(here())

##############
### Model selection -- Some quadratic (based on F-test, BIC, AIC)
##############

# Note F-test, BIC, and AIC all suggest the same quadratic terms. We use BIC table to select the quadratic fitted model. 
# Note that the positive quadratic term will not be considered.

remove_positive <- cur.model.compare.BIC %>% 
  inner_join(cur.model.2.sum) %>%
  filter(order == 2, term == "I(density^2)" & estimate > 0) %>%
  select(focal, competitor, decay, treatment)

ind_row <- (cur.model.compare.BIC$focal == remove_positive$focal)&
  (cur.model.compare.BIC$competitor == remove_positive$competitor)&
  (cur.model.compare.BIC$decay == remove_positive$decay)&
  (cur.model.compare.BIC$treatment == remove_positive$treatment)

cur.model.compare.BIC$order[ind_row] <- 1
  

cur.model.sum <- cur.model.compare.BIC %>% 
  filter(order == 2) %>% 
  select(focal, competitor, decay, treatment) %>%
  inner_join(cur.model.2.sum)

cur.model.sum <- cur.model.compare.BIC %>% 
  filter(order == 1) %>% 
  select(focal, competitor, decay, treatment) %>%
  inner_join(cur.model.1.sum) %>%
  rbind(cur.model.sum)


coef <- cur.model.sum %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "a",
    term == "density" ~ "b1",
    term == "I(density^2)" ~ "b2",
    TRUE ~ "error"
  )) %>%
  select(1:6) %>%
  group_by(focal, competitor, decay, treatment) %>%
  group_modify(~ spread(.x, key = term, value = estimate)) %>%
  mutate(competition = if_else(focal == competitor, 'conspecific', 'heterospecific'),
         inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  mutate(b2 = ifelse(is.na(b2), 0, b2)) %>%
  ungroup()


# save(coef, file = "Data\\CleanData\\density_coef.rdata")

### Plotting
strip_colors <- c("ENGERO x ENGERO" = "#1984c5", "MACHZU x MACHZU" = "#c23728", "MACHZU x ENGERO" = "#e1a692", "ENGERO x MACHZU" = "#a7d5ed")

cur.pred <- with(cur.data, expand.grid(
  focal = unique(focal),
  competitor = unique(competitor),
  decay = unique(decay),
  treatment = unique(treatment),
  density = seq(0, 4.5, length.out = 1000)
)) %>%
  left_join(coef) %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead")),
         `focal x comp` = factor(`focal x comp`, levels = c("ENGERO x ENGERO", "MACHZU x MACHZU", "MACHZU x ENGERO", "ENGERO x MACHZU"))) %>% 
  mutate(performance = a + b1*density + b2*density^2)

cur.data %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead")),
         `focal x comp` = factor(`focal x comp`, levels = c("ENGERO x ENGERO", "MACHZU x MACHZU", "MACHZU x ENGERO", "ENGERO x MACHZU"))) %>%
  ggplot(aes(x = density, y = performance, color = inoculum)) +
  geom_point(size = 3) +
  geom_line(data = cur.pred,
            mapping = aes(x = density, y = performance, color = inoculum), size = 2) +
  facet_wrap2(~ `focal x comp`, nrow = 2) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  scale_color_manual(values = c("#548235", "#B3D89C", "#A17B01", "#FFCF37")) +
  scale_y_continuous(limits = c(-1,6)) +
  labs(#title = "performance ~ competitor density",
       x = "number of competitors",
       y = "focal species growth")

ggsave(filename = "Data\\Figures\\density-competition_model-fit.pdf", width = 7, height = 7)
ggsave(filename = "Data\\Figures\\density-competition_model-fit.jpg", width = 7, height = 7)


########################
### 2. Invasion analysis
########################

monoculture <- coef %>% 
  filter(competition == "conspecific") %>%
  mutate(N_star = case_when(
    b2 > 0 ~ Inf,
    b2 < 0 ~ 1/(2*b2)*(-b1 + sqrt(b1^2 - 4*b2*a)),
    b1 < 0 ~ -a/b1,
    TRUE ~ Inf
  )) %>%
  select(competitor, treatment, decay, N_star) # focal = competitor

invasion <- coef %>% 
  filter(competition == "heterospecific") %>%
  left_join(monoculture) %>%
  mutate(IGR = case_when(
    b2 == 0 ~ a + b1*N_star,
    TRUE ~ a + b1*N_star + b2*N_star^2
  )) %>%
  select(focal, competitor, treatment, decay, IGR)


comp_outcome <- invasion %>%
  group_by(treatment, decay) %>%
  group_modify(~ .x %>% 
                 mutate(outcome = case_when(
                   sum(.x$IGR > 0) == 2 ~ "coexist",
                   .x$IGR[.x$focal == "MACHZU"] > 0 ~ "M win",
                   .x$IGR[.x$focal == "ENGERO"] > 0 ~ "E win",
                   TRUE ~ "priority effect"
                 )) %>% 
                 select(outcome) %>%
                 unique
  )


##############
### 3. Wrap coef as a function
##############

COEF <- function(cur.data, cur.model.compare.obs){
  
  cur.model.1.sum <- cur.data %>% 
    group_by(focal, competitor, decay, treatment) %>%
    group_modify(~ tidy(lm(performance ~ density, .x)))
  
  cur.model.2.sum <- cur.data %>% 
    group_by(focal, competitor, decay, treatment) %>%
    group_modify(~ tidy(lm(performance ~ density + I(density^2), .x)))
  
  cur.model.sum <- cur.model.compare.obs %>% 
    filter(order == 2) %>% 
    select(focal, competitor, decay, treatment) %>%
    inner_join(cur.model.2.sum)
  
  cur.model.sum <- cur.model.compare.obs %>% 
    filter(order == 1) %>% 
    select(focal, competitor, decay, treatment) %>%
    inner_join(cur.model.1.sum) %>%
    rbind(cur.model.sum)
  
  
  coef <- cur.model.sum %>%
    mutate(term = case_when(
      term == "(Intercept)" ~ "a",
      term == "density" ~ "b1",
      term == "I(density^2)" ~ "b2",
      TRUE ~ "error"
    )) %>%
    select(1:6) %>%
    group_by(focal, competitor, decay, treatment) %>%
    group_modify(~ spread(.x, key = term, value = estimate)) %>%
    mutate(competition = if_else(focal == competitor, 'conspecific', 'heterospecific'),
           inoculum = paste0(treatment, "+", decay)) %>%
    mutate(b2 = ifelse(is.na(b2), 0, b2)) %>%
    ungroup()
  
  
  return(coef)
}

##########################################################################
##########################################################################
### Bootstrapping
##########################################################################
##########################################################################

########################
### 1. Bootstrapping
########################

boot.replicates <- 10000
boot.metrics <- NULL
boot.metrics.unnest <- NULL

boot.metrics <- cur.data %>%
  mutate(strata = paste0(focal, competitor, inoculum)) %>%
  # Stratified bootstrap
  # (the warning is a over-cautious rule of thumb from the 'rsample' library;
  # it can be ignored)
  bootstraps(times = boot.replicates, strata = strata, pool = 0, apparent = TRUE) %>%
  # Calculate all metrics for the outcomes
  mutate(metrics = map(
    splits, function(split) COEF(analysis(split), 
                                 cur.model.compare.obs = cur.model.compare.BIC)))

# Clean this up a bit so it's easier to work with for further analysis
boot.metrics.unnest <- boot.metrics %>%
  # Unnest: restore one column per metric and make four rows for each bootstrap 
  # replicate, one per treatment
  unnest(metrics) %>% 
  select(-c(splits))

#save(boot.metrics.unnest, file = "Data\\Bootstrap\\density-competition_boot-metrics-unnest.rdata")




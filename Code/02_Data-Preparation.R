library(tidyverse)
library(readxl)

##########################################################################
##########################################################################
### Data preparation
##########################################################################
##########################################################################

##############
### 1. Load data
##############


raw_dat <- read_xlsx(path = "Data\\RawData\\harvest-dry weight.xlsx") 

dat <- raw_dat %>%
  select(1:14, Remove_density)
dat[is.na(dat)] <- 0


initial_raw_dat <- read_xlsx(path = "Data\\RawData\\initial biomass.xlsx")

##############
### 2. Partition the broken roots weight into focal and competitors
##############

n_row <- nrow(dat)
for(ind in 1:n_row){
  if(dat$Below_E_broke[ind] > 0){  # if it has broken roots
    focal <- dat$Above_focal[ind]
    comp <- dat$Above_comp[ind]
    
    w_f <- focal/(focal + comp)
    w_c <- comp/(focal + comp)
    
    dat$Below_E_focal[ind] <- dat$Below_E_focal[ind] + w_f*dat$Below_E_broke[ind]
    dat$Below_E_comp[ind] <- dat$Below_E_comp[ind] + w_c*dat$Below_E_broke[ind]
  }
  if(dat$Below_M_broke[ind] > 0){  # if it has broken roots
    focal <- dat$Above_focal[ind]
    comp <- dat$Above_comp[ind]
    
    w_f <- focal/(focal + comp)
    w_c <- comp/(focal + comp)
    
    dat$Below_M_focal[ind] <- dat$Below_M_focal[ind] + w_f*dat$Below_M_broke[ind]
    dat$Below_M_comp[ind] <- dat$Below_M_comp[ind] + w_c*dat$Below_M_broke[ind]
  }
}

##############
### 3. Mutate and summarize variable for analysis
##############


dat <- dat %>% 
  mutate(competition = case_when(
    focal == competitor ~ "conspecific",
    TRUE ~ "heterospecific"
  )) %>%
  mutate(Below_focal = case_when(
    focal == "MACHZU" ~ Below_M_focal,
    focal == "ENGERO" ~ Below_E_focal,
    TRUE ~ -1
  )) %>%
  mutate(Below_comp = case_when(
    competitor == "MACHZU" ~ Below_M_comp,
    competitor == "ENGERO" ~ Below_E_comp,
    TRUE ~ -1
  )) %>%
  mutate(tot_focal = Above_focal + Below_focal,
         tot_comp = Above_comp + Below_comp) 

dat <- dat %>%
  mutate(inoculum = paste0(treatment, "+", decay)) %>% 
  select(ID, focal, competitor, density, competition, inoculum, decay, treatment,
         tot_focal, tot_comp, Above_focal, Above_comp, Below_focal, Below_comp, Remove_density)


initial_dat <- initial_raw_dat %>% 
  select(-Tag) %>%
  group_by(Species) %>% 
  summarise_all(mean) %>%
  mutate(above = above_dry_bm_g,
         below = below_dry_bm_g,
         tot = above_dry_bm_g + below_dry_bm_g) %>%
  select(Species, above, below, tot)

dat <- dat %>%
  mutate(Above_focal_initial = case_when(
    focal == "MACHZU" ~ initial_dat$above[2],
    focal == "ENGERO" ~ initial_dat$above[1],
    TRUE ~ -1
  )) %>%
  mutate(Below_focal_initial = case_when(
    focal == "MACHZU" ~ initial_dat$below[2],
    focal == "ENGERO" ~ initial_dat$below[1],
    TRUE ~ -1
  )) %>%
  mutate(tot_focal_initial = case_when(
    focal == "MACHZU" ~ initial_dat$tot[2],
    focal == "ENGERO" ~ initial_dat$tot[1],
    TRUE ~ -1
  )) %>%
  mutate(Above_comp_initial = case_when(
    competitor == "MACHZU" ~ initial_dat$above[2],
    competitor == "ENGERO" ~ initial_dat$above[1],
    TRUE ~ -1
  )) %>%
  mutate(Below_comp_initial = case_when(
    competitor == "MACHZU" ~ initial_dat$below[2],
    competitor == "ENGERO" ~ initial_dat$below[1],
    TRUE ~ -1
  )) %>%
  mutate(tot_comp_initial = case_when(
    competitor == "MACHZU" ~ initial_dat$tot[2],
    competitor == "ENGERO" ~ initial_dat$tot[1],
    TRUE ~ -1
  ))

dat <- dat %>% 
  left_join(moss) %>%
  ungroup()

##############
### 4. Save prepared data
##############

# write.csv(dat, file = "Data\\CleanData\\performance.csv")
# save(dat, file = "Data\\CleanData\\performance.rdata")

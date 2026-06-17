library(tidyverse)
library(gt)
library(here)

########################
### 1. Load data
########################

setwd(here())
setwd("Data\\Bootstrap")
boot_a <- read.csv(file = "biomass-competition_boot_a.csv", row.names = 1)
boot_b1 <- read.csv(file = "biomass-competition_boot_b1.csv", row.names = 1)
boot_b2 <- read.csv(file = "biomass-competition_boot_b2.csv", row.names = 1)
boot_comp <- read.csv(file = "biomass-competition_boot_comp_outcome.csv", row.names = 1)
boot_E_persist <- read.csv(file = "biomass-competition_boot_E_persist.csv", row.names = 1) %>%
  select(-label)
boot_M_persist <- read.csv(file = "biomass-competition_boot_M_persist.csv", row.names = 1) %>%
  select(-label)

setwd(here())
setwd("Data\\Model Selection")
AIC <- read.csv(file = "biomass-competition_AIC.csv") %>% select(-order)
BIC <- read.csv(file = "biomass-competition_BIC.csv") %>% select(-order)
AIC_d <- read.csv(file = "density-competition_AIC.csv") %>% select(-order)
BIC_d <- read.csv(file = "density-competition_BIC.csv") %>% select(-order)
F.1_0 <- read.csv(file = "biomass-competition_F-test.1_0.csv")
F.2_1 <- read.csv(file = "biomass-competition_F-test.2_1.csv")
F.2_0 <- read.csv(file = "biomass-competition_F-test.2_0.csv")
setwd(here())

########################
### 2. Build tables
########################

# boot_a
tb <- boot_a %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  rename(`focal x comp` = focal.x.comp) %>%
  select(inoculum, `focal x comp`, estimate, q.05, q.95) %>%
  gt

tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_a.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_a.png")

# boot_b
tb1 <- boot_b1 %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  rename(`focal x comp` = focal.x.comp) %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  select(inoculum, `focal x comp`, estimate, q.05, q.95)

tb2 <- boot_b2 %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  rename(`focal x comp` = focal.x.comp) %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  filter(across(where(is.numeric), function(x) x != 0)) %>%
  select(inoculum, `focal x comp`, estimate, q.05, q.95)

tb <- rbind(tb1, tb2) %>%
  gt %>%
  tab_row_group(
    label = "2nd order",
    rows = 17
  )%>%
  tab_row_group(
    label = "1st order",
    rows = 1:16
  ) 
  

tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_b.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_b.png")


# boot_comp
tb <- boot_comp %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  gt

tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_comp.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_comp.png")

# persistence
tb <- boot_E_persist %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead")),
         IGR = factor(IGR, levels = c("IGR > 0", "IGR < 0"))) %>%
  as.tibble %>%
  arrange(inoculum, IGR) %>%
  rename("ENGERO" = "ratio") %>%
  left_join(boot_M_persist, by = c("inoculum", "IGR")) %>%
  rename("MACHZU" = "ratio") %>%
  gt %>%
  tab_spanner(
    label = "  ",
    columns = 1:2
  ) %>%
  tab_spanner(
    label = "ratio",
    columns = 3:4
  )
tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_persist.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_boot_persist.png")

# AIC and BIC
tb <- AIC %>% left_join(BIC) %>%
  select(focal, competitor, decay, treatment, AIC.1_ = AIC.1, AIC.2_ = AIC.2, BIC.1_ = BIC.1, BIC.2_ = BIC.2) %>%
  left_join(AIC_d) %>%
  left_join(BIC_d) %>%
  mutate(inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  select(inoculum, `focal x comp`, AIC.1_, AIC.2_, BIC.1_, BIC.2_, AIC.1, AIC.2, BIC.1, BIC.2) %>%
  gt %>%
  tab_spanner(
    label = "biomass",
    columns = 3:6
  ) %>% 
  tab_spanner(
    label = "density",
    columns = 7:10
  )

tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_AIC_BIC.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_AIC_BIC.png")


#F.2_0
tb <- F.2_0 %>% 
  mutate(inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  select(`focal x comp`, inoculum, 5:10) %>%
  rename("P value" = "Pr..F.") %>%
  mutate(across(where(is.numeric), round, 5)) %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  select(inoculum, `focal x comp`, Res.Df, RSS, Df, Sum.of.Sq, F, `P value`) %>%
  gt

tb
gtsave(tb, file = "Data\\Tables\\biomass-competition_F-test.2_0.tex")
gtsave(tb, file = "Data\\Tables\\biomass-competition_F-test.2_0.png")



library(tidyverse)
library(gt)
library(here)

########################
### 1. Load data
########################

setwd(here())
setwd("Data\\Bootstrap")
boot_a <- read.csv(file = "density-competition_boot_a.csv", row.names = 1)
boot_b1 <- read.csv(file = "density-competition_boot_b1.csv", row.names = 1)
boot_b2 <- read.csv(file = "density-competition_boot_b2.csv", row.names = 1)
boot_comp <- read.csv(file = "density-competition_boot_comp_outcome.csv", row.names = 1)
boot_E_persist <- read.csv(file = "density-competition_boot_E_persist.csv", row.names = 1) %>%
  select(-label)
boot_M_persist <- read.csv(file = "density-competition_boot_M_persist.csv", row.names = 1) %>%
  select(-label)

setwd(here())
setwd("Data\\Model Selection")
AIC <- read.csv(file = "density-competition_AIC.csv") %>% select(-order)
BIC <- read.csv(file = "density-competition_BIC.csv") %>% select(-order)
F.1_0 <- read.csv(file = "density-competition_F-test_1_0.csv")
F.2_1 <- read.csv(file = "density-competition_F-test_2_1.csv")
F.2_0 <- read.csv(file = "density-competition_F-test_2_0.csv")
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
gtsave(tb, file = "Data\\Tables\\density-competition_boot_a.tex")
gtsave(tb, file = "Data\\Tables\\density-competition_boot_a.png")

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
gtsave(tb, file = "Data\\Tables\\density-competition_boot_b.tex")
gtsave(tb, file = "Data\\Tables\\density-competition_boot_b.png")


# boot_comp
tb <- boot_comp %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  as_tibble %>%
  arrange(inoculum) %>%
  gt

tb
gtsave(tb, file = "Data\\Tables\\density-competition_boot_comp.tex")
gtsave(tb, file = "Data\\Tables\\density-competition_boot_comp.png")

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
gtsave(tb, file = "Data\\Tables\\density-competition_boot_persist.tex")
gtsave(tb, file = "Data\\Tables\\density-competition_boot_persist.png")

# AIC and BIC
# See 07_tables_Density.R


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
gtsave(tb, file = "Data\\Tables\\density-competition_F-test.2_0.tex")
gtsave(tb, file = "Data\\Tables\\density-competition_F-test.2_0.png")



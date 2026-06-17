library(tidyverse)
library(patchwork)
library(here)

########################
### 1. Load bootstrapping data
########################

setwd(here())
load("Data\\Bootstrap\\biomass-competition-08_boot-metrics-unnest.rdata")
load("Data\\CleanData\\biomass-08_coef.rdata")

########################
### 2. Coef distribution
########################

coef.a <- coef %>%
  select(`focal x comp`, inoculum, estimate = a)

boot.a <- boot.metrics.unnest %>%
  mutate(inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  select(id, `focal x comp`, inoculum, a)

summary_boot.a <- boot.a %>% 
  select(-id) %>% 
  group_by(`focal x comp`, inoculum) %>%
  summarise(q.05 = quantile(a, probs = 0.05), 
            q.95 = quantile(a, probs = 0.95))%>%
  left_join(coef.a)


coef.b1 <- coef %>%
  select(`focal x comp`, inoculum, estimate = b1)

boot.b1 <- boot.metrics.unnest %>%
  mutate(inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  select(id, `focal x comp`, inoculum, b1)

summary_boot.b1 <- boot.b1 %>% 
  select(-id) %>% 
  group_by(`focal x comp`, inoculum) %>%
  summarise(q.05 = quantile(b1, probs = 0.05), 
            q.95 = quantile(b1, probs = 0.95)) %>%
  left_join(coef.b1)

coef.b2 <- coef %>%
  select(`focal x comp`, inoculum, estimate = b2)

boot.b2 <- boot.metrics.unnest %>%
  mutate(inoculum = paste0(treatment, "+", decay),
         `focal x comp` = paste0(focal, " x ", competitor)) %>%
  select(id, `focal x comp`, inoculum, b2)

summary_boot.b2 <- boot.b2 %>% 
  select(-id) %>% 
  group_by(`focal x comp`, inoculum) %>%
  summarise(q.05 = quantile(b2, probs = 0.05), 
            q.95 = quantile(b2, probs = 0.95))%>%
  left_join(coef.b2)



### Plotting
p1 <- summary_boot.a %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  ggplot(aes(x = inoculum, y = estimate, fill = `focal x comp`)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = q.05, ymax = q.95), 
                width = 0.4, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#1984c5", "#a7d5ed", "#e1a692", "#c23728")) +
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  labs(y = "growth without competitor")

left_axis_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.0, label = " ", angle = 90, size = 6) +
  theme_void()

left_axis_label + p1+ plot_layout(ncol = 2, widths = c(0.1, 1))


# ggsave(filename = "Data\\Figures\\biomass-competition-08_intrinsic-growth.jpg", width = 10, height = 7)
# ggsave(filename = "Data\\Figures\\biomass-competition-08_intrinsic-growth.pdf", width = 10, height = 7)


p1 <- summary_boot.b1 %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  ggplot(aes(x = inoculum, y = estimate, fill = `focal x comp`)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = q.05, ymax = q.95), 
                width = 0.4, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("#1984c5", "#a7d5ed", "#e1a692", "#c23728")) +
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "none") + 
  labs(#title = "Competition strength distribution affected by tot_comp dependent competition",
    y = "linear")

p2 <- summary_boot.b2 %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized+living", "sterilized+living", "unsterilized+dead" , "sterilized+dead"))) %>%
  ggplot(aes(x = inoculum, y = estimate, fill = `focal x comp`)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = q.05, ymax = q.95), 
                width = 0.4, position = position_dodge(0.9)) +
  scale_fill_manual(values=c("#1984c5", "#a7d5ed", "#e1a692", "#c23728")) +
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") + 
  labs(#title = "Competition strength distribution affected by tot_comp dependent competition",
    y = "quadratic") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))


left_axis_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.0, label = "competitor effect", angle = 90, size = 8) +
  theme_void()

combined_plot <- p1 + p2 + 
  plot_layout(nrow = 2, 
              heights = c(5, 1))

left_axis_label + combined_plot + plot_layout(ncol = 2, widths = c(0.1, 1))


# ggsave(filename = "Data\\Figures\\biomass-competition-08_comp-coef.jpg", width = 10, height = 7)
# ggsave(filename = "Data\\Figures\\biomass-competition-08_comp-coef.pdf", width = 10, height = 7)


########################
### 3. Invasion analysis
########################

boot.monoculture <- boot.metrics.unnest %>% 
  filter(competition == "conspecific") %>%
  mutate(N_star = case_when(
    b2 > 0 ~ 99999,
    b2 < 0 ~ 1/(2*b2)*(-b1 + sqrt(b1^2 - 4*b2*a)),
    b1 < 0 ~ -a/b1,
    TRUE ~ 99999
  )) %>%
  select(id, competitor, treatment, decay, N_star) # focal = competitor

boot.invasion <- boot.metrics.unnest %>% 
  filter(competition == "heterospecific") %>%
  left_join(boot.monoculture) %>%
  mutate(IGR = a + b1*N_star + b2*N_star^2) %>%
  select(id, focal, competitor, treatment, decay, IGR)

boot.comp_outcome <- boot.invasion %>%
  group_by(treatment, decay, id) %>%
  group_modify(~ .x %>% 
                 mutate(outcome = case_when(
                   sum(.x$IGR > 0) == 2 ~ "coexist",
                   .x$IGR[.x$focal == "MACHZU"] > 0 ~ "MACHZU wins",
                   .x$IGR[.x$focal == "ENGERO"] > 0 ~ "ENGERO wins",
                   TRUE ~ "priority effect"
                 )) %>% 
                 select(outcome) %>%
                 unique
  ) %>%
  ungroup

boot.apparent.comp_outcome <- boot.comp_outcome %>%
  filter(id == "Apparent") %>%
  mutate(inoculum = paste0(treatment, "\n+", decay)) %>%
  select(inoculum, outcome) %>%
  mutate(label = "\u25B2")

boot.clean.comp_outcome <- boot.comp_outcome %>%
  filter(id != "Apparent") %>%
  mutate(inoculum = paste0(treatment, "\n+", decay))

boot.summary.comp_outcome <- boot.clean.comp_outcome %>% 
  select(inoculum, outcome) %>%
  group_by(inoculum) %>%
  count(outcome) %>% 
  mutate(ratio = n/10000) 

### Plotting
boot.summary.comp_outcome %>% 
  left_join(boot.apparent.comp_outcome) %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized\n+living", "sterilized\n+living", "unsterilized\n+dead" , "sterilized\n+dead"))) %>%
  ggplot(aes(x = inoculum, y = ratio, fill = outcome)) +
  geom_bar(stat = "identity", color="black", position=position_dodge()) +
  geom_text(aes(label = label), vjust = 1.6, color = "black",
            position = position_dodge(0.9), size = 8) +
  annotate("text", x = 1.8, y = 1, size = 7,
           label = "\u25B2 denotes the estimated outcome") +
  scale_fill_manual(values = c("#a4a2a8", "#1984c5", "#c23728", "#e2e2e2")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = "")) +
  theme_bw() + 
  theme(text = element_text(size = 25),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.85, 0.87)) + 
  labs(#title = "competitive outcome ratio affected by biomass dependent competition", 
    y = "proportion of bootstrap replicates (%)")


# ggsave(filename = "Data\\Figures\\biomass-competition-08_comp-outcome.jpg", width = 10, height = 7)
# ggsave(filename = "Data\\Figures\\biomass-competition-08_comp-outcome.pdf", width = 10, height = 7, device = cairo_pdf)

boot.apparent.E_persist <- boot.apparent.comp_outcome %>%
  mutate(outcome = case_when(
    outcome %in% c("ENGERO wins", "coexist") ~ "IGR > 0",
    TRUE ~ "IGR < 0"
  )) %>%
  rename("IGR" = "outcome")

boot.apparent.M_persist <- boot.apparent.comp_outcome %>%
  mutate(outcome = case_when(
    outcome %in% c("MACHZU wins", "coexist") ~ "IGR > 0",
    TRUE ~ "IGR < 0"
  )) %>%
  rename("IGR" = "outcome")

E_persist <- boot.summary.comp_outcome %>%
  filter(outcome %in% c("ENGERO wins", "coexist")) %>% 
  select(-outcome, -n) %>%
  group_by(inoculum) %>%
  summarise_all(sum) %>%
  mutate(`IGR < 0` = 1 - ratio) %>%
  rename(`IGR > 0` = "ratio") %>%
  gather(IGR, ratio, 2:3) %>%
  left_join(boot.apparent.E_persist)

M_persist <- boot.summary.comp_outcome %>%
  filter(outcome %in% c("MACHZU wins", "coexist")) %>% 
  select(-outcome, -n) %>%
  group_by(inoculum) %>%
  summarise_all(sum) %>%
  mutate(`IGR < 0` = 1 - ratio) %>%
  rename(`IGR > 0` = "ratio") %>%
  gather(IGR, ratio, 2:3) %>%
  left_join(boot.apparent.M_persist)

p1 <- E_persist %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized\n+living", "sterilized\n+living", "unsterilized\n+dead" , "sterilized\n+dead")),
         IGR = factor(IGR, levels = c("IGR > 0", "IGR < 0"))) %>%
  ggplot(aes(x = inoculum, y = ratio, fill = IGR)) +
  geom_bar(stat = "identity", color="black", position=position_dodge(), width = 0.5) +
  geom_text(aes(label = label), vjust = 1.4, color = "black",
            position = position_dodge(0.5), size = 8) +
  scale_fill_manual(values = c("#1984c5", "#e2e2e2")) +
  theme_bw() + 
  theme(text = element_text(size = 25),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank()) + 
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1, suffix = "")) +
  labs(fill = "ENGERO's IGR")

p2 <- M_persist %>%
  mutate(inoculum = factor(inoculum, levels = c("unsterilized\n+living", "sterilized\n+living", "unsterilized\n+dead" , "sterilized\n+dead")),
         IGR = factor(IGR, levels = c("IGR > 0", "IGR < 0"))) %>%
  ggplot(aes(x = inoculum, y = ratio, fill = IGR)) +
  geom_bar(stat = "identity", color="black", position=position_dodge(), width = 0.5) +
  geom_text(aes(label = label), vjust = 1.4, color = "black",
            position = position_dodge(0.5), size = 8) +
  annotate("text", x = 1.8, y = 0.98, size = 6,
           label = "\u25B2 denotes the estimated outcome") +
  scale_fill_manual(values = c("#c23728", "#e2e2e2")) +
  theme_bw() + 
  theme(text = element_text(size = 25),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()) + 
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent_format(accuracy = 1, suffix = "")) +
  labs(fill = "MACHZU's IGR")

left_axis_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.0, label = "proportion of bootstrap replicates (%)", angle = 90, size = 8.5) +
  theme_void()

combined_plot <- p1 + p2 + 
  plot_layout(nrow = 2, 
              heights = c(1, 1))

left_axis_label + combined_plot + plot_layout(ncol = 2, widths = c(0.1, 1))


# ggsave(filename = "Data\\Figures\\biomass-competition-08_persist.jpg", width = 10, height = 7)
# ggsave(filename = "Data\\Figures\\biomass-competition-08_persist.pdf", width = 10, height = 7, device = cairo_pdf)


########################
### 4. Write out bootstrapping data
########################


# setwd("Data\\Bootstrap")
# write.csv(summary_boot.a, file = "biomass-competition_boot_a.csv")
# write.csv(summary_boot.b1, file = "biomass-competition_boot_b1.csv")
# write.csv(summary_boot.b2, file = "biomass-competition_boot_b2.csv")
# write.csv(boot.summary.comp_outcome, file = "biomass-competition_boot_comp_outcome.csv")
# write.csv(E_persist, file = "biomass-competition_boot_E_persist.csv")
# write.csv(M_persist, file = "biomass-competition_boot_M_persist.csv")


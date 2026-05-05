## -----------------------------------------------------------------------------
## Title: Analyze CPB egg predation rates in the control and pheromone-treated plots
##
## Author: Gen-Chang Hsu
##
## Date: 2026-04-29
##
## Description:
## 1. Organize the raw data
## 2. Analyze egg predation rates in the control and pheromone-treated plots (shrunk eggs included)
## 3. Analyze egg predation rates in the control and pheromone-treated plots (shrunk eggs excluded)
## 4. Visualize the data
##
## -----------------------------------------------------------------------------
set.seed(123)


# Libraries --------------------------------------------------------------------
library(tidyverse)
library(readxl)
library(glmmTMB)
library(DHARMa)
library(performance)
library(lmtest)
library(car)
library(sjPlot)
library(emmeans)
library(multcomp)
library(grid)
library(ggplotify)
library(ggeffects)
library(officer)


# Functions for extracting LMM, LMM, GLM, and GLMM results ---------------------
extract_model_results_lm <- function(models) {
  
  model_dfs <- map(models, function(x) {
    
    df_residual <- summary(x)$df[2]
    
    x %>% 
      Anova() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Predictor") %>% 
      rename(P = `Pr(>F)`) %>% 
      mutate(n = nobs(x), .before = "Predictor") %>% 
      mutate(`F` = if_else(`F value` > 0.05, round(`F value`, 1), round(`F value`, 2)),
             P = case_when(P > 0.01 ~ as.character(round(P, 2)),
                           P < 0.01 & P > 0.001 ~ as.character(round(P, 3)),
                           P < 0.001 ~ "< 0.001")) %>% 
      filter(Predictor != "Residuals") %>% 
      mutate(Df_residual = df_residual, .after = Df) %>% 
      dplyr::select(n, Predictor, `F`, Df, Df_residual, P)
  })
  
  model_results <- model_dfs %>% 
    bind_rows(.id = "Model_response")
  
  return(model_results)
}

extract_model_results_lmm <- function(models) {
  
  model_dfs <- map(models, function(x) {
    x %>% 
      Anova(., test.statistic = "F") %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Predictor") %>% 
      rename(P = `Pr(>F)`) %>% 
      mutate(n = nobs(x), .before = "Predictor") %>% 
      mutate(`F` = if_else(`F` > 0.05, round(`F`, 1), round(`F`, 2)),
             Df.res = round(Df.res, 1),
             P = case_when(P > 0.01 ~ as.character(round(P, 2)),
                           P < 0.01 & P > 0.001 ~ as.character(round(P, 3)),
                           P < 0.001 ~ "< 0.001")) %>% 
      dplyr::select(n, Predictor, `F`, Df, Df_residual = Df.res, P)
  })
  
  model_results <- model_dfs %>% 
    bind_rows(.id = "Model_response")
  
  return(model_results)
}

extract_model_results_glm <- function(models) {
  
  model_dfs <- map(models, function(x) {
    x %>% 
      Anova() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Predictor") %>% 
      rename(P = `Pr(>Chisq)`,
             Chisq = `LR Chisq`) %>% 
      mutate(n = x$modelInfo$nobs, .before = "Predictor") %>% 
      mutate(Chisq = if_else(Chisq > 0.05, round(Chisq, 1), round(Chisq, 2)),
             P = case_when(P > 0.01 ~ as.character(round(P, 2)),
                           P < 0.01 & P > 0.001 ~ as.character(round(P, 3)),
                           P < 0.001 ~ "< 0.001"))
  })
  
  model_results <- model_dfs %>% 
    bind_rows(.id = "Model_response")
  
  return(model_results)
}

extract_model_results_glmm <- function(models) {
  
  model_dfs <- map(models, function(x) {
    x %>% 
      Anova() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "Predictor") %>% 
      rename(P = `Pr(>Chisq)`) %>% 
      mutate(n = x$modelInfo$nobs, .before = "Predictor") %>% 
      mutate(Chisq = if_else(Chisq > 0.05, round(Chisq, 1), round(Chisq, 2)),
             P = case_when(P > 0.01 ~ as.character(round(P, 2)),
                           P < 0.01 & P > 0.001 ~ as.character(round(P, 3)),
                           P < 0.001 ~ "< 0.001"))
  })
  
  model_results <- model_dfs %>% 
    bind_rows(.id = "Model_response")
  
  return(model_results)
}


# ggplot theme -----------------------------------------------------------------
my_ggtheme <- 
  theme(# axis
    axis.text.x = element_text(size = 14, color = "black", margin = margin(t = 3)),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, margin = margin(t = 10)),
    axis.title.y = element_text(size = 16, margin = margin(r = 8)),
    axis.ticks.length.x = unit(0.18, "cm"),
    axis.ticks.length.y = unit(0.15, "cm"),
    
    # plot
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(size = 16),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    plot.background = element_rect(colour = "transparent"),
    
    # panel
    panel.background = element_rect(fill = "transparent"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    # legend
    legend.position = "right",
    legend.key.spacing.x = unit(0.2, "cm"),
    legend.key.spacing.y = unit(0.2, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.key.size = unit(0.5, "line"),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, hjust = 0.5),
    legend.box.just = "center",
    legend.justification = c(0.5, 0.5),
    legend.background = element_rect(fill = "transparent", linewidth = 0.25, linetype = "solid", colour = "black"),
    
    # facet strip
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 13, hjust = 0.5)
  )


# Import files -----------------------------------------------------------------
sentinel_prey_raw <- readxl::read_xlsx("./01_Data_Raw/Sentinel_Prey.xlsx", sheet = "Sentinel_Prey", na = "NA")


############################### Code starts here ###############################

# 1. Organize the raw data -----------------------------------------------------
sentinel_prey_clean <- sentinel_prey_raw %>% 
  mutate(Trial_id = as.factor(Trial_id)) %>% 
  mutate(N_control_diff = N_prey_control_start - N_prey_control_end,
         N_predation_diff = N_prey_predation_start - N_prey_predation_end,
         N_predation_diff_corrected = N_predation_diff - N_control_diff * (N_prey_predation_start / N_prey_control_start)) %>% 
  mutate(N_predation_diff_corrected = if_else(N_predation_diff_corrected < 0, 0, round(N_predation_diff_corrected))) %>% 
  drop_na(N_control_diff, N_predation_diff) %>%
  filter(is.na(Notes) | !str_detect(Notes, "hatched|disappeared")) 

### Exclude the observations with shrunk eggs
sentinel_prey_clean_shrunk_eggs_excluded <- sentinel_prey_raw %>% 
  mutate(Trial_id = as.factor(Trial_id)) %>% 
  mutate(N_control_diff = N_prey_control_start - N_prey_control_end,
         N_predation_diff = N_prey_predation_start - N_prey_predation_end,
         N_predation_diff_corrected = N_predation_diff - N_control_diff * (N_prey_predation_start / N_prey_control_start)) %>%  
  mutate(N_predation_diff_corrected = if_else(N_predation_diff_corrected < 0, 0, N_predation_diff_corrected)) %>% 
  drop_na(N_control_diff, N_predation_diff) %>%
  filter(is.na(Notes) | !str_detect(Notes, "hatched|disappeared")) %>% 
  filter(is.na(Notes) | !str_detect(Notes, "shrank"))


# 2. Analyze egg predation rates in the control and pheromone-treated plots (shrunk eggs included) ----
### (1) Test overdispersion
predation_mortality_poisson <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                         data = sentinel_prey_clean,
                                         family = "poisson",
                                         na.action = na.omit)

predation_mortality_nb <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                  data = sentinel_prey_clean,
                                  family = "nbinom2",
                                  na.action = na.omit)

lrtest(predation_mortality_poisson, predation_mortality_nb)  # overdispersion is significant
AIC(predation_mortality_poisson, predation_mortality_nb)  # negative bimonial is better

### (2) Test zero inflation
predation_mortality_zi_nb <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                  data = sentinel_prey_clean,
                                  ziformula = ~ Pheromone_treatment,
                                  family = "nbinom2",
                                  na.action = na.omit)

testZeroInflation(predation_mortality_nb)  # zero inflation not significant
lrtest(predation_mortality_nb, predation_mortality_zi_nb)  # zero inflation not significant
AIC(predation_mortality_nb, predation_mortality_zi_nb)  # model without zero inflation is better

### (3) Model diagnostics
plot(simulateResiduals(predation_mortality_nb))  # no pattern

### (4) Model summary
summary(predation_mortality_nb)
Anova(predation_mortality_nb, type = 2)

### (5) emmeans
emmeans_predation_mortality <- emmeans(predation_mortality_nb, "Pheromone_treatment", type = "response", infer = c(T, T))
pairs(emmeans_predation_mortality)
cld(emmeans_predation_mortality, Letters = letters)

### (6) Model visualization
plot_model(predation_mortality_nb, 
           type = "pred", 
           terms = c("Pheromone_treatment"))
  

# 3. Analyze egg predation rates in the control and pheromone-treated plots (shrunk eggs excluded) ----
### (1) Test overdispersion
predation_mortality_poisson_shrunk_eggs_excluded <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                       data = sentinel_prey_clean_shrunk_eggs_excluded,
                                       family = "poisson",
                                       na.action = na.omit)

predation_mortality_nb_shrunk_eggs_excluded <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                  data = sentinel_prey_clean_shrunk_eggs_excluded,
                                  family = "nbinom2",
                                  na.action = na.omit)

lrtest(predation_mortality_poisson_shrunk_eggs_excluded, predation_mortality_nb_shrunk_eggs_excluded)  # overdispersion is significant
AIC(predation_mortality_poisson_shrunk_eggs_excluded, predation_mortality_nb_shrunk_eggs_excluded)  # negative binomial model is better

### (2) Test zero inflation
predation_mortality_zi_nb_shrunk_eggs_excluded <- glmmTMB(N_predation_diff_corrected ~ Pheromone_treatment + Trial_id,
                                     data = sentinel_prey_clean_shrunk_eggs_excluded,
                                     ziformula = ~ 1,
                                     family = "nbinom2",
                                     na.action = na.omit)

testZeroInflation(predation_mortality_nb_shrunk_eggs_excluded)  # zero inflation not significant
lrtest(predation_mortality_nb_shrunk_eggs_excluded, predation_mortality_zi_nb_shrunk_eggs_excluded)  # zero inflation not significant
AIC(predation_mortality_nb_shrunk_eggs_excluded, predation_mortality_zi_nb_shrunk_eggs_excluded)  # model without zero inflation is better

### (3) Model diagnostics
plot(simulateResiduals(predation_mortality_nb_shrunk_eggs_excluded))  # no pattern

### (4) Model summary
summary(predation_mortality_nb_shrunk_eggs_excluded)
Anova(predation_mortality_nb_shrunk_eggs_excluded, type = 2)

### (5) emmeans
emmeans_predation_mortality_shrunk_eggs_excluded <- emmeans(predation_mortality_nb_shrunk_eggs_excluded, "Pheromone_treatment", type = "response", infer = c(T, T))
pairs(emmeans_predation_mortality_shrunk_eggs_excluded)
cld(emmeans_predation_mortality_shrunk_eggs_excluded, Letters = letters)

### (6) Model visualization
plot_model(predation_mortality_nb_shrunk_eggs_excluded, 
           type = "pred", 
           terms = c("Pheromone_treatment"))


# 3. Visualize the data --------------------------------------------------------
### (1) Egg predation rates in the control and pheromone plots by week (shrunk eggs included)
sentinel_prey_clean_summary_by_treatment_week <- sentinel_prey_clean %>% 
  group_by(Trial_id, Date_start, Pheromone_treatment) %>% 
  summarise(Mean_n_predation_diff_corrected = mean(N_predation_diff_corrected, na.rm = T),
            N = n(),
            SE = sd(N_predation_diff_corrected, na.rm = T)/sqrt(N))

ggplot(data = sentinel_prey_clean_summary_by_treatment_week) + 
  geom_point(aes(x = Date_start, y = Mean_n_predation_diff_corrected, color = Pheromone_treatment, group = Pheromone_treatment),
             position = position_dodge2(width = 2e5), size = 1.5) + 
  geom_path(aes(x = Date_start, y = Mean_n_predation_diff_corrected, color = Pheromone_treatment, group = Pheromone_treatment),
            position = position_dodge2(width = 2e5)) + 
  geom_errorbar(aes(x = Date_start, ymin = Mean_n_predation_diff_corrected - SE, ymax = Mean_n_predation_diff_corrected + SE, color = Pheromone_treatment, group = Pheromone_treatment),
                position = position_dodge2(width = 0.75e5, padding = 0, preserve = "single"), width = 2e5) +
  labs(x = NULL, y = "Predation rates (eggs/3 days)") + 
  scale_x_datetime(breaks = seq(min(sentinel_prey_clean_summary_by_treatment_week$Date_start), max(sentinel_prey_clean_summary_by_treatment_week$Date_start), "7 days"), date_labels = "%B %d") + 
  scale_color_manual(values = c("#e66101", "#5e3c99"), labels = c("Control", "Pheromone")) +
  scale_y_continuous(limits = c(-0.5, 41), expand = c(0, 0)) +
  guides(color = guide_legend(byrow = T)) + 
  my_ggtheme + 
  theme(legend.position = c(0.25, 0.86),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.key.width = unit(0.3, "in"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("./03_Outputs/Figures/Predation_Rates_by_Treatment_and_Week_Shrunk_Eggs_Included.tiff", width = 6, height = 4, dpi = 600, device = "tiff")

### (2) Egg predation rates by week pooled across the control and pheromone treatment (shrunk eggs included)
sentinel_prey_clean_summary_by_week <- sentinel_prey_clean %>% 
  group_by(Trial_id, Date_start) %>% 
  summarise(Mean_n_predation_diff_corrected = mean(N_predation_diff_corrected, na.rm = T),
            N = n(),
            SE = sd(N_predation_diff_corrected, na.rm = T)/sqrt(N))

ggplot(data = sentinel_prey_clean_summary_by_week) + 
  geom_point(aes(x = Date_start, y = Mean_n_predation_diff_corrected), size = 1.5) + 
  geom_path(aes(x = Date_start, y = Mean_n_predation_diff_corrected)) + 
  geom_errorbar(aes(x = Date_start, ymin = Mean_n_predation_diff_corrected - SE, ymax = Mean_n_predation_diff_corrected + SE), width = 1e5) +
  labs(x = NULL, y = "Predation rates (eggs/3 days)") + 
  scale_x_datetime(breaks = seq(min(sentinel_prey_clean_summary_by_treatment_week$Date_start), max(sentinel_prey_clean_summary_by_treatment_week$Date_start), "7 days"), date_labels = "%B %d") + 
  scale_y_continuous(limits = c(-0.5, 31), expand = c(0, 0)) +
  my_ggtheme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("./03_Outputs/Figures/Predation_Rates_by_Week_Shrunk_Eggs_Included.tiff", width = 6, height = 4, dpi = 600, device = "tiff")

### (3) Egg predation rates pooled across weeks (shrunk eggs included)
sentinel_prey_clean_summary_by_treatment <- sentinel_prey_clean %>% 
  group_by(Pheromone_treatment) %>% 
  summarise(Mean_n_predation_diff_corrected = mean(N_predation_diff_corrected, na.rm = T),
            N = n(),
            SE = sd(N_predation_diff_corrected, na.rm = T)/sqrt(N))

ggplot(data = sentinel_prey_clean_summary_by_treatment) + 
  geom_col(aes(x = Pheromone_treatment, y = Mean_n_predation_diff_corrected, fill = Pheromone_treatment), show.legend = F, width = 0.5) + 
  geom_errorbar(aes(x = Pheromone_treatment, ymin = Mean_n_predation_diff_corrected - SE, ymax = Mean_n_predation_diff_corrected + SE),
                width = 0.2, linewidth = 0.6) +
  labs(x = "Plot treatment", y = "Predation rates (eggs/3 days)") + 
  scale_y_continuous(limits = c(0, 15.2), expand = c(0, 0)) +
  scale_fill_manual(values = c("#e66101", "#5e3c99"), labels = c("Control", "Pheromone")) +
  my_ggtheme + 
  annotate(geom = "text", x = 1, y = 2.8, label = "a") + 
  annotate(geom = "text", x = 2, y = 14, label = "b")

ggsave("./03_Outputs/Figures/Predation_Rates_by_Treatment_Shrunk_Eggs_Included.tiff", width = 4, height = 4.5, dpi = 600, device = "tiff")








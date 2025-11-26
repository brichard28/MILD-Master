# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

####################################################
##    Hit Rates    ##
####################################################

lead_hit_rates <- read.csv("/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER_Lead_Hit_Rates.csv")
lag_hit_rates <- read.csv("/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER_Lag_Hit_Rates.csv")

# Remove unneeded columns, put in long format
lead_hit_rates$OriginalVariableNames <- array(0:39)
colnames(lead_hit_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
lead_hit_rates <- pivot_longer(lead_hit_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                          names_to = c("Spatialization"), values_to = "HitRate")

# Remove unneeded columns, put in long format
lag_hit_rates$OriginalVariableNames <- array(0:39)
colnames(lag_hit_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
lag_hit_rates <- pivot_longer(lag_hit_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                          names_to = c("Spatialization"), values_to = "HitRate")

lead_hit_rates$WordPosition <- "Lead"
lag_hit_rates$WordPosition <- "Lag"

# Merge all data frames on the common identifier columns
hit_rates <- rbind(lead_hit_rates, lag_hit_rates)

# Organize Factors
to.factor <- c('S','Spatialization','WordPosition')
hit_rates[, to.factor] <- lapply(hit_rates[, to.factor], as.factor)






####################################################
##    FA Rates    ##
####################################################

lead_FA_rates <- read.csv("/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER_Lead_FA_Rates.csv")
lag_FA_rates <- read.csv("/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/MILD-MASTER_Lag_FA_Rates.csv")

# Remove unneeded columns, put in long format
lead_FA_rates$OriginalVariableNames <- array(0:39)
colnames(lead_FA_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
lead_FA_rates <- pivot_longer(lead_FA_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                               names_to = c("Spatialization"), values_to = "FARate")

# Remove unneeded columns, put in long format
lag_FA_rates$OriginalVariableNames <- array(0:39)
colnames(lag_FA_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
lag_FA_rates <- pivot_longer(lag_FA_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                              names_to = c("Spatialization"), values_to = "FARate")

lead_FA_rates$WordPosition <- "Lead"
lag_FA_rates$WordPosition <- "Lag"



# Merge all data frames on the common identifier columns
FA_rates <- rbind(lead_FA_rates, lag_FA_rates)

# Organize Factors
to.factor <- c('S','Spatialization','WordPosition')
FA_rates[, to.factor] <- lapply(FA_rates[, to.factor], as.factor)




hit_rates$Spatialization <- factor(hit_rates$Spatialization,
                                   levels = c("ITD5","ITD15","ILD5","ILD15"))
FA_rates$Spatialization <- factor(FA_rates$Spatialization,
                                  levels = c("ITD5","ITD15","ILD5","ILD15"))

# --- Combine all plots with shared legend ---
# ------------- build interdigitated x positions (by SpatializationGroup) -------------
library(stringr)

# create SpatializationGroup in the per-measure dataframes
hit_rates <- hit_rates %>%
  mutate(SpatializationGroup = case_when(
    str_detect(Spatialization, "ITD") ~ "ITD",
    str_detect(Spatialization, "ILD") ~ "ILD"
  ))
hit_rates <- hit_rates %>%
  mutate(MagnitudeGroup = if_else(str_detect(Spatialization, "15"),"15","5")
  )

FA_rates <- FA_rates %>%
  mutate(SpatializationGroup = case_when(
    str_detect(Spatialization, "ITD") ~ "ITD",
    str_detect(Spatialization, "ILD") ~ "ILD"
  ))

FA_rates <- FA_rates %>%
  mutate(MagnitudeGroup = if_else(str_detect(Spatialization, "15"),"15","5")
  )

# explicit group ordering (change if you want ILD first)
group_levels <- c("ITD", "ILD")
mag_levels <- c("5","15")

# spatialization levels (keeps the order you already set)
spat_levels <- levels(hit_rates$Spatialization)

# build an ordered mapping: for each SpatializationGroup -> (all Leads across that group, then all Lags)
order_list <- list()
for (g in group_levels) {
  spats_in_group <- spat_levels[str_detect(spat_levels, g)]
  for (m in mag_levels) {
    mags_in_group <-  spat_levels[str_detect(mag_levels, m)]
    for (wp in c("Lead", "Lag")) {
      for (s in spats_in_group) {
        order_list[[length(order_list) + 1]] <- tibble(
          Spatialization = s,
          WordPosition = wp,
          SpatializationGroup = g,
          MagnitudeGroup = m
        )
      }
    }
  }
}
order_df <- bind_rows(order_list) %>%
  mutate(x_pos = row_number())
# --- Combine all plots with shared legend ---
# ------------- build explicit x positions in desired order -------------

# Define the explicit desired order
desired_order <- tibble(
  Spatialization = rep(c("ITD5", "ITD15", "ILD5", "ILD15"), times = 2),
  WordPosition   = rep(c("Lead", "Lag"), each = 4)
) %>%
  mutate(x_pos = row_number())

# Determine group (ITD vs ILD)
desired_order <- desired_order %>%
  mutate(SpatializationGroup = ifelse(grepl("ITD", Spatialization), "ITD", "ILD"))

desired_order <- desired_order %>%
  mutate(MagnitudeGroup = ifelse(grepl("15", Spatialization), "15", "5"))

# Join this ordering into your datasets
hit_rates <- hit_rates %>%
  left_join(desired_order, by = c("Spatialization", "WordPosition", "SpatializationGroup", "MagnitudeGroup"))

FA_rates <- FA_rates %>%
  left_join(desired_order, by = c("Spatialization", "WordPosition", "SpatializationGroup", "MagnitudeGroup"))

# ---------------------------
# Prepare long-format data (use the new x_pos)
# ---------------------------
long_data <- bind_rows(
  hit_rates %>%
    select(S, Spatialization, WordPosition, SpatializationGroup, MagnitudeGroup, x_pos, Value = HitRate) %>%
    mutate(Measure = "Hit Rate"),
  
  FA_rates %>%
    select(S, Spatialization, WordPosition, SpatializationGroup, MagnitudeGroup, x_pos, Value = FARate) %>%
    mutate(Measure = "FA Rate")
)

# Summary statistics
summary_data <- long_data %>%
  group_by(Spatialization, WordPosition, x_pos, Measure, SpatializationGroup, MagnitudeGroup) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    sem_value  = sd(Value, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

# ---------------------------
# Build x-axis labels (to match your order)
# ---------------------------
Spatialization_positions <- desired_order %>%
  group_by(Spatialization) %>%
  summarise(mean_x = mean(x_pos), .groups = "drop")

Spatialization_labels_wrapped <- Spatialization_positions$Spatialization

# ---------------------------
# Custom y-limits per Measure (re-add so geom_blank() has ymin/ymax)
# ---------------------------
y_limits <- tibble(
  Measure = levels(factor(long_data$Measure)),
  ymin = c(0, 0),
  ymax = c(1, 1)
)



summary_data <- summary_data %>% left_join(y_limits, by = "Measure")
long_data    <- long_data %>% left_join(y_limits, by = "Measure")

long_data$Measure <- factor(long_data$Measure,
                            levels = c("Hit Rate", "FA Rate"))
summary_data$Measure <- factor(summary_data$Measure,
                               levels = c("Hit Rate", "FA Rate"))

# ------------------------------------------------------------------------------------


p <- ggplot(long_data, aes(x = x_pos, y = Value, shape = WordPosition, fill = MagnitudeGroup)) +
  geom_line(aes(group = interaction(S, WordPosition, SpatializationGroup)), alpha = 0.3) +
  geom_point(aes(group = interaction(S, WordPosition)), alpha = 0.3) +
  geom_errorbar(data = summary_data,
                aes(x = x_pos, ymin = mean_value - sem_value, ymax = mean_value + sem_value),
                width = 0.2, size = 1.2, inherit.aes = FALSE) +
  geom_point(data = summary_data,
             aes(x = x_pos, y = mean_value, fill = MagnitudeGroup, shape = WordPosition),
            size = 3, stroke = 1.2, inherit.aes = FALSE) +
  geom_blank(data = long_data, aes(y = ymin)) +
  geom_blank(data = long_data, aes(y = ymax)) +
  scale_shape_manual(values = c("Lead" = 21, "Lag" = 24)) +
  scale_fill_manual(values = c("5" = "black", "15" = "white")) + 
  scale_x_continuous(breaks = Spatialization_positions$mean_x,
                     labels = Spatialization_labels_wrapped) +
  labs(x = "Spatialization", y = "") +
  facet_grid(. ~ Measure, scales = "free_y", switch = "y") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y.left = element_text(margin = margin(l = 10), size = 14, face = "bold"),
    axis.text.y.left = element_text(margin = margin(r = 5), size = 12),
    axis.text.x = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    plot.margin = unit(c(5,5,5,40), "pt"),
    strip.placement = "outside",
    panel.grid.major.x = element_blank(),  # << remove vertical lines
    panel.grid.minor.x = element_blank()   # << remove vertical minor lines too
  )
p

ggsave("/Users/benrichardson/Documents/GitHub/MILD-Master/PAPER FIGURES/behavior_raw.svg", p, device="svg", width = 11, height = 6, units = "in", bg = "transparent")









## Refactor for stat analysis
library(emmeans)


# LMEM for Hit Rates
model_hitrate <- mixed(HitRate ~ Spatialization*WordPosition + (1|S),
                       data = hit_rates,
                       control = lmerControl(optimizer = "bobyqa"),
                       method = 'LRT')
# significant effect of Spatialization post hoc
emm_spatialization <- emmeans(model_hitrate, ~ Spatialization)
pairs(emm_spatialization, adjust = "bonferroni")

# Significant effect of Word position post hoc 
emm_wp <- emmeans(model_hitrate, ~ WordPosition)
pairs(emm_wp, adjust = "bonferroni")


# LMEM for False Alarm Rates
model_FArate <- mixed(FARate ~ Spatialization*WordPosition + (1|S),
                       data = FA_rates,
                       control = lmerControl(optimizer = "bobyqa"),
                       method = 'LRT')

# Significant effect of Spatialization post hoc
emm_fa <- emmeans(model_FArate, ~ Spatialization)
pairs(emm_fa, adjust = "bonferroni")



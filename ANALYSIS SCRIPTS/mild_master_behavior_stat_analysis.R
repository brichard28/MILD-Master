# Author: Benjamin Richardson
# Uses information from srm_nirs_eeg_analyze_behavior.m


library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)

# SummarySE Function ####
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


####################################################
##    Hit Rates    ##
####################################################

lead_hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_Lead_Hit_Rates.csv")
lag_hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_Lag_Hit_Rates.csv")

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

lead_FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_Lead_FA_Rates.csv")
lag_FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_Lag_FA_Rates.csv")

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

FA_rates <- FA_rates %>%
  mutate(SpatializationGroup = case_when(
    str_detect(Spatialization, "ITD") ~ "ITD",
    str_detect(Spatialization, "ILD") ~ "ILD"
  ))

# explicit group ordering (change if you want ILD first)
group_levels <- c("ITD", "ILD")

# spatialization levels (keeps the order you already set)
spat_levels <- levels(hit_rates$Spatialization)

# build an ordered mapping: for each SpatializationGroup -> (all Leads across that group, then all Lags)
order_list <- list()
for (g in group_levels) {
  spats_in_group <- spat_levels[str_detect(spat_levels, g)]
  for (wp in c("Lead", "Lag")) {
    for (s in spats_in_group) {
      order_list[[length(order_list) + 1]] <- tibble(
        Spatialization = s,
        WordPosition = wp,
        SpatializationGroup = g
      )
    }
  }
}
order_df <- bind_rows(order_list) %>%
  mutate(x_pos = row_number())

# join the new x_pos into your per-measure data frames
hit_rates <- hit_rates %>%
  left_join(order_df, by = c("Spatialization", "WordPosition", "SpatializationGroup"))

FA_rates <- FA_rates %>%
  left_join(order_df, by = c("Spatialization", "WordPosition", "SpatializationGroup"))

# ---------------------------
# Prepare long-format data (use the new x_pos)
# ---------------------------
long_data <- bind_rows(
  hit_rates %>%
    select(S, Spatialization, WordPosition, SpatializationGroup, x_pos, Value = HitRate) %>%
    mutate(Measure = "Hit Rate"),
  
  FA_rates %>%
    select(S, Spatialization, WordPosition, SpatializationGroup, x_pos, Value = FARate) %>%
    mutate(Measure = "FA Rate")
)

# Summary statistics
summary_data <- long_data %>%
  group_by(Spatialization, WordPosition, x_pos, Measure) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    sem_value  = sd(Value, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

# ---------------------------
# Build breaks and labels so each Spatialization label sits between its Lead/Lag points
# ---------------------------
Spatialization_positions <- order_df %>%
  group_by(Spatialization) %>%
  summarise(mean_x = x_pos, .groups = "drop") %>%
  arrange(mean_x)

Spatialization_labels_wrapped <- str_wrap(Spatialization_positions$Spatialization, width = 1)

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


p <- ggplot(long_data, aes(x = x_pos, y = Value, color = WordPosition)) +
  geom_line(aes(group = interaction(S, WordPosition, SpatializationGroup)), alpha = 0.3) +
  geom_point(aes(group = interaction(S, WordPosition)), alpha = 0.3) +
  geom_point(data = summary_data,
             aes(x = x_pos, y = mean_value, color = WordPosition, fill = WordPosition),
             size = 3, shape = 21, stroke = 1.2, inherit.aes = FALSE) +
  geom_errorbar(data = summary_data,
                aes(x = x_pos, ymin = mean_value - sem_value, ymax = mean_value + sem_value, color = WordPosition),
                width = 0.2, size = 1.2, inherit.aes = FALSE) +
  geom_blank(data = long_data, aes(y = ymin)) +
  geom_blank(data = long_data, aes(y = ymax)) +
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










## Refactor for stat analysis
library(emmeans)

# For hit_rates
hit_rates <- hit_rates %>%
  mutate(CueType = str_extract(Spatialization, "[A-Z]+"),   # extract letters
         CueMagnitude = str_extract(Spatialization, "[0-9]+") %>% as.numeric())  # extract digits
hit_rates$CueMagnitude <- factor(hit_rates$CueMagnitude)

# For FA_rates
FA_rates <- FA_rates %>%
  mutate(CueType = str_extract(Spatialization, "[A-Z]+"),
         CueMagnitude = str_extract(Spatialization, "[0-9]+") %>% as.numeric())
FA_rates$CueMagnitude <- factor(FA_rates$CueMagnitude)

# LMEM for Hit Rates
model_hitrate <- mixed(HitRate ~ CueType*CueMagnitude*WordPosition + (1|S),
                       data = hit_rates,
                       control = lmerControl(optimizer = "bobyqa"),
                       method = 'LRT')
# significant effect of cue magnitude post hoc
emm_cue <- emmeans(model_hitrate, ~ CueMagnitude)
pairs(emm_cue, adjust = "bonferroni")

# Significant effect of Word position post hoc 
emm_wp <- emmeans(model_hitrate, ~ WordPosition)
pairs(emm_wp, adjust = "bonferroni")


# LMEM for False Alarm Rates
model_FArate <- mixed(FARate ~ CueType*CueMagnitude*WordPosition + (1|S),
                       data = FA_rates,
                       control = lmerControl(optimizer = "bobyqa"),
                       method = 'LRT')

# Significant interaction between cue type and cue magnitude
# Get estimated marginal means
emm_fa <- emmeans(model_FArate, ~ CueMagnitude | CueType)

# Pairwise comparisons of CueMagnitude within each CueType
pairs(emm_fa, adjust = "bonferroni")



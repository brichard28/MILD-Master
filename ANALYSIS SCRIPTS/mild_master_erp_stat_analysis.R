# Author: Benjamin Richardson
# plot and conduct statistics for MILD-Master p1,n1, and p3 components experiment 1

library(tidyverse)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
#library(rhdf5)
library(janitor)

# Read the .mat file
p1_data <- read.csv('/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/all_subs_p1.csv')
n1_data <- read.csv('/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/all_subs_n1.csv')
p2_data <- read.csv('/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/all_subs_p2.csv')
p3_data <- read.csv('/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/all_subs_p3.csv')

names(p1_data)[names(p1_data) == "Lead_Amplitude"] <- "Lead_p1"
names(n1_data)[names(n1_data) == "Lead_Amplitude"] <- "Lead_n1"
names(p2_data)[names(p2_data) == "Lead_Amplitude"] <- "Lead_p2"
names(p3_data)[names(p3_data) == "Lead_Amplitude"] <- "Lead_p3"

names(p1_data)[names(p1_data) == "Lag_Amplitude"] <- "Lag_p1"
names(n1_data)[names(n1_data) == "Lag_Amplitude"] <- "Lag_n1"
names(p2_data)[names(p2_data) == "Lag_Amplitude"] <- "Lag_p2"
names(p3_data)[names(p3_data) == "Lag_Amplitude"] <- "Lag_p3"


p1_data <- p1_data %>%
  mutate(Lag_Stream = ifelse(Lead_Stream == "Target", "Masker","Target"))
p2_data <- p2_data %>%
  mutate(Lag_Stream = ifelse(Lead_Stream == "Target", "Masker","Target"))
n1_data <- n1_data %>%
  mutate(Lag_Stream = ifelse(Lead_Stream == "Target", "Masker","Target"))
p3_data <- p3_data %>%
  mutate(Lag_Stream = ifelse(Lead_Stream == "Target", "Masker","Target"))

# Merge all data frames on the common identifier columns
all_data <- Reduce(function(x, y) merge(x, y, by = c("S", "Condition", "Lead_Stream", "Lag_Stream", "Lead_Word", "Lag_Word", "Electrode")), list(p1_data, n1_data, p2_data, p3_data))

frontocentral_electrodes <- c("Fz","FC1","FC2","C3","Cz","C4","CP1","CP2")
parietooccipital_electrodes <- c("P3","Pz","PO3","O1","Oz","O2","PO4","P4")

# Define ERP variables to pivot
erp_vars <- c("p1", "n1", "p2", "p3")

# For frontocentral electrodes: calculate mean(n1 - p1) per wordposition, subject, lead word, and lag word
lead_frontocentral_summary <- all_data %>%
  filter(Electrode %in% frontocentral_electrodes) %>%
  mutate(p1n1 = Lead_p1 - Lead_n1) %>%
  group_by(S, Condition, Lead_Word, Lag_Word, Lead_Stream, Lag_Stream) %>%
  summarise(mean_p1n1 = mean(p1n1, na.rm = TRUE)) %>%
  ungroup()

lag_frontocentral_summary <- all_data %>%
  filter(Electrode %in% frontocentral_electrodes) %>%
  mutate(p1n1 = Lag_p1 - Lag_n1) %>%
  group_by(S, Condition, Lead_Word, Lag_Word, Lead_Stream, Lag_Stream) %>%
  summarise(mean_p1n1 = mean(p1n1, na.rm = TRUE)) %>%
  ungroup()

#names(lead_frontocentral_summary)[names(lead_frontocentral_summary) == "Lead_Stream"] <- "Stream"
#names(lag_frontocentral_summary)[names(lag_frontocentral_summary) == "Lag_Stream"] <- "Stream"


# For parietooccipital electrodes: calculate mean p3 per subject and group
lead_parietooccipital_summary <- all_data %>%
  filter(Electrode %in% parietooccipital_electrodes) %>%
  mutate(p3 = Lead_p3) %>%
  group_by(S, Condition, Lead_Word, Lag_Word, Lead_Stream, Lag_Stream) %>%
  summarise(mean_p3 = mean(p3, na.rm = TRUE)) %>%
  ungroup()

lag_parietooccipital_summary <- all_data %>%
  filter(Electrode %in% parietooccipital_electrodes) %>%
  mutate(p3 = Lag_p3) %>%
  group_by(S, Condition, Lead_Word, Lag_Word, Lead_Stream, Lag_Stream) %>%
  summarise(mean_p3 = mean(p3, na.rm = TRUE)) %>%
  ungroup()

#names(lead_parietooccipital_summary)[names(lead_parietooccipital_summary) == "Lead_Stream"] <- "Stream"
#names(lag_parietooccipital_summary)[names(lag_parietooccipital_summary) == "Lag_Stream"] <- "Stream"










## P1-N1 PLOTTING
library(dplyr)
library(ggplot2)

# Add a column to indicate whether the row is Lead or Lag
lead_frontocentral_summary <- lead_frontocentral_summary %>%
  mutate(Position = "Lead")

lag_frontocentral_summary <- lag_frontocentral_summary %>%
  mutate(Position = "Lag")

# Combine
all_summary_p1n1 <- bind_rows(lead_frontocentral_summary, lag_frontocentral_summary)

summary_data_p1n1 <- all_summary_p1n1 %>%
  group_by(Condition, Lead_Stream, Lead_Word, Lag_Word, Position) %>%
  summarise(
    mean_value = mean(mean_p1n1, na.rm = TRUE),
    sem_value  = sd(mean_p1n1, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"),
         # set factor order for x-axis
         WordPair = factor(WordPair, levels = c("bash_non-bash", "non-bash_bash", "non-bash_non-bash")),
         # relabel Lead_Stream facets
         Lead_Stream = factor(Lead_Stream, levels = c("Target", "Masker")),
         Condition = factor(Condition, levels = c("small_itd","large_itd","small_ild","large_ild")),
         Position = factor(Position,level = c("Lead","Lag")))

labels_for_column <- c("Target" = "Target Stream Leads", "Masker" = "Masker Stream Leads")


pos <- position_dodge(width = 0.4)
p1n1_figure <- ggplot(summary_data_p1n1,
       aes(x = WordPair,
           y = mean_value,
           color = WordPair,
           fill = Position,
           group = Position)) +
  geom_point(position = pos, size = 3, shape = 21, stroke = 1.2) +
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value),
                width = 0.4, position = pos) +
  facet_grid(Condition ~ Lead_Stream, labeller = labeller(Lead_Stream = labels_for_column)) +
  labs(x = "Lead-Lag Word Pair", y = "Mean P1-N1", color = "WordPair", fill = "Position") +
  scale_color_manual(values = c("bash_non-bash" = "blue", "non-bash_bash"="orange","non-bash_non-bash"="green")) +
  scale_fill_manual(values = c("Lead" = "black", "Lag" = "white")) +
  scale_x_discrete(expand = expansion(mult = 0.4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("/Users/benrichardson/Documents/GitHub/MILD-Master/CASUAL FIGURES/p1n1_results.svg", p1n1_figure,device="svg", width = 10, height = 7, units = "in", bg = "transparent")





## P300 PLOTTING
library(dplyr)
library(ggplot2)

# Add a column to indicate whether the row is Lead or Lag
lead_parietooccipital_summary <- lead_parietooccipital_summary %>%
  mutate(Position = "Lead")

lag_parietooccipital_summary <- lag_parietooccipital_summary %>%
  mutate(Position = "Lag")

# Combine
all_summary_p3 <- bind_rows(lead_parietooccipital_summary, lag_parietooccipital_summary)

summary_data_p3 <- all_summary_p3 %>%
  group_by(Condition, Lead_Stream, Lead_Word, Lag_Word, Position) %>%
  summarise(
    mean_value = mean(mean_p3, na.rm = TRUE),
    sem_value  = sd(mean_p3, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"),
         # set factor order for x-axis
         WordPair = factor(WordPair, levels = c("bash_non-bash", "non-bash_bash", "non-bash_non-bash")),
         # relabel Stream facets
         Lead_Stream = factor(Lead_Stream, levels = c("Target", "Masker")),
         Condition = factor(Condition, levels = c("small_itd","large_itd","small_ild","large_ild")),
         Position = factor(Position,level = c("Lead","Lag")))

labels_for_column <- c("Target" = "Target Stream Leads", "Masker" = "Masker Stream Leads")

pos <- position_dodge(width = 0.4)
p3_figure <- ggplot(summary_data_p3,
       aes(x = WordPair,
           y = mean_value,
           color = WordPair,
           fill = Position,
           group = Position)) +
  geom_point(position = pos, size = 3, shape = 21, stroke = 1.2) +
  geom_errorbar(aes(ymin = mean_value - sem_value, ymax = mean_value + sem_value),
                width = 0.4, position = pos) +
  facet_grid(Condition ~ Lead_Stream, labeller = labeller(Lead_Stream = labels_for_column)) +
  labs(x = "Lead-Lag Word Pair", y = "Mean P300", color = "WordPair", fill = "Position") +
  scale_color_manual(values = c("bash_non-bash" = "blue", "non-bash_bash"="orange","non-bash_non-bash"="green")) +
  scale_fill_manual(values = c("Lead" = "black", "Lag" = "white")) +
  scale_x_discrete(expand = expansion(mult = 0.4)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("/Users/benrichardson/Documents/GitHub/MILD-Master/CASUAL FIGURES/p3_results.svg", p3_figure, device="svg", width = 10, height = 7, units = "in", bg = "transparent")














## STATISTICAL ANALYSIS
library(emmeans)
# Lead p1n1
lead_p1n1_for_stats <- lead_frontocentral_summary %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"))

model_lead_p1n1 <- mixed(mean_p1n1 ~ Condition*WordPair*Lead_Stream + (1|S),
                         data = lead_p1n1_for_stats,
                         control = lmerControl(optimizer = "bobyqa"),
                         method = 'LRT')
# Significant effect of word pair
emm_wordpair <- emmeans(model_lead_p1n1, ~ WordPair)
pairs(emm_wordpair, adjust = "bonferroni")

# lag p1n1
lag_p1n1_for_stats <- lag_frontocentral_summary %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"))

model_lag_p1n1 <- mixed(mean_p1n1 ~ Condition*WordPair*Lead_Stream + (1|S),
                         data = lag_p1n1_for_stats,
                         control = lmerControl(optimizer = "bobyqa"),
                         method = 'LRT')
# Significant effect of Word Pair
emm_wordpair <- emmeans(model_lag_p1n1, ~ WordPair)
pairs(emm_wordpair, adjust = "bonferroni")

# Significant effect of Stream
emm_stream <- emmeans(model_lag_p1n1, ~ Lead_Stream)
pairs(emm_stream, adjust = "bonferroni")



# lead p300
lead_p3_for_stats <- lead_parietooccipital_summary %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"))

model_lead_p3 <- mixed(mean_p3 ~ Condition*WordPair*Lead_Stream + (1|S),
                         data = lead_p3_for_stats,
                         control = lmerControl(optimizer = "bobyqa"),
                         method = 'LRT')

# Significant interaction between word pair and stream
emm_wordpair_stream <- emmeans(model_lead_p3, ~ Lead_Stream | WordPair)
pairs(emm_wordpair_stream,adjust="bonferroni")


# lag p300
lag_p3_for_stats <- lag_parietooccipital_summary %>%
  # create combined label for Lead + Lag words
  mutate(WordPair = paste(Lead_Word, Lag_Word, sep = "_"))

model_lag_p3 <- mixed(mean_p3 ~ Condition*WordPair*Lead_Stream + (1|S),
                       data = lag_p3_for_stats,
                       control = lmerControl(optimizer = "bobyqa"),
                       method = 'LRT')

# Significant interaction between word pair and stream
emm_wordpair_stream <- emmeans(model_lag_p3, ~ Lead_Stream | WordPair)
pairs(emm_wordpair_stream,adjust="bonferroni")





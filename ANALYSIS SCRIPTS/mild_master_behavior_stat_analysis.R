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

hit_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_Hit_Rates.csv")

# Remove unneeded columns, put in long format
hit_rates$OriginalVariableNames <- array(0:39)
colnames(hit_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
hit_rates <- pivot_longer(hit_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                          names_to = c("Spatialization"), values_to = "HitRate")

# Organize Factors
to.factor <- c('S','Spatialization')
hit_rates[, to.factor] <- lapply(hit_rates[, to.factor], as.factor)

# Summary Statistics
# hit_rates %>% group_by(Condition, Masker) %>% get_summary_stats(HitRate, type = "mean_sd")
# 
# # Boxplot
# bxp <- ggboxplot(hit_rates, x = "Condition", y = "HitRate", color = "Masker", palette = "jco")
# bxp

# Check for normality, remove outliers
# hit_rates %>% group_by(Condition, Masker) %>% shapiro_test(HitRate)

#hit_rates_no_outliers <- hit_rates %>%
#  group_by(Condition, Masker) %>%
#  identify_outliers("HitRate") %>%
#  filter(!is.outlier)



# Create a QQ plot
# ggqqplot(hit_rates, "HitRate", ggtheme = theme_bw()) + facet_grid(Condition ~ Masker, labeller = "label_both")


model_hitrate <- mixed(HitRate ~ Spatialization + (1|S),
                       data= hit_rates, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

model_hitrate


# Compare All spatializations to each other
# ITD5 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD5")
posthoc_hitrate_itd5_speech <- lmer(HitRate ~ Spatialization + (1|S),
                                     data= hit_rates, 
                                     control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd5_speech)

# ITD15 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ITD15")
posthoc_hitrate_itd15_speech <- lmer(HitRate ~ Spatialization + (1|S),
                                    data= hit_rates, 
                                    control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_itd15_speech)

# ILD5 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD5")
posthoc_hitrate_ild5_speech <- lmer(HitRate ~ Spatialization + (1|S),
                                    data= hit_rates, 
                                    control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild5_speech)

# ILD15 as reference
hit_rates$Spatialization <- relevel(hit_rates$Spatialization, "ILD15")
posthoc_hitrate_ild15_speech <- lmer(HitRate ~ Spatialization + (1|S),
                                    data= hit_rates, 
                                    control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_hitrate_ild15_speech)


####################################################
##    FA Rates    ##
####################################################

FA_rates <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_FA_Rates.csv")

FA_rates$OriginalVariableNames <- array(0:39)
colnames(FA_rates) <- c("S","ITD5","ITD15","ILD5","ILD15")
FA_rates <- pivot_longer(FA_rates, cols=c("ITD5","ITD15","ILD5","ILD15"),
                         names_to = c("Spatialization"), values_to = "FARate")

# Organize Factors
to.factor <- c('S','Spatialization')
FA_rates[, to.factor] <- lapply(FA_rates[, to.factor], as.factor)

# Create a QQ plot
# ggqqplot(FA_rates, "FARate", ggtheme = theme_bw()) + facet_grid(Condition ~ Masker, labeller = "label_both")


model_FArate <- mixed(FARate ~ Spatialization + (1|S),
                      data= FA_rates, 
                      control = lmerControl(optimizer = "bobyqa"), method = 'LRT')

model_FArate


# Compare All spatializations to each other
# ITD5 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ITD5")
posthoc_FArate_itd5_speech <- lmer(FARate ~ Spatialization + (1|S),
                                   data= FA_rates, 
                                   control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_FArate_itd5_speech)

# ITD15 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ITD15")
posthoc_FArate_itd15_speech <- lmer(FARate ~ Spatialization + (1|S),
                                    data= FA_rates, 
                                    control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_FArate_itd15_speech)

# ILD5 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ILD5")
posthoc_FArate_ild5_speech <- lmer(FARate ~ Spatialization + (1|S),
                                   data= FA_rates, 
                                   control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_FArate_ild5_speech)

# ILD15 as reference
FA_rates$Spatialization <- relevel(FA_rates$Spatialization, "ILD15")
posthoc_FArate_ild15_speech <- lmer(FARate ~ Spatialization + (1|S),
                                    data= FA_rates, 
                                    control = lmerControl(optimizer = "bobyqa"))
summary(posthoc_FArate_ild15_speech)


####################################################
##    D primes    ##
####################################################

d_primes <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\MILD-MASTER_d_primes.csv")

# Remove unneeded columns, put in long format
d_primes$OriginalVariableNames <- array(0:39)
colnames(d_primes) <- c("S","ITD5","ITD15","ILD5","ILD15")
d_primes <- pivot_longer(d_primes, cols=c("ITD5","ITD15","ILD5","ILD15"),
                         names_to = c("Spatialization"), values_to = "DPrime")

# Organize Factors
to.factor <- c('S','Spatialization')
d_primes[, to.factor] <- lapply(d_primes[, to.factor], as.factor)

# LMEM
model_dprime <- mixed(DPrime ~ Spatialization + (1|S),data= d_primes,control = lmerControl(optimizer = "bobyqa"))

model_dprime


# Post hocs
# ITD5 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ITD5")
posthoc_dprime_itd5 <- lmer(DPrime ~ Spatialization + (1|S),
                             data= d_primes, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_itd5)


# ITD15 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ITD15")
posthoc_dprime_itd15 <- lmer(DPrime ~ Spatialization + (1|S),
                            data= d_primes, 
                            control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_itd15)

# ILD5 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ILD5")
posthoc_dprime_ild5 <- lmer(DPrime ~ Spatialization + (1|S),
                            data= d_primes, 
                            control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_ild5)


# ITD15 as reference
d_primes$Spatialization <- relevel(d_primes$Spatialization, "ILD15")
posthoc_dprime_ild15 <- lmer(DPrime ~ Spatialization + (1|S),
                             data= d_primes, 
                             control = lmerControl(optimizer = "bobyqa"))

summary(posthoc_dprime_ild15)






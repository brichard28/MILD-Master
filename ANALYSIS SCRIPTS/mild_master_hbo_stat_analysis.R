library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(car)
require(gridExtra)

#### Data Preparation ###########
# Load in Data
speech_masker_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\ANALYSIS SCRIPTS\\all_subjects_mean_during_stim_speech_masker.csv")

colnames(speech_masker_data) <- c("S","Channel","ITD50","ITD500","ILD5","ILD5Mag")
speech_masker_data <- pivot_longer(speech_masker_data, cols=c("ITD50","ITD500","ILD5","ILD5Mag"), names_to = "Spatialization", values_to = "MeanHbO")


# Combine speech and noise
all_data<-speech_masker_data


# Organize Factors
to.factor <- c('S','Spatialization')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)




all_data_cleaned <- na.omit(all_data)


### Summary SE function ########
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)
  
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


# Check for normality, remove outliers
all_data_cleaned %>% group_by(Spatialization) %>% shapiro_test(MeanHbO)




##### LMEM Model #############
model <- mixed(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                       data= all_data, 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model

# Pairwise Comparisons (treatment coding)

# ITD50 as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "ITD50")
lmer_itd50 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                              data= all_data,
                              control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd50)

# ITD500 as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "ITD500")
lmer_itd500 <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                   data= all_data,
                   control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd500)


# ILD5deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "ILD5")
lmer_ild5deg <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                    data= all_data,
                    control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild5deg)

# ILD5degMag as reference

all_data$Spatialization <- relevel(all_data$Spatialization, "ILD5Mag")
lmer_ild5degMag <- lmer(MeanHbO ~ Spatialization + (1|S) + (1|Channel),
                     data= all_data,
                     control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild5degMag)


##### PFC Plot ##########
pfc_se_data_speech <- summarySE(all_data_cleaned_pfc_speech, measurevar="MeanHbO", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_speech <- summarySE(pfc_se_data_speech, measurevar="MeanHbO", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotspeech <- ggplot(pfc_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Speech Masker PFC") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
#geom_signif(comparisons = list(c("ITD50","ITD500")), y_position = 0.09, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD70n")), y_position = 0.10, tip_length = 0, color="black", annotation = c("***"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.11, tip_length = 0, color="black", annotation = c("***"), textsize = 5)


pfc_se_data_noise <- summarySE(all_data_cleaned_pfc_noise, measurevar="MeanHbO", groupvars=c("S","Masker","Spatialization"), na.rm = TRUE)
pfc_se_data_noise <- summarySE(pfc_se_data_noise, measurevar="MeanHbO", groupvars=c("Masker","Spatialization"), na.rm = TRUE)
plotnoise <- ggplot(pfc_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) + 
  ggtitle("Noise Masker PFC") +
  labs(x="",y="") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))




######### STG Plot ###############

stg_se_data_speech <- summarySE(all_data_cleaned_stg_speech, measurevar="MeanHbO", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_speech <- summarySE(stg_se_data_speech, measurevar="MeanHbO", groupvars=c("Spatialization"),  na.rm=TRUE)
plotspeech <- ggplot(stg_se_data_speech, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Speech Masker STG") +
  labs(x="",y="Mean \u0394HbO (\u03BCM)", parse=TRUE) +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none")
#geom_signif(comparisons = list(c("ITD50","ILD10")), y_position = 0.090, tip_length = 0, color="black", annotation = c("**"), textsize = 5) +
#geom_signif(comparisons = list(c("ITD500","ILD10")), y_position = 0.080, tip_length = 0, color="black", annotation =c("*"), textsize = 5) +
#geom_signif(comparisons = list(c("ILD70n","ILD10")), y_position = 0.070, tip_length = 0, color="black", annotation =c("**"), textsize = 5) 


stg_se_data_noise <- summarySE(all_data_cleaned_stg_noise, measurevar="MeanHbO", groupvars=c("S","Spatialization"), na.rm=TRUE)
stg_se_data_noise <- summarySE(stg_se_data_noise, measurevar="MeanHbO", groupvars=c("Spatialization"),  na.rm=TRUE)
plotnoise <- ggplot(stg_se_data_noise, aes(x=factor(Spatialization, level=c('ITD50', 'ITD500', 'ILD70n', 'ILD10')), y=MeanHbO, color = Spatialization)) + 
  scale_color_manual(values = c("ITD50" = "red","ITD500" =  "blue","ILD70n" = "magenta","ILD10" = "green")) +
  geom_errorbar(aes(ymin=MeanHbO-se, ymax=MeanHbO+se, color=Spatialization), width=.1, position=position_dodge(width=0.5)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.5)) +  
  ggtitle("Noise Masker STG") +
  labs(x="", y="") +
  ylim(0,0.12) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_blank()) +
  scale_x_discrete(labels=c("ITD50" = "Small\nITD", "ITD500" = "Large\nITD","ILD70n" = "Natural\nILD","ILD10" = "Broadband\nILD")) +
  theme(legend.position="none") 

grid.arrange(plotspeech,plotnoise, ncol=2, widths = c(1,0.9))



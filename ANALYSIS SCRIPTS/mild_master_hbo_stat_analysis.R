library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(plyr)
library(car)
require(gridExtra)

#### Data Preparation ###########
# Load in Data
all_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\group_results.csv")
colnames(all_data)[which(names(all_data) == "Condition")] <- "Spatialization"

# Organize Factors
to.factor <- c('Spatialization','ID','ch_name')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)


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
all_data %>% group_by(Spatialization) %>% shapiro_test(mean_hbo)









##### LMEM Model HbO #############
model_hbo <- mixed(mean_hbo ~ Spatialization + (1|ID) + (1|ch_name),
                       data= subset(all_data, Chroma == "hbo"), 
                       control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_hbo



# Pairwise Comparisons (treatment coding)

# ITD5deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=5_az=0")
lmer_itd5deg <- lmer(mean_hbo ~ Spatialization + (1|ID) + (1|ch_name),
                   data= subset(all_data, Chroma == "hbo"), 
                              control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd5deg)

# ITD15deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=15_az=0")
lmer_itd15deg <- lmer(mean_hbo ~ Spatialization + (1|ID) + (1|ch_name),
                    data= subset(all_data, Chroma == "hbo"), 
                   control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd15deg)


# ILD5deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=0_az=5")
lmer_ild5deg <- lmer(mean_hbo ~ Spatialization + (1|ID) + (1|ch_name),
                     data= subset(all_data, Chroma == "hbo"), 
                    control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild5deg)

# ILD15deg as reference

all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=0_az=15")
lmer_ild15deg <- lmer(mean_hbo ~ Spatialization + (1|ID) + (1|ch_name),
                        data= subset(all_data, Chroma == "hbo"),
                     control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild15deg)


##### Mean HbO Plot ##########
mean_hbo_data <- summarySE(subset(all_data,Chroma == "hbo"), measurevar="mean_hbo", groupvars=c("ch_name","Spatialization"))
mean_hbo_data$Spatialization <- ordered(mean_hbo_data$Spatialization, levels = c("az_itd=5_az=0","az_itd=15_az=0","az_itd=0_az=5","az_itd=0_az=15"))
plothbo <- ggplot(data = mean_hbo_data, aes(x = Spatialization, y = mean_hbo,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(x=Spatialization, ymin=mean_hbo-se, ymax=mean_hbo+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(x = Spatialization, y = mean_hbo, color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean \u0394HbO (\u03BCM)") +
  ylim(-0.15,0.15) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")
grid.arrange(plothbo, ncol=1, widths = c(1))



##### LMEM Model Beta #############
model_beta <- mixed(theta ~ Spatialization + (1|ID) + (1| ch_name),
                    data= subset(all_data, Chroma == "hbo"), 
                    control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_beta

# ITD5deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=5_az=0")
lmer_itd5deg_beta <- lmer(theta ~ Spatialization + (1|ID) + (1|ch_name),
                     data= subset(all_data, Chroma == "hbo"), 
                     control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd5deg_beta)

# ITD15deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=15_az=0")
lmer_itd15deg_beta <- lmer(theta ~ Spatialization + (1|ID) + (1|ch_name),
                          data= subset(all_data, Chroma == "hbo"), 
                          control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_itd15deg_beta)

# ILD5deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=0_az=5")
lmer_ild5deg_beta <- lmer(theta ~ Spatialization + (1|ID) + (1|ch_name),
                          data= subset(all_data, Chroma == "hbo"), 
                          control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild5deg_beta)

# ILD15deg as reference
all_data$Spatialization <- relevel(all_data$Spatialization, "az_itd=0_az=15")
lmer_ild15deg_beta <- lmer(theta ~ Spatialization + (1|ID) + (1|ch_name),
                          data= subset(all_data, Chroma == "hbo"), 
                          control = lmerControl(optimizer = "bobyqa"))#
summary(lmer_ild15deg_beta)

##### Beta Plot ##########
beta_data <- summarySE(subset(all_data,Chroma == "hbo"), measurevar="theta", groupvars=c("Spatialization","ch_name"), na.rm = TRUE)
beta_data$Spatialization <- ordered(beta_data$Spatialization, levels = c("az_itd=5_az=0","az_itd=15_az=0","az_itd=0_az=5","az_itd=0_az=15"))
plot_beta <- ggplot(data = beta_data, aes(x = Spatialization, y = theta,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=theta-se, ymax=theta+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean Beta") +
  ylim(-0.1,0.2) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")

grid.arrange(plot_beta, ncol=1, widths = c(1))


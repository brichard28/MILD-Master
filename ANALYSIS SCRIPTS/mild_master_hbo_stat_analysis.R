library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(plyr)
library(car)
require(gridExtra)
library(emmeans)

#### Data Preparation ###########
# Load in Data
all_data <- read.csv("/Users/benrichardson/Documents/GitHub/MILD-Master/RESULTS DATA/group_results_glm_dur_6.csv")
colnames(all_data)[which(names(all_data) == "Condition")] <- "Spatialization"

# Organize Factors
to.factor <- c('Spatialization','ID','ch_name')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)

all_data_hbo <- subset(all_data,Chroma == "hbo")
all_data_hbr <- subset(all_data,Chroma == "hbr")

# Define ROIs
ch_names_unique <- unique(all_data_hbo$ch_name)

roi_left_hbo <- list("S8_D17 hbo","S8_D16 hbo","S9_D17 hbo","S9_D16 hbo","S9_D15 hbo","S10_D16 hbo","S16_D16 hbo","S16_D17 hbo",
                 "S16_D22 hbo","S16_D23 hbo","S17_D15 hbo","S17_D16 hbo","S17_D17 hbo","S17_D21 hbo",
                 "S17_D22 hbo","S17_D23 hbo","S18_D15 hbo","S18_D16 hbo","S18_D21 hbo","S18_D22 hbo","S18_D23 hbo",
                 "S1_D7 hbo","S1_D8 hbo","S2_D6 hbo","S2_D7 hbo","S2_D8 hbo","S3_D5 hbo",
                  "S3_D6 hbo","S3_D7 hbo","S7_D7 hbo","S7_D8 hbo","S7_D17 hbo","S7_D18 hbo",
                  "S8_D6 hbo","S8_D7 hbo","S8_D8 hbo","S8_D18 hbo","S9_D5 hbo","S9_D6 hbo","S9_D7 hbo",
                 "S10_D5 hbo","S10_D6 hbo","S10_D14 hbo","S10_D15 hbo","S10_D21 hbo","S15_D17 hbo",
                 "S15_D18 hbo","S15_D22 hbo","S16_D18 hbo","S19_D14 hbo","S19_D15 hbo","S19_D21 hbo") 

roi_right_hbo <- list("S11_D11 hbo","S12_D10 hbo","S12_D11 hbo","S12_D12 hbo","S13_D10 hbo","S13_D11 hbo","S21_D11 hbo",
              "S21_D12 hbo","S21_D19 hbo","S21_D20 hbo","S22_D10 hbo","S22_D11 hbo","S22_D12 hbo","S22_D19 hbo",
             "S22_D20 hbo","S23_D10 hbo","S23_D11 hbo","S23_D19 hbo","S4_D2 hbo","S4_D3 hbo","S4_D4 hbo","S5_D1 hbo","S5_D2 hbo","S5_D3 hbo",
              "S6_D1 hbo","S6_D2 hbo","S11_D3 hbo","S11_D4 hbo","S11_D12 hbo","S11_D13 hbo",
              "S11_D20 hbo","S12_D2 hbo","S12_D3 hbo","S12_D4 hbo","S13_D1 hbo","S13_D2 hbo",
              "S13_D3 hbo","S13_D9 hbo","S14_D1 hbo","S14_D2 hbo","S14_D9 hbo","S14_D10 hbo",
              "S20_D12 hbo","S20_D13 hbo","S20_D20 hbo","S23_D9 hbo","S24_D9 hbo","S24_D10 hbo","S24_D19 hbo") # right DLPFC
roi_left_hbr <- list("S8_D17 hbr","S8_D16 hbr","S9_D17 hbr","S9_D16 hbr","S9_D15 hbr","S10_D16 hbr","S16_D16 hbr","S16_D17 hbr",
                     "S16_D22 hbr","S16_D23 hbr","S17_D15 hbr","S17_D16 hbr","S17_D17 hbr","S17_D21 hbr",
                     "S17_D22 hbr","S17_D23 hbr","S18_D15 hbr","S18_D16 hbr","S18_D21 hbr","S18_D22 hbr","S18_D23 hbr",
                     "S1_D7 hbr","S1_D8 hbr","S2_D6 hbr","S2_D7 hbr","S2_D8 hbr","S3_D5 hbr",
                     "S3_D6 hbr","S3_D7 hbr","S7_D7 hbr","S7_D8 hbr","S7_D17 hbr","S7_D18 hbr",
                     "S8_D6 hbr","S8_D7 hbr","S8_D8 hbr","S8_D18 hbr","S9_D5 hbr","S9_D6 hbr","S9_D7 hbr",
                     "S10_D5 hbr","S10_D6 hbr","S10_D14 hbr","S10_D15 hbr","S10_D21 hbr","S15_D17 hbr",
                     "S15_D18 hbr","S15_D22 hbr","S16_D18 hbr","S19_D14 hbr","S19_D15 hbr","S19_D21 hbr") 

roi_right_hbr <- list("S11_D11 hbr","S12_D10 hbr","S12_D11 hbr","S12_D12 hbr","S13_D10 hbr","S13_D11 hbr","S21_D11 hbr",
                      "S21_D12 hbr","S21_D19 hbr","S21_D20 hbr","S22_D10 hbr","S22_D11 hbr","S22_D12 hbr","S22_D19 hbr",
                      "S22_D20 hbr","S23_D10 hbr","S23_D11 hbr","S23_D19 hbr","S4_D2 hbr","S4_D3 hbr","S4_D4 hbr","S5_D1 hbr","S5_D2 hbr","S5_D3 hbr",
                      "S6_D1 hbr","S6_D2 hbr","S11_D3 hbr","S11_D4 hbr","S11_D12 hbr","S11_D13 hbr",
                      "S11_D20 hbr","S12_D2 hbr","S12_D3 hbr","S12_D4 hbr","S13_D1 hbr","S13_D2 hbr",
                      "S13_D3 hbr","S13_D9 hbr","S14_D1 hbr","S14_D2 hbr","S14_D9 hbr","S14_D10 hbr",
                      "S20_D12 hbr","S20_D13 hbr","S20_D20 hbr","S23_D9 hbr","S24_D9 hbr","S24_D10 hbr","S24_D19 hbr") # right DLPFC

all_data_hbo$Roi<- "NA"
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_left_hbo)] <- "Left"
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_right_hbo)] <- "Right"

all_data_hbr$Roi <- "NA"
all_data_hbr$Roi[which(all_data_hbr$ch_name %in% roi_left_hbr)] <- "Left"
all_data_hbr$Roi[which(all_data_hbr$ch_name %in% roi_right_hbr)] <- "Right"


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






### Plot data in each ROI #######
hbo_beta_data <- summarySE(all_data_hbo, measurevar="theta", groupvars=c("Spatialization","ch_name","Roi"), na.rm = TRUE)
hbo_beta_data$Spatialization <- ordered(hbo_beta_data$Spatialization, levels = c("az_itd=5_az=0","az_itd=15_az=0","az_itd=0_az=5","az_itd=0_az=15"))

hbr_beta_data <- summarySE(all_data_hbr, measurevar="theta", groupvars=c("Spatialization","ch_name","Roi"), na.rm = TRUE)
hbr_beta_data$Spatialization <- ordered(hbr_beta_data$Spatialization, levels = c("az_itd=5_az=0","az_itd=15_az=0","az_itd=0_az=5","az_itd=0_az=15"))

plot_roi_hbo_left <- ggplot(data = subset(hbo_beta_data, Roi == "Left"), aes(x = Spatialization, y = theta,group = Spatialization)) +
  geom_violin(fill = "red", color = "white") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black", size = 1.5) +
  labs(x="",y="Mean Beta") +
  ylim(-0.055,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "Small\nITDs", "az_itd=15_az=0" = "Large\nITDs","az_itd=0_az=5" = "Small\nILDs","az_itd=0_az=15" = "Large\nILDs")) +
  theme(legend.position="none")

plot_roi_hbo_right <- ggplot(data = subset(hbo_beta_data, Roi == "Right"), aes(x = Spatialization, y = theta,group = Spatialization)) +
  geom_violin(fill = "red", color = "white") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black", size = 1.5) +
  labs(x="",y="Mean Beta") +
  ylim(-0.055,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "Small\nITDs", "az_itd=15_az=0" = "Large\nITDs","az_itd=0_az=5" = "Small\nILDs","az_itd=0_az=15" = "Large\nILDs")) +
  theme(legend.position="none")

plot_roi_hbr_left <- ggplot(data = subset(hbr_beta_data, Roi == "Left"), aes(x = Spatialization, y = theta,group = Spatialization)) +
  geom_violin(fill = "blue", color = "white") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black", size = 1.5) +
  labs(x="",y="Mean Beta") +
  ylim(-0.055,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "Small\nITDs", "az_itd=15_az=0" = "Large\nITDs","az_itd=0_az=5" = "Small\nILDs","az_itd=0_az=15" = "Large\nILDs")) +
  theme(legend.position="none")

plot_roi_hbr_right <- ggplot(data = subset(hbr_beta_data, Roi == "Right"), aes(x = Spatialization, y = theta,group = Spatialization)) +
  geom_violin(fill = "blue", color = "white") +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black", size = 1.5) +
  labs(x="",y="Mean Beta") +
  ylim(-0.055,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "Small\nITDs", "az_itd=15_az=0" = "Large\nITDs","az_itd=0_az=5" = "Small\nILDs","az_itd=0_az=15" = "Large\nILDs")) +
  theme(legend.position="none")


fnirs_raw_beta_plot <- grid.arrange(plot_roi_hbo_left,plot_roi_hbo_right,plot_roi_hbr_left,plot_roi_hbr_right, ncol=2, nrow=2, widths = c(1,1))

ggsave("/Users/benrichardson/Documents/GitHub/MILD-Master/PAPER FIGURES/fnirs_beta_raw.svg", plot = fnirs_raw_beta_plot, width = 8, height = 6, units = "in")


## LMEM Tests beta

# HbO
# Left Hemisphere
model_roi_left_beta_hbo <- mixed(theta ~ Spatialization + (1|ID) + (1|ch_name),
                          data= subset(all_data_hbo, Roi == "Left"), 
                          control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_roi_left_beta_hbo

# Significant effect of Spatialization post hoc
emm_spatialization <- emmeans(model_roi_left_beta_hbo$full_model, ~ Spatialization)
pairs(emm_spatialization, adjust = "bonferroni")


# Right Hemisphere
model_roi_right_beta_hbo <- mixed(theta ~ Spatialization + (1|ID) + (1|ch_name),
                             data= subset(all_data_hbo, Roi == "Right"), 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_roi_right_beta_hbo

# Significant effect of Spatialization post hoc
emm_spatialization <- emmeans(model_roi_right_beta_hbo$full_model, ~ Spatialization)
pairs(emm_spatialization, adjust = "bonferroni")





# HbR
model_roi_left_beta_hbr <- mixed(theta ~ Spatialization + (1|ID) + (1|ch_name),
                             data= subset(all_data_hbr, Roi == "Left"), 
                             control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_roi_left_beta_hbr


model_roi_right_beta_hbr <- mixed(theta ~ Spatialization + (1|ID) + (1|ch_name),
                                 data= subset(all_data_hbr, Roi == "Right"), 
                                 control = lmerControl(optimizer = "bobyqa"), method = 'LRT')
model_roi_right_beta_hbr

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
all_data <- read.csv("C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\group_results_glm_dur_11.csv")
colnames(all_data)[which(names(all_data) == "Condition")] <- "Spatialization"

# Organize Factors
to.factor <- c('Spatialization','ID','ch_name')
all_data[, to.factor] <- lapply(all_data[, to.factor], as.factor)

all_data_hbo <- subset(all_data,Chroma == "hbo")
# Define ROIs
ch_names_unique <- unique(all_data_hbo$ch_name)
roi_1 <- list("S8_D17 hbo","S8_D16 hbo","S9_D17 hbo","S9_D16 hbo","S9_D15 hbo","S10_D16 hbo","S16_D16 hbo","S16_D17 hbo",
              "S16_D22 hbo","S16_D23 hbo","S17_D15 hbo","S17_D16 hbo","S17_D17 hbo","S17_D21 hbo",
              "S17_D22 hbo","S17_D23 hbo","S18_D15 hbo","S18_D16 hbo","S18_D21 hbo","S18_D22 hbo","S18_D23 hbo") # left IFG

roi_2 <- list("S11_D11 hbo","S12_D10 hbo","S12_D11 hbo","S12_D12 hbo","S13_D10 hbo","S13_D11 hbo","S21_D11 hbo",
              "S21_D12 hbo","S21_D19 hbo","S21_D20 hbo","S22_D10 hbo","S22_D11 hbo","S22_D12 hbo","S22_D19 hbo",
              "S22_D20 hbo","S23_D10 hbo","S23_D11 hbo","S23_D19 hbo")  # right IFG

roi_3 <- list("S1_D7 hbo","S1_D8 hbo","S2_D6 hbo","S2_D7 hbo","S2_D8 hbo","S3_D5 hbo",
              "S3_D6 hbo","S3_D7 hbo","S7_D7 hbo","S7_D8 hbo","S7_D17 hbo","S7_D18 hbo",
              "S8_D6 hbo","S8_D7 hbo","S8_D8 hbo","S8_D18 hbo","S9_D5 hbo","S9_D6 hbo","S9_D7 hbo",
              "S10_D5 hbo","S10_D6 hbo","S10_D14 hbo","S10_D15 hbo","S10_D21 hbo","S15_D17 hbo",
              "S15_D18 hbo","S15_D22 hbo","S16_D18 hbo","S19_D14 hbo","S19_D15 hbo","S19_D21 hbo") # left DLPFC

roi_4 <- list("S4_D2 hbo","S4_D3 hbo","S4_D4 hbo","S5_D1 hbo","S5_D2 hbo","S5_D3 hbo",
              "S6_D1 hbo","S6_D2 hbo","S11_D3 hbo","S11_D4 hbo","S11_D12 hbo","S11_D13 hbo",
              "S11_D20 hbo","S12_D2 hbo","S12_D3 hbo","S12_D4 hbo","S13_D1 hbo","S13_D2 hbo",
              "S13_D3 hbo","S13_D9 hbo","S14_D1 hbo","S14_D2 hbo","S14_D9 hbo","S14_D10 hbo",
              "S20_D12 hbo","S20_D13 hbo","S20_D20 hbo","S23_D9 hbo","S24_D9 hbo","S24_D10 hbo",
              "S24_D19 hbo") # right DLPFC

all_data_hbo$Roi<- "NA"
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_1)] <- 1
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_2)] <- 2
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_3)] <- 3
all_data_hbo$Roi[which(all_data_hbo$ch_name %in% roi_4)] <- 4


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


### Plot data in each ROI ###
beta_data <- summarySE(all_data_hbo, measurevar="theta", groupvars=c("Spatialization","ch_name","Roi"), na.rm = TRUE)
beta_data$Spatialization <- ordered(beta_data$Spatialization, levels = c("az_itd=5_az=0","az_itd=15_az=0","az_itd=0_az=5","az_itd=0_az=15"))
plot_roi_1 <- ggplot(data = subset(beta_data, Roi == 1), aes(x = Spatialization, y = theta,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=theta-se, ymax=theta+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean Beta") +
  ylim(-0.05,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")

plot_roi_2 <- ggplot(data = subset(beta_data, Roi == 2), aes(x = Spatialization, y = theta,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=theta-se, ymax=theta+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean Beta") +
  ylim(-0.05,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")

plot_roi_3 <- ggplot(data = subset(beta_data, Roi == 3), aes(x = Spatialization, y = theta,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=theta-se, ymax=theta+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean Beta") +
  ylim(-0.05,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")

plot_roi_4 <- ggplot(data = subset(beta_data, Roi == 4), aes(x = Spatialization, y = theta,group = ch_name)) +
  geom_line(size=1, position=position_dodge(width=0.3)) + 
  geom_errorbar(aes(ymin=theta-se, ymax=theta+se, color=Spatialization), width=.1, position=position_dodge(width=0.3)) +
  geom_point(aes(color = Spatialization),size = 4, position=position_dodge(width=0.3)) + 
  labs(x="",y="Mean Beta") +
  ylim(-0.05,0.07) +
  theme_bw() +
  theme(plot.title = element_text(size = 18), axis.title=element_text(size=18), axis.text.x= element_text(size=12), axis.text.y= element_text(size=12)) +
  scale_x_discrete(labels=c("az_itd=5_az=0" = "5 deg\nITDs", "az_itd=15_az=0" = "15deg\nITDs","az_itd=0_az=5" = "5deg\nILDs","az_itd=0_az=15" = "15deg\nILDs")) +
  theme(legend.position="none")

grid.arrange(plot_roi_1,plot_roi_2,plot_roi_3,plot_roi_4, ncol=2, widths = c(1,1))



##### Wilcoxon signed-rank test between spatialization, within each ROI ##
beta_data <- summarySE(all_data_hbo, measurevar="theta", groupvars=c("Spatialization","Roi"), na.rm = TRUE)
roi_1_data <- subset(all_data_hbo, Roi == 1)
foo <- pairwise.wilcox.test(roi_1_data$theta, roi_1_data$Spatialization, p.adjust.method="bonferroni")
foo
broom::tidy(foo)

wilcox.test(roi_1_theta, mu = 0, alternative = "two.sided")
# follow up tests for ROI 1
result <- wilcox.test(before, after, paired = TRUE)
print(result)

wilcox.test(subset(beta_data, Roi == 2)$theta, mu = 0, alternative = "two.sided")
wilcox.test(subset(beta_data, Roi == 3)$theta, mu = 0, alternative = "two.sided")
wilcox.test(subset(beta_data, Roi == 4)$theta, mu = 0, alternative = "two.sided")


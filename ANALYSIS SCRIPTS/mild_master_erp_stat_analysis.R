# Author: Benjamin Richardson
# plot and conduct statistics for MILD-Master p1,n1, and p3 components experiment 1

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(rstatix)
library(afex)
library(dplyr)
library(rhdf5)
# Read the .mat file
library(R.matlab)
p1_data <- read.csv('C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\all_subs_p1.csv')
n1_data <- read.csv('C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\all_subs_n1.csv')
p3_data <- read.csv('C:\\Users\\benri\\Documents\\GitHub\\MILD-Master\\RESULTS DATA\\all_subs_p3.csv')

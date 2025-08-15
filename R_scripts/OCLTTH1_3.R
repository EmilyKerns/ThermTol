# Testing OCLTTH Predictions 1-3 at the start & end of the experiment
# Figure 3, Supplemental Figures 1 & 2


# Load Packages -----------------------------------------------------------

library(car)
library(ggplot2)
library(MASS)
library(rstatix)
library(plyr)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(FSA)
library(dunn.test)
library(emmeans)
library(lmtest)
library(patchwork)
library(AER)
library(viridisLite)
library(viridis)
library(modelsummary)
library(readxl)
library(permuco)
library(pgirmess)


# Load data ---------------------------------------------------------------


ThermTol <- read_xlsx("R:/ThermTol_MultStress/Sasser_Ch3/Ch3/ThermalToleranceProject_SuspectDataRemoved.xlsx")
head(ThermTol)

ThermTol$Ecotype <- as.factor(ThermTol$Ecotype)
ThermTol$Fibrosis <- as.factor(ThermTol$Fibrosis)
ThermTol$Trial <- as.factor(ThermTol$Trial)

EcotypeMortality <- read_xlsx("R:/ThermTol_MultStress/RawData/EcotypeMortality.xlsx")
head(EcotypeMortality)

LakeData <- read_xlsx("R:/ThermTol_MultStress/RawData/depth_profiles.xlsx")
LakeData <- subset(LakeData, Lake == "Finger" | Lake == "Watson" | Lake == "Spirit" | Lake == "Wik")

LakeData <- LakeData %>% mutate(
  Lake = case_when(
    Lake == "Finger" ~ "FG",
    Lake == "Spirit" ~ "SL",
    Lake == "Watson" ~ "WT",
    Lake == "Wik" ~ "WK",
    TRUE ~ "Yer missing something"  # A default value if none of the conditions match
  ))

colnames(LakeData)[3] <- "DepthM"
colnames(LakeData)[7] <- "SA"
colnames(LakeData)[1] <- "Pop"
head(LakeData)

LakeData <- LakeData %>% 
  group_by(Pop) %>% 
  summarise(MinDepth = max(DepthM, na.rm = TRUE), 
            MaxDepth = min(DepthM, na.rm = TRUE),
            MedDepth = median(DepthM, na.rm = TRUE),
            MeanDepth = mean(DepthM, na.rm = TRUE),
            SA = first(SA),
            .groups = 'drop')

ThermTol <- merge.data.frame(ThermTol, LakeData, by = "Pop")
head(ThermTol)


# Relabel data ------------------------------------------------------------


ThermTol$Trial <- factor(ThermTol$Trial, levels = c("S", "F"), 
                         labels = c("Start", "End"))

ThermTolTAF <- ThermTol 
ThermTolTAF$Temp <- as.factor(ThermTolTAF$Temp)

EcotypeMortalityTAF <- EcotypeMortality 
EcotypeMortalityTAF$Temp <- as.factor(EcotypeMortalityTAF$Temp)


# Clean up health metrics -------------------------------------------------


hist(ThermTol$SpleenMass)
ThermTol$SpleenMass <- log(ThermTol$SpleenMass)
hist(ThermTol$SpleenMass)

hist(ThermTol$LiverMass)
ThermTol$LiverMass <- log(ThermTol$LiverMass)
hist(ThermTol$LiverMass)

hist(ThermTol$GonadMass)
ThermTol$GonadMass <- log(ThermTol$GonadMass)
hist(ThermTol$GonadMass)

hist(ThermTol$Mass)
ThermTol$LogMass <- log(ThermTol$Mass)
hist(ThermTol$Mass)

hist(ThermTol$SL)
ThermTol$LogSL <- log(ThermTol$SL)
hist(ThermTol$LogSL)

ThermTol$BodyCond <- ((ThermTol$Mass/(ThermTol$SL)^3)*100)
hist(ThermTol$BodyCond)

hist(ThermTol$HSI)
ThermTol$HSI <- (ThermTol$LiverMass/ThermTol$LogMass)
hist(ThermTol$HSI)

hist(ThermTol$SSI)
ThermTol$SSI <- (ThermTol$SpleenMass/ThermTol$LogMass)
hist(ThermTol$SSI)

hist(ThermTol$GSI)
ThermTol$GSI <- (ThermTol$GonadMass/ThermTol$LogMass)
hist(ThermTol$GSI)

## Subset data
Start <-subset(ThermTol, Trial== "Start" )
End <-subset(ThermTol, Trial == "End" )
StartTAF <-subset(ThermTolTAF, Trial== "Start" )
EndTAF <-subset(ThermTolTAF, Trial == "End" )
Limno <-subset(ThermTol, Ecotype== "Limnetic" )
Benno <-subset(ThermTol, Ecotype== "Benthic" )


# Test the OCLTTH ---------------------------------------------------------


### Test the OCLTT hypothesis. Assumptions: 1) As temp increases, MMR will stay steady. 2) As temp increases, RMR will increase. 3) As temp increases, AS will decrease.


# Start of experiment ---------------------------------------


# Assumption 1 ------------------------------------------------------------

hist(Start$MMR)
MMRL <- glm(MMR ~ Temp*Ecotype, data = Start)
summary(MMRL)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              1.49     295.50   0.005    0.996    
# Temp                    84.66      13.46   6.291 2.27e-09 ***
# EcotypeLimnetic        529.19     418.76   1.264    0.208    
# Temp:EcotypeLimnetic   -19.97      19.05  -1.049    0.296    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 126726.1)
# 
# Null deviance: 31583845  on 186  degrees of freedom
# Residual deviance: 23190881  on 183  degrees of freedom
# (25 observations deleted due to missingness)
# AIC: 2733.8
MMRL <- lm(MMR ~ Temp*Ecotype, data = Start)
summary(MMRL)
# Multiple R-squared:  0.2657,	Adjusted R-squared:  0.2537 
# F-statistic: 22.08 on 3 and 183 DF,  p-value: 3.013e-12

MMRQ <- glm(MMR ~ I(Temp^2) + Temp*Ecotype, data = Start)
summary(MMRQ)
# AIC: 2735.4

MMRpoint <- ggplot(Start, aes(x=as.factor(Temp), y=MMR, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), 
              method = "glm", 
              formula = y ~ x, 
              size = .5,
              se = TRUE, 
              alpha = .2)   +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
MMRpoint

# By lake attribute instead of ecotype
MMRMD <- glm(MMR ~ Temp*MaxDepth, data = Start)
summary(MMRMD)
# AIC: 2732.5

MMRMDQ <- glm(MMR ~ I(Temp^2) + Temp*MaxDepth, data = Start)
summary(MMRMDQ)
# AIC: 2734.2

MMRAvgD <- glm(MMR ~ Temp*MeanDepth, data = Start)
summary(MMRAvgD)
# AIC: 2732.2

MMRAvgDQ <- glm(MMR ~ I(Temp^2) + Temp*MeanDepth, data = Start)
summary(MMRAvgDQ)
# AIC: 2733.9

MMRSA <- glm(MMR ~ Temp*SA, data = Start)
summary(MMRSA)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) -1587.166    855.877  -1.854   0.0653 .  
#   Temp          163.665     38.830   4.215 3.92e-05 ***
#   SA            190.222     83.476   2.279   0.0238 *  
#   Temp:SA        -9.162      3.829  -2.393   0.0177 *  
# 
# (Dispersion parameter for gaussian family taken to be 125368)
# 
# Null deviance: 31583845  on 186  degrees of freedom
# Residual deviance: 22942346  on 183  degrees of freedom
# AIC: 2731.8
MMRSAL <- lm(MMR ~ Temp*SA, data = Start)
summary(MMRSAL)
# Multiple R-squared:  0.2736,	Adjusted R-squared:  0.2617 
# F-statistic: 22.98 on 3 and 183 DF,  p-value: 1.14e-12

MMRSA <- glm(MMR ~ I(Temp^2) + Temp*SA, data = Start)
summary(MMRSA)
# AIC: 2731.9

MMRSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=MMR, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), 
              method = "glm", 
              formula = y ~ x, 
              size = 1,
              se = TRUE, 
              alpha = .1)   +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
MMRSApoint

# Assumption 2 ------------------------------------------------------------

Start$RMR <- as.numeric(Start$RMR)

hist(Start$RMR)
RMRL <- glm(RMR~Temp*Ecotype, data = Start)
summary(RMRL)
# AIC: 2519.5


RMRQ <- glm(RMR~ I(Temp^2) + Temp*Ecotype, data = Start)
summary(RMRQ)
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)  
# (Intercept)          1682.020   1038.256   1.620   0.1070  
# I(Temp^2)               4.904      2.211   2.218   0.0278 *
# Temp                 -147.953     96.293  -1.536   0.1262  
# EcotypeLimnetic      -519.283    234.608  -2.213   0.0281 *
# Temp:EcotypeLimnetic   24.166     10.667   2.266   0.0247 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 39423.41)
# 
# Null deviance: 16146914  on 186  degrees of freedom
# Residual deviance:  7175061  on 182  degrees of freedom
# (25 observations deleted due to missingness)
# AIC: 2516.5
RMRQ <- lm(RMR~ I(Temp^2) + Temp*Ecotype, data = Start)
summary(RMRQ)
# Multiple R-squared:  0.5556,	Adjusted R-squared:  0.5459 
# F-statistic: 56.89 on 4 and 182 DF,  p-value: < 2.2e-16

RMRpoint <- ggplot(Start, aes(x=as.factor(Temp), y=RMR, color=Ecotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 0.5), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual(values=c('black','blue')) +
  labs(x = expression("Temperature"*degree*C), y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr))) +
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = .5,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
RMRpoint

# By lake attribute instead of ecotype
RMRMD <- glm(RMR ~ Temp*MaxDepth, data = Start)
summary(RMRMD)
# AIC: 2517.3

RMRMDQ <- glm(RMR ~ I(Temp^2) + Temp*MaxDepth, data = Start)
summary(RMRMDQ)
# AIC: 2513.9

RMRAvgD <- glm(RMR ~ Temp*MeanDepth, data = Start)
summary(RMRAvgD)
# AIC: 2516.5

RMRAvgDQ <- glm(RMR ~ I(Temp^2) + Temp*MeanDepth, data = Start)
summary(RMRAvgDQ)
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)    2091.728   1040.914   2.010  0.04596 * 
# I(Temp^2)         5.061      2.196   2.304  0.02235 * 
# Temp           -170.615     95.767  -1.782  0.07649 . 
# MeanDepth        84.792     31.718   2.673  0.00819 **
# Temp:MeanDepth   -3.976      1.419  -2.803  0.00561 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 38725.82)
# 
# Null deviance: 16146914  on 186  degrees of freedom
# Residual deviance:  7048099  on 182  degrees of freedom
# (25 observations deleted due to missingness)
# AIC: 2513.1
RMRAvgDQ <- lm(RMR ~ I(Temp^2) + Temp*MeanDepth, data = Start)
summary(RMRAvgDQ)
# Multiple R-squared:  0.5635,	Adjusted R-squared:  0.5539 
# F-statistic: 58.74 on 4 and 182 DF,  p-value: < 2.2e-16

RMRSA <- glm(RMR ~ Temp*SA, data = Start)
summary(RMRSA)
# AIC: 2523.1

RMRSA <- glm(RMR ~ I(Temp^2) + Temp*SA, data = Start)
summary(RMRSA)
# AIC: 2517.3

RMRSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=RMR, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = 1,se = TRUE,  alpha = .1)   +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
RMRSApoint


# Assumption 3 ------------------------------------------------------------

hist(Start$AS)

ASL <- glm(AS ~ Temp*Ecotype, data = Start)
summary(ASL)
# AIC: 2708.7

ASQ <- glm(AS ~ I(Temp^2) +Temp*Ecotype, data = Start)
summary(ASQ)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          -2881.683   1725.610  -1.670  0.09665 . 
# I(Temp^2)               -7.495      3.675  -2.040  0.04284 * 
# Temp                   345.105    160.042   2.156  0.03237 * 
# EcotypeLimnetic       1022.575    389.925   2.622  0.00947 **
# Temp:EcotypeLimnetic   -43.002     17.728  -2.426  0.01626 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 108900.6)
# 
# Null deviance: 21393842  on 186  degrees of freedom
# Residual deviance: 19819918  on 182  degrees of freedom
# (25 observations deleted due to missingness)
# AIC: 2706.5

ASQ <- lm(AS ~ I(Temp^2) +Temp*Ecotype, data = Start)
summary(ASQ)
# Multiple R-squared:  0.07357,	Adjusted R-squared:  0.05321 
# F-statistic: 3.613 on 4 and 182 DF,  p-value: 0.007348

ASpoint <- ggplot(Start, aes(x=as.factor(Temp), y=AS, color=Ecotype)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 0.5), alpha= 0.1) +
  theme_classic() +
  scale_color_manual(values=c('black','blue')) +
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr))) +
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = .5,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
ASpoint

# By lake attribute instead of ecotype
ASMD <- glm(AS ~ Temp*MaxDepth, data = Start)
summary(ASMD)
# AIC: 2705.2

ASMDQ <- glm(AS ~ I(Temp^2) + Temp*MaxDepth, data = Start)
summary(ASMDQ)
# AIC: 2703.1

ASAvgD <- glm(AS ~ Temp*MeanDepth, data = Start)
summary(ASAvgD)
# AIC: 2704.5

ASAvgDQ <- glm(AS ~ I(Temp^2) + Temp*MeanDepth, data = Start)
summary(ASAvgDQ)
# AIC: 2702.6

ASSA <- glm(AS ~ Temp*SA, data = Start)
summary(ASSA)
# AIC: 2707.3

ASSA <- glm(AS ~ I(Temp^2) + Temp*SA, data = Start)
summary(ASSA)
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) -6729.001   1995.911  -3.371 0.000913 ***
#   I(Temp^2)     -11.742      3.674  -3.196 0.001645 ** 
#   Temp          627.455    169.944   3.692 0.000294 ***
#   SA            247.104     77.236   3.199 0.001625 ** 
#   Temp:SA       -12.372      3.553  -3.482 0.000623 ***
#   --
# (Dispersion parameter for gaussian family taken to be 104704.1)
# 
# Null deviance: 21393842  on 186  degrees of freedom
# Residual deviance: 19056145  on 182  degrees of freedom
# AIC: 2699.1
ASSA <- lm(AS ~ I(Temp^2) + Temp*SA, data = Start)
summary(ASSA)
# Multiple R-squared:  0.1093,	Adjusted R-squared:  0.08969 
# F-statistic: 5.582 on 4 and 182 DF,  p-value: 0.0002925

ASSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=AS, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = 1,se = TRUE,  alpha = .1)   +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
ASSApoint

# End of Experiment -------------------------------------------------------


# Assumption 1 ------------------------------------------------------------

MMREndL <- glm(MMR ~ Temp*Ecotype + Sex, data = End)
summary(MMREndL)
# AIC: 1049.4

MMREndQ <- glm(MMR ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(MMREndQ)
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          13225.46    4223.10   3.132  0.00264 **
# I(Temp^2)               26.39       8.97   2.942  0.00456 **
# Temp                 -1125.75     393.12  -2.864  0.00568 **
# EcotypeLimnetic       1064.58    1171.45   0.909  0.36693   
# SexM                    30.33     108.54   0.279  0.78082   
# Temp:EcotypeLimnetic   -38.83      58.58  -0.663  0.50978   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 190899.9)
# 
# Null deviance: 15999694  on 68  degrees of freedom
# Residual deviance: 12026691  on 63  degrees of freedom
# (105 observations deleted due to missingness)
# AIC: 1042.5
MMRQ<- lm(MMR ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(MMRQ)
# Multiple R-squared:  0.2483,	Adjusted R-squared:  0.1887 
# F-statistic: 4.162 on 5 and 63 DF,  p-value: 0.002488

hist(End$MMR)
MMREnd <- ggplot(End, aes(x=factor(Temp), y=MMR, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = .5,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
MMREnd

MMRSex <- ggplot(End, aes(x=factor(Temp), y=MMR, color=Sex)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 2, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = 1,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
MMRSex


# Assumption 2 ------------------------------------------------------------


End$RMR <- as.numeric(End$RMR)
hist(End$RMR)
shapiro.test(End$RMR)
RMREndL <- glm(RMR~Temp*Ecotype + Sex, data = End)
summary(RMREndL)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -292.319    185.073  -1.579 0.119156    
# Temp                   46.064      8.367   5.505 6.97e-07 ***
# EcotypeLimnetic      1364.198    371.917   3.668 0.000499 ***
# SexM                  -88.253     41.688  -2.117 0.038153 *  
# Temp:EcotypeLimnetic  -67.094     18.275  -3.671 0.000494 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 28649.92)
# 
# Null deviance: 2916751  on 68  degrees of freedom
# Residual deviance: 1833595  on 64  degrees of freedom
# (105 observations deleted due to missingness)
# AIC: 910.76

RMREndQ <- glm(RMR ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(RMREndQ)

# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          4589.363   1528.364   3.003  0.00383 **
# I(Temp^2)              10.436      3.246   3.215  0.00206 **
# Temp                 -410.611    142.273  -2.886  0.00534 **
# EcotypeLimnetic       583.220    423.953   1.376  0.17380   
# SexM                  -71.758     39.281  -1.827  0.07247 . 
# Temp:EcotypeLimnetic  -26.694     21.199  -1.259  0.21260   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 25003.25)
# 
# Null deviance: 2916751  on 68  degrees of freedom
# Residual deviance: 1575205  on 63  degrees of freedom
# (105 observations deleted due to missingness)
# AIC: 902.28
RMREndQ <- lm(RMR ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(RMREndQ)
# Multiple R-squared:  0.4599,	Adjusted R-squared:  0.4171 
# F-statistic: 10.73 on 5 and 63 DF,  p-value: 1.731e-07

RMREnd <- ggplot(End, aes(x=factor(Temp), y=RMR, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = .5,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
RMREnd

RMRSex <- ggplot(End, aes(x=factor(Temp), y=RMR, color=Sex)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 2, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = 1,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
RMRSex

# Assumption 3 ------------------------------------------------------------


hist(End$AS)
shapiro.test(End$AS)
ASEndL <- glm(AS ~ Temp*Ecotype + Sex, data = End)
summary(ASEndL)
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)           1172.41     373.74   3.137  0.00258 **
# Temp                   -16.92      16.90  -1.002  0.32033   
# EcotypeLimnetic       1675.42     751.05   2.231  0.02921 * 
# SexM                    76.87      84.19   0.913  0.36462   
# Temp:EcotypeLimnetic   -73.91      36.90  -2.003  0.04945 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 116835.4)
# 
# Null deviance: 9655681  on 68  degrees of freedom
# Residual deviance: 7477465  on 64  degrees of freedom
# (105 observations deleted due to missingness)
# AIC: 1007.8

ASEndQ <- glm(AS ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(ASEndQ)
# Coefficients:
#                      Estimate Std. Error t value Pr(>|t|)   
# (Intercept)          8636.095   3192.610   2.705  0.00877 **
# I(Temp^2)              15.956      6.781   2.353  0.02176 * 
# Temp                 -715.141    297.196  -2.406  0.01906 * 
# EcotypeLimnetic       481.364    885.598   0.544  0.58867   
# SexM                  102.089     82.055   1.244  0.21805   
# Temp:EcotypeLimnetic  -12.139     44.282  -0.274  0.78488   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for gaussian family taken to be 109102.5)
# 
# Null deviance: 9655681  on 68  degrees of freedom
# Residual deviance: 6873456  on 63  degrees of freedom
# (105 observations deleted due to missingness)
# AIC: 1003.9
ASEndQ <- lm(AS ~ I(Temp^2) +Temp*Ecotype + Sex, data = End)
summary(ASEndQ)
# Multiple R-squared:  0.2881,	Adjusted R-squared:  0.2316 
# F-statistic:   5.1 on 5 and 63 DF,  p-value: 0.0005463

ASEnd <- ggplot(End, aes(x=factor(Temp), y=AS, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = .5,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
ASEnd

ASSex <- ggplot(End, aes(x=factor(Temp), y=AS, color=Sex)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 2, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), size = 1,se = TRUE,  alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
ASSex


#pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/Figure3.pdf",
#    width = 6.5, 
#    height = 7)

MMRpoint + RMRpoint + ASpoint + MMREnd + RMREnd + ASEnd

#dev.off()

pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/SuppFigure2.pdf",
     width = 6.6, 
     height = 4)

MMRSex + RMRSex + ASSex

 dev.off()

pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/SuppFigure1.pdf",
     width = 6.5, 
    height = 4)

MMRSApoint + RMRSApoint + ASSApoint

dev.off()


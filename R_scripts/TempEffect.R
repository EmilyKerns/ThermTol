# Direct effect of temperature on phenotypes
# Figure 5, Supplemental Figure 4

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

# Direct effect of temperature --------------------------------------------


model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(BodyCond))
BodyCondTemp <- lme(BodyCond ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = model_data)
summary(BodyCondTemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC       BIC   logLik
# -2375.453 -2353.585 1194.726
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 9.070103e-05 0.0001748276
# 
# Fixed effects:  BodyCond ~ Temp * Ecotype + Sex 
#                             Value    Std.Error  DF   t-value p-value
# (Intercept)           0.0013195307 0.0002556926 152  5.160614  0.0000
# Temp                 -0.0000090925 0.0000116955 152 -0.777436  0.4381
# EcotypeLimnetic      -0.0009604431 0.0003974127 152 -2.416740  0.0168
# SexM                  0.0002204429 0.0000275546 152  8.000231  0.0000
# Temp:EcotypeLimnetic  0.0000389244 0.0000184629 152  2.108251  0.0366
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.990                     
# EcotypeLimnetic      -0.543  0.552              
# SexM                 -0.007 -0.041 -0.156       
# Temp:EcotypeLimnetic  0.543 -0.563 -0.992  0.152
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.6442347 -0.6359798  0.0169038  0.6566455  2.7124727 
# 
# Number of Observations: 173
# Number of Groups: 17 

BCPLOT <-  ggplot(End, aes(x=Temp, y=BodyCond, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm",linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(Body~Condition)) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
BCPLOT

BCPLOTsex <-  ggplot(End, aes(x=Temp, y=BodyCond, color=Sex)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm",linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(Body~Condition)) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
BCPLOTsex

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(HSI))
HSITemp <- lme(HSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = model_data)
summary(HSITemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 634.4166 656.2005 -310.2083
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0001156835 1.453642
# 
# Fixed effects:  HSI ~ Temp * Ecotype + Sex 
#                         Value Std.Error  DF    t-value p-value
# (Intercept)          -3.224274 1.1477211 150 -2.8092837  0.0056
# Temp                  0.054763 0.0523158 150  1.0467729  0.2969
# EcotypeLimnetic      -0.570144 2.0117925 150 -0.2834009  0.7773
# SexM                 -0.028812 0.2266585 150 -0.1271186  0.8990
# Temp:EcotypeLimnetic  0.012043 0.0940031 150  0.1281136  0.8982
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.987                     
# EcotypeLimnetic      -0.559  0.570              
# SexM                 -0.059 -0.035 -0.170       
# Temp:EcotypeLimnetic  0.541 -0.562 -0.994  0.172
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -5.3959168 -0.2029609  0.2594957  0.5557691  1.0581642 
# 
# Number of Observations: 171
# Number of Groups: 17 

HSIPLOT <-  ggplot(End, aes(x=Temp, y=HSI, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(HSI)) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
HSIPLOT

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(GSI))
GSITemp <- lme(GSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = model_data)
summary(GSITemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 369.3534 385.4818 -177.6767
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:    0.860918 2.258998
# 
# Fixed effects:  GSI ~ Temp * Ecotype + Sex 
#                           Value Std.Error DF    t-value p-value
# (Intercept)          -5.076972  4.121638 59 -1.2317851  0.2229
# Temp                  0.027140  0.172473 59  0.1573573  0.8755
# EcotypeLimnetic      -1.539633  6.262843 59 -0.2458361  0.8067
# SexM                 -0.625562  1.240305 59 -0.5043613  0.6159
# Temp:EcotypeLimnetic  0.046735  0.292158 59  0.1599646  0.8735
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.952                     
# EcotypeLimnetic      -0.555  0.584              
# SexM                 -0.365  0.085  0.022       
# Temp:EcotypeLimnetic  0.533 -0.569 -0.994 -0.024
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -5.4667617 -0.2996948  0.2028616  0.4772907  1.3183379 
# 
# Number of Observations: 79
# Number of Groups: 16 

GSIPLOT <-  ggplot(End, aes(x=Temp, y=GSI, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm",linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(GSI)) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
GSIPLOT

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(SSI))
SSITemp <- lme(SSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = model_data)
summary(SSITemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 802.7696 824.3394 -394.3848
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0003088419 2.593232
# 
# Fixed effects:  SSI ~ Temp * Ecotype + Sex 
#                         Value Std.Error  DF   t-value p-value
# (Intercept)          -6.708033  2.086790 145 -3.214522  0.0016
# Temp                  0.031285  0.095704 145  0.326895  0.7442
# EcotypeLimnetic      -0.854209  3.585981 145 -0.238208  0.8121
# SexM                  0.090342  0.410345 145  0.220162  0.8261
# Temp:EcotypeLimnetic  0.018493  0.167869 145  0.110163  0.9124
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.987                     
# EcotypeLimnetic      -0.578  0.588              
# SexM                 -0.023 -0.073 -0.179       
# Temp:EcotypeLimnetic  0.559 -0.580 -0.994  0.183
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -4.3944015 -0.2115759  0.2476515  0.5727505  1.1157639 
# 
# Number of Observations: 166
# Number of Groups: 17

SSIPLOT <-  ggplot(End, aes(x=Temp, y=SSI, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(SSI)) +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
SSIPLOT

End <- End %>% 
  mutate(
    Fibrosis = case_when(
      Fibrosis == "Y" ~ 1,
      Fibrosis == "N" ~ 0,
      TRUE ~ NA_real_))

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(Fibrosis))
FibPresTemp <- lme(Fibrosis ~ Temp*Ecotype, random = ~ 1 | CCD, data = model_data)
summary(FibPresTemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC     BIC    logLik
# 272.4936 291.273 -130.2468
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.1164672 0.4829438
# 
# Fixed effects:  Fibrosis ~ Temp * Ecotype 
#                           Value Std.Error  DF    t-value p-value
# (Intercept)           0.2008747 0.4845974 153  0.4145187  0.6791
# Temp                  0.0082111 0.0220678 153  0.3720845  0.7103
# EcotypeLimnetic       0.3753644 0.8160757 153  0.4599627  0.6462
# Temp:EcotypeLimnetic -0.0133266 0.0376667 153 -0.3538026  0.7240
# Correlation: 
#   (Intr) Temp   EctypL
# Temp                 -0.992              
# EcotypeLimnetic      -0.566  0.564       
# Temp:EcotypeLimnetic  0.557 -0.564 -0.993
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.1178497 -0.7998728 -0.5687094  1.0732522  1.5095200 
# 
# Number of Observations: 173
# Number of Groups: 17 

FibPresPLOT <- ggplot(End, aes(x=Temp, y=Fibrosis, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 1, jitter.height = 0.01, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual( values=c('black','blue'))+
  stat_smooth(method="glm", method.args = list(family = "binomial"), se = TRUE, linewidth = 0.00001, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(Fibrosis~Presence))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
FibPresPLOT

model_data <- End %>%
  filter(!is.na(FibrosisScore), !is.na(Temp), !is.na(Ecotype), !is.na(CCD))

FibSevTemp <- lme(FibrosisScore ~ Temp*Ecotype, random = ~ 1 | CCD, data = model_data)
summary(FibSevTemp)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 529.7531 548.5679 -258.8766
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.4735568 0.9990041
# 
# Fixed effects:  FibrosisScore ~ Temp * Ecotype 
#                           Value Std.Error  DF    t-value p-value
# (Intercept)          -0.2776892 1.3806629 154 -0.2011275  0.8409
# Temp                  0.0527283 0.0629501 154  0.8376210  0.4035
# EcotypeLimnetic       0.4893031 2.1592904 154  0.2266036  0.8210
# Temp:EcotypeLimnetic -0.0290716 0.1000862 154 -0.2904656  0.7719
# Correlation: 
#   (Intr) Temp   EctypL
# Temp                 -0.991              
# EcotypeLimnetic      -0.550  0.552       
# Temp:EcotypeLimnetic  0.548 -0.560 -0.992
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.2323188 -0.6137974 -0.4280177  0.4099177  2.4679690 
# 
# Number of Observations: 174
# Number of Groups: 17 
PermTest(FibSevTemp, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lme(obj = FibSevTemp, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#   p.value
# (Intercept)   0.8840
# Temp          0.3891
# Ecotype       0.8092
# Temp:Ecotype  0.7572

FibScorePLOT <- ggplot(End, aes(x=Temp, y=as.numeric(as.character(FibrosisScore)), color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression("Temperature"*degree*C), y = expression(Fibrosis~Score))+
  geom_smooth(method = "glm", method.args = list(family = poisson(link = "log")), 
              se = TRUE, linewidth = 0.00001, alpha = .2)+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") + 
  ylim(0,NA) + 
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
FibScorePLOT

# pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/Figure5.pdf",
#     width = 6.5, 
#     height = 6.5)

BCPLOT + HSIPLOT + GSIPLOT + SSIPLOT + FibPresPLOT + FibScorePLOT + plot_layout(ncol = 2)

# dev.off()


# Lake Data as a Predictor of Thermal Tolerance ---------------------------


model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(BodyCond))
BCMean <- lme(BodyCond ~ Temp*MeanDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(BCMean)
# AIC       BIC  logLik
# -2361.641 -2339.773 1187.82

BCMax <- lme(BodyCond ~ Temp*MaxDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(BCMax)
# AIC       BIC   logLik
# -2359.057 -2337.189 1186.529

BCSA <- lme(BodyCond ~ Temp*SA + Sex, random = ~ 1 | CCD, data = model_data)
summary(BCSA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC       BIC   logLik
# -2368.627 -2346.759 1191.313
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 8.441752e-05 0.0001754398
# 
# Fixed effects:  BodyCond ~ Temp * SA + Sex 
#                     Value    Std.Error  DF   t-value p-value
# (Intercept)  0.0003290755 0.0007257222 152  0.453446  0.6509
# Temp         0.0000212379 0.0000336319 152  0.631483  0.5287
# SA           0.0000614681 0.0000656713 152  0.935996  0.3508
# SexM         0.0002186539 0.0000274436 152  7.967399  0.0000
# Temp:SA     -0.0000015492 0.0000030541 152 -0.507235  0.6127
# Correlation: 
#   (Intr) Temp   SA     SexM  
# Temp    -0.992                     
# SA      -0.959  0.955              
# SexM    -0.052  0.022  0.020       
# Temp:SA  0.950 -0.961 -0.992 -0.006
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.84383233 -0.57302455  0.04526031  0.58924805  2.49418273 
# 
# Number of Observations: 173
# Number of Groups: 17 

BCLake <-  ggplot(model_data, aes(x=as.numeric(as.character(Temp)), y=BodyCond, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26)) +
  geom_smooth(method = "lm",linewidth = 0.5, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(Body~Condition)) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BCLake

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(HSI))
HSIMean <- lme(HSI ~ Temp*MeanDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(HSIMean)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC     BIC  logLik
# 642.14 663.924 -314.07
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 7.910203e-05 1.452549
# 
# Fixed effects:  HSI ~ Temp * MeanDepth + Sex 
# Value Std.Error  DF    t-value p-value
# (Intercept)    -3.194969 1.9083083 150 -1.6742416  0.0962
# Temp            0.061676 0.0874156 150  0.7055441  0.4816
# MeanDepth       0.036596 0.2744492 150  0.1333433  0.8941
# SexM           -0.041903 0.2253515 150 -0.1859433  0.8527
# Temp:MeanDepth  0.000436 0.0127048 150  0.0343571  0.9726
# Correlation: 
#   (Intr) Temp   MnDpth SexM  
# Temp           -0.991                     
# MeanDepth       0.868 -0.870              
# SexM            0.022 -0.082  0.134       
# Temp:MeanDepth -0.854  0.870 -0.994 -0.138
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -5.4165694 -0.2352132  0.2601810  0.5370446  1.0211898 
# 
# Number of Observations: 171
# Number of Groups: 17 

HSIMax <- lme(HSI ~ Temp*MaxDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(HSIMax)
# AIC      BIC   logLik
# 645.5819 667.3659 -315.791

HSISA <- lme(HSI ~ Temp*SA + Sex, random = ~ 1 | CCD, data = model_data)
summary(HSISA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 642.9418 664.7257 -314.4709
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0001043084 1.460725
# 
# Fixed effects:  HSI ~ Temp * SA + Sex 
#                 Value Std.Error  DF    t-value p-value
# (Intercept) -5.312928  3.568833 150 -1.4887020  0.1387
# Temp         0.138457  0.163714 150  0.8457223  0.3991
# SA           0.166680  0.342692 150  0.4863843  0.6274
# SexM        -0.027664  0.225799 150 -0.1225154  0.9027
# Temp:SA     -0.007069  0.015817 150 -0.4469029  0.6556
# Correlation: 
#   (Intr) Temp   SA     SexM  
# Temp    -0.991                     
# SA      -0.964  0.960              
# SexM    -0.016 -0.027 -0.038       
# Temp:SA  0.952 -0.964 -0.992  0.051
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -5.2632010 -0.1771511  0.2643272  0.5515514  1.1813434 
# 
# Number of Observations: 171
# Number of Groups: 17

HSILake <-  ggplot(model_data, aes(x=as.numeric(as.character(Temp)), y=HSI, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26)) +
  geom_smooth(method = "lm",linewidth = 0.5, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(HSI)) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSILake

model_data <- End %>%
  filter(!is.na(Temp), !is.na(GSI))
GSIMean <- lme(GSI ~ Temp*MeanDepth, random = ~ 1 | CCD, data = model_data)
summary(GSIMean)
# AIC      BIC    logLik
# 378.0534 391.9584 -183.0267

GSIMax <- lme(GSI ~ Temp*MeanDepth, random = ~ 1 | CCD, data = model_data)
summary(GSIMax)
# AIC      BIC    logLik
# 378.0534 391.9584 -183.0267

GSISA <- lme(GSI ~ Temp*SA, random = ~ 1 | CCD, data = model_data)
summary(GSISA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 367.3801 381.2851 -177.6901
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0002181751 2.192104
# 
# Fixed effects:  GSI ~ Temp * SA 
#                   Value Std.Error DF    t-value p-value
# (Intercept) -17.282221  8.889054 60 -1.9442137  0.0566
# Temp          0.401061  0.401598 60  0.9986628  0.3220
# SA            0.847291  0.846872 60  1.0004944  0.3211
# Temp:SA      -0.024219  0.038834 60 -0.6236500  0.5352
# Correlation: 
#   (Intr) Temp   SA    
# Temp    -0.994              
# SA      -0.969  0.969       
# Temp:SA  0.954 -0.967 -0.993
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -5.35315231 -0.31213363  0.06005205  0.54658547  1.63900076 
# 
# Number of Observations: 79
# Number of Groups: 16 

GSILake <-  ggplot(model_data, aes(x=as.numeric(as.character(Temp)), y=GSI, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26)) +
  geom_smooth(method = "lm",linewidth = 0.5, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(GSI)) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
GSILake

model_data <- End %>%
  filter(!is.na(Temp), !is.na(Sex), !is.na(SSI))
SSIMean <- lme(SSI ~ Temp*MeanDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(SSIMean)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 810.9003 832.4702 -398.4502
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0002344411 2.594852
# 
# Fixed effects:  SSI ~ Temp * MeanDepth + Sex 
#                     Value Std.Error  DF    t-value p-value
# (Intercept)    -5.814187  3.501077 145 -1.6606854  0.0989
# Temp           -0.001861  0.161255 145 -0.0115384  0.9908
# MeanDepth       0.204072  0.496515 145  0.4110092  0.6817
# SexM            0.090739  0.409627 145  0.2215153  0.8250
# Temp:MeanDepth -0.006891  0.023036 145 -0.2991403  0.7653
# Correlation: 
#   (Intr) Temp   MnDpth SexM  
# Temp           -0.990                     
# MeanDepth       0.874 -0.876              
# SexM            0.060 -0.122  0.156       
# Temp:MeanDepth -0.861  0.876 -0.994 -0.163
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -4.4254886 -0.2276604  0.2326439  0.5722820  1.1315176 
# 
# Number of Observations: 166
# Number of Groups: 17 

SSIMax <- lme(SSI ~ Temp*MaxDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(SSIMax)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 814.1715 835.7413 -400.0857
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 0.0002593883 2.594962
# 
# Fixed effects:  SSI ~ Temp * MaxDepth + Sex 
#                   Value Std.Error  DF    t-value p-value
# (Intercept)   -6.129610 3.1520642 145 -1.9446337  0.0538
# Temp           0.010968 0.1452172 145  0.0755249  0.9399
# MaxDepth       0.075846 0.2191559 145  0.3460819  0.7298
# SexM           0.086638 0.4103175 145  0.2111476  0.8331
# Temp:MaxDepth -0.002366 0.0101961 145 -0.2320693  0.8168
# Correlation: 
#   (Intr) Temp   MxDpth SexM  
# Temp          -0.990                     
# MaxDepth       0.842 -0.846              
# SexM           0.053 -0.122  0.165       
# Temp:MaxDepth -0.827  0.844 -0.993 -0.172
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -4.4172940 -0.2265183  0.2321170  0.5682290  1.1102899 
# 
# Number of Observations: 166
# Number of Groups: 17 

SSILake <-  ggplot(model_data, aes(x=as.numeric(as.character(Temp)), y=SSI, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26)) +
  geom_smooth(method = "lm",linewidth = 0.5, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(SSI)) +
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSILake

End <-subset(ThermTol, Trial == "End" ) 

model_data <- End %>%
  filter(!is.na(Sex), !is.na(Temp), !is.na(Ecotype), !is.na(Fibrosis))

model_data <- model_data %>%
  mutate(Fibrosis = case_when(Fibrosis == "Y" ~ 1, 
                              Fibrosis == "N" ~ 0,
                              TRUE ~ NA_real_))

FibPresMean <- lme(Fibrosis ~ Temp*MeanDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(FibPresMean)
# AIC      BIC    logLik
# 281.9748 303.8426 -133.9874

FibPresMax <- lme(Fibrosis ~ Temp*MaxDepth + Sex, random = ~ 1 | CCD, data = model_data)
summary(FibPresMax)
# AIC      BIC    logLik
# 285.4615 307.3292 -135.7307

FibPresSA <- lme(Fibrosis ~ Temp*SA + Sex, random = ~ 1 | CCD, data = model_data)
summary(FibPresSA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 273.1676 295.0353 -129.5838
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 2.133788e-05 0.4755856
# 
# Fixed effects:  Fibrosis ~ Temp * SA + Sex 
#                 Value Std.Error  DF    t-value p-value
# (Intercept)  1.0636793 1.1582533 152  0.9183478  0.3599
# Temp        -0.0573361 0.0531162 152 -1.0794469  0.2821
# SA          -0.0732832 0.1106602 152 -0.6622361  0.5088
# SexM         0.0750982 0.0732016 152  1.0259092  0.3066
# Temp:SA      0.0059414 0.0051046 152  1.1639341  0.2463
# Correlation: 
#   (Intr) Temp   SA     SexM  
# Temp    -0.991                     
# SA      -0.965  0.961              
# SexM    -0.017 -0.026 -0.037       
# Temp:SA  0.953 -0.965 -0.992  0.051
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6174799 -0.7522687 -0.4848007  0.9874082  1.7274313 
# 
# Number of Observations: 173
# Number of Groups: 17 

FibPresLake <- ggplot(model_data, aes(x=Temp, y=Fibrosis, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 1, jitter.height = 0.01, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual( values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  stat_smooth(method="glm", method.args = list(family = "binomial"), se = TRUE, linewidth = 0.5, alpha = .2) +
  labs(x = expression("Temperature"*degree*C), y = expression(Fibrosis~Presence))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
FibPresLake

model_data <- End %>%
  filter(!is.na(FibrosisScore), !is.na(Temp), !is.na(Ecotype), !is.na(Sex))
FibSevMean <- lme(FibrosisScore ~ Temp * MeanDepth+Sex, random = ~ 1 | CCD, data = model_data)
summary(FibSevMean)
# AIC      BIC    logLik
# 536.555 558.4227 -261.2775

FibSevMax <- lme(FibrosisScore ~ Temp * MaxDepth+Sex, random = ~ 1 | CCD, data = model_data)
summary(FibSevMax)
# AIC      BIC   logLik
# 540.03 561.8978 -263.015

FibSevSA <- lme(FibrosisScore ~ Temp * SA+Sex, random = ~ 1 | CCD, data = model_data)
summary(FibSevSA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 527.4687 549.3365 -256.7344
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:  0.06610736 1.012198
# 
# Fixed effects:  FibrosisScore ~ Temp * SA + Sex 
#                 Value Std.Error  DF   t-value p-value
# (Intercept)  3.170230 2.5230182 152  1.256523  0.2109
# Temp        -0.158757 0.1156461 152 -1.372785  0.1718
# SA          -0.296696 0.2408417 152 -1.231913  0.2199
# SexM         0.130456 0.1559405 152  0.836575  0.4041
# Temp:SA      0.018095 0.0110990 152  1.630336  0.1051
# Correlation: 
#   (Intr) Temp   SA     SexM  
# Temp    -0.991                     
# SA      -0.964  0.960              
# SexM    -0.019 -0.023 -0.033       
# Temp:SA  0.953 -0.965 -0.992  0.047
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.5500476 -0.6313380 -0.4448665  0.4511303  2.5598707 
# 
# Number of Observations: 173
# Number of Groups: 17 

FibScoreLake <- ggplot(model_data, aes(x=Temp, y=as.numeric(as.character(FibrosisScore)), color=Pop)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 1), alpha= 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  labs(x = expression("Temperature"*degree*C), y = expression(Fibrosis~Score))+
  geom_smooth(method = "glm", method.args = list(family = poisson(link = "log")), 
              se = TRUE, linewidth = 0.5, alpha = .2)+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") + 
  ylim(0,NA) + 
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26))
FibScoreLake

# pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/SuppFigure4.pdf",
#     width = 6.5, 
#     height = 6.5)

BCLake + HSILake + GSILake + SSILake + FibPresLake + FibScoreLake + plot_layout(ncol = 2)

# dev.off()

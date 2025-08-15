# OCLTTH Prediction 4
# Figure 4

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

# Assumption 4 ------------------------------------------------------------



###### Does AS influence performance metrics? Body Condition, HSI, GSI, SSI, fibrosis? Use tank density as a random effect. CCD = cumulative competitor days (i.e. density)

complete_data <- End[complete.cases(End[c("BodyCond", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$BodyCond)
shapiro.test(complete_data$BodyCond)
BodyCond <- lme(BodyCond ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(BodyCond)
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -897.2912 -882.1791 455.6456
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)    Residual
# StdDev: 3.301216e-05 0.000141199
# 
# Fixed effects:  BodyCond ~ AS * Ecotype + Sex 
#                             Value    Std.Error DF   t-value p-value
# (Intercept)         0.0009379616 7.229853e-05 61 12.973453  0.0000
# AS                  0.0000001124 7.977000e-08 61  1.409029  0.1639
# EcotypeLimnetic     0.0001823549 1.077706e-04  3  1.692066  0.1892
# SexM                0.0002122160 3.553555e-05 61  5.971934  0.0000
# AS:EcotypeLimnetic -0.0000002318 1.028500e-07 61 -2.253585  0.0278
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.875                     
# EcotypeLimnetic    -0.662  0.625              
# SexM               -0.056 -0.232 -0.125       
# AS:EcotypeLimnetic  0.679 -0.773 -0.898  0.167
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.30762962 -0.62546682  0.07094662  0.60470225  2.16414487 
# 
# Number of Observations: 69
# Number of Groups: 5 


BodyCond <- ggplot(End, aes(x=AS, y=BodyCond, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCond

BodyCondSex <- ggplot(End, aes(x=AS, y=BodyCond, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondSex

complete_data <- End[complete.cases(End[c("HSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$HSI)
HSI <- lme(HSI ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(HSI)
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 271.8368 286.8388 -128.9184
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:   0.4058095 1.336893
# 
# Fixed effects:  HSI ~ AS * Ecotype + Sex 
#                         Value Std.Error DF   t-value p-value
# (Intercept)        -3.242425 0.7078876 60 -4.580423  0.0000
# AS                  0.001252 0.0007677 60  1.630512  0.1082
# EcotypeLimnetic     1.927603 1.0603472  3  1.817898  0.1667
# SexM               -0.136506 0.3427397 60 -0.398278  0.6918
# AS:EcotypeLimnetic -0.002073 0.0009899 60 -2.094388  0.0405
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.857                     
# EcotypeLimnetic    -0.660  0.613              
# SexM               -0.043 -0.245 -0.137       
# AS:EcotypeLimnetic  0.665 -0.772 -0.877  0.177
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -5.3253351 -0.1914553  0.2365350  0.5063263  1.2870953 
# 
# Number of Observations: 68
# Number of Groups: 5 

HSI <- ggplot(End, aes(x=AS, y=HSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSI

HSISex <- ggplot(End, aes(x=AS, y=HSI, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISex

complete_data <- End[complete.cases(End[c("GSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$GSI)
GSI <- lme(GSI ~ AS*Ecotype, random = ~ 1 | CCD, data = complete_data)
summary(GSI)
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC     BIC    logLik
# 192.3121 200.916 -90.15606
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:    1.463822 2.376081
# 
# Fixed effects:  GSI ~ AS * Ecotype 
#                         Value Std.Error DF   t-value p-value
# (Intercept)        -9.739645 2.0492718 28 -4.752735  0.0001
# AS                  0.003567 0.0019454 28  1.833361  0.0774
# EcotypeLimnetic     5.519751 2.7138031  3  2.033954  0.1348
# AS:EcotypeLimnetic -0.004660 0.0023113 28 -2.016092  0.0535
# Correlation: 
#   (Intr) AS     EctypL
# AS                 -0.867              
# EcotypeLimnetic    -0.755  0.654       
# AS:EcotypeLimnetic  0.729 -0.842 -0.812
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -4.1417877 -0.1727776  0.1234606  0.3718742  2.0041465 
# 
# Number of Observations: 35
# Number of Groups: 5 

GSI <- ggplot(End, aes(x=AS, y=GSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(GSI))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
GSI

complete_data <- End[complete.cases(End[c("SSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$SSI)
SSI <- lme(SSI ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(SSI)
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 335.8517 350.8537 -160.9259
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:   0.6680251 2.222535
# 
# Fixed effects:  SSI ~ AS * Ecotype + Sex 
# Value Std.Error DF   t-value p-value
# (Intercept)        -7.978655 1.1828206 60 -6.745448  0.0000
# AS                  0.001634 0.0012892 60  1.267436  0.2099
# EcotypeLimnetic     3.226141 1.7630140  3  1.829901  0.1647
# SexM                0.214204 0.5640434 60  0.379765  0.7055
# AS:EcotypeLimnetic -0.003225 0.0016559 60 -1.947473  0.0562
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.861                     
# EcotypeLimnetic    -0.662  0.614              
# SexM               -0.051 -0.227 -0.129       
# AS:EcotypeLimnetic  0.670 -0.777 -0.878  0.171
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -4.5808120 -0.3490345  0.2252840  0.5424333  1.3709496 
# 
# Number of Observations: 68
# Number of Groups: 5 


SSI <- ggplot(End, aes(x=AS, y=SSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.00001, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSI

End <- End %>% 
  mutate(
    Fibrosis = case_when(
      Fibrosis == "Y" ~ 1,
      Fibrosis == "N" ~ 0,
      TRUE ~ NA_real_))

complete_data <- End[complete.cases(End[c("Fibrosis", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_)

FibPres <- lme(Fibrosis ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(FibPres)
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 135.627 150.7392 -60.81352
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 1.713519e-05 0.4564633
# 
# Fixed effects:  Fibrosis ~ AS * Ecotype + Sex 
#                         Value Std.Error DF    t-value p-value
# (Intercept)         0.5052759 0.2173742 61  2.3244517  0.0234
# AS                 -0.0003740 0.0002454 61 -1.5241566  0.1326
# EcotypeLimnetic     0.2274843 0.3207299  3  0.7092704  0.5293
# SexM                0.1005876 0.1129805 61  0.8903096  0.3768
# AS:EcotypeLimnetic  0.0000441 0.0003161 61  0.1393865  0.8896
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.900                     
# EcotypeLimnetic    -0.663  0.641              
# SexM               -0.098 -0.201 -0.085       
# AS:EcotypeLimnetic  0.702 -0.769 -0.930  0.122
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.5656663 -0.6615424 -0.4354992  1.0561615  2.0287494 
# 
# Number of Observations: 69
# Number of Groups: 5 

FibPresA <- lme(Fibrosis ~ AS + Ecotype, random = ~ 1 | CCD, data = complete_data)
summary(FibPresA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 115.5901 126.5384 -52.79505
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 1.643394e-05 0.4522724
# 
# Fixed effects:  Fibrosis ~ AS + Ecotype 
#                       Value  Std.Error DF   t-value p-value
# (Intercept)      0.5194058 0.14800459 63  3.509390  0.0008
# AS              -0.0003242 0.00015305 63 -2.118233  0.0381
# EcotypeLimnetic  0.2610335 0.11659363  3  2.238831  0.1111
# Correlation: 
#   (Intr) AS    
# AS              -0.879       
# EcotypeLimnetic -0.017 -0.309
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.4677392 -0.7359494 -0.4457871  1.0714861  1.9459711 
# 
# Number of Observations: 69
# Number of Groups: 5 

FibPres <- ggplot(End, aes(x=AS, y=Fibrosis, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual( values=c('black','blue'))+
  stat_smooth(method="glm", method.args = list(family = "binomial"), se = TRUE, linewidth = 0.00001, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Fibrosis~Presence))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
FibPres


model_data <- End %>%
  filter(!is.na(FibrosisScore), !is.na(AS), !is.na(Ecotype), !is.na(CCD))

FibSev <- lme(FibrosisScore ~ AS*Ecotype, random = ~ 1 | CCD, data = model_data)
summary(FibSev)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 215.8521 228.8985 -101.9261
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 2.69025e-05 0.8693328
# 
# Fixed effects:  FibrosisScore ~ AS * Ecotype 
#                         Value Std.Error DF    t-value p-value
# (Intercept)         0.8148639 0.4119806 62  1.9779181  0.0524
# AS                 -0.0004997 0.0004577 62 -1.0916606  0.2792
# EcotypeLimnetic     0.4760682 0.6086261  3  0.7822014  0.4912
# AS:EcotypeLimnetic -0.0000316 0.0005975 62 -0.0529401  0.9579
# Correlation: 
#   (Intr) AS     EctypL
# AS                 -0.944              
# EcotypeLimnetic    -0.677  0.639       
# AS:EcotypeLimnetic  0.723 -0.766 -0.930
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.2651211 -0.6153933 -0.3765244  0.3641491  2.9708135 
# 
# Number of Observations: 69
# Number of Groups: 5 
PermTest(FibSev, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lme(obj = FibSev, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#   p.value
# (Intercept)  0.2141
# AS           0.2857
# Ecotype      0.4288
# AS:Ecotype   0.9560

FibSevA <- lme(FibrosisScore ~ AS + Ecotype, random = ~ 1 | CCD, data = model_data)
summary(FibSevA)
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 200.8396 211.7878 -95.41978
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 2.551703e-05 0.8627404
# 
# Fixed effects:  FibrosisScore ~ AS + Ecotype 
#                       Value  Std.Error DF   t-value p-value
# (Intercept)      1.8306393 0.28232885 63  6.484067  0.0000
# AS              -0.0005182 0.00029195 63 -1.775046  0.0807
# EcotypeLimnetic  0.4461114 0.22241031  3  2.005803  0.1385
# Correlation: 
#   (Intr) AS    
# AS              -0.879       
# EcotypeLimnetic -0.017 -0.309
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.2637985 -0.6252069 -0.3766920  0.3660521  2.9923458 
# 
# Number of Observations: 69
# Number of Groups: 5
PermTest(FibSevA, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lme(obj = FibSevA, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#   p.value
# (Intercept)  0.0711
# AS           0.0852
# Ecotype      0.0369

End$FibrosisScore <- factor(End$FibrosisScore,
                            levels = c("0", "1", "2", "3", "4"),
                            ordered = TRUE)

FibScore <- ggplot(End, aes(x=AS, y=as.numeric(as.character(FibrosisScore)), color=Ecotype)) +
  geom_point(size = .5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Fibrosis~Score))+
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.00001, alpha = .2)+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
FibScore

# pdf(file = "R:/ThermTol_MultStress/RespManuscriptFigs/Figure4.pdf",
#     width = 6, 
#     height = 6.5)

BodyCond + HSI + GSI + SSI + FibPres + FibScore + plot_layout(ncol = 2)

# dev.off()

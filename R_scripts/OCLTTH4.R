# OCLTTH Prediction 4
# Figure 4

# Load Packages -----------------------------------------------------------
{
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
library(nlme)
}

# Load data ---------------------------------------------------------------

{
ThermTol <- read_xlsx("/home/ekerns/ThermTol/RawData/12125_ThermTol_OutliersRemoved.xlsx")
head(ThermTol)

ThermTol_CCD <- read_xlsx("/home/ekerns/ThermTol/RawData/ThermalToleranceProject_SuspectDataRemoved.xlsx")
head(ThermTol_CCD)

ThermTol$Ecotype <- as.factor(ThermTol$Ecotype)
ThermTol$Fibrosis <- as.factor(ThermTol$Fibrosis)
ThermTol$Trial <- as.factor(ThermTol$Trial)

ThermTol$CCD <- ThermTol_CCD$CCD[match(ThermTol$FishID, ThermTol_CCD$FishID)]

EcotypeMortality <- read_xlsx("/home/ekerns/ThermTol/RawData/EcotypeMortality.xlsx")
head(EcotypeMortality)

LakeData <- read_xlsx("/home/ekerns/ThermTol/RawData/depth_profiles.xlsx")
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
colnames(LakeData)[1] <- "Pop"
head(LakeData)

LakeData <- LakeData %>% 
  group_by(Pop) %>% 
  mutate(MinDepth = max(DepthM, na.rm = TRUE), 
            MaxDepth = min(DepthM, na.rm = TRUE),
            MedDepth = median(DepthM, na.rm = TRUE),
            MeanDepth = mean(DepthM, na.rm = TRUE),
            SA = first(SA))
LakeData

LakeSummary <- LakeData %>%
  ungroup() %>%
  dplyr::select(Pop, DepthM, SA) %>%
  dplyr::group_by(Pop) %>%
  dplyr::summarise(
    SA        = first(SA),
    MinDepth  = max(DepthM, na.rm = TRUE),
    MaxDepth  = min(DepthM, na.rm = TRUE),
    MedDepth  = median(DepthM, na.rm = TRUE),
    MeanDepth = mean(DepthM, na.rm = TRUE),
    .groups   = "drop"
  )


ThermTol <- left_join(ThermTol, LakeSummary, by = "Pop")
head(ThermTol)
}

# Relabel data ------------------------------------------------------------


{
ThermTol$Trial <- factor(ThermTol$Trial, levels = c("S", "F"), 
                         labels = c("Start", "End"))

ThermTolTAF <- ThermTol 
ThermTolTAF$Temp <- as.factor(ThermTolTAF$Temp)

EcotypeMortalityTAF <- EcotypeMortality 
EcotypeMortalityTAF$Temp <- as.factor(EcotypeMortalityTAF$Temp)
}

# Clean up data -------------------------------------------------
{
# Fulton's Body condition = (Mass / SL^3)100
ThermTol$BodyCond <- ((ThermTol$Mass/(ThermTol$SL)^3)*100)
hist(ThermTol$BodyCond)


## Subset data
Start <-subset(ThermTol, Trial== "Start" )
End <-subset(ThermTol, Trial == "End" )
StartTAF <-subset(ThermTolTAF, Trial== "Start" )
EndTAF <-subset(ThermTolTAF, Trial == "End" )
Limno <-subset(ThermTol, Ecotype== "Limnetic" )
Benno <-subset(ThermTol, Ecotype== "Benthic" )
}


# Test Assumption 4 of the OCLTTH ---------------------------------------------------------


## Start of experiment -------------------------------------------------------


###### Does AS influence performance metrics? Body Condition?
{
  complete_data <- Start[complete.cases(Start[c("BodyCond", "AS", "Ecotype")]), ]
  nrow(complete_data)
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- glm(BodyCond ~ AS*Ecotype, data = complete_data)
  summary(BodyCond)
}
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.119e-03  7.364e-05  15.195   <2e-16 ***
# AS                 6.594e-08  6.241e-08   1.057   0.2922    
# EcotypeBenthic     2.520e-04  1.116e-04   2.259   0.0251 *  
# AS:EcotypeBenthic -1.985e-07  1.005e-07  -1.975   0.0498 *  
# 
# Null deviance: 7.7167e-06  on 178  degrees of freedom
# Residual deviance: 7.4709e-06  on 175  degrees of freedom
# AIC: -2523.6
emtrends(BodyCond, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE  df  lower.CL upper.CL
# Limnetic  6.59e-08 6.24e-08 175 -5.72e-08 1.89e-07
# Benthic  -1.33e-07 7.88e-08 175 -2.88e-07 2.29e-08
plot(BodyCond)

BodyCondP <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black'))+
  labs(x = expression(Post~Acclimation~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Post~Acclimation~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")+
  annotate("text", x = 1800, y = .0017, label = "Ecotype: p<0.05\nAS*Ecotype: p<0.05", hjust = 1, vjust = -1, size = 2, fontface = "bold")
BodyCondP


{
  complete_data <- Start[complete.cases(Start[c("BodyCond", "Temp", "Ecotype")]), ]
  nrow(complete_data)
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- glm(BodyCond ~ Temp*Ecotype, data = complete_data)
  summary(BodyCond)
}
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          1.438e-03  1.770e-04   8.123 7.87e-14 ***
# Temp                -1.121e-05  8.063e-06  -1.390   0.1663    
# EcotypeBenthic      -3.884e-04  2.475e-04  -1.569   0.1184    
# Temp:EcotypeBenthic  1.973e-05  1.127e-05   1.751   0.0817 .  
# 
# Null deviance: 7.7167e-06  on 178  degrees of freedom
# Residual deviance: 7.5063e-06  on 175  degrees of freedom
# AIC: -2522.7
emtrends(BodyCond, "Ecotype", var = "Temp")
# Ecotype  Temp.trend       SE  df  lower.CL upper.CL
# Limnetic  -1.12e-05 8.06e-06 175 -2.71e-05 4.71e-06
# Benthic    8.52e-06 7.87e-06 175 -7.01e-06 2.40e-05
plot(BodyCond)

BodyCond_Temps <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black'))+
  labs(x = "Temperature (\u00B0C)", 
       y = "Post Acclimation Body Condition")+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
BodyCond_Temps

## End of experiment -------------------------------------------------------



###### Does AS influence performance metrics? Body Condition, HSI, GSI, SSI, fibrosis? Use tank density as a random effect. CCD = cumulative competitor days (i.e. density)
{
complete_data <- End[complete.cases(End[c("BodyCond", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$BodyCond)
shapiro.test(complete_data$BodyCond)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
BodyCond <- lme(BodyCond ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(BodyCond)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC   logLik
# -881.901 -866.899 447.9505
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 3.327158e-05 0.0001422991
# 
# Fixed effects:  BodyCond ~ AS * Ecotype + Sex 
#                           Value    Std.Error DF   t-value p-value
# (Intercept)        0.0011259256 9.411005e-05 60 11.963924  0.0000
# AS                -0.0000001250 8.114000e-08 60 -1.540872  0.1286
# EcotypeBenthic    -0.0001879070 1.182277e-04  3 -1.589365  0.2102
# SexM               0.0002117219 3.605266e-05 60  5.872574  0.0000
# AS:EcotypeBenthic  0.0000002377 1.149300e-07 60  2.068130  0.0429
# Correlation: 
#   (Intr) AS     EctypB SexM  
# AS                -0.902                     
# EcotypeBenthic    -0.788  0.716              
# SexM              -0.244  0.051  0.159       
# AS:EcotypeBenthic  0.677 -0.714 -0.915 -0.200
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.28791591 -0.64447760  0.02848734  0.60486784  2.14653689 
# 
# Number of Observations: 68
# Number of Groups: 5 
emtrends(BodyCond, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Limnetic -1.25e-07 8.11e-08 60 -2.87e-07 3.73e-08
# Benthic   1.13e-07 8.04e-08 60 -4.82e-08 2.74e-07
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(BodyCond)

BodyCond <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")+
  annotate("text", x = 1800, y = .0016, label = "AS*Ecotype: p<0.05", hjust = 1, vjust = -1, size = 2, fontface = "bold")
BodyCond

BodyCondSex <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondSex


{
complete_data <- End[complete.cases(End[c("HSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$HSI)
HSI <- lme(HSI ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(HSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -304.6516 -289.7617 159.3258
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)   Residual
# StdDev:  0.00729815 0.01294831
# 
# Fixed effects:  HSI ~ AS * Ecotype + Sex 
#                           Value   Std.Error DF   t-value p-value
# (Intercept)         0.04074918 0.007849398 59  5.191377  0.0000
# AS                  0.00002152 0.000007675 59  2.803981  0.0068
# EcotypeLimnetic     0.03771250 0.012746138  3  2.958740  0.0596
# SexM               -0.01107452 0.003385640 59 -3.271025  0.0018
# AS:EcotypeLimnetic -0.00003312 0.000011002 59 -3.010148  0.0038
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.769                     
# EcotypeLimnetic    -0.613  0.523              
# SexM               -0.017 -0.267 -0.176       
# AS:EcotypeLimnetic  0.536 -0.711 -0.810  0.235
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.90838381 -0.49877342 -0.09451304  0.54376505  2.28686737 
# 
# Number of Observations: 67
# Number of Groups: 5 
emtrends(HSI, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Benthic   2.15e-05 7.68e-06 59  6.16e-06 3.69e-05
# Limnetic -1.16e-05 7.74e-06 59 -2.71e-05 3.89e-06
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(HSI)

HSI <- ggplot(complete_data, aes(x=AS, y=HSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")+
  annotate("text", x = 1800, y = .12, label = "AS: p<0.01\nAS*Ecotype: p<0.01", hjust = 1, vjust = 1, size = 2, fontface = "bold")
HSI

HSISex <- ggplot(complete_data, aes(x=AS, y=HSI, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISex


{
complete_data <- End[complete.cases(End[c("GSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$GSI)
GSI <- lme(GSI ~ AS*Ecotype, random = ~ 1 | CCD, data = complete_data)
summary(GSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -258.0829 -249.6757 135.0415
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)    Residual
# StdDev: 8.352405e-08 0.001499047
# 
# Fixed effects:  GSI ~ AS * Ecotype 
#                             Value    Std.Error DF    t-value p-value
# (Intercept)         0.0020993899 0.0009521355 27  2.2049278  0.0362
# AS                 -0.0000001459 0.0000010339 27 -0.1410665  0.8889
# EcotypeLimnetic    -0.0009548305 0.0014016259  3 -0.6812306  0.5446
# AS:EcotypeLimnetic  0.0000017514 0.0000013985 27  1.2523923  0.2212
# Correlation: 
#   (Intr) AS     EctypL
# AS                 -0.932              
# EcotypeLimnetic    -0.679  0.633       
# AS:EcotypeLimnetic  0.689 -0.739 -0.926
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.5571296 -0.6278084 -0.4046428  0.3139224  2.3490563 
# 
# Number of Observations: 34
# Number of Groups: 5 
emtrends(GSI, "Ecotype", var = "AS")
#  Ecotype   AS.trend       SE df  lower.CL upper.CL
# Benthic  -1.46e-07 1.03e-06 27 -2.27e-06 1.98e-06
# Limnetic  1.61e-06 9.42e-07 27 -3.27e-07 3.54e-06
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(GSI)

GSI <- ggplot(complete_data, aes(x=AS, y=GSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(GSI))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") 
GSI

{
complete_data <- End[complete.cases(End[c("SSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$SSI)
SSI <- lme(SSI ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(SSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC  logLik
# -718.5819 -703.692 366.291
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 2.771884e-08 0.0004769213
# 
# Fixed effects:  SSI ~ AS * Ecotype + Sex 
#                             Value    Std.Error DF   t-value p-value
# (Intercept)         0.0008424243 0.0002296699 59  3.667978  0.0005
# AS                  0.0000001752 0.0000002610 59  0.671228  0.5047
# EcotypeLimnetic     0.0004615705 0.0003678703  3  1.254710  0.2984
# SexM               -0.0000331179 0.0001189926 59 -0.278319  0.7817
# AS:EcotypeLimnetic  0.0000001237 0.0000003695 59  0.334914  0.7389
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.902                     
# EcotypeLimnetic    -0.605  0.597              
# SexM               -0.108 -0.185 -0.114       
# AS:EcotypeLimnetic  0.636 -0.709 -0.941  0.144
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.4439895 -0.8142956 -0.1826019  0.4147996  2.8065060 
# 
# Number of Observations: 67
# Number of Groups: 5 
emtrends(SSI, "Ecotype", var = "AS")
# Ecotype  AS.trend       SE df  lower.CL upper.CL
# Benthic  1.75e-07 2.61e-07 59 -3.47e-07 6.97e-07
# Limnetic 2.99e-07 2.61e-07 59 -2.23e-07 8.21e-07
plot(SSI)

SSI <- ggplot(complete_data, aes(x=AS, y=SSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSI

{
End <- End %>% 
  mutate(
    Fibrosis = case_when(
      Fibrosis == "Y" ~ 1,
      Fibrosis == "N" ~ 0,
      TRUE ~ NA_real_))

complete_data <- End[complete.cases(End[c("Fibrosis", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data) #68

FibPresA <- lme(AS ~ Fibrosis + Ecotype, random = ~ 1 | CCD, data = complete_data)
summary(FibPresA)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 946.8562 957.7281 -468.4281
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:    159.3972 291.8172
# 
# Fixed effects:  AS ~ Fibrosis + Ecotype 
#                     Value Std.Error DF   t-value p-value
# (Intercept)      892.0897 104.62278 62  8.526725  0.0000
# Fibrosis        -164.8308  77.66278 62 -2.122392  0.0378
# EcotypeLimnetic  213.2639 163.32801  3  1.305740  0.2827
# Correlation: 
#   (Intr) Fibrss
# Fibrosis        -0.183       
# EcotypeLimnetic -0.602 -0.094
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.1114594 -0.7062971  0.1603954  0.6499488  2.0581975 
# 
# Number of Observations: 68
# Number of Groups: 5 
emmeans(FibPresA, ~ Fibrosis, at = list(Fibrosis = c(0,1)))
# Fibrosis emmean   SE df lower.CL upper.CL
# 0    999 85.6  3      726     1271
# 1    834 95.9  3      529     1139
# 
# Results are averaged over the levels of: Ecotype 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(FibPresA)

anno <- data.frame(x1 = c(1), x2 = c(2), 
                   y1 = c(2400), y2 = c(2500), 
                   xstar = c(1.5), ystar = c(2560),
                   lab = c("p<0.05"),
                   Ecotype = c("Benthic"))
anno


FibPres <- ggplot(complete_data, aes(x = factor(Fibrosis), y = AS, color = Ecotype)) +
  geom_boxplot(size = 0.5, width = 0.5, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'blue')) +
  scale_x_discrete(labels = c("Absent", "Present")) +
  labs(x = expression(Fibrosis), y = expression(Aerobic~Scope~(mgO[2]/kg/hr))) +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = lab), size = 2, fontface = "bold") +
  geom_segment(data = anno, aes(x = x1, xend = x1, y = y1, yend = y2), colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, y = y1, yend = y2), colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), colour = "black") +
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
FibPres

{
model_data <- End %>%
  filter(!is.na(FibrosisScore), !is.na(AS), !is.na(Ecotype), !is.na(CCD))

FibSevA <- lme(AS ~ FibrosisScore + Ecotype, random = ~ 1 | CCD, data = model_data)
summary(FibSevA)
}
# Linear mixed-effects model fit by REML
# Data: model_data 
# AIC      BIC    logLik
# 950.5021 961.3741 -470.2511
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:    149.3007 298.1272
# 
# Fixed effects:  AS ~ FibrosisScore + Ecotype 
#                     Value Std.Error DF   t-value p-value
# (Intercept)     874.5421  99.49253 62  8.790027  0.0000
# FibrosisScore   -59.3478  41.95387 62 -1.414595  0.1622
# EcotypeLimnetic 202.0076 155.86288  3  1.296060  0.2857
# Correlation: 
#   (Intr) FbrssS
# FibrosisScore   -0.164       
# EcotypeLimnetic -0.606 -0.095
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -2.0147924 -0.6506353  0.1397220  0.7102469  2.0908604 
# 
# Number of Observations: 68
# Number of Groups: 5 
PermTest(FibSevA, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lme(obj = FibSevA, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#               p.value
# (Intercept)    0.9926
# FibrosisScore  0.1665
# Ecotype        0.1671
emmeans(FibSevA, ~ FibrosisScore, at = list(FibrosisScore = c(0,1, 2, 3)))
# FibrosisScore emmean    SE df lower.CL upper.CL
# 0    976  81.1  3      717     1234
# 1    916  79.7  3      663     1170
# 2    857  98.2  3      544     1169
# 3    798 128.0  3      389     1206
# 
# Results are averaged over the levels of: Ecotype 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(FibSevA)


FibScore <- ggplot(model_data, aes(x=as.factor(FibrosisScore), y=AS, color=Ecotype)) +
  geom_boxplot(size = 0.5, width = 0.5)+
  geom_jitter(aes(color = Ecotype), position = position_jitterdodge(jitter.width = 0.15, dodge.width = .75), size = 1, alpha = 0.2)+
  theme_classic() +
  scale_color_manual(values=c('black','blue'))+
  labs(x = expression(Fibrosis~Score), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
  annotate("text", x = -Inf, y = Inf, label = "F", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
FibScore

#pdf(file = "/home/ekerns/ThermTol/Figures/Figure4.pdf",
#    width = 8,
#    height = 8)

BodyCond + HSI + GSI + SSI + FibPres + FibScore + plot_layout(ncol = 2)

#dev.off()


# Scaled Mass Index -------------------------------------------------------



# based on Peig & Green, New perspectives for estimating body condition from mass/length data:
#     the scaled mass index as an alternative method"
#     Oikos 118: 1883-1891, 2009
# script from https://apansharing.blogspot.com/2018/05/an-r-function-olsrobust-caled-mass-index.html




## Format data -------------------------------------------------------------

{
data <- read_xlsx("/home/ekerns/ThermTol/RawData/12125_ThermTol_OutliersRemoved.xlsx")
Ben_18 <- filter(data, Ecotype == "Benthic" & Temp == "18")
Lim_18 <- filter(data, Ecotype == "Limnetic" & Temp == "18")
Ben_all <- filter(data, Ecotype == "Benthic")
Lim_all <- filter(data, Ecotype == "Limnetic")
}

# create function
scaledMassIndex <-
  function(x, y, x.0 = mean(x)) {
    require(smatr)
    require(magrittr)
    require(MASS)
    require(data.table)
    logM.ols <- lm(log(y) ~ log(x))
    logM.rob <- rlm(log(y) ~ log(x), method = "M")
    b.msa.ols <- coef(sma(log(y) ~ log(x)))[2]
    b.msa.rob <- coef(sma(log(y) ~ log(x), robust = T))[2]
    SMI.ols <- y * (x.0 / x) ^ b.msa.ols
    SMI.rob <- y * (x.0 / x) ^ b.msa.rob
    res <- data.frame(SMI.ols, SMI.rob, x, y)
    pred.DT <-
      data.table(x = seq(min(x), max(x), length = 100)) %>%
      .[, y.ols := predict(logM.ols, newdata = .) %>% exp] %>%
      .[, y.rob := predict(logM.rob, newdata = .) %>% exp]
    attr(res, "b.msa") <- c(ols = b.msa.ols, rob = b.msa.rob)
    attr(res, "pred")  <- pred.DT 
    return(res)
  }




## calculate scaled mass index for Benthic Control fish --------------------


{
mean(Ben_18$SL) # 61.65132
median(Ben_18$SL) # 61.185

Ben.dt.SMI <-
  scaledMassIndex(Ben_18$SL, Ben_18$Mass)
summary(Ben.dt.SMI)
#       SMI.ols         SMI.rob            x               y        
# Min.   :1.584   Min.   :1.501   Min.   :51.60   Min.   :1.432  
# 1st Qu.:2.429   1st Qu.:2.477   1st Qu.:59.13   1st Qu.:2.361  
# Median :2.732   Median :2.757   Median :61.19   Median :2.811  
# Mean   :2.772   Mean   :2.781   Mean   :61.65   Mean   :2.780  
# 3rd Qu.:3.114   3rd Qu.:3.177   3rd Qu.:63.71   3rd Qu.:3.270  
# Max.   :3.661   Max.   :3.741   Max.   :74.16   Max.   :4.343 


g1 <-
  ggplot(Ben_18, aes(SL, Mass)) +
  geom_point() +
  geom_line(aes(x, y, color = Method),
            data = attr(Ben.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g1

g2 <-
  ggplot(Ben_18, aes(log10(SL), log10(Mass))) +
  geom_point() +
  geom_line(
    aes(log10(x), log10(y), color = Method),
    data = attr(Ben.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
  ) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g2

g3 <-
  ggplot(reshape2::melt(Ben.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
         aes(x, value)) +
  geom_point(aes(color = Method)) +
  xlab("Body length") +
  ylab("SMI (at body length = 61)") +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g3

SMI_Ben <- ggarrange(g1,
          g2,
          g3,
          ncol = 1,
          nrow = 3,
          align = "hv")
SMI_Ben
}


## calculate scaled mass index for Limnetic Control fish --------------------
{
mean(Lim_18$SL) #61.49156
median(Lim_18$SL) # 61.575

Lim.dt.SMI <-
  scaledMassIndex(Lim_18$SL, Lim_18$Mass, x.0 = 61)
summary(Lim.dt.SMI)
#         SMI.ols         SMI.rob            x               y        
# Min.   :1.722   Min.   :1.719   Min.   :54.71   Min.   :1.703  
# 1st Qu.:2.506   1st Qu.:2.502   1st Qu.:57.84   1st Qu.:2.338  
# Median :2.765   Median :2.768   Median :61.58   Median :2.874  
# Mean   :2.736   Mean   :2.736   Mean   :61.49   Mean   :2.820  
# 3rd Qu.:3.110   3rd Qu.:3.108   3rd Qu.:64.35   3rd Qu.:3.072  
# Max.   :3.617   Max.   :3.623   Max.   :72.76   Max.   :4.084 


g1 <-
  ggplot(Lim_18, aes(SL, Mass)) +
  geom_point() +
  geom_line(aes(x, y, color = Method),
            data = attr(Lim.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g1

g2 <-
  ggplot(Lim_18, aes(log10(SL), log10(Mass))) +
  geom_point() +
  geom_line(
    aes(log10(x), log10(y), color = Method),
    data = attr(Lim.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
  ) +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g2

g3 <-
  ggplot(reshape2::melt(Lim.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
         aes(x, value)) +
  geom_point(aes(color = Method)) +
  xlab("Body length") +
  ylab("SMI (at body length = 61)") +
  scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
g3

SMI_Lim <- ggarrange(g1,
          g2,
          g3,
          ncol = 1,
          nrow = 3,
          align = "hv")
SMI_Lim
}

SMI_Ben + SMI_Lim

Ben.dt.SMI


## calculate scaled mass index for all benthic fish --------------------
{
  mean(Ben_all$SL) # 61.8378
  median(Ben_all$SL) # 61.77
  
  Ben.all.dt.SMI <-
    scaledMassIndex(Ben_all$SL, Ben_all$Mass, x.0 = 61)
  summary(Ben.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.560   Min.   :1.563   Min.   :51.09   Min.   :1.432  
  # 1st Qu.:2.377   1st Qu.:2.379   1st Qu.:59.38   1st Qu.:2.523  
  # Median :2.698   Median :2.699   Median :61.77   Median :2.774  
  # Mean   :2.737   Mean   :2.738   Mean   :61.84   Mean   :2.861  
  # 3rd Qu.:3.063   3rd Qu.:3.064   3rd Qu.:63.89   3rd Qu.:3.231  
  # Max.   :4.306   Max.   :4.299   Max.   :74.16   Max.   :4.765 
  

  g1 <-
    ggplot(Ben_all, aes(x = SL, y = Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(Ben.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))+
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  
  Ben_less <- filter(Ben_all, !FishID == "FG20S_3")
  
  g2 <-ggplot(Ben_less, aes(x = log10(SL), y = log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(Ben.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(Ben.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI at mean body length of benthic fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_Ben <- ggarrange(g1,
                       g2,
                       g3,
                       ncol = 1,
                       nrow = 3,
                       align = "hv")
  SMI_Ben
}


## calculate scaled mass index for all limetic fish --------------------
{
  mean(Lim_all$SL) # 61.43267
  median(Lim_all$SL) # 61.72
  
  Lim.all.dt.SMI <-
    scaledMassIndex(Lim_all$SL, Lim_all$Mass)
  summary(Lim.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.733   Min.   :1.732   Min.   :49.49   Min.   :1.238  
  # 1st Qu.:2.390   1st Qu.:2.389   1st Qu.:58.78   1st Qu.:2.336  
  # Median :2.696   Median :2.696   Median :61.72   Median :2.735  
  # Mean   :2.730   Mean   :2.730   Mean   :61.43   Mean   :2.741  
  # 3rd Qu.:3.005   3rd Qu.:3.005   3rd Qu.:63.64   3rd Qu.:3.116  
  # Max.   :3.928   Max.   :3.929   Max.   :72.76   Max.   :4.214 
  
  
  g1 <-
    ggplot(Lim_all, aes(SL, Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(Lim.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)")) +
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  
  g2 <-
    ggplot(Lim_all, aes(log10(SL), log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(Lim.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(Lim.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI (at mean body length of limnetic fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_Lim <- ggarrange(g1,
                       g2,
                       g3,
                       ncol = 1,
                       nrow = 3,
                       align = "hv")
  SMI_Lim
}

SMI_Ben + SMI_Lim


## The effect of AS & temp on SMI ---------------------------------


{
colnames(Ben.all.dt.SMI)[3] <- "SL"
colnames(Ben.all.dt.SMI)[4] <- "Mass"

colnames(Lim.all.dt.SMI)[3] <- "SL"
colnames(Lim.all.dt.SMI)[4] <- "Mass"

Ben_all$SMI.ols <- Ben.all.dt.SMI$SMI.ols
Ben_all$SMI.rob <- Ben.all.dt.SMI$SMI.rob
Lim_all$SMI.ols <- Lim.all.dt.SMI$SMI.ols
Lim_all$SMI.rob <- Lim.all.dt.SMI$SMI.rob

data <- bind_rows(Ben_all, Lim_all)

data <- left_join(data, ThermTol[, c("CCD", "FishID")], by = "FishID")

table(data$Ecotype)
is.na(data$SMI.rob)

SMI_end <- filter(data, Trial == "F")
SMI_start <- filter(data, Trial == "S")

data$Temp <- as.factor(data$Temp)

SMI_AS <- ggplot(data, aes(x=AS, y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'blue')) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = "SMI") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_AS
}

## Experiment end --------------------------------------------------------

{
complete_data <- SMI_end[complete.cases(SMI_end[c("SMI.rob", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)

hist(complete_data$SMI.rob)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
SMI.rob <- lme(SMI.rob ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(SMI.rob)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 90.68031 105.6823 -38.34015
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:  0.09153585 0.318888
# 
# Fixed effects:  SMI.rob ~ AS * Ecotype + Sex 
#                       Value  Std.Error DF   t-value p-value
# (Intercept)        2.5516259 0.21631199 60 11.796045  0.0000
# AS                -0.0002582 0.00018382 60 -1.404805  0.1652
# EcotypeBenthic    -0.4607935 0.27208001  3 -1.693596  0.1889
# SexM               0.5324593 0.08117568 60  6.559346  0.0000
# AS:EcotypeBenthic  0.0005091 0.00026086 60  1.951724  0.0556
# Correlation: 
#   (Intr) AS     EctypB SexM  
# AS                -0.890                     
# EcotypeBenthic    -0.788  0.706              
# SexM              -0.246  0.060  0.168       
# AS:EcotypeBenthic  0.669 -0.715 -0.902 -0.211
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.30995669 -0.61755591 -0.06000501  0.61872598  2.18767892 
# 
# Number of Observations: 68
# Number of Groups: 5 

slopes <- emtrends(SMI.rob, "Ecotype", var = "AS")
slopes
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Limnetic -0.000258 0.000184 60 -0.000626 0.000109
# Benthic   0.000251 0.000182 60 -0.000114 0.000616
plot(SMI.rob)

SMI_AS <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = "SMI") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  annotate("text", x = 1800, y = 3.5, label = "AS*Ecotype: p=0.0556", hjust = 1, vjust = -1, size = 2, fontface = "bold")
SMI_AS

SMI_Sex <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  scale_color_manual(values = c('black', 'blue')) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = "SMI") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_Sex

{
SMI.rob <- lme(SMI.rob ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(SMI.rob)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 72.04941 87.05136 -29.02471
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:  0.07443952 0.3245581
# 
# Fixed effects:  SMI.rob ~ Temp * Ecotype + Sex 
#                           Value Std.Error DF   t-value p-value
# (Intercept)          3.0350213 0.8322531 62  3.646753  0.0005
# Temp                -0.0381461 0.0411870  1 -0.926168  0.5244
# EcotypeBenthic      -0.2358912 0.9420745  1 -0.250396  0.8438
# SexM                 0.5541494 0.0805218 62  6.881982  0.0000
# Temp:EcotypeBenthic  0.0148746 0.0460541  1  0.322980  0.8011
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.994                     
# EcotypeBenthic      -0.872  0.870              
# SexM                -0.175  0.126  0.091       
# Temp:EcotypeBenthic  0.886 -0.892 -0.993 -0.098
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.52284098 -0.59754621  0.02250419  0.59780116  1.96512310 
# 
# Number of Observations: 68
# Number of Groups: 5 
emtrends(SMI.rob, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE df lower.CL upper.CL
# Limnetic    -0.0381 0.0412  1   -0.561    0.485
# Benthic     -0.0233 0.0208  1   -0.287    0.241
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95
plot(SMI.rob)

complete_data$Temp <- as.factor(complete_data$Temp)

SMI_temp <- ggplot(complete_data, aes(x=Temp, y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'blue')) +
  labs(x = "Temp", y = "SMI") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_temp

SMI_temp_sex <- ggplot(complete_data, aes(x=Temp, y=SMI.rob, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('black', 'blue')) +
  labs(x = "Temp", y = "SMI") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_temp_sex


## Experiment start --------------------------------------------------------
{
complete_data <- SMI_start[complete.cases(SMI_start[c("SMI.rob", "AS", "Ecotype")]), ]
nrow(complete_data)

hist(complete_data$SMI.rob)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
SMI.rob <- glm(SMI.rob ~ AS*Ecotype, data = complete_data)
summary(SMI.rob)
}
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.5390374  0.1786698  14.211   <2e-16 ***
# AS                 0.0002127  0.0001514   1.405   0.1618    
# EcotypeBenthic     0.5682670  0.2706757   2.099   0.0372 *  
# AS:EcotypeBenthic -0.0005115  0.0002439  -2.098   0.0374 *  
# 
# Null deviance: 45.105  on 178  degrees of freedom
# Residual deviance: 43.974  on 175  degrees of freedom
# AIC: 266.7
emtrends(SMI.rob, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE  df  lower.CL upper.CL
# Limnetic  0.000213 0.000151 175 -8.61e-05 5.12e-04
# Benthic  -0.000299 0.000191 175 -6.76e-04 7.85e-05
plot(SMI.rob)

SMI_AS <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression(Post~Acclimation~Aerobic~Scope~(mgO[2]/kg/hr)), y = "Post Acclimation SMI") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  annotate("text", x = 1800, y = 4, label = "Ecotype: p<0.05\nAS*Ecotype: p<0.05", hjust = 1, vjust = -1, size = 2, fontface = "bold")
SMI_AS

BodyCondP + SMI_AS

{
complete_data <- SMI_start[complete.cases(SMI_start[c("SMI.rob", "Temp", "Ecotype")]), ]
nrow(complete_data)

hist(complete_data$SMI.rob)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
SMI.rob <- glm(SMI.rob ~ Temp*Ecotype, data = complete_data)
summary(SMI.rob)
}
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          3.45045    0.42923   8.039  1.3e-13 ***
# Temp                -0.03081    0.01956  -1.575   0.1170    
# EcotypeBenthic      -1.10387    0.60035  -1.839   0.0677 .  
# Temp:EcotypeBenthic  0.05164    0.02732   1.890   0.0604 .  
# 
# Null deviance: 45.105  on 178  degrees of freedom
# Residual deviance: 44.157  on 175  degrees of freedom
# AIC: 267.45
emtrends(SMI.rob, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE  df lower.CL upper.CL
# Limnetic    -0.0308 0.0196 175  -0.0694  0.00779
# Benthic      0.0208 0.0191 175  -0.0168  0.05849
plot(SMI.rob)

SMI_Temp <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression("Temperature (\u00B0C)"), y = "Post Acclimation SMI") +
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_Temp

BodyCondP + SMI_AS + BodyCond_Temps + SMI_Temp


# Effect of Temperature Directly on Condition -----------------------------

###### Does Temp influence performance metrics? Body Condition, HSI, GSI, SSI, fibrosis? Use tank density Temp a random effect. CCD = cumulative competitor days (i.e. density)
{
  complete_data <- End[complete.cases(End[c("BodyCond", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- lme(BodyCond ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(BodyCond)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -899.9137 -884.9117 456.9568
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 1.210147e-07 0.0001455592
# 
# Fixed effects:  BodyCond ~ Temp * Ecotype + Sex 
#                             Value    Std.Error DF   t-value p-value
# (Intercept)          0.0012687630 0.0002881005 62  4.403890  0.0000
# Temp                -0.0000138406 0.0000142073  1 -0.974188  0.5083
# EcotypeBenthic      -0.0000096421 0.0003239913  1 -0.029760  0.9811
# SexM                 0.0002195571 0.0000359801 62  6.102178  0.0000
# Temp:EcotypeBenthic  0.0000032576 0.0000158869  1  0.205052  0.8712
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.993                     
# EcotypeBenthic      -0.871  0.870              
# SexM                -0.225  0.163  0.120       
# Temp:EcotypeBenthic  0.884 -0.892 -0.993 -0.129
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.55321828 -0.70319956  0.09819431  0.66987611  1.97981652 
# 
# Number of Observations: 68
# Number of Groups: 5 
emtrends(BodyCond, "Ecotype", var = "Temp")
# Ecotype  Temp.trend       SE df  lower.CL upper.CL
# Limnetic  -1.38e-05 1.42e-05  1 -0.000194 1.67e-04
# Benthic   -1.06e-05 7.20e-06  1 -0.000102 8.08e-05
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(BodyCond)

BodyCond <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = "Temperature (\u00B0C)", y = "Body Condition")+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
BodyCond

BodyCondSex <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "right", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = "Temperature (\u00B0C)", y = "Body Condition")+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
BodyCondSex


{
  complete_data <- End[complete.cases(End[c("HSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$HSI)
  HSI <- lme(HSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(HSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -322.1056 -307.2157 168.0528
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)   Residual
# StdDev: 5.471947e-07 0.01365051
# 
# Fixed effects:  HSI ~ Temp * Ecotype + Sex 
#                           Value   Std.Error DF   t-value p-value
# (Intercept)           0.03863660 0.014932504 61  2.587416  0.0121
# Temp                  0.00084609 0.000674824  1  1.253789  0.4286
# EcotypeLimnetic      -0.07187344 0.031020703  1 -2.316951  0.2594
# SexM                 -0.00755736 0.003421042 61 -2.209082  0.0309
# Temp:EcotypeLimnetic  0.00404941 0.001517579  1  2.668338  0.2283
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.981                     
# EcotypeLimnetic      -0.444  0.464              
# SexM                 -0.165  0.037 -0.150       
# Temp:EcotypeLimnetic  0.408 -0.438 -0.993  0.157
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.03678937 -0.76192305  0.02431427  0.62039089  2.65341720 
# 
# Number of Observations: 67
# Number of Groups: 5 

emtrends(HSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend       SE df lower.CL upper.CL
# Benthic    0.000846 0.000675  1 -0.00773  0.00942
# Limnetic   0.004896 0.001360  1 -0.01244  0.02223
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(HSI)

HSI <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
HSI

HSISex <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISex


{
  complete_data <- End[complete.cases(End[c("GSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$GSI)
  GSI <- lme(GSI ~ Temp*Ecotype, random = ~ 1 | CCD, data = complete_data)
  summary(GSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC  logLik
# -277.6099 -269.2028 144.805
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)    Residual
# StdDev: 1.103456e-07 0.001518139
# 
# Fixed effects:  GSI ~ Temp * Ecotype 
#                             Value   Std.Error DF    t-value p-value
# (Intercept)           0.001455470 0.002839581 29  0.5125651  0.6121
# Temp                  0.000023353 0.000126883  1  0.1840479  0.8841
# EcotypeLimnetic       0.006955088 0.004856772  1  1.4320392  0.3881
# Temp:EcotypeLimnetic -0.000311165 0.000236880  1 -1.3135961  0.4142
# Correlation: 
#   (Intr) Temp   EctypL
# Temp                 -0.992              
# EcotypeLimnetic      -0.585  0.580       
# Temp:EcotypeLimnetic  0.532 -0.536 -0.992
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.4568847 -0.6731828 -0.1260429  0.3244948  2.2499836 
# 
# Number of Observations: 34
# Number of Groups: 5 
emtrends(GSI, "Ecotype", var = "Temp")
#  Ecotype  Temp.trend       SE df lower.CL upper.CL
# Benthic    2.34e-05 0.000127  1 -0.00159  0.00164
# Limnetic  -2.88e-04 0.000200  1 -0.00283  0.00225
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(GSI)

GSI <- ggplot(complete_data, aes(x=as.factor(Temp), y=GSI, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(GSI))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") 
GSI

{
  complete_data <- End[complete.cases(End[c("SSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$SSI)
  SSI <- lme(SSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC     BIC logLik
# -737.3199 -722.43 375.66
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)     Residual
# StdDev: 3.748528e-08 0.0004795104
# 
# Fixed effects:  SSI ~ Temp * Ecotype + Sex 
#                               Value    Std.Error DF    t-value p-value
# (Intercept)           0.0010874308 0.0005333518 61  2.0388621  0.0458
# Temp                 -0.0000043032 0.0000240878  1 -0.1786449  0.8875
# EcotypeLimnetic       0.0014879918 0.0010740497  1  1.3854032  0.3980
# SexM                 -0.0000411473 0.0001191858 61 -0.3452368  0.7311
# Temp:EcotypeLimnetic -0.0000437391 0.0000525819  1 -0.8318282  0.5583
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.982                     
# EcotypeLimnetic      -0.469  0.484              
# SexM                 -0.140  0.017 -0.131       
# Temp:EcotypeLimnetic  0.429 -0.456 -0.993  0.138
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.5068370 -0.7830405 -0.1569127  0.5279036  2.8039601 
# 
# Number of Observations: 67
# Number of Groups: 5 
emtrends(SSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend       SE df  lower.CL upper.CL
# Benthic    -4.3e-06 2.41e-05  1 -0.000310 0.000302
# Limnetic   -4.8e-05 4.68e-05  1 -0.000643 0.000547
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
plot(SSI)

SSI <- ggplot(complete_data, aes(x=as.factor(Temp), y=GSI, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSI


BodyCond + HSI + GSI + SSI


# Change BC by ecotype ----------------------------------------------------

head(ThermTol)

summary <- ThermTol %>%
  dplyr::group_by(Trial, Ecotype) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(BodyCond, na.rm = TRUE),
    median = median(BodyCond, na.rm = TRUE),
    sd = sd(BodyCond, na.rm = TRUE),
    se = (sd(BodyCond, na.rm = TRUE))/(sqrt(n())),
    min = min(BodyCond, na.rm = TRUE),
    max = max(BodyCond, na.rm = TRUE)
  )
summary
# Trial Ecotype      n    mean  median       sd        se      min     max
# Start Benthic     86 0.00124 0.00120 0.000227 0.0000245 0.000721 0.00183
# Start Limnetic    93 0.00119 0.00118 0.000188 0.0000195 0.000818 0.00168
# End   Benthic     41 0.00115 0.00117 0.000187 0.0000292 0.000848 0.00155
# End   Limnetic    27 0.00111 0.00108 0.000183 0.0000351 0.000766 0.00139

BenS<- filter(Benno, Trial == "Start") 
BenF <- filter(Benno, Trial == "End")

t.test(x = BenF$BodyCond, y = BenS$BodyCond, paired = FALSE)
# Welch Two Sample t-test
# 
# data:  BenF$BodyCond and BenS$BodyCond
# t = -2.1337, df = 94.168, p-value = 0.03547
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.572244e-04 -5.657283e-06
# sample estimates:
#   mean of x   mean of y 
# 0.001153592 0.001235033 

LimS <- filter(Limno, Trial == "Start")
LimF <- filter(Limno, Trial == "End")

t.test(x = LimF$BodyCond, y = LimS$BodyCond, paired = FALSE)
# Welch Two Sample t-test
# 
# data:  LimF$BodyCond and LimS$BodyCond
# t = -2.1562, df = 43.238, p-value = 0.03667
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -1.676289e-04 -5.618064e-06
# sample estimates:
#   mean of x   mean of y 
# 0.001106821 0.001193444 

anno <- data.frame(x1 = c(.8), x2 = c(1.2), x3 = c(1.82), x4 = c(2.18),
                   y1 = c(.00185), y2 = c(.0019),  y3 = c(.00185), y4 = c(.0019),
                   xstar = c(1), ystar = c(.00195), xstar2 = c(2), ystar2 = c(.00195),
                   Trial = c("p<0.05"))
anno

ggplot(ThermTol, aes(x = Ecotype, y = BodyCond, color = Trial)) +
  geom_boxplot(size = 0.5, width = 0.5, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Trial), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), 
              size = 1, alpha = 0.2) +
  theme_classic()+
  theme(legend.position = "right") +
  scale_color_manual(values = c('cadetblue', 'blue4', "black")) +
  labs(x = "Ecotype", y = "Body Condition") +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = Trial), size = 3, fontface = "bold") +
  geom_segment(data = anno, aes(x = x1, xend = x1, y = y1, yend = y2), colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, y = y1, yend = y2), colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), colour = "black") +
  geom_text(data = anno, aes(x = xstar2, y = ystar2, label = Trial), size = 3, fontface = "bold") +
  geom_segment(data = anno, aes(x = x3, xend = x3, y = y3, yend = y4), colour = "black") +
  geom_segment(data = anno, aes(x = x4, xend = x4, y = y3, yend = y4), colour = "black") +
  geom_segment(data = anno, aes(x = x3, xend = x4, y = y4, yend = y4), colour = "black")

mean(BenF$BodyCond) - mean(BenS$BodyCond) # -8.144086e-05 avg lost from starting BC
mean(LimF$BodyCond) - mean(LimS$BodyCond) # -8.662349e-05 avg lost from starting BC

(mean(BenF$BodyCond) - mean(BenS$BodyCond))/mean(BenS$BodyCond) # -0.06594227 avg proportion lost from starting BC for benthic fish

(mean(LimF$BodyCond) - mean(LimS$BodyCond))/mean(LimS$BodyCond) # -0.07258275 avg proportion lost from starting BC for limnetic fish

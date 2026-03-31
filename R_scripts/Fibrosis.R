# Analyze effect of/on fibrosis

# Three hypotheses:
# 1) With respect to immunity, we expected that higher temperatures would weaken fibrosis responses and increase the incidence of cestode infections due to higher parasite growth rates at warmer temperatures.
#    a) Model: Fibrosis ~ Temperature * Ecotype
# 2) When fibrosis did occur, we expected it to drive reduced AS (via higher SMR or lower MMR) and body condition.
#    b) Model: AS/SMR/MMR ~ Fibrosis + Ecotype
#    c) Model: BC ~ Fibrosis + Ecotype
# 3) Finally, we expected heat-induced decreases in fibrosis to be most pronounced in cestode-immune limnetic fish.
#    d) See model a: Fibrosis ~ Temperature * Ecotype
# 


# Hypotheses 1 & 3: Fibrosis ~ Temp * Ecotype -----------------------------

# Does temperature affect propensity to fibrose? ---------------------------------------
# Do limnetic fish show the most decline in fibrosis with temperature?



### No effect of temperature on fibrosis severity or propensity for either ecotype


complete_data <- End[complete.cases(End[c("Fibrosis", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data) #69

complete_data$Temp <- as.numeric(complete_data$Temp)
FibPresA <- glmer(Fibrosis ~ Temp * Ecotype + (1 | CCD), 
                  data = complete_data, 
                  family = binomial(link = "logit"))
plot(FibPresA)
summary(FibPresA)
#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: binomial  ( logit )
# Formula: Fibrosis ~ Temp * Ecotype + (1 | CCD)
# Data: complete_data
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 93.4     104.6     -41.7      83.4        64 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -0.9354 -0.6338 -0.5718  1.0690  1.9387 
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# CCD    (Intercept) 0        0       
# Number of obs: 69, groups:  CCD, 5
# 
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)
# (Intercept)          -2.25083    2.54022  -0.886    0.376
# Temp                  0.05149    0.11494   0.448    0.654
# EcotypeLimnetic       3.63142    4.59746   0.790    0.430
# Temp:EcotypeLimnetic -0.13561    0.22420  -0.605    0.545
# 
# Correlation of Fixed Effects:
#   (Intr) Temp   EctypL
# Temp        -0.990              
# EcotypLmntc -0.553  0.547       
# Tmp:EctypLm  0.507 -0.513 -0.992
# optimizer (Nelder_Mead) convergence code: 0 (OK)
# boundary (singular) fit: see help('isSingular')


# complete_data <- End[complete.cases(End[c("Ecotype", "CCD", "Sex", "Fibrosis", "Temp")]), ]
# complete_data <- complete_data %>%
#   mutate(
#     Temp_num = as.numeric(as.character(Temp)),
#     Fibrosis_num = ifelse(Fibrosis == "Y", 1, 0)
#   )



plot_data <- complete_data %>%
  mutate(Fibrosis_num = ifelse(Fibrosis == "Y", 1, 0),
         Temp_num = as.numeric(as.character(Temp)))

FibPres <- ggplot(plot_data, aes(x = Temp_num, y = Fibrosis_num, color = Ecotype)) +
  geom_jitter(height = 0.03, alpha = 0.4, size = 1.5) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  labs(x = "Temperature (°C)", y = "Probability of Fibrosis", color = "Ecotype", tag = "E") +
  scale_x_continuous(breaks = c(18, 20, 22, 24, 26)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_color_manual(values = c("Benthic" = "black", "Limnetic" = "blue")) +
  theme_classic()+
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98),
        legend.position = "none")
FibPres

complete_data <- End2[complete.cases(End2[c("Ecotype", "CCD", "Sex", "FibrosisScore", "Temp")]), ]
class(complete_data$Fibrosis)
class(complete_data$FibrosisScore)
class(complete_data$Temp)
complete_data$Temp <- as.numeric(complete_data$Temp)
complete_data$Fibrosis <- factor(complete_data$Fibrosis, levels = c("0", "1"))
complete_data$FibrosisScore <- factor(complete_data$FibrosisScore, levels = c("0", "1", "2", "3", "4"))
FibSevA <- clmm(FibrosisScore ~ Temp * Ecotype + (1 | CCD), 
                data = complete_data)
summary(FibSevA)
# Cumulative Link Mixed Model fitted with the Laplace approximation
# 
# formula: FibrosisScore ~ Temp * Ecotype + (1 | CCD)
# data:    complete_data
# 
# link  threshold nobs logLik  AIC    niter     max.grad cond.H 
# logit flexible  173  -194.77 405.55 738(1193) 6.73e-04 2.5e+10
# 
# Random effects:
#   Groups Name        Variance Std.Dev.
# CCD    (Intercept) 0.000878 0.02963 
# Number of groups:  CCD 17 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# Temp                  0.03816    0.09604   0.397    0.691
# EcotypeLimnetic       1.84465    7.83332   0.235    0.814
# Temp:EcotypeLimnetic -0.06859    0.28020  -0.245    0.807
# 
# Threshold coefficients:
#   Estimate Std. Error z value
# 0|1    1.370      1.613   0.849
# 1|2    2.324      1.676   1.386
# 2|3    3.176      1.737   1.829
# 3|4    6.163      2.058   2.994
FibSevA_2 <- clmm(FibrosisScore ~ Temp + Ecotype + (1 | CCD), 
                data = complete_data)
FibSevA_null <- clmm(FibrosisScore ~ 1 + (1 | CCD), data = complete_data)
anova(FibSevA_null, FibSevA)
# Likelihood ratio tests of cumulative link models:
#   
# formula:                                   link: threshold:
# FibSevA_null FibrosisScore ~ 1 + (1 | CCD)              logit flexible  
# FibSevA      FibrosisScore ~ Temp * Ecotype + (1 | CCD) logit flexible  
# 
#               no.par    AIC  logLik LR.stat df Pr(>Chisq)
# FibSevA_null      5 400.45 -195.23                      
# FibSevA           8 405.55 -194.77  0.9036  3     0.8246
anova(FibSevA_null, FibSevA_2)
# Likelihood ratio tests of cumulative link models:
#   
# formula:                                   link: threshold:
# FibSevA_null FibrosisScore ~ 1 + (1 | CCD)              logit flexible  
# FibSevA_2    FibrosisScore ~ Temp + Ecotype + (1 | CCD) logit flexible  
# 
#               no.par    AIC  logLik LR.stat df Pr(>Chisq)
# FibSevA_null      5 400.45 -195.23                      
# FibSevA_2         7 403.80 -194.90  0.6456  2     0.7241

complete_data$Temp <- factor(complete_data$Temp, levels = c("18", "20", "22", "24", "26"))

FibSev <- ggplot(complete_data, aes(x = Temp, y = FibrosisScore, color = Ecotype, group = Ecotype)) +
  geom_jitter(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.2, dodge.width = .5, jitter.height = .2)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar", width = 0.2, linewidth = 1,
               position = position_dodge(width = 0.5)) +
  labs(x = "Temperature", y = "Fibrosis Score", tag = "F") +
  scale_color_manual(values = c("Benthic" = "black", "Limnetic" = "blue")) +
  theme_classic()+
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98))
FibSev




# Hypothesis 2: Resp ~ Fibrosis + Ecotype -------------------------------------------

### Both fibrosis presence and severity affect SMR, MMR, and AS

## Assess the effects of fibrosis on MMR ----------------------------

###### Best model ######
# Add fibrosis to the best model for MMR
MMREndQ.fib <- lm(MMR ~ I(Temp^2) + Temp+ Ecotype + Sex + Fibrosis, data = End)
plot(MMREndQ.fib)
AIC(MMREndQ.fib) #1033.849
summary(MMREndQ.fib)
# Call:
#   lm(formula = MMR ~ I(Temp^2) + Temp + Ecotype + Sex + Fibrosis, 
#      data = End)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -944.07 -255.32   33.62  244.24 1305.26 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     15246.772   3088.403   4.937 6.14e-06 ***
# I(Temp^2)          30.419      6.785   4.483 3.17e-05 ***
# Temp            -1305.376    292.022  -4.470 3.33e-05 ***
# EcotypeLimnetic   356.246    110.204   3.233  0.00195 ** 
# SexM               63.501    100.217   0.634  0.52862    
# FibrosisY        -324.169    108.315  -2.993  0.00394 ** 
# 
# Residual standard error: 410.2 on 63 degrees of freedom
# Multiple R-squared:  0.3373,	Adjusted R-squared:  0.2847 
# F-statistic: 6.413 on 5 and 63 DF,  p-value: 7.104e-05

# what is the effect of fibrosis on MMR at the end of the experiment?
emmeans(MMREndQ.fib,~Fibrosis)
# Fibrosis emmean    SE df lower.CL upper.CL
# N          1461  81.9 63     1298     1625
# Y          1137 104.0 63      929     1346
# 
# Results are averaged over the levels of: Ecotype, Sex 
# Confidence level used: 0.95 

(1461-1137)/1461 #0.2217659 == 22% decrease in MMR in fish that fibrosed

# Did one ecotype fibrose more often than another?
#Use full dataset. Some fish may have been excluded earlier, but their fibrosis data still addresses ecotype differences
Fib.data <- read_xlsx('/home/ekerns/ThermTol/RawData/ThermalToleranceProject_OG.xlsx')
Fib.data <- Fib.data[!is.na(Fib.data$Fibrosis),]

Fib.data <- Fib.data %>%
  mutate(
    Fibrosis = case_when(
      Fibrosis == "Y" ~ "1",
      Fibrosis == "N" ~ "0",
      TRUE ~ "Yer missing something"  # A default value if none of the conditions match
    )
  ) 
Fib.data$Fibrosis <- as.numeric(Fib.data$Fibrosis)

#Did Fibrosis rates increase with temp and vary by ecotype?
Fib.ecotype.temp <- glm(Fibrosis ~ Ecotype + Temp,family = "binomial", data = Fib.data)
summary(Fib.ecotype.temp)
#                   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)     -0.638784   1.323204  -0.483   0.6293  
# EcotypeLimnetic  0.536734   0.314301   1.708   0.0877 .
# Temp             0.002344   0.060031   0.039   0.9689  
# 
# Null deviance: 236.02  on 173  degrees of freedom
# Residual deviance: 233.05  on 171  degrees of freedom
# AIC: 239.05

emmeans(Fib.ecotype.temp,~Ecotype)
# Ecotype   emmean    SE  df asymp.LCL asymp.UCL
# Benthic  -0.5885 0.212 Inf    -1.003    -0.174
# Limnetic -0.0518 0.231 Inf    -0.504     0.400
# 
# Results are given on the logit (not the response) scale. 
# Confidence level used: 0.95

exp(-0.5885) / (1 + exp(-0.5885)) #0.3569791 = Benthic fibrosis rate
exp(-0.0518) / (1 + exp(-0.0518)) #0.4870529 = Limnetic fibrosis rate

effectsize::standardize_parameters(Fib.ecotype.temp, exp = TRUE)
# # Standardization method: refit
# 
# Parameter          | Std_Odds_Ratio |       95% CI
# --------------------------------------------------
# (Intercept)        |           0.56 | [0.36, 0.83]
# Ecotype [Limnetic] |           1.71 | [0.93, 3.18]
# Temp               |           1.01 | [0.74, 1.37]
#
# - Response is unstandardized.



MMREndQ.fibS <- lm(MMR ~ I(Temp^2) + Temp+ Ecotype + Sex + FibrosisScore, data = End)
plot(MMREndQ.fibS)
AIC(MMREndQ.fibS) #1037.56
summary(MMREndQ.fibS)

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     14371.697   3182.580   4.516 2.83e-05 ***
# I(Temp^2)          28.604      6.991   4.091 0.000124 ***
# Temp            -1226.065    301.022  -4.073 0.000132 ***
# EcotypeLimnetic   336.360    112.736   2.984 0.004049 ** 
# SexM               63.818    103.095   0.619 0.538136    
# FibrosisScore    -134.977     59.248  -2.278 0.026117 *  
# 
# Residual standard error: 421.4 on 63 degrees of freedom
# Multiple R-squared:  0.3007,	Adjusted R-squared:  0.2452 
# F-statistic: 5.418 on 5 and 63 DF,  p-value: 0.0003305
PermTest(MMREndQ.fibS, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lm(obj = MMREndQ.fibS, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#   p.value
# I(Temp^2)      0.7265
# Temp           0.0000
# Ecotype        0.0099
# Sex            0.6725
# FibrosisScore  0.0256




## Assess the effects of fibrosis on SMR ----------------------------

###### Best model ######
# Add fibrosis to the best model for MMR
RMREndQ.a <- lm(RMR ~ I(Temp^2) + Temp + Ecotype + Sex + Fibrosis, data = End)
plot(RMREndQ.a)
AIC(RMREndQ.a) #895.0609
summary(RMREndQ.a)
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     5892.511   1129.686   5.216 2.17e-06 ***
# I(Temp^2)         13.040      2.482   5.254 1.87e-06 ***
# Temp            -527.506    106.817  -4.938 6.10e-06 ***
# EcotypeLimnetic   75.230     40.311   1.866  0.06666 .  
# SexM             -55.269     36.658  -1.508  0.13663    
# FibrosisY       -116.943     39.620  -2.952  0.00444 ** 
# 
# Residual standard error: 150.1 on 63 degrees of freedom
# Multiple R-squared:  0.5136,	Adjusted R-squared:  0.475 
# F-statistic: 13.31 on 5 and 63 DF,  p-value: 7.479e-09

# what is the effect of fibrosis on MMR at the end of the experiment?
emmeans(RMREndQ.a,~Fibrosis)
# Fibrosis emmean   SE df lower.CL upper.CL
# N           574 30.0 63      514      634
# Y           457 38.2 63      380      533
# 
# Results are averaged over the levels of: Ecotype, Sex 
# Confidence level used: 0.95 

(574-457)/574 #0.2038328 == 20% decrease in SMR in fish that fibrosed



SMREndQ.fibS <- lm(RMR ~ I(Temp^2) + Temp + Ecotype + Sex + FibrosisScore, data = End)
plot(SMREndQ.fibS)
AIC(SMREndQ.fibS) #896.1145
summary(SMREndQ.fibS)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -331.73  -73.10   -9.51   96.32  421.43 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     5531.523   1141.939   4.844 8.63e-06 ***
# I(Temp^2)         12.288      2.509   4.898 7.07e-06 ***
# Temp            -494.506    108.010  -4.578 2.26e-05 ***
# EcotypeLimnetic   71.367     40.451   1.764  0.08253 .  
# SexM             -53.698     36.992  -1.452  0.15157    
# FibrosisScore    -58.704     21.259  -2.761  0.00753 ** 
# 
# Residual standard error: 151.2 on 63 degrees of freedom
# Multiple R-squared:  0.5061,	Adjusted R-squared:  0.4669 
# F-statistic: 12.91 on 5 and 63 DF,  p-value: 1.185e-08
PermTest(SMREndQ.fibS, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lm(obj = SMREndQ.fibS, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#               p.value
# I(Temp^2)      0.0000
# Temp           0.0000
# Ecotype        0.1574
# Sex            0.0959
# FibrosisScore  0.0071


## Assess the effects of fibrosis on AS ----------------------------

# Add fibrosis to the best model
ASEndQ.fib <- lm(AS ~ I(Temp^2) + Temp+  Ecotype + Fibrosis, data = End)
plot(ASEndQ.fib)
AIC(ASEndQ.fib) # 998.106
summary(ASEndQ.fib)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -629.39 -218.70   14.45  181.09 1066.24 
# 
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      9391.38    2399.38   3.914 0.000223 ***
# I(Temp^2)          17.20       5.27   3.263 0.001770 ** 
# Temp             -772.61     226.86  -3.406 0.001144 ** 
# EcotypeLimnetic   271.49      85.39   3.179 0.002276 ** 
# Fibrosis         -198.96      83.98  -2.369 0.020860 *  
# 
# Residual standard error: 318.7 on 64 degrees of freedom
# Multiple R-squared:  0.3266,	Adjusted R-squared:  0.2845 
# F-statistic:  7.76 on 4 and 64 DF,  p-value: 3.659e-05

emmeans(ASEndQ.fib,~Fibrosis)
# Fibrosis emmean   SE df lower.CL upper.CL
# 0    892 63.6 64      765     1019
# 1    693 80.7 64      532      854
# 
# Results are averaged over the levels of: Ecotype 
# Confidence level used: 0.95 

(892-693)/693 #0.2871573, ~29% decrease with fibrosis


ASEndQ.fibS <- lm(AS ~ I(Temp^2) + Temp+  Ecotype + FibrosisScore, data = End)
plot(ASEndQ.fibS)
AIC(ASEndQ.fibS) #1001.381
summary(ASEndQ.fibS)
#Residuals:
# Min      1Q  Median      3Q     Max 
# -602.20 -202.89   28.03  174.99 1107.67 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     8909.006   2464.451   3.615 0.000592 ***
# I(Temp^2)         16.205      5.414   2.993 0.003922 ** 
# Temp            -729.349    233.135  -3.128 0.002645 ** 
# EcotypeLimnetic  255.347     87.067   2.933 0.004657 ** 
# FibrosisScore    -70.626     45.726  -1.545 0.127390    
# 
# Residual standard error: 326.4 on 64 degrees of freedom
# Multiple R-squared:  0.2939,	Adjusted R-squared:  0.2497 
# F-statistic: 6.659 on 4 and 64 DF,  p-value: 0.0001518
PermTest(ASEndQ.fibS, B = 10000)
# Monte-Carlo test
# 
# Call: 
#   PermTest.lm(obj = ASEndQ.fibS, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#               p.value
# I(Temp^2)      0.0034
# Temp           0.0065
# Ecotype        0.0081
# FibrosisScore  0.1284


## Make figure -------------------------------------------------------------

End$Fibrosis <- as.factor(End$Fibrosis)
End$Fibrosis

FibLabels <- c("Absent", "Present")

ASfib <- ggplot(End, aes(x=Fibrosis, y=AS, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98))+
  scale_color_manual(values = c('black', "blue")) +
  geom_segment(aes(x = 1, xend = 2, y = max(End$AS, na.rm=TRUE) * 1.05, 
                   yend = max(End$AS, na.rm=TRUE) * 1.05), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 1, xend = 1, y = max(End$AS, na.rm=TRUE) * 1.05, 
                   yend = max(End$AS, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 2, xend = 2, y = max(End$AS, na.rm=TRUE) * 1.05, 
                   yend = max(End$AS, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(End$AS, na.rm=TRUE) * 1.07, label = "*", size = 6) +
  scale_x_discrete(labels= FibLabels) +
  labs(x = "Fibrosis", y = expression(Aerobic~Scope~(mgO[2]/kg/hr)), tag = "C")
ASfib

SMRfib <- ggplot(End, aes(x=Fibrosis, y=RMR, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(legend.position = "none", plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98)) +
  scale_color_manual(values = c('black', "blue")) +
  geom_segment(aes(x = 1, xend = 2, y = max(End$RMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$RMR, na.rm=TRUE) * 1.05), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 1, xend = 1, y = max(End$RMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$RMR, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 2, xend = 2, y = max(End$RMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$RMR, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(End$RMR, na.rm=TRUE) * 1.07, label = "*", size = 6) +
  scale_x_discrete(labels= FibLabels) +
  labs(x = "Fibrosis", y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)), tag = "B")
SMRfib

MMRfib <- ggplot(End, aes(x=Fibrosis, y=MMR, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(legend.position = "none", plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98))+
  scale_color_manual(values = c('black', "blue")) +
  geom_segment(aes(x = 1, xend = 2, y = max(End$MMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$MMR, na.rm=TRUE) * 1.05), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 1, xend = 1, y = max(End$MMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$MMR, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  geom_segment(aes(x = 2, xend = 2, y = max(End$MMR, na.rm=TRUE) * 1.05, 
                   yend = max(End$MMR, na.rm=TRUE) * 1.03), 
               color = "black", inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = max(End$MMR, na.rm=TRUE) * 1.07, 
           label = "*", size = 6) +  
  scale_x_discrete(labels= FibLabels) +
  labs(x = "Fibrosis", y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)), tag = "A")
MMRfib

MMRFibSev <- ggplot(End, aes(x = FibrosisScore, y = MMR, color = Ecotype, group = Ecotype)) +
  geom_jitter(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.2, dodge.width = .5, jitter.height = .2)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar", width = 0.2, linewidth = 1,
               position = position_dodge(width = 0.5)) +
  labs(x = "Fibrosis Score", y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)), tag = "A") +
  scale_color_manual(values = c("Benthic" = "black", "Limnetic" = "blue")) +
  theme_classic()+
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98), 
        legend.position = "none")
MMRFibSev

SMRFibSev <- ggplot(End, aes(x = FibrosisScore, y = RMR, color = Ecotype, group = Ecotype)) +
  geom_jitter(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.2, dodge.width = .5, jitter.height = .2)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar", width = 0.2, linewidth = 1,
               position = position_dodge(width = 0.5)) +
  labs(x = "Fibrosis Score", y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)), tag = "B") +
  scale_color_manual(values = c("Benthic" = "black", "Limnetic" = "blue")) +
  theme_classic()+
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98),
        legend.position = "none")
SMRFibSev

ASFibSev <- ggplot(End, aes(x = FibrosisScore, y = AS, color = Ecotype, group = Ecotype)) +
  geom_jitter(alpha = 0.4, position = position_jitterdodge(jitter.width = 0.2, dodge.width = .5, jitter.height = .2)) +
  stat_summary(fun = mean, geom = "point", size = 4, shape = 18,
               position = position_dodge(width = 0.5)) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar", width = 0.2, linewidth = 1,
               position = position_dodge(width = 0.5)) +
  labs(x = "Fibrosis Score", y = expression(Aerobic~Scope~(mgO[2]/kg/hr)), tag = "C") +
  scale_color_manual(values = c("Benthic" = "black", "Limnetic" = "blue")) +
  theme_classic()+
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98))
ASFibSev

pdf(file = "NewFig4.pdf",
    width = 6.5, 
    height = 4)

MMRfib+ SMRfib + ASfib 

dev.off()

pdf(file = "NewSuppFig4.pdf",
    width = 6.5, 
    height = 4)

MMRFibSev + SMRFibSev + ASFibSev

dev.off()
http://144.92.58.154:8787/graphics/00173552-221d-4853-a3b3-2f65d9f85915.png
# Hypothesis 2: BC ~ Fibrosis ---------------------------------------------

BCEnd.fib.q.i <- lm(BodyCond ~ I(Temp^2)*Ecotype*Fibrosis + Temp*Ecotype*Fibrosis + Sex, data = End)
plot(BCEnd.fib.q.i)
AIC(BCEnd.fib.q.i) #-57.79971
summary(BCEnd.fib.q.i)

BCEnd.fib.q.a <- lm(BodyCond ~ I(Temp^2) + Temp + Ecotype + Sex + Fibrosis, data = End)
plot(BCEnd.fib.q.a)
AIC(BCEnd.fib.q.a) #-64.60953
summary(BCEnd.fib.q.a)

BCEnd.fib.i <- lm(BodyCond ~ Temp * Ecotype * Fibrosis + Sex, data = End)
plot(BCEnd.fib.i)
AIC(BCEnd.fib.i) #-58.4282
summary(BCEnd.fib.i)


BCEnd.fib.a <- lm(BodyCond ~ Temp + Ecotype + Sex + Fibrosis, data = End)
plot(BCEnd.fib.a)
AIC(BCEnd.fib.a) # -65.77794
summary(BCEnd.fib.a)
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      1.271789   0.141718   8.974 6.28e-13 ***
# Temp            -0.010519   0.006308  -1.668    0.100    
# EcotypeLimnetic -0.049741   0.037519  -1.326    0.190    
# SexM             0.220083   0.034923   6.302 3.07e-08 ***
# FibrosisY       -0.058816   0.037743  -1.558    0.124    
# 
# Residual standard error: 0.143 on 64 degrees of freedom
# Multiple R-squared:  0.433,	Adjusted R-squared:  0.3976 
# F-statistic: 12.22 on 4 and 64 DF,  p-value: 1.936e-07

BC.fib.plot <- ggplot(End, aes(x=Fibrosis, y=BodyCond, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98))+
  scale_color_manual(values = c('black', "blue")) +
  scale_x_discrete(labels= FibLabels) +
  labs(x = "Fibrosis", y = "Body Condition", tag = "A")
BC.fib.plot

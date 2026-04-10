# OCLTTH Prediction 4 - Effect of AS on condition metrics
# Figure 4
# Effect of temperature directly on condition metrics


# Test Assumption 4 of the OCLTTH ---------------------------------------------------------


## Start of experiment -------------------------------------------------------

### Fulton's K Body Condition -------------------------------------------------------


#### Does AS influence Body Condition (by ecotype)? ------------------


{
  complete_data <- Start[complete.cases(Start[c("BodyCond", "AS", "Ecotype")]), ]
  nrow(complete_data) # 183
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- glm(BodyCond ~ AS*Ecotype, data = complete_data)
  summary(BodyCond)
}
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        1.129e+00  6.694e-02  16.869   <2e-16 ***
# AS                 5.803e-05  5.731e-05   1.012   0.3127    
# EcotypeBenthic     2.435e-01  1.050e-01   2.319   0.0215 *  
# AS:EcotypeBenthic -1.911e-04  9.429e-05  -2.026   0.0442 *  
# 
# Null deviance: 7.9080  on 182  degrees of freedom
# Residual deviance: 7.6525  on 179  degrees of freedom
# AIC: -51.593
slopes <- emtrends(BodyCond, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE  df  lower.CL upper.CL
# Limnetic  0.000058 5.73e-05 179 -5.51e-05 1.71e-04
# Benthic  -0.000133 7.49e-05 179 -2.81e-04 1.47e-05
pairs(slopes)
# contrast           estimate       SE  df t.ratio p.value
# Limnetic - Benthic 0.000191 9.43e-05 179   2.026  0.0442

plot(BodyCond)

BodyCondP <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm"), axis.title = element_text(size = 10)) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black'))+
  labs(x = expression(Experiment~Start~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Experiment~Start~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")+
  annotate("text", x = 2200, y = .85, label = "Ecotype: p<0.05\nAS*Ecotype: p<0.05", hjust = 1, vjust = 1, size = 2, fontface = "bold")
BodyCondP

#### Does temperature directly influence Body Condition (by ecotype)? ------------------

{
  # temperature is not treated as a factor, not dropping outliers or mortalities from during respirometry measurements
  complete_data <- Start2[complete.cases(Start2[c("BodyCond", "Temp", "Ecotype")]), ]
  nrow(complete_data) # 208
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- glm(BodyCond ~ Temp*Ecotype, data = complete_data)
  summary(BodyCond)
}
#                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          1.339301   0.174826   7.661 7.32e-13 ***
# Temp                -0.006805   0.007893  -0.862    0.390    
# EcotypeBenthic      -0.334886   0.242344  -1.382    0.169    
# Temp:EcotypeBenthic  0.017784   0.010934   1.627    0.105    
# 
# (Dispersion parameter for gaussian family taken to be 0.04808878)
# 
# Null deviance: 10.1115  on 207  degrees of freedom
# Residual deviance:  9.8101  on 204  degrees of freedom
# AIC: -34.979
slopes <- emtrends(BodyCond, "Ecotype", var = "Temp")
# Ecotype  Temp.trend      SE  df lower.CL upper.CL
# Limnetic   -0.00681 0.00789 204 -0.02237  0.00876
# Benthic     0.01098 0.00757 204 -0.00394  0.02590
pairs(slopes)
# contrast           estimate     SE  df t.ratio p.value
# Limnetic - Benthic  -0.0178 0.0109 204  -1.627  0.1054
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
       y = "Experiment Start Body Condition")+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCond_Temps



#### Does AS influence Fulton's K Body Condition (by population)? ------------------
{
  complete_data <- Start[complete.cases(Start[c("BodyCond", "AS", "Pop")]), ]
  nrow(complete_data) #183
  
  hist(complete_data$BodyCond)
  BodyCond <- glm(BodyCond ~ AS*Pop, data = complete_data)
  summary(BodyCond)
}
#                   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  1.7903241  0.1831811   9.774  < 2e-16 ***
#   AS          -0.0005071  0.0001814  -2.795 0.005766 ** 
#   PopSL       -0.6739603  0.1995958  -3.377 0.000904 ***
#   PopWK       -0.7562845  0.2200382  -3.437 0.000735 ***
#   PopWT       -0.5312568  0.2034846  -2.611 0.009817 ** 
#   AS:PopSL     0.0006125  0.0001950   3.140 0.001981 ** 
#   AS:PopWK     0.0005775  0.0002046   2.823 0.005314 ** 
#   AS:PopWT     0.0004619  0.0001983   2.330 0.020965 *  
# 
# (Dispersion parameter for gaussian family taken to be 0.04026403)
# 
# Null deviance: 7.9080  on 182  degrees of freedom
# Residual deviance: 7.0462  on 175  degrees of freedom
# AIC: -58.699
slopes <- emtrends(BodyCond, "Pop", var = "AS")
# Pop  AS.trend       SE  df  lower.CL  upper.CL
# FG  -5.07e-04 1.81e-04 175 -8.65e-04 -0.000149
# SL   1.05e-04 7.16e-05 175 -3.59e-05  0.000247
# WK   7.04e-05 9.46e-05 175 -1.16e-04  0.000257
# WT  -4.52e-05 8.00e-05 175 -2.03e-04  0.000113
pairs(slopes)
# contrast  estimate       SE  df t.ratio p.value
# FG - SL  -0.000612 0.000195 175  -3.140  0.0106
# FG - WK  -0.000577 0.000205 175  -2.823  0.0270
# FG - WT  -0.000462 0.000198 175  -2.330  0.0953
# SL - WK   0.000035 0.000119 175   0.295  0.9910
# SL - WT   0.000151 0.000107 175   1.403  0.4989
# WK - WT   0.000116 0.000124 175   0.934  0.7869
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 
plot(BodyCond)

BodyCondPops <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Pop)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = expression(Experiment~Start~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Experiment~Start~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondPops

#### Does temperature influence Fulton's K Body Condition (by population)? ------------------
{
  # temp is not a factor, using all fish
  complete_data <- Start2[complete.cases(Start2[c("BodyCond", "Temp", "Pop")]), ]
  nrow(complete_data) # 208
  
  hist(complete_data$BodyCond)
  BodyCond <- glm(BodyCond ~ Temp*Pop, data = complete_data)
  summary(BodyCond)
}
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  0.30235    0.27688   1.092 0.276144    
#   Temp         0.04648    0.01252   3.713 0.000265 ***
#   PopSL        1.06532    0.33987   3.135 0.001980 ** 
#   PopWK        0.82217    0.41607   1.976 0.049525 *  
#   PopWT        1.03718    0.33833   3.066 0.002472 ** 
#   Temp:PopSL  -0.05282    0.01540  -3.430 0.000733 ***
#   Temp:PopWK  -0.04750    0.01862  -2.551 0.011475 *  
#   Temp:PopWT  -0.05245    0.01527  -3.434 0.000722 ***
# 
# (Dispersion parameter for gaussian family taken to be 0.04321741)
# 
# Null deviance: 10.1115  on 207  degrees of freedom
# Residual deviance:  8.6435  on 200  degrees of freedom
# AIC: -53.314


slopes <- emtrends(BodyCond, "Pop", var = "Temp")
# Pop Temp.trend      SE  df lower.CL upper.CL
# FG     0.04648 0.01250 200   0.0218   0.0712
# SL    -0.00634 0.00897 200  -0.0240   0.0113
# WK    -0.00103 0.01380 200  -0.0282   0.0262
# WT    -0.00598 0.00875 200  -0.0232   0.0113
pairs(slopes)
# contrast  estimate     SE  df t.ratio p.value
# FG - SL   0.052815 0.0154 200   3.430  0.0041
# FG - WK   0.047502 0.0186 200   2.551  0.0553
# FG - WT   0.052454 0.0153 200   3.434  0.0040
# SL - WK  -0.005313 0.0164 200  -0.323  0.9883
# SL - WT  -0.000361 0.0125 200  -0.029  1.0000
# WK - WT   0.004952 0.0163 200   0.303  0.9903
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 

plot(BodyCond)

BodyCond_TempsPop <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = "Temperature (\u00B0C)", 
       y = "Experiment Start Body Condition")+
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCond_TempsPop



## End of experiment -------------------------------------------------------



### Fulton's K Body Condition -------------------------------------------------------


#### Does AS influence Body Condition (by ecotype)? ------------------

{
complete_data <- End[complete.cases(End[c("BodyCond", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data) #69

hist(complete_data$BodyCond)
shapiro.test(complete_data$BodyCond)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
BodyCond <- lme(BodyCond ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(BodyCond)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC   logLik
# -13.09856 2.013624 13.54928
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:   0.0330121 0.141199
# 
# Fixed effects:  BodyCond ~ AS * Ecotype + Sex 
#                         Value  Std.Error DF   t-value p-value
# (Intercept)        1.1203165 0.08080022 61 13.865266  0.0000
# AS                -0.0001194 0.00006529 61 -1.828459  0.0724
# EcotypeBenthic    -0.1823549 0.10777054  3 -1.692066  0.1892
# SexM               0.2122160 0.03553555 61  5.971934  0.0000
# AS:EcotypeBenthic  0.0002318 0.00010285 61  2.253585  0.0278
# Correlation: 
#   (Intr) AS     EctypB SexM  
# AS                -0.868                     
# EcotypeBenthic    -0.742  0.651              
# SexM              -0.216 -0.020  0.125       
# AS:EcotypeBenthic  0.590 -0.631 -0.898 -0.167
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.30762961 -0.62546694  0.07094707  0.60470216  2.16414472 
# 
# Number of Observations: 69
# Number of Groups: 5 
slopes <- emtrends(BodyCond, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Limnetic -0.000119 6.53e-05 61 -2.50e-04 1.12e-05
# Benthic   0.000112 7.98e-05 61 -4.71e-05 2.72e-04
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95  
pairs(slopes)
# contrast            estimate       SE df t.ratio p.value
# Limnetic - Benthic -0.000232 0.000103 61  -2.254  0.0278
plot(BodyCond)

BodyCondEnd <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm"), axis.title = element_text(size = 10)) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Experiment~End~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Experiment~End~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")+
  annotate("text", x = 2200, y = .8, label = "AS*Ecotype: p<0.05", hjust = 1, vjust = 1, size = 2, fontface = "bold")
BodyCondEnd

BodyCondSexEnd <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondSexEnd

#### Does temperature directly influence Body Condition (by ecotype)? ------------------

{
  complete_data <- End2[complete.cases(End2[c("BodyCond", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #173
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  BodyCond <- lme(BodyCond ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(BodyCond)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -54.44707 -32.57933 34.22354
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:  0.09069957 0.1748278
# 
# Fixed effects:  BodyCond ~ Temp * Ecotype + Sex 
#                         Value Std.Error  DF   t-value p-value
# (Intercept)          0.3590912 0.3362246 152  1.068010  0.2872
# Temp                 0.0298317 0.0153191 152  1.947359  0.0533
# EcotypeBenthic       0.9604373 0.3974103 152  2.416740  0.0168
# SexM                 0.2204428 0.0275546 152  8.000222  0.0000
# Temp:EcotypeBenthic -0.0389242 0.0184628 152 -2.108253  0.0366
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.992                     
# EcotypeBenthic      -0.769  0.774              
# SexM                -0.189  0.152  0.156       
# Temp:EcotypeBenthic  0.760 -0.776 -0.992 -0.152
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.64423004 -0.63598217  0.01690286  0.65664318  2.71247241 
# 
# Number of Observations: 173
# Number of Groups: 17 
slopes <- emtrends(BodyCond, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE  df  lower.CL upper.CL
# Limnetic    0.02983 0.0153 152 -0.000434   0.0601
# Benthic    -0.00909 0.0117 152 -0.032199   0.0140
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate     SE  df t.ratio p.value
# Limnetic - Benthic   0.0389 0.0185 152   2.108  0.0366
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(BodyCond)

BodyCondTempEnd <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Ecotype)) +
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
BodyCondTempEnd

BodyCondTempSex <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Sex)) +
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
BodyCondTempSex

#### Does AS influence Fulton's K Body Condition (by population)? ------------------

{
  complete_data <- End[complete.cases(End[c("BodyCond", "AS", "Pop", "CCD", "Sex")]), ]
  nrow(complete_data) #69
  
  hist(complete_data$BodyCond)
  shapiro.test(complete_data$BodyCond)
  BodyCond <- lme(BodyCond ~ AS*Pop + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(BodyCond)
}

slopes <- emtrends(BodyCond, "Pop", var = "AS")

pairs(slopes)

plot(BodyCond)

BodyCondPopEnd <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Pop)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm"), axis.title = element_text(size = 10)) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue', "grey", "lightblue2"))+
  labs(x = expression(Experiment~End~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Experiment~End~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
  # annotate("text", x = 2200, y = .8, label = "AS*Ecotype: p<0.05", hjust = 1, vjust = 1, size = 2, fontface = "bold")
BodyCondPopEnd

BodyCondSexEnd <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none")+
  scale_color_manual( values=c('black','blue'))+
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondSexEnd


{
complete_data <- End[complete.cases(End[c("HSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
hist(complete_data$HSI)
HSI <- lme(HSI ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
summary(HSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 273.1591 288.1611 -129.5796
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:   0.7689908 1.323656
# 
# Fixed effects:  HSI ~ AS * Ecotype + Sex 
#                       Value Std.Error DF   t-value p-value
# (Intercept)        6.942132 0.9324836 60  7.444777  0.0000
# AS                -0.000242 0.0006410 60 -0.377795  0.7069
# EcotypeBenthic    -2.875747 1.2332186  3 -2.331904  0.1020
# SexM              -1.032460 0.3441984 60 -2.999606  0.0039
# AS:EcotypeBenthic  0.002356 0.0010146 60  2.321919  0.0236
# Correlation: 
#   (Intr) AS     EctypB SexM  
# AS                -0.738                     
# EcotypeBenthic    -0.754  0.558              
# SexM              -0.209  0.008  0.148       
# AS:EcotypeBenthic  0.509 -0.633 -0.773 -0.211
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.98606855 -0.51243502 -0.03741655  0.44481887  2.28886095 
# 
# Number of Observations: 68
# Number of Groups: 5
slopes  <- emtrends(HSI, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Limnetic -0.000242 0.000641 60 -0.001524  0.00104
# Benthic   0.002114 0.000785 60  0.000543  0.00368
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate      SE df t.ratio p.value
# Limnetic - Benthic -0.00236 0.00101 60  -2.322  0.0236
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
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
  annotate("text", x = 2200, y = 2, label = "AS: p<0.01\nAS*Ecotype: p<0.01", hjust = 1, vjust = 1, size = 2, fontface = "bold")
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
# AIC      BIC    logLik
# 16.99931 25.60324 -2.499656
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev: 2.146464e-06 0.1478172
# 
# Fixed effects:  GSI ~ AS * Ecotype 
#                           Value  Std.Error DF    t-value p-value
# (Intercept)         0.20993899 0.09388760 28  2.2360674  0.0335
# AS                 -0.00001459 0.00010195 28 -0.1430588  0.8873
# EcotypeLimnetic    -0.11584372 0.12757515  3 -0.9080429  0.4308
# AS:EcotypeLimnetic  0.00019803 0.00012428 28  1.5933680  0.1223
# Correlation: 
#   (Intr) AS     EctypL
# AS                 -0.932              
# EcotypeLimnetic    -0.736  0.686       
# AS:EcotypeLimnetic  0.765 -0.820 -0.913
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.6542484 -0.6016640 -0.3972967  0.3082361  2.3822313 
# 
# Number of Observations: 35
# Number of Groups: 5 
emtrends(GSI, "Ecotype", var = "AS")
#  Ecotype   AS.trend       SE df  lower.CL upper.CL
# Benthic  -1.46e-05 1.02e-04 28 -2.23e-04 0.000194
# Limnetic  1.83e-04 7.11e-05 28  3.78e-05 0.000329
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
# AIC       BIC   logLik
# -149.7611 -134.7592 81.88057
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)   Residual
# StdDev: 1.196499e-06 0.04790254
# 
# Fixed effects:  SSI ~ AS * Ecotype + Sex 
#                           Value  Std.Error DF   t-value p-value
# (Intercept)         0.08456832 0.02306683 60  3.666231  0.0005
# AS                  0.00001816 0.00002621 60  0.692669  0.4912
# EcotypeLimnetic     0.06470861 0.03379753  3  1.914596  0.1514
# SexM               -0.00487584 0.01188527 60 -0.410242  0.6831
# AS:EcotypeLimnetic -0.00000736 0.00003354 60 -0.219519  0.8270
# Correlation: 
#   (Intr) AS     EctypL SexM  
# AS                 -0.903                     
# EcotypeLimnetic    -0.666  0.644              
# SexM               -0.107 -0.184 -0.078       
# AS:EcotypeLimnetic  0.709 -0.775 -0.929  0.110
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.3125291 -0.8208460 -0.1773247  0.3707781  2.8113168 
# 
# Number of Observations: 68
# Number of Groups: 5 
emtrends(SSI, "Ecotype", var = "AS")
# Ecotype  AS.trend       SE df  lower.CL upper.CL
# Benthic  1.82e-05 2.62e-05 60 -3.43e-05 7.06e-05
# Limnetic 1.08e-05 2.12e-05 60 -3.16e-05 5.32e-05
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
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
# 973.477 984.4252 -481.7385
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:     183.234 319.8836
# 
# Fixed effects:  AS ~ Fibrosis + Ecotype 
#                     Value Std.Error DF   t-value p-value
# (Intercept)      901.2891  119.0035 63  7.573635  0.0000
# Fibrosis        -202.0651   84.4671 63 -2.392234  0.0197
# EcotypeLimnetic  258.9081  185.5241  3  1.395550  0.2572
# Correlation: 
#   (Intr) Fibrss
# Fibrosis        -0.175       
# EcotypeLimnetic -0.608 -0.082
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -1.98736249 -0.66688853  0.07106841  0.59264074  3.34260167 
# 
# Number of Observations: 69
# Number of Groups: 5 
emmeans(FibPresA, ~ Fibrosis, at = list(Fibrosis = c(0,1)))
# Fibrosis emmean    SE df lower.CL upper.CL
# 0   1031  96.7  3      723     1339
# 1    829 108.0  3      485     1173
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
# 977.5992 988.5475 -483.7996
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:    171.2687 327.9803
# 
# Fixed effects:  AS ~ FibrosisScore + Ecotype 
#                     Value Std.Error DF   t-value p-value
# (Intercept)     880.6867 112.97386 63  7.795491  0.0000
# FibrosisScore   -75.0849  45.93739 63 -1.634505  0.1071
# EcotypeLimnetic 247.6439 176.62582  3  1.402082  0.2554
# Correlation: 
#   (Intr) FbrssS
# FibrosisScore   -0.158       
# EcotypeLimnetic -0.610 -0.086
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8881286 -0.5707999  0.1316754  0.5519239  3.4042313 
# 
# Number of Observations: 69
# Number of Groups: 5 
PermTest(FibSevA, B = 10000)
# Monte-Carlo test
# 
# Call: 
# PermTest.lme(obj = FibSevA, B = 10000)
# 
# Based on 10000 replicates
# Simulated p-value:
#   p.value
# (Intercept)    0.9952
# FibrosisScore  0.1074
# Ecotype        0.1332
emmeans(FibSevA, ~ FibrosisScore, at = list(FibrosisScore = c(0,1, 2, 3)))
# FibrosisScore emmean    SE df lower.CL upper.CL
# 0   1005  91.6  3      713     1296
# 1    929  90.4  3      642     1217
# 2    854 110.0  3      503     1205
# 3    779 143.0  3      325     1234
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

pdf(file = "NewFigure5.pdf",
   width = 8,
   height = 8)

BodyCond + HSI + GSI + SSI + plot_layout(ncol = 2)

dev.off()


# Scaled Mass Index -------------------------------------------------------



# based on Peig & Green, New perspectives for estimating body condition from mass/length data:
#     the scaled mass index as an alternative method"
#     Oikos 118: 1883-1891, 2009
# script from https://apansharing.blogspot.com/2018/05/an-r-function-olsrobust-caled-mass-index.html




## Format data -------------------------------------------------------------

{
data <- read_xlsx("12125_ThermTol_OutliersRemoved.xlsx")
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
#                         Value  Std.Error DF   t-value p-value
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
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
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
  complete_data <- End2[complete.cases(End2[c("BodyCond", "Temp", "Ecotype", "CCD", "Sex")]), ]
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
# -54.44707 -32.57933 34.22354
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:  0.09069957 0.1748278
# 
# Fixed effects:  BodyCond ~ Temp * Ecotype + Sex 
#                         Value Std.Error  DF   t-value p-value
# (Intercept)          0.3590912 0.3362246 152  1.068010  0.2872
# Temp                 0.0298317 0.0153191 152  1.947359  0.0533
# EcotypeBenthic       0.9604373 0.3974103 152  2.416740  0.0168
# SexM                 0.2204428 0.0275546 152  8.000222  0.0000
# Temp:EcotypeBenthic -0.0389242 0.0184628 152 -2.108253  0.0366
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.992                     
# EcotypeBenthic      -0.769  0.774              
# SexM                -0.189  0.152  0.156       
# Temp:EcotypeBenthic  0.760 -0.776 -0.992 -0.152
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.64423004 -0.63598217  0.01690286  0.65664318  2.71247241 
# 
# Number of Observations: 173
# Number of Groups: 17 
slopes <- emtrends(BodyCond, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE  df  lower.CL upper.CL
# Limnetic    0.02983 0.0153 152 -0.000434   0.0601
# Benthic    -0.00909 0.0117 152 -0.032199   0.0140
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate     SE  df t.ratio p.value
# Limnetic - Benthic   0.0389 0.0185 152   2.108  0.0366
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
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
  complete_data <- End2[complete.cases(End2[c("HSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$HSI)
  HSI <- lme(HSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(HSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 640.1442 661.9281 -313.0721
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:   0.8593782 1.399091
# 
# Fixed effects:  HSI ~ Temp * Ecotype + Sex 
#                         Value Std.Error  DF    t-value p-value
# (Intercept)           4.862056  2.273721 150  2.1383695  0.0341
# Temp                  0.039566  0.104260 150  0.3794964  0.7049
# EcotypeLimnetic      -2.629109  3.505910 150 -0.7499080  0.4545
# SexM                 -0.479119  0.222071 150 -2.1575049  0.0326
# Temp:EcotypeLimnetic  0.103714  0.164244 150  0.6314668  0.5287
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.990                     
# EcotypeLimnetic      -0.526  0.537              
# SexM                 -0.003 -0.039 -0.158       
# Temp:EcotypeLimnetic  0.528 -0.549 -0.993  0.153
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -3.07361304 -0.66933307 -0.09814904  0.69956598  3.07538879 
# 
# Number of Observations: 171
# Number of Groups: 17 

slopes <- emtrends(HSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend    SE  df lower.CL upper.CL
# Benthic      0.0396 0.104 150   -0.166    0.246
# Limnetic     0.1433 0.138 150   -0.129    0.416
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
  complete_data <- End2[complete.cases(End2[c("GSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$GSI)
  GSI <- lme(GSI ~ Temp*Ecotype, random = ~ 1 | CCD, data = complete_data)
  summary(GSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -34.93656 -21.03163 23.46828
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.1016333 0.1438734
# 
# Fixed effects:  GSI ~ Temp * Ecotype 
#                             Value Std.Error DF    t-value p-value
# (Intercept)           0.28649545 0.3167809 60  0.9043961  0.3694
# Temp                 -0.00049958 0.0143268 60 -0.0348704  0.9723
# EcotypeLimnetic      -0.28888097 0.5001082 60 -0.5776369  0.5657
# Temp:EcotypeLimnetic  0.01584557 0.0234527 60  0.6756401  0.5019
# Correlation: 
#   (Intr) Temp   EctypL
# Temp                 -0.992              
# EcotypeLimnetic      -0.561  0.561       
# Temp:EcotypeLimnetic  0.547 -0.557 -0.993
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.7148280 -0.6691286 -0.1769502  0.5501603  2.2241670 
# 
# Number of Observations: 79
# Number of Groups: 16 
emtrends(GSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE df lower.CL upper.CL
# Benthic     -0.0005 0.0143 60  -0.0292   0.0282
# Limnetic     0.0153 0.0195 60  -0.0237   0.0544
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
  complete_data <- End2[complete.cases(End2[c("SSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data)
  
  hist(complete_data$SSI)
  SSI <- lme(SSI ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -344.5849 -323.0151 179.2925
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)   Residual
# StdDev: 0.008375683 0.07319862
# 
# Fixed effects:  SSI ~ Temp * Ecotype + Sex 
#                           Value  Std.Error  DF    t-value p-value
# (Intercept)           0.16028067 0.06287392 145  2.5492392  0.0118
# Temp                 -0.00148126 0.00287941 145 -0.5144322  0.6077
# EcotypeLimnetic       0.02442189 0.10931135 145  0.2234158  0.8235
# SexM                 -0.02590979 0.01161944 145 -2.2298666  0.0273
# Temp:EcotypeLimnetic  0.00042133 0.00509331 145  0.0827231  0.9342
# Correlation: 
#   (Intr) Temp   EctypL SexM  
# Temp                 -0.988                     
# EcotypeLimnetic      -0.566  0.576              
# SexM                 -0.019 -0.070 -0.175       
# Temp:EcotypeLimnetic  0.551 -0.571 -0.994  0.178
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.2695439 -0.6049009 -0.1958677  0.3002341  7.7166816 
# 
# Number of Observations: 166
# Number of Groups: 17 
emtrends(SSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend      SE  df lower.CL upper.CL
# Benthic    -0.00148 0.00288 145 -0.00717  0.00421
# Limnetic   -0.00106 0.00418 145 -0.00933  0.00721
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


BodyCond + HSI + GSI + SSI + FibPres + FibSev + plot_layout(ncol = 2)


# Change BC by ecotype ----------------------------------------------------

head(ThermTol)

summary <- ThermTol %>%
  dplyr::group_by(Trial, Ecotype, Temp) %>%
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
# Trial Ecotype  Temp      n  mean median    sd     se   min   max
# Start Benthic  18       23  1.17   1.16 0.215 0.0448 0.721  1.52
# Start Benthic  20       10  1.26   1.16 0.300 0.0949 0.833  1.67
# Start Benthic  22       21  1.28   1.26 0.240 0.0523 0.924  1.83
# Start Benthic  24       21  1.27   1.22 0.230 0.0502 0.844  1.69
# Start Benthic  26       13  1.20   1.14 0.161 0.0447 0.980  1.43
# Start Limnetic 18       18  1.24   1.22 0.176 0.0414 0.938  1.57
# Start Limnetic 20       22  1.19   1.23 0.202 0.0430 0.818  1.63
# Start Limnetic 22       21  1.18   1.10 0.183 0.0399 0.942  1.57
# Start Limnetic 24       21  1.25   1.20 0.198 0.0431 0.923  1.68
# Start Limnetic 26       13  1.07   1.06 0.107 0.0298 0.896  1.22
# End   Benthic  18       15  1.20   1.16 0.227 0.0585 0.848  1.55
# End   Benthic  22       15  1.14   1.20 0.159 0.0410 0.874  1.36
# End   Benthic  26       11  1.10   1.17 0.164 0.0494 0.853  1.29
# End   Limnetic 18       15  1.16   1.16 0.167 0.0432 0.770  1.39
# End   Limnetic 22       13  1.05   1.06 0.182 0.0504 0.766  1.38

BenS<- filter(Benno, Trial == "Start") 
BenF <- filter(Benno, Trial == "End")

t.test(x = BenF$BodyCond, y = BenS$BodyCond, paired = FALSE)
# Welch Two Sample t-test
# 
# data:  BenF$BodyCond and BenS$BodyCond
# t = -2.1228, df = 94.165, p-value = 0.03639
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.156498285 -0.005231633
# sample estimates:
#   mean of x mean of y 
# 1.153592  1.234457 

LimS <- filter(Limno, Trial == "Start")
LimF <- filter(Limno, Trial == "End")

t.test(x = LimF$BodyCond, y = LimS$BodyCond, paired = FALSE)
# Welch Two Sample t-test
# 
# data:  LimF$BodyCond and LimS$BodyCond
# t = -2.2521, df = 45.628, p-value = 0.02918
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.165960558 -0.009288293
# sample estimates:
#   mean of x mean of y 
# 1.105866  1.193490 

anno <- data.frame(x1 = c(.8), x2 = c(1.2), x3 = c(1.82), x4 = c(2.18),
                   y1 = c(1.85), y2 = c(1.9),  y3 = c(1.85), y4 = c(1.9),
                   xstar = c(1), ystar = c(1.95), xstar2 = c(2), ystar2 = c(1.95),
                   lab = c("p<0.05"))
anno

ThermTol$Trial <- factor(ThermTol$Trial, levels = c("Start", "End"))
ThermTol$Temp <- factor(ThermTol$Temp, levels = c("18", "20", "22", "24", "26"))

ChangeBC <- ggplot(ThermTol, aes(x = Ecotype, y = BodyCond, color = Trial)) +
  geom_boxplot(size = 0.5, width = 0.4, position = position_dodge(0.75)) +
  geom_jitter(aes(color = Trial, shape = Temp, group = Trial), 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.75), 
              size = 2, alpha = 0.4) +
  theme_classic(base_size = 12) +  
  theme(plot.tag = element_text(size = 15, face = "bold"),
        plot.tag.position = c(0.02, .98),
        legend.position = "right", 
        legend.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text("Timepoint")) +
  scale_color_manual(name = "Timepoint", labels = c("Start", "End"), values = c('cadetblue', "blue4")) +
  labs(x = "Ecotype", y = "Body Condition", tag = "C") +
  geom_text(data = anno, aes(x = xstar, y = ystar, label = lab), 
            size = 3, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = anno, aes(x = xstar2, y = ystar2, label = lab), size = 3, fontface = "bold", inherit.aes = FALSE) +  
  geom_segment(data = anno, aes(x = x1, xend = x1, y = y1, yend = y2), 
               colour = "black", inherit.aes = FALSE) +
  geom_segment(data = anno, aes(x = x2, xend = x2, y = y1, yend = y2), 
               colour = "black", inherit.aes = FALSE) +
  geom_segment(data = anno, aes(x = x1, xend = x2, y = y2, yend = y2), 
               colour = "black", inherit.aes = FALSE) +
  geom_segment(data = anno, aes(x = x3, xend = x3, y = y3, yend = y4), 
               colour = "black", inherit.aes = FALSE) +
  geom_segment(data = anno, aes(x = x4, xend = x4, y = y3, yend = y4), 
               colour = "black", inherit.aes = FALSE) +
  geom_segment(data = anno, aes(x = x3, xend = x4, y = y4, yend = y4), 
               colour = "black", inherit.aes = FALSE)
ChangeBC

pdf(file = "Figure2.pdf",
    width = 8.5,
    height = 4)

p1 + p2 + ChangeBC

dev.off()

KapMei + ChangeBC

mean(BenF$BodyCond) - mean(BenS$BodyCond) # -0.08086496 avg lost from starting BC
mean(LimF$BodyCond) - mean(LimS$BodyCond) # -0.08762443 avg lost from starting BC

(mean(BenF$BodyCond) - mean(BenS$BodyCond))/mean(BenS$BodyCond) #-0.06550651 avg proportion lost from starting BC

(mean(LimF$BodyCond) - mean(LimS$BodyCond))/mean(LimS$BodyCond) # -0.07341865 avg proportion lost from starting BC


# find 95% CI of avg proportion lost between start and end BC
mean(BenS$BodyCond) - 0.156498285 #1.077959
(1.077959 - mean(BenS$BodyCond))/mean(BenS$BodyCond) # -0.1267747

mean(BenS$BodyCond) - 0.005231633 #1.229225
(1.229225 - mean(BenS$BodyCond))/mean(BenS$BodyCond) # -0.004238225


mean(LimS$BodyCond) - 0.165960558 #1.027529
(1.027529 - mean(LimS$BodyCond))/mean(LimS$BodyCond) # -0.1390552

mean(LimS$BodyCond) - 0.009288293 #1.184202
(1.184202 - mean(LimS$BodyCond))/mean(LimS$BodyCond) # -0.007782182

mean(LimS$BodyCond)


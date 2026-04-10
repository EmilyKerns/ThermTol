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
              size = 1, alpha = 0.5) +
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

BodyCondPopsTemp <- ggplot(complete_data, aes(x=AS, y=BodyCond, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = expression(Experiment~Start~Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(Experiment~Start~Body~Condition))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCondPopsTemp

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
  annotate("text", x = -Inf, y = Inf, label = "A", 
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
  scale_color_manual( values=c('blue','black'))+
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
  labs(x = "Temperature (\u00B0C)", y = "Experiment End Body Condition")+
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
  labs(x = "Temperature (\u00B0C)", y = "Experiment End Body Condition")+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
BodyCondTempSex


#### Does temperature influence Fulton's K Body Condition (by population)? ------------------

## Only 2 populations were used to measure MMR, SMR, and AS at the end of the experiment, so repeating Condition ~ AS*Ecotype would be the exact same as Condition ~ AS*Population. However, condition metrics were collected on all populations, so we can assess the affect of temperature directly for each condition metric.

{
  # temp is not a factor, using all fish
  complete_data <- End2[complete.cases(End2[c("BodyCond", "Temp", "Pop", "Sex", "CCD")]), ]
  nrow(complete_data) # 173
  
  hist(complete_data$BodyCond)
  BodyCond <- lme(BodyCond ~ Temp*Pop + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(BodyCond)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -37.74783 -3.649303 29.87392
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.0579623 0.1750518
# 
# Fixed effects:  BodyCond ~ Temp * Pop + Sex 
#                   Value Std.Error  DF   t-value p-value
# (Intercept)  1.2672916 0.3151394 151  4.021368  0.0001
# Temp        -0.0025560 0.0144733 151 -0.176600  0.8601
# PopSL       -0.5943366 0.4574763  13 -1.299164  0.2165
# PopWK       -1.1480315 0.5683912 151 -2.019791  0.0452
# PopWT       -0.1092517 0.4107143  13 -0.266004  0.7944
# SexM         0.2245847 0.0275847 151  8.141633  0.0000
# Temp:PopSL   0.0195091 0.0209292 151  0.932144  0.3528
# Temp:PopWK   0.0428453 0.0267196 151  1.603518  0.1109
# Temp:PopWT  -0.0029648 0.0187418  13 -0.158189  0.8767
# Correlation: 
#   (Intr) Temp   PopSL  PopWK  PopWT  SexM   Tm:PSL Tm:PWK
# Temp       -0.990                                                 
# PopSL      -0.689  0.686                                          
# PopWK      -0.428  0.437  0.305                                   
# PopWT      -0.767  0.761  0.532  0.330                            
# SexM       -0.001 -0.032 -0.124 -0.082 -0.024                     
# Temp:PopSL  0.685 -0.695 -0.992 -0.310 -0.529  0.125              
# Temp:PopWK  0.426 -0.443 -0.302 -0.994 -0.329  0.072  0.314       
# Temp:PopWT  0.765 -0.772 -0.529 -0.336 -0.991  0.017  0.536  0.341
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.98413447 -0.59369922  0.05188157  0.63370583  2.36571293 
# 
# Number of Observations: 173
# Number of Groups: 17 

slopes <- emtrends(BodyCond, "Pop", var = "Temp")
# Pop Temp.trend     SE  df lower.CL upper.CL
# FG    -0.00256 0.0145 151 -0.03115   0.0260
# SL     0.01695 0.0151 151 -0.01279   0.0467
# WK     0.04029 0.0241 151 -0.00733   0.0879
# WT    -0.00552 0.0119  13 -0.03126   0.0202
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast estimate     SE  df t.ratio p.value
# FG - SL  -0.01951 0.0209 151  -0.932  0.7876
# FG - WK  -0.04285 0.0267 151  -1.604  0.3797
# FG - WT   0.00296 0.0187  13   0.158  0.9985
# SL - WK  -0.02334 0.0283 151  -0.825  0.8428
# SL - WT   0.02247 0.0192  13   1.170  0.6552
# WK - WT   0.04581 0.0269  13   1.703  0.3607
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 4 estimates 
# 
# P value adjustment: tukey method for comparing a family of 4 estimates 

plot(BodyCond)

BodyCond_TempsPopEnd <- ggplot(complete_data, aes(x=as.factor(Temp), y=BodyCond, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = "Temperature (\u00B0C)", 
       y = "Experiment End Body Condition")+
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
BodyCond_TempsPopEnd


### Hepatosomatic Index (HSI) -------------------------------------------------------


#### Does AS influence HSI (by ecotype)? ------------------


{
  complete_data <- End[complete.cases(End[c("HSI", "AS", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #68
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

HSIP <- ggplot(complete_data, aes(x=AS, y=HSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm"), legend.position = "none") +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")+
  annotate("text", x = 2200, y = 2.55, label = "AS: p<0.01\nAS*Ecotype: p<0.05", hjust = 1, vjust = 1, size = 2, fontface = "bold")
HSIP

HSISex <- ggplot(complete_data, aes(x=AS, y=HSI, color=Sex)) +
  geom_point(size = 0.5) +
  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISex



#### Does temperature directly influence HSI (by ecotype)? ------------------


{
  complete_data <- End2[complete.cases(End2[c("HSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #171
  
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

HSIPTemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Ecotype)) +
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
HSIPTemp

HSISexTemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISexTemp


#### Does temperature influence HSI (by population)? ------------------


{
  complete_data <- End2[complete.cases(End2[c("HSI", "Temp", "Pop", "CCD", "Sex")]), ]
  nrow(complete_data) #171
  
  hist(complete_data$HSI)
  HSI <- lme(HSI ~ Temp*Pop + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(HSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC     BIC    logLik
# 633.9324 667.896 -305.9662
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:  0.08601202 1.412512
# 
# Fixed effects:  HSI ~ Temp * Pop + Sex 
#                 Value Std.Error  DF    t-value p-value
# (Intercept)  3.814271  2.029539 149  1.8793780  0.0621
# Temp         0.080025  0.093252 149  0.8581580  0.3922
# PopSL       -3.223457  2.766063  13 -1.1653592  0.2648
# PopWK        0.521977  4.014431 149  0.1300253  0.8967
# PopWT       -0.519815  2.451060  13 -0.2120778  0.8353
# SexM        -0.449602  0.221937 149 -2.0258081  0.0446
# Temp:PopSL   0.194958  0.128935 149  1.5120581  0.1326
# Temp:PopWK  -0.083787  0.185552 149 -0.4515534  0.6522
# Temp:PopWT   0.027634  0.112226  13  0.2462326  0.8093
# Correlation: 
#   (Intr) Temp   PopSL  PopWK  PopWT  SexM   Tm:PSL Tm:PWK
# Temp       -0.990                                                 
# PopSL      -0.731  0.731                                          
# PopWK      -0.497  0.497  0.378                                   
# PopWT      -0.828  0.821  0.611  0.414                            
# SexM       -0.015 -0.028 -0.143 -0.090 -0.022                     
# Temp:PopSL  0.714 -0.727 -0.993 -0.370 -0.596  0.141              
# Temp:PopWK  0.490 -0.498 -0.372 -0.995 -0.409  0.082  0.370       
# Temp:PopWT  0.823 -0.831 -0.606 -0.412 -0.991  0.012  0.602  0.413
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -3.0049643 -0.6948291 -0.0492015  0.7101337  2.8667158 
# 
# Number of Observations: 171
# Number of Groups: 17 

slopes <- emtrends(HSI, "Pop", var = "Temp")
#  Pop Temp.trend     SE  df lower.CL upper.CL
# FG     0.08002 0.0933 149  -0.1042    0.264
# SL     0.27498 0.0886 149   0.0999    0.450
# WK    -0.00376 0.1610 149  -0.3216    0.314
# WT     0.10766 0.0625  13  -0.0273    0.243
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast estimate    SE  df t.ratio p.value
# FG - SL   -0.1950 0.129 149  -1.512  0.4329
# FG - WK    0.0838 0.186 149   0.452  0.9692
# FG - WT   -0.0276 0.112  13  -0.246  0.9945
# SL - WK    0.2787 0.183 149   1.527  0.4241
# SL - WT    0.1673 0.109  13   1.541  0.4430
# WK - WT   -0.1114 0.173  13  -0.645  0.9153
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 4 estimates 
plot(HSI)

HSI_TempsPopEnd <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = "Temperature (\u00B0C)", 
       y = "HSI")+
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSI_TempsPopEnd


HSISexTemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=HSI, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(HSI))+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
HSISexTemp


### Gonadosomatic Index (GSI - males only) -------------------------------------------------------


#### Does AS influence GSI (by ecotype)? ------------------

{
complete_data <- End[complete.cases(End[c("GSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data) #35

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
slopes <- emtrends(GSI, "Ecotype", var = "AS")
#  Ecotype   AS.trend       SE df  lower.CL upper.CL
# Benthic  -1.46e-05 1.02e-04 28 -2.23e-04 0.000194
# Limnetic  1.83e-04 7.11e-05 28  3.78e-05 0.000329
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast            estimate       SE df t.ratio p.value
# Benthic - Limnetic -0.000198 0.000124 28  -1.593  0.1223
# 
# Degrees-of-freedom method: containment 
plot(GSI)

GSIP <- ggplot(complete_data, aes(x=AS, y=GSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(GSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") 
GSIP


#### Does temperature directly influence GSI (by ecotype)? ------------------

{
  complete_data <- End2[complete.cases(End2[c("GSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #79
  
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
#                           Value Std.Error DF    t-value p-value
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
slopes <- emtrends(GSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE df lower.CL upper.CL
# Benthic     -0.0005 0.0143 60  -0.0292   0.0282
# Limnetic     0.0153 0.0195 60  -0.0237   0.0544
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate     SE df t.ratio p.value
# Benthic - Limnetic  -0.0158 0.0235 60  -0.676  0.5019
# 
# Degrees-of-freedom method: containment 
plot(GSI)

GSITemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=GSI, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(GSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold") 
GSITemp


#### Does temperature influence GSI (by population)? ------------------

{
  complete_data <- End2[complete.cases(End2[c("GSI", "Temp", "Pop", "CCD", "Sex")]), ]
  nrow(complete_data) #79
  
  hist(complete_data$GSI)
  GSI <- lme(GSI ~ Temp*Pop, random = ~ 1 | CCD, data = complete_data)
  summary(GSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -25.22716 -2.600356 22.61358
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev: 3.942687e-06  0.14322
# 
# Fixed effects:  GSI ~ Temp * Pop 
#                 Value Std.Error DF    t-value p-value
# (Intercept)  0.1185953 0.3619314 59  0.3276734  0.7443
# Temp         0.0133994 0.0164719 59  0.8134716  0.4192
# PopSL        0.6295068 0.4546168 12  1.3846979  0.1914
# PopWK       -0.7062464 0.6278783 59 -1.1248142  0.2652
# PopWT        0.1171760 0.4329074 12  0.2706721  0.7912
# Temp:PopSL  -0.0346260 0.0213428 59 -1.6223705  0.1101
# Temp:PopWK   0.0322775 0.0293071 59  1.1013554  0.2752
# Temp:PopWT  -0.0159617 0.0195064 12 -0.8182822  0.4291
# Correlation: 
#   (Intr) Temp   PopSL  PopWK  PopWT  Tm:PSL Tm:PWK
# Temp       -0.992                                          
# PopSL      -0.796  0.790                                   
# PopWK      -0.576  0.572  0.459                            
# PopWT      -0.836  0.829  0.666  0.482                     
# Temp:PopSL  0.766 -0.772 -0.992 -0.441 -0.640              
# Temp:PopWK  0.558 -0.562 -0.444 -0.995 -0.466  0.434       
# Temp:PopWT  0.838 -0.844 -0.667 -0.483 -0.993  0.652  0.475
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.8447488 -0.5288671 -0.2106220  0.4602073  2.6441268 
# 
# Number of Observations: 79
# Number of Groups: 16 

slopes <- emtrends(GSI, "Pop", var = "Temp")
#  Pop Temp.trend     SE df lower.CL upper.CL
# FG     0.01340 0.0165 59 -0.01956  0.04636
# SL    -0.02123 0.0136 59 -0.04838  0.00593
# WK     0.04568 0.0242 59 -0.00283  0.09418
# WT    -0.00256 0.0104 12 -0.02533  0.02020
# 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast estimate     SE df t.ratio p.value
# FG - SL    0.0346 0.0213 59   1.622  0.3742
# FG - WK   -0.0323 0.0293 59  -1.101  0.6901
# FG - WT    0.0160 0.0195 12   0.818  0.8448
# SL - WK   -0.0669 0.0278 59  -2.408  0.0866
# SL - WT   -0.0187 0.0171 12  -1.090  0.7021
# WK - WT    0.0482 0.0264 12   1.828  0.3080
# 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 4 estimates 
plot(GSI)

GSI_TempsPopEnd <- ggplot(complete_data, aes(x=as.factor(Temp), y=GSI, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = "Temperature (\u00B0C)", 
       y = "GSI")+
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
GSI_TempsPopEnd


### Spleenosomatic Index (SSI) -------------------------------------------------------


#### Does AS influence SSI (by ecotype)? ------------------


{
complete_data <- End[complete.cases(End[c("SSI", "AS", "Ecotype", "CCD", "Sex")]), ]
nrow(complete_data) #68

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
slopes <- emtrends(SSI, "Ecotype", var = "AS")
# Ecotype  AS.trend       SE df  lower.CL upper.CL
# Benthic  1.82e-05 2.62e-05 60 -3.43e-05 7.06e-05
# Limnetic 1.08e-05 2.12e-05 60 -3.16e-05 5.32e-05
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate       SE df t.ratio p.value
# Benthic - Limnetic 7.36e-06 3.35e-05 60   0.220  0.8270
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(SSI)

SSIP <- ggplot(complete_data, aes(x=AS, y=SSI, color=Ecotype)) +
  geom_point(size = 0.5) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(method = "lm", linewidth = 0.5, alpha = .2) +
  labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSIP


#### Does temperature directly influence SSI (by ecotype)? ------------------

{
  complete_data <- End2[complete.cases(End2[c("SSI", "Temp", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #166
  
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
slopes <- emtrends(SSI, "Ecotype", var = "Temp")
# Ecotype  Temp.trend      SE  df lower.CL upper.CL
# Benthic    -0.00148 0.00288 145 -0.00717  0.00421
# Limnetic   -0.00106 0.00418 145 -0.00933  0.00721
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast            estimate      SE  df t.ratio p.value
# Benthic - Limnetic -0.000421 0.00509 145  -0.083  0.9342
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(SSI)

SSIPTemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=SSI, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "D", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSIPTemp

BodyCond + HSI + GSI + SSI + FibPres + FibSev + plot_layout(ncol = 2)



#### Does temperature influence SSI (by population)? ------------------


{
  complete_data <- End2[complete.cases(End2[c("SSI", "Temp", "Pop", "CCD", "Sex")]), ]
  nrow(complete_data) #166
  
  hist(complete_data$SSI)
  SSI <- lme(SSI ~ Temp*Pop + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SSI)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC       BIC   logLik
# -314.7798 -281.1611 168.3899
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)   Residual
# StdDev: 0.001288735 0.07291942
# 
# Fixed effects:  SSI ~ Temp * Pop + Sex 
#                 Value  Std.Error  DF    t-value p-value
# (Intercept)  0.06259904 0.10367811 144  0.6037826  0.5469
# Temp         0.00252710 0.00476625 144  0.5302066  0.5968
# PopSL        0.15069725 0.13901831  13  1.0840101  0.2981
# PopWK       -0.03295722 0.20617053 144 -0.1598542  0.8732
# PopWT        0.14711403 0.12584579  13  1.1690024  0.2634
# SexM        -0.02475745 0.01166171 144 -2.1229687  0.0355
# Temp:PopSL  -0.00440381 0.00649284 144 -0.6782565  0.4987
# Temp:PopWK   0.00235221 0.00952280 144  0.2470079  0.8053
# Temp:PopWT  -0.00609516 0.00577293  13 -1.0558168  0.3103
# Correlation: 
#   (Intr) Temp   PopSL  PopWK  PopWT  SexM   Tm:PSL Tm:PWK
# Temp       -0.990                                                 
# PopSL      -0.743  0.742                                          
# PopWK      -0.501  0.500  0.387                                   
# PopWT      -0.824  0.816  0.614  0.413                            
# SexM       -0.016 -0.026 -0.135 -0.092  0.003                     
# Temp:PopSL  0.725 -0.737 -0.992 -0.377 -0.599  0.133              
# Temp:PopWK  0.494 -0.502 -0.380 -0.994 -0.408  0.083  0.378       
# Temp:PopWT  0.818 -0.825 -0.607 -0.409 -0.991 -0.017  0.603  0.411
# 
# Standardized Within-Group Residuals:
#   Min         Q1        Med         Q3        Max 
# -1.3322100 -0.5223266 -0.1805554  0.2760204  7.6159276 
# 
# Number of Observations: 166
# Number of Groups: 17 

slopes <- emtrends(SSI, "Pop", var = "Temp")
#  Pop Temp.trend      SE  df lower.CL upper.CL
# FG     0.00253 0.00477 144 -0.00689  0.01195
# SL    -0.00188 0.00439 144 -0.01055  0.00680
# WK     0.00488 0.00824 144 -0.01140  0.02116
# WT    -0.00357 0.00327  13 -0.01062  0.00349
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast estimate      SE  df t.ratio p.value
# FG - SL   0.00440 0.00649 144   0.678  0.9052
# FG - WK  -0.00235 0.00952 144  -0.247  0.9947
# FG - WT   0.00610 0.00577  13   1.056  0.7210
# SL - WK  -0.00676 0.00928 144  -0.728  0.8858
# SL - WT   0.00169 0.00550  13   0.308  0.9894
# WK - WT   0.00845 0.00888  13   0.952  0.7783
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# P value adjustment: tukey method for comparing a family of 4 estimates 
plot(SSI)

SSI_TempsPopEnd <- ggplot(complete_data, aes(x=as.factor(Temp), y=SSI, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('blue','black', "grey", "lightblue2"))+
  labs(x = "Temperature (\u00B0C)", 
       y = "SSI")+
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "E", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSI_TempsPopEnd


SSISexTemp <- ggplot(complete_data, aes(x=as.factor(Temp), y=SSI, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +  theme_classic() +
  scale_color_manual( values=c('black','blue'))+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  labs(x = "Temperature (\u00B0C)", y = expression(SSI))+
  annotate("text", x = -Inf, y = Inf, label = "E", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SSISexTemp


# Make Figures ------------------------------------------------------------


pdf(file = "NewFigure5.pdf",
    width = 8,
    height = 8)

BodyCondP + BodyCondEnd + HSIP + GSIP + SSIP + plot_layout(ncol = 2)

dev.off()

pdf(file = "NewSuppFig5.pdf",
    width = 6,
    height = 4)

BodyCondSexEnd + HSISex

dev.off()


# Scaled Mass Index -------------------------------------------------------



# based on Peig & Green, New perspectives for estimating body condition from mass/length data:
#     the scaled mass index as an alternative method"
#     Oikos 118: 1883-1891, 2009
# script from https://apansharing.blogspot.com/2018/05/an-r-function-olsrobust-caled-mass-index.html




## Format data -------------------------------------------------------------

{
data <- ThermTol2
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
mean(Ben_18$SL) # 62.01739
median(Ben_18$SL) # 61.565

Ben.dt.SMI <-
  scaledMassIndex(Ben_18$SL, Ben_18$Mass)
summary(Ben.dt.SMI)
#       SMI.ols         SMI.rob            x               y        
# Min.   :1.609   Min.   :1.553   Min.   :51.60   Min.   :1.432  
# 1st Qu.:2.428   1st Qu.:2.478   1st Qu.:59.61   1st Qu.:2.393  
# Median :2.767   Median :2.793   Median :61.56   Median :2.847  
# Mean   :2.830   Mean   :2.836   Mean   :62.02   Mean   :2.827  
# 3rd Qu.:3.185   3rd Qu.:3.235   3rd Qu.:64.24   3rd Qu.:3.349  
# Max.   :4.132   Max.   :4.181   Max.   :74.16   Max.   :4.343  


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
mean(Lim_18$SL) #61.60083
median(Lim_18$SL) # 61.575

Lim.dt.SMI <-
  scaledMassIndex(Lim_18$SL, Lim_18$Mass)
summary(Lim.dt.SMI)
#         SMI.ols         SMI.rob            x               y        
# Min.   :1.772   Min.   :1.766   Min.   :54.71   Min.   :1.449  
# 1st Qu.:2.526   1st Qu.:2.525   1st Qu.:57.48   1st Qu.:2.331  
# Median :2.831   Median :2.841   Median :61.58   Median :2.874  
# Mean   :2.776   Mean   :2.777   Mean   :61.60   Mean   :2.789  
# 3rd Qu.:3.048   3rd Qu.:3.052   3rd Qu.:64.38   3rd Qu.:3.125  
# Max.   :3.756   Max.   :3.772   Max.   :72.90   Max.   :4.084  


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
  ylab("Scaled Mass Index of Limnetic Control Fish") +
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



## calculate scaled mass index for all benthic fish --------------------
{
  mean(Ben_all$SL, na.rm = TRUE) # 62.59357
  median(Ben_all$SL, na.rm = TRUE) # 62.14
  
  Ben_all <- Ben_all[!is.na(Ben_all$Mass),]
  Ben_all <- Ben_all[!is.na(Ben_all$SL),]
  nrow(Ben_all) #205
  
  Ben.all.dt.SMI <-
    scaledMassIndex(Ben_all$SL, Ben_all$Mass, x.0 = 61)
  summary(Ben.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.458   Min.   :1.464   Min.   :51.09   Min.   :1.266  
  # 1st Qu.:2.321   1st Qu.:2.360   1st Qu.:59.93   1st Qu.:2.590  
  # Median :2.680   Median :2.695   Median :62.13   Median :2.931  
  # Mean   :2.753   Mean   :2.762   Mean   :62.50   Mean   :3.001  
  # 3rd Qu.:3.193   3rd Qu.:3.177   3rd Qu.:65.16   3rd Qu.:3.381  
  # Max.   :4.443   Max.   :4.332   Max.   :75.50   Max.   :5.521  
  

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
  mean(Lim_all$SL) # 61.54148
  median(Lim_all$SL) # 61.77
  
  Lim_all <- Lim_all[!is.na(Lim_all$Mass),]
  Lim_all <- Lim_all[!is.na(Lim_all$SL),]
  nrow(Lim_all) #181
  
  Lim.all.dt.SMI <-
    scaledMassIndex(Lim_all$SL, Lim_all$Mass)
  summary(Lim.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.680   Min.   :1.674   Min.   :49.49   Min.   :1.238  
  # 1st Qu.:2.379   1st Qu.:2.379   1st Qu.:59.00   1st Qu.:2.336  
  # Median :2.709   Median :2.705   Median :61.74   Median :2.716  
  # Mean   :2.717   Mean   :2.715   Mean   :61.54   Mean   :2.732  
  # 3rd Qu.:3.035   3rd Qu.:3.027   3rd Qu.:63.85   3rd Qu.:3.107  
  # Max.   :4.013   Max.   :3.987   Max.   :72.90   Max.   :4.970  
  
  
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



# SMI by population -------------------------------------------------------


## Format data -------------------------------------------------------------

{
  WK_18 <- filter(data, Pop == "WK" & Temp == "18")
  SL_18 <- filter(data, Pop == "SL" & Temp == "18")
  WT_18 <- filter(data, Pop == "WT" & Temp == "18")
  FG_18 <- filter(data, Pop == "FG" & Temp == "18")
  WK_all <- filter(data, Pop == "WK")
  SL_all <- filter(data, Pop == "SL")
  WT_all <- filter(data, Pop == "WT")
  FG_all <- filter(data, Pop == "FG")
}


## calculate scaled mass index for all WK fish --------------------
{
  mean(WK_all$SL) # 62.13075
  median(WK_all$SL) # 62.62
  
  WK.all.dt.SMI <-
    scaledMassIndex(WK_all$SL, WK_all$Mass, x.0 = 61)
  summary(WK.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.643   Min.   :1.571   Min.   :49.49   Min.   :1.238  
  # 1st Qu.:2.243   1st Qu.:2.176   1st Qu.:59.43   1st Qu.:2.240  
  # Median :2.385   Median :2.391   Median :62.62   Median :2.708  
  # Mean   :2.462   Mean   :2.455   Mean   :62.13   Mean   :2.693  
  # 3rd Qu.:2.730   3rd Qu.:2.726   3rd Qu.:64.45   3rd Qu.:3.170  
  # Max.   :3.478   Max.   :3.558   Max.   :72.90   Max.   :4.596  
  
  
  g1 <-
    ggplot(WK_all, aes(x = SL, y = Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(WK.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))+
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  

  g2 <-ggplot(WK_all, aes(x = log10(SL), y = log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(WK.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(WK.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI at mean body length of WK fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_WK_all <- ggarrange(g1,
                       g2,
                       g3,
                       ncol = 1,
                       nrow = 3,
                       align = "hv")
  SMI_WK_all
}


## calculate scaled mass index for all SL fish --------------------
{
  mean(SL_all$SL) # 61.29938
  median(SL_all$SL) # 61.31
  
  SL_all <- SL_all[!is.na(SL_all$Mass),]
  SL_all <- SL_all[!is.na(SL_all$SL),]
  nrow(SL_all) #128
  
  SL.all.dt.SMI <-
    scaledMassIndex(SL_all$SL, SL_all$Mass)
  summary(SL.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  #  Min.   :1.654   Min.   :1.643   Min.   :53.00   Min.   :1.432  
  # 1st Qu.:2.390   1st Qu.:2.403   1st Qu.:58.95   1st Qu.:2.392  
  # Median :2.725   Median :2.725   Median :61.27   Median :2.720  
  # Mean   :2.743   Mean   :2.740   Mean   :61.29   Mean   :2.747  
  # 3rd Qu.:3.075   3rd Qu.:3.067   3rd Qu.:63.31   3rd Qu.:3.072  
  # Max.   :3.948   Max.   :3.905   Max.   :71.15   Max.   :4.970  
  
  
  g1 <-
    ggplot(SL_all, aes(SL, Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(SL.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)")) +
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  
  g2 <-
    ggplot(SL_all, aes(log10(SL), log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(SL.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(SL.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI (at mean body length of SL fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_SL <- ggarrange(g1,
                       g2,
                       g3,
                       ncol = 1,
                       nrow = 3,
                       align = "hv")
  SMI_SL
}

SMI_WK_all + SMI_SL



## calculate scaled mass index for all WT fish --------------------
{
  mean(WT_all$SL, na.rm = TRUE) # 62.45688
  median(WT_all$SL, na.rm = TRUE) # 62.13
  
  WT_all <- WT_all[!is.na(WT_all$Mass),]
  WT_all <- WT_all[!is.na(WT_all$SL),]
  nrow(WT_all) #139
  
  WT.all.dt.SMI <-
    scaledMassIndex(WT_all$SL, WT_all$Mass)
  summary(WT.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.958   Min.   :1.890   Min.   :51.09   Min.   :1.266  
  # 1st Qu.:2.528   1st Qu.:2.540   1st Qu.:59.80   1st Qu.:2.554  
  # Median :2.809   Median :2.798   Median :62.11   Median :2.836  
  # Mean   :2.861   Mean   :2.857   Mean   :62.32   Mean   :2.871  
  # 3rd Qu.:3.146   3rd Qu.:3.121   3rd Qu.:64.44   3rd Qu.:3.163  
  # Max.   :4.130   Max.   :4.155   Max.   :75.50   Max.   :5.165   
  
  
  g1 <-
    ggplot(WT_all, aes(SL, Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(WT.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)")) +
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  
  g2 <-
    ggplot(WT_all, aes(log10(SL), log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(WT.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(WT.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI (at mean body length of WT fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_WT <- ggarrange(g1,
                      g2,
                      g3,
                      ncol = 1,
                      nrow = 3,
                      align = "hv")
  SMI_WT
}

SMI_WK_all + SMI_SL + SMI_WT





## calculate scaled mass index for all WT fish --------------------
{
  mean(FG_all$SL) # 62.88561
  median(FG_all$SL) # 62.265
  
  FG.all.dt.SMI <-
    scaledMassIndex(FG_all$SL, FG_all$Mass)
  summary(FG.all.dt.SMI)
  #       SMI.ols         SMI.rob            x               y        
  # Min.   :1.599   Min.   :1.628   Min.   :52.58   Min.   :1.505  
  # 1st Qu.:2.576   1st Qu.:2.649   1st Qu.:60.11   1st Qu.:2.771  
  # Median :3.370   Median :3.416   Median :62.27   Median :3.352  
  # Mean   :3.313   Mean   :3.301   Mean   :62.89   Mean   :3.273  
  # 3rd Qu.:4.098   3rd Qu.:4.028   3rd Qu.:65.75   3rd Qu.:3.708  
  # Max.   :5.231   Max.   :5.005   Max.   :74.16   Max.   :5.521  
  
  
  g1 <-
    ggplot(FG_all, aes(SL, Mass)) +
    geom_point() +
    geom_line(aes(x, y, color = Method),
              data = attr(FG.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)")) +
    geom_text(aes(label=FishID), size=3, hjust=0, vjust=0)
  g1
  
  g2 <-
    ggplot(FG_all, aes(log10(SL), log10(Mass))) +
    geom_point() +
    geom_line(
      aes(log10(x), log10(y), color = Method),
      data = attr(FG.all.dt.SMI, "pred") %>% melt(., "x", c("y.ols", "y.rob"), "Method", "y")
    ) +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g2
  
  g3 <-
    ggplot(reshape2::melt(FG.all.dt.SMI, "x", c("SMI.ols", "SMI.rob"), "Method"),
           aes(x, value)) +
    geom_point(aes(color = Method)) +
    xlab("Body length") +
    ylab("SMI (at mean body length of FG fish)") +
    scale_colour_discrete(labels = c("OLS", "Robust (M-estimation)"))
  g3
  
  SMI_FG <- ggarrange(g1,
                      g2,
                      g3,
                      ncol = 1,
                      nrow = 3,
                      align = "hv")
  SMI_FG
}

SMI_WK_all + SMI_SL + SMI_WT + SMI_FG



# The effect of AS & temp on SMI ---------------------------------


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

todrop.1 <- which(data$FishID %in% c("WK24S_7","SL26S_6","SL26S_9","SL26S_10","SL26S_11","FG26S_2","WT_24S_15","WT22S_2","WT24S_15","WK18S_4", "WT20S_10", "FG22S_6", "FG22S_1"))

data <- data[-todrop.1,]

table(data$Ecotype)
# Benthic Limnetic 
#  199      175 
is.na(data$SMI.rob)

SMI_end <- filter(data, Trial == "End")
table(SMI_end$Ecotype)
# Benthic Limnetic 
# 100       78 
table(SMI_end$Pop)
# FG SL WK WT 
# 31 57 21 69  

SMI_start <- filter(data, Trial == "Start")
table(SMI_start$Ecotype)
# Benthic Limnetic 
# 99       97 
table(SMI_start$Pop)
# FG SL WK WT 
# 32 67 30 67 

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
  SMI_end<- left_join(SMI_end, Density[ , c("FishID", "CCD")], by = c("FishID"))
  complete_data <- SMI_end[complete.cases(SMI_end[c("SMI.rob", "AS", "Ecotype", "CCD", "Sex")]), ]
  nrow(complete_data) #69

  hist(complete_data$SMI.rob)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  SMI.rob <- lme(SMI.rob ~ AS*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SMI.rob)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 91.29819 106.4104 -38.64909
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept) Residual
# StdDev:  0.09278112 0.317726
# 
# Fixed effects:  SMI.rob ~ AS * Ecotype + Sex 
#                       Value  Std.Error DF   t-value p-value
# (Intercept)        2.5189728 0.18785705 61 13.408987  0.0000
# AS                -0.0002159 0.00014872 61 -1.451981  0.1516
# EcotypeBenthic    -0.4372802 0.25032769  3 -1.746831  0.1790
# SexM               0.5475745 0.08034519 61  6.815274  0.0000
# AS:EcotypeBenthic  0.0004644 0.00023466 61  1.979159  0.0523
# Correlation: 
#   (Intr) AS     EctypB SexM  
# AS                -0.851                     
# EcotypeBenthic    -0.744  0.639              
# SexM              -0.218 -0.010  0.134       
# AS:EcotypeBenthic  0.580 -0.632 -0.881 -0.180
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.32357527 -0.62516838 -0.08451066  0.59830049  2.19903290 
# 
# Number of Observations: 69
# Number of Groups: 5 

slopes <- emtrends(SMI.rob, "Ecotype", var = "AS")
slopes
# Ecotype   AS.trend       SE df  lower.CL upper.CL
# Limnetic -0.000216 0.000149 61 -0.000513 8.14e-05
# Benthic   0.000248 0.000182 61 -0.000115 6.12e-04
#
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast            estimate       SE df t.ratio p.value
# Limnetic - Benthic -0.000464 0.000235 61  -1.979  0.0523
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(SMI.rob)

SMI_AS_end <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Ecotype)) +
  geom_point(size = .5) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression(Experiment~End~Aerobic~Scope~(mgO[2]/kg/hr)), y = "Experiment End Scaled Mass Index") +
  annotate("text", x = -Inf, y = Inf, label = "B", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  annotate("text", x = 2100, y = 3.5, label = "AS*Ecotype: p=0.0523", hjust = 1, vjust = -1, size = 2, fontface = "bold")
SMI_AS_end

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
  # Temp is numeric, including all fish
  complete_data <- SMI_end[complete.cases(SMI_end[c("SMI.rob", "Temp", "Ecotype", "CCD", "Sex")]), ]
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  
  SMI.rob <- lme(SMI.rob ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SMI.rob)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC   logLik
# 218.3479 240.2157 -102.174
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.2296614 0.3910287
# 
# Fixed effects:  SMI.rob ~ Temp * Ecotype + Sex 
#                           Value Std.Error  DF   t-value p-value
# (Intercept)          0.6317519 0.7940777 152  0.795579  0.4275
# Temp                 0.0763124 0.0361857 152  2.108909  0.0366
# EcotypeBenthic       2.3054698 0.9359716 152  2.463183  0.0149
# SexM                 0.5856956 0.0617089 152  9.491268  0.0000
# Temp:EcotypeBenthic -0.0972636 0.0436503 152 -2.228247  0.0273
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.992                     
# EcotypeBenthic      -0.757  0.764              
# SexM                -0.184  0.149  0.155       
# Temp:EcotypeBenthic  0.746 -0.764 -0.992 -0.150
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.63142629 -0.61385562  0.08125319  0.60384327  2.57979464 
# 
# Number of Observations: 173
# Number of Groups: 17 
slopes <- emtrends(SMI.rob, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE df lower.CL upper.CL
# Limnetic     0.0763 0.0362 152  0.00482    0.148
# Benthic     -0.0210 0.0283 152 -0.07686    0.035
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95
pairs(slopes)
# contrast           estimate     SE  df t.ratio p.value
# Limnetic - Benthic   0.0973 0.0437 152   2.228  0.0273
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(SMI.rob)


SMI_temp_end <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 1, alpha = 0.5) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression("Temperature (\u00B0C)"), y = "Experiment End Scaled Mass Index") +
  annotate("text", x = 5, y = 1, label = "Temp: p<0.05\nEcotype: p<0.05\nTemp*Ecotype: p<0.05", hjust = 1, vjust = 0, size = 2, fontface = "bold")+
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_temp_end

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
nrow(complete_data) #183

hist(complete_data$SMI.rob)
complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
SMI.rob <- glm(SMI.rob ~ AS*Ecotype, data = complete_data)
summary(SMI.rob)
}
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        2.5664237  0.1657590  15.483   <2e-16 ***
# AS                 0.0002100  0.0001419   1.480   0.1407    
# EcotypeBenthic     0.5485111  0.2599541   2.110   0.0362 *  
# AS:EcotypeBenthic -0.0005142  0.0002335  -2.202   0.0289 *  
# 
# (Dispersion parameter for gaussian family taken to be 0.2621471)
# 
# Null deviance: 48.204  on 182  degrees of freedom
# Residual deviance: 46.924  on 179  degrees of freedom
# AIC: 280.28

slopes <- emtrends(SMI.rob, "Ecotype", var = "AS")
# Ecotype   AS.trend       SE  df  lower.CL upper.CL
# Limnetic  0.000210 0.000142 179 -7.01e-05 4.90e-04
# Benthic  -0.000304 0.000185 179 -6.70e-04 6.17e-05
pairs(slopes)
# contrast           estimate       SE  df t.ratio p.value
# Limnetic - Benthic 0.000514 0.000233 179   2.202  0.0289

plot(SMI.rob)

SMI_AS_start <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Ecotype)) +
  geom_point(size = .5) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression(Experiment~Start~Aerobic~Scope~(mgO[2]/kg/hr)), y = "Experiment Start Scaled Mass Index") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold") +
  annotate("text", x = 1800, y = 4, label = "Ecotype: p<0.05\nAS*Ecotype: p<0.05", hjust = 1, vjust = -1, size = 2, fontface = "bold")
SMI_AS_start

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
# (Intercept)          3.50207    0.42847   8.173 3.97e-14 ***
# Temp                -0.03228    0.01948  -1.657   0.0991 .  
# EcotypeBenthic      -1.22035    0.58303  -2.093   0.0377 *  
# Temp:EcotypeBenthic  0.05507    0.02642   2.084   0.0385 *  
# 
# (Dispersion parameter for gaussian family taken to be 0.2600501)
# 
# Null deviance: 51.079  on 195  degrees of freedom
# Residual deviance: 49.930  on 192  degrees of freedom
# AIC: 298.19
slopes <- emtrends(SMI.rob, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE  df lower.CL upper.CL
# Limnetic    -0.0323 0.0195 192  -0.0707  0.00614
# Benthic      0.0228 0.0179 192  -0.0124  0.05800
pairs(slopes)
# contrast           estimate     SE  df t.ratio p.value
# Limnetic - Benthic  -0.0551 0.0264 192  -2.084  0.0385
plot(SMI.rob)

SMI_Temp <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 1, alpha = 0.5) +
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black')) +
  labs(x = expression("Temperature (\u00B0C)"), y = "Experiment Start Scaled Mass Index") +
  annotate("text", x = 5, y = 1.2, label = "Ecotype: p<0.05\nTemp*Ecotype: p<0.05", hjust = 1, vjust = -1, size = 2, fontface = "bold")+
  annotate("text", x = -Inf, y = Inf, label = "C", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_Temp

BodyCondP + SMI_AS + BodyCond_Temps + SMI_Temp


## The effect of AS & temp on SMI by Pop ---------------------------------


{
  colnames(WK.all.dt.SMI)[3] <- "SL"
  colnames(WK.all.dt.SMI)[4] <- "Mass"
  
  colnames(SL.all.dt.SMI)[3] <- "SL"
  colnames(SL.all.dt.SMI)[4] <- "Mass"
  
  colnames(WT.all.dt.SMI)[3] <- "SL"
  colnames(WT.all.dt.SMI)[4] <- "Mass"
  
  colnames(FG.all.dt.SMI)[3] <- "SL"
  colnames(FG.all.dt.SMI)[4] <- "Mass"
  
  WK_all$SMI.ols <- WK.all.dt.SMI$SMI.ols
  WK_all$SMI.rob <- WK.all.dt.SMI$SMI.rob
  SL_all$SMI.ols <- SL.all.dt.SMI$SMI.ols
  SL_all$SMI.rob <- SL.all.dt.SMI$SMI.rob
  WT_all$SMI.ols <- WT.all.dt.SMI$SMI.ols
  WT_all$SMI.rob <- WT.all.dt.SMI$SMI.rob
  FG_all$SMI.ols <- FG.all.dt.SMI$SMI.ols
  FG_all$SMI.rob <- FG.all.dt.SMI$SMI.rob
  
  data <- bind_rows(WK_all, SL_all, WT_all, FG_all)
  
  data <- left_join(data, ThermTol2[, c("CCD", "FishID")], by = "FishID")
  
  table(data$Pop)
  is.na(data$SMI.rob)
  
  SMI_end <- left_join(data, End2[, c("CCD", "FishID")], by = "FishID")
  SMI_end <- filter(SMI_end, Trial == "F")
  SMI_start <- filter(data, Trial == "S")
  
  data$Temp <- as.factor(data$Temp)
  
  SMI_AS <- ggplot(data, aes(x=AS, y=SMI.rob, color=Pop)) +
    geom_jitter(aes(color = Pop, group = Pop), 
                position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
                size = 3, alpha = 0.2) +
    geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
    theme_classic() +
    #theme(legend.position = "none") +
    scale_color_manual(values = c('black', 'blue', "grey", "lightblue")) +
    labs(x = expression(Aerobic~Scope~(mgO[2]/kg/hr)), y = "Scaled Mass Index") +
    annotate("text", x = -Inf, y = Inf, label = "A", 
             hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
  SMI_AS
}


## Experiment start --------------------------------------------------------
{
  complete_data <- SMI_start[complete.cases(SMI_start[c("SMI.rob", "AS", "Pop")]), ]
  nrow(complete_data)
  
  hist(complete_data$SMI.rob)
  SMI.rob <- glm(SMI.rob ~ AS*Pop, data = complete_data)
  summary(SMI.rob)
}
#                 Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)  5.1163804  0.4933076  10.372  < 2e-16 ***
#   AS          -0.0019956  0.0005009  -3.984 0.000100 ***
#   PopSL       -2.7494746  0.5415038  -5.077 9.92e-07 ***
#   PopWK       -2.7766731  0.5752119  -4.827 3.05e-06 ***
#   PopWT       -2.2468155  0.5406177  -4.156 5.11e-05 ***
#   AS:PopSL     0.0024199  0.0005389   4.491 1.30e-05 ***
#   AS:PopWK     0.0021334  0.0005510   3.872 0.000153 ***
#   AS:PopWT     0.0019764  0.0005402   3.659 0.000337 ***
# 
# (Dispersion parameter for gaussian family taken to be 0.2371074)
# 
# Null deviance: 51.827  on 178  degrees of freedom
# Residual deviance: 40.545  on 171  degrees of freedom
# AIC: 260.17
emtrends(SMI.rob, "Pop", var = "AS")
# Pop  AS.trend       SE  df  lower.CL  upper.CL
# FG  -2.00e-03 0.000501 171 -2.98e-03 -0.001007
# SL   4.24e-04 0.000199 171  3.23e-05  0.000816
# WK   1.38e-04 0.000229 171 -3.15e-04  0.000591
# WT  -1.92e-05 0.000202 171 -4.18e-04  0.000380
plot(SMI.rob)

SMI_AS <- ggplot(complete_data, aes(x=AS, y=SMI.rob, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 1, alpha = 0.5) +
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual(values = c('blue', 'black', "grey", "lightblue3")) +
  labs(x = expression(Experiment~Start~Aerobic~Scope~(mgO[2]/kg/hr)), y = "Experiment Start Scaled Mass Index") +
  annotate("text", x = -Inf, y = Inf, label = "A", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_AS

BodyCondP + SMI_AS

{
  complete_data <- SMI_start[complete.cases(SMI_start[c("SMI.rob", "Temp", "Pop")]), ]
  nrow(complete_data)
  
  hist(complete_data$SMI.rob)
  SMI.rob <- glm(SMI.rob ~ Temp*Pop, data = complete_data)
  summary(SMI.rob)
}
#               Estimate Std. Error t value Pr(>|t|)    
#   (Intercept) -1.21121    0.88696  -1.366 0.173865    
#   Temp         0.21078    0.04223   4.991 1.47e-06 ***
#   PopSL        4.70565    1.01685   4.628 7.28e-06 ***
#   PopWK        3.75897    1.16530   3.226 0.001506 ** 
#   PopWT        4.38132    1.00114   4.376 2.09e-05 ***
#   Temp:PopSL  -0.24194    0.04809  -5.031 1.23e-06 ***
#   Temp:PopWK  -0.21249    0.05385  -3.946 0.000116 ***
#   Temp:PopWT  -0.22522    0.04705  -4.787 3.64e-06 ***
# 
# (Dispersion parameter for gaussian family taken to be 0.2294179)
# 
# Null deviance: 51.827  on 178  degrees of freedom
# Residual deviance: 39.230  on 171  degrees of freedom
# AIC: 254.27
emtrends(SMI.rob, "Pop", var = "Temp")
# Pop Temp.trend     SE  df lower.CL upper.CL
# FG     0.21078 0.0422 171   0.1274   0.2941
# SL    -0.03116 0.0230 171  -0.0766   0.0142
# WK    -0.00172 0.0334 171  -0.0677   0.0642
# WT    -0.01445 0.0207 171  -0.0554   0.0265
plot(SMI.rob)

SMI_Temp_Pop <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Pop)) +
  geom_jitter(aes(color = Pop, group = Pop), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  geom_smooth(aes(group = Pop), method = "lm", linewidth = .5, alpha = 0.2) +
  theme_classic() +
  #theme(legend.position = "none") +
  scale_color_manual(values = c('blue', 'black', "grey", "lightblue2")) +
  labs(x = expression("Temperature (\u00B0C)"), y = "Experiment Start Scaled Mass Index") +
  annotate("text", x = -Inf, y = Inf, label = "D", 
           hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
SMI_Temp_Pop


# Make Figures with SMI ---------------------------------------------------

pdf(file = "NewSuppFig7.pdf",
    width = 8,
    height = 7)

SMI_AS_start + SMI_AS_end + SMI_Temp + SMI_temp_end

dev.off()

# Effect of Temperature Directly on Condition -----------------------------


{
  complete_data <- SMI_end[complete.cases(SMI_end[c("SMI.rob", "Ecotype", "Pop")]), ]
  nrow(complete_data)
  
  hist(complete_data$SMI.rob)
  shapiro.test(complete_data$SMI.rob)
  complete_data$Ecotype <- relevel(factor(complete_data$Ecotype), "Limnetic")
  SMI.rob <- lme(SMI.rob ~ Temp*Ecotype + Sex, random = ~ 1 | CCD, data = complete_data)
  summary(SMI.rob)
}
# Linear mixed-effects model fit by REML
# Data: complete_data 
# AIC      BIC    logLik
# 74.47807 89.48001 -30.23903
# 
# Random effects:
#   Formula: ~1 | CCD
# (Intercept)  Residual
# StdDev:   0.0693352 0.3310659
# 
# Fixed effects:  SMI.rob ~ Temp * Ecotype + Sex 
#                           Value Std.Error DF   t-value p-value
# (Intercept)          2.9709152 0.8199320 62  3.623368  0.0006
# Temp                -0.0371838 0.0405619  1 -0.916716  0.5276
# EcotypeBenthic      -0.0682592 0.9274850  1 -0.073596  0.9532
# SexM                 0.5556066 0.0821052 62  6.767006  0.0000
# Temp:EcotypeBenthic  0.0130177 0.0453555  1  0.287015  0.8221
# Correlation: 
#   (Intr) Temp   EctypB SexM  
# Temp                -0.994                     
# EcotypeBenthic      -0.872  0.870              
# SexM                -0.181  0.131  0.094       
# Temp:EcotypeBenthic  0.886 -0.892 -0.993 -0.102
# 
# Standardized Within-Group Residuals:
#   Min          Q1         Med          Q3         Max 
# -2.55350839 -0.59407624  0.03858452  0.59393676  1.98425941 
# 
# Number of Observations: 68
# Number of Groups: 5 
slopes <- emtrends(SMI.rob, "Ecotype", var = "Temp")
# Ecotype  Temp.trend     SE df lower.CL upper.CL
# Limnetic    -0.0372 0.0406  1   -0.553    0.478
# Benthic     -0.0242 0.0205  1   -0.284    0.236
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
# Confidence level used: 0.95 
pairs(slopes)
# contrast           estimate     SE df t.ratio p.value
# Limnetic - Benthic   -0.013 0.0454  1  -0.287  0.8221
# 
# Results are averaged over the levels of: Sex 
# Degrees-of-freedom method: containment 
plot(SMI.rob)

SMI.rob <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Ecotype)) +
  geom_jitter(aes(color = Ecotype, group = Ecotype), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = "Temperature (\u00B0C)", y = "Scaled Mass Index")+
  geom_smooth(aes(group = Ecotype), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
SMI.rob

SMI.robSex <- ggplot(complete_data, aes(x=as.factor(Temp), y=SMI.rob, color=Sex)) +
  geom_jitter(aes(color = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5), 
              size = 3, alpha = 0.2) +
  theme_classic() +
  theme(legend.position = "right", plot.margin = margin(1.1, .4, .4, .4, "cm")) +
  coord_cartesian(clip = 'off') +
  scale_color_manual( values=c('black','blue'))+
  labs(x = "Temperature (\u00B0C)", y = "Scaled Mass Index")+
  geom_smooth(aes(group = Sex), method = "lm", linewidth = .5, alpha = 0.2) +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = -1, size = 5, fontface = "bold")
SMI.robSex

BodyCond + SMI.rob

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

p1 + p2

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


# Testing OCLTTH Predictions 1-3 at the start & end of the experiment
# Figure 3, Supplemental Figures 1 & 2
#By Lake attribute

# Test the OCLTTH ---------------------------------------------------------


### Test the OCLTT hypothesis. Assumptions: 1) As temp increases, MMR will plateau and/or decrease. 2) As temp increases, RMR will increase. 3) At high temp  AS will decrease.


# Start of experiment ---------------------------------------


# Assumption 1 ------------------------------------------------------------
###################
# By lake attribute instead of ecotype. Mean depth and surface area most likely affect natural temp variation

## Mean depth --------------------------------------------------------------


#no effect of depth
MMRAvgD.i <- lm(MMR ~ Temp*MeanDepth, data = Start)
AIC(MMRAvgD.i) # 2690.805
summary(MMRAvgD.i)

#### Best mean depth model ####
MMRAvgD.a <- lm(MMR ~ Temp + MeanDepth, data = Start) 
AIC(MMRAvgD.a) # 2690.244
plot(MMRAvgD.a)
summary(MMRAvgD.a)
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   78.261    230.308   0.340   0.7344    
# Temp          78.783     10.087   7.810 4.55e-13 ***
# MeanDepth    -13.305      7.554  -1.761   0.0799 .  
# 
# Residual standard error: 371.6 on 180 degrees of freedom
# Multiple R-squared:  0.2597,	Adjusted R-squared:  0.2515 
# F-statistic: 31.57 on 2 and 180 DF,  p-value: 1.77e-12

MMRAvgDQ.i.q <- lm(MMR ~ I(Temp^2)*MeanDepth + Temp*MeanDepth, data = Start)
AIC(MMRAvgDQ.i.q) # 2694.479
summary(MMRAvgDQ.i.q)

MMRAvgDQ.q.a <- lm(MMR ~ I(Temp^2) + Temp + MeanDepth, data = Start)
AIC(MMRAvgDQ.q.a) # 2692.233
summary(MMRAvgDQ.q.a)

## Surface Area --------------------------------------------------------------


#### Best SA model ####
MMRSA.i <- lm(MMR ~ Temp*SA, data = Start)
AIC(MMRSA.i) # 2690.669
plot(MMRSA.i)
summary(MMRSA.i)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1444.027    906.019  -1.594 0.112742    
# Temp          157.224     41.144   3.821 0.000183 ***
# SA            170.203     88.861   1.915 0.057037 .  
# Temp:SA        -8.267      4.086  -2.023 0.044547 *  
# 
# Residual standard error: 371.1 on 179 degrees of freedom
# Multiple R-squared:  0.266,	Adjusted R-squared:  0.2537 
# F-statistic: 21.63 on 3 and 179 DF,  p-value: 5.343e-12

MMRSA.a <- lm(MMR ~ Temp+SA, data = Start)
AIC(MMRSA.a) # 2692.806
summary(MMRSA.i)

#evidence of quadratic trends, but AIC not significantly improved over linear model
MMRSA.i.q <- lm(MMR ~ I(Temp^2)*SA + Temp*SA, data = Start)
AIC(MMRSA.i.q) # 2689.817
summary(MMRSA.i.q)


MMRSA.a.q <- lm(MMR ~ I(Temp^2)+ Temp + SA, data = Start)
summary(MMRSA.a.q)
# AIC: 2694.8

MMRSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=MMR, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)), color = "Population (Surface Area)")+
  scale_color_manual(labels = c("Finger (13.8)", "Spirit (12.2)", "Wik (8.0)", "Watson (7.5)"), values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = 1)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = 1))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), 
              method = "lm", 
              formula = y ~ x, 
              size = 1,
              se = TRUE, 
              alpha = .1)   +
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
MMRSApoint

# Assumption 2 ------------------------------------------------------------
Start$RMR <- as.numeric(Start$RMR)


## Mean depth --------------------------------------------------------------



RMRAvgD.i <- lm(RMR ~ Temp*MeanDepth, data = Start)
AIC(RMRAvgD.i) # 2448.416
summary(RMRAvgD.i)

RMRAvgD.a <- lm(RMR ~ Temp*MeanDepth, data = Start)
AIC(RMRAvgD.a) # 2448.416
summary(RMRAvgD.a)

#### Best mean depth model ####
RMRAvgD.q.i <- lm(RMR ~ I(Temp^2)*MeanDepth + Temp*MeanDepth, data = Start)
summary(RMRAvgD.q.i)
plot(RMRAvgD.q.i)
AIC(RMRAvgD.q.i) #2438.554
#                       Estimate Std. Error t value Pr(>|t|)   
# (Intercept)         5508.9948  2107.0577   2.615  0.00970 **
# I(Temp^2)             12.5867     4.4648   2.819  0.00536 **
# MeanDepth            426.1368   264.8567   1.609  0.10941   
# Temp                -495.3290   195.3893  -2.535  0.01211 * 
# I(Temp^2):MeanDepth    0.7394     0.5558   1.330  0.18512   
# MeanDepth:Temp       -36.1191    24.4262  -1.479  0.14100   
# 
# Residual standard error: 185.3 on 177 degrees of freedom
# Multiple R-squared:  0.6089,	Adjusted R-squared:  0.5979 
# F-statistic: 55.12 on 5 and 177 DF,  p-value: < 2.2e-16

#additive quadratic. No effect of Mean Depth.
RMRAvg.q.a <- lm(RMR ~ I(Temp^2) + Temp + MeanDepth, data = Start)
summary(RMRAvg.q.a)
AIC(RMRAvg.q.a) #2443.845


## Surface Area --------------------------------------------------------------

#No interaction
RMRSA.i <- lm(RMR ~ Temp*SA, data = Start)
summary(RMRSA.i)
AIC(RMRSA.i) #2454.943

RMRSA.a <- lm(RMR ~ Temp+SA, data = Start)
summary(RMRSA.a)
AIC(RMRSA.a) # 2453.101

#Significant quadratic interactions #slightly better than mean depth#both better predictors than ecotype
#driven by trend in FG, NO 26C values for FG
RMRSA.q.i <- lm(RMR ~ I(Temp^2)*SA + Temp*SA, data = Start)
summary(RMRSA.q.i)
AIC(RMRSA.q.i) #2440.675

#### Best SA model ####
#Additive ecotype, quad temp
RMRSA.q.a <- lm(RMR ~ I(Temp^2) + Temp + SA, data = Start)
summary(RMRSA.q.a)
AIC(RMRSA.q.a) #2441.215
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2667.636    997.708   2.674 0.008195 ** 
# I(Temp^2)      8.040      2.140   3.756 0.000233 ***
# Temp        -268.361     93.253  -2.878 0.004492 ** 
# SA            11.937      5.581   2.139 0.033813 *  
# 
# Residual standard error: 187.7 on 179 degrees of freedom
# Multiple R-squared:  0.5944,	Adjusted R-squared:  0.5877 
# F-statistic: 87.46 on 3 and 179 DF,  p-value: < 2.2e-16

RMRSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=RMR, color=Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  #theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)), color = "Population (Depth)") +
         scale_color_manual(labels = c("Finger (-7.3)", "Spirit (-9.6)", "Wik (-12.2)", "Watson (-2.6)"), values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
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


## Mean depth --------------------------------------------------------------

#### Best mean depth model ####
ASAvgD.i <- lm(AS ~ Temp*MeanDepth, data = Start)
AIC(ASAvgD.i)
plot(ASAvgD.i)
summary(ASAvgD.i)
# AIC: 2649.287
#                 Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     -82.453    426.278  -0.193  0.84685   
# Temp             48.892     19.058   2.565  0.01112 * 
# MeanDepth      -168.217     53.654  -3.135  0.00201 **
# Temp:MeanDepth    7.163      2.402   2.982  0.00326 **
# 
# Residual standard error: 331.4 on 179 degrees of freedom
# Multiple R-squared:  0.05746,	Adjusted R-squared:  0.04166 
# F-statistic: 3.637 on 3 and 179 DF,  p-value: 0.01396

ASAvgD.a <- lm(AS ~ Temp+MeanDepth, data = Start)
AIC(ASAvgD.a) #2656.159
summary(ASAvgD.a)

ASAvgD.q.i <- lm(AS ~ I(Temp^2)*MeanDepth + Temp*MeanDepth, data = Start)
summary(ASAvgD.q.i)
AIC(ASAvgD.q.i)
# AIC: 2650.109

ASAvgD.q.a <- lm(AS ~ I(Temp^2) + Temp + MeanDepth, data = Start)
summary(ASAvgD.q.a)
AIC(ASAvgD.q.a) # 2654.492


## Surface Area --------------------------------------------------------------

ASSA.i <- lm(AS ~ Temp*SA, data = Start)
summary(ASSA.i)
AIC(ASSA.i) # 2650.898

ASSA.a <- lm(AS ~ Temp+SA, data = Start)
summary(ASSA.a)
AIC(ASSA.a) # 2655.106

#### Best SA model ####
ASSA.q.i <- lm(AS ~ I(Temp^2)*SA + Temp*SA, data = Start)
summary(ASSA.q.i)
plot(ASSA.q.i)
AIC(ASSA.q.i) # 2634.064
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -29154.358   6941.585  -4.200 4.22e-05 ***
# I(Temp^2)       -60.905     14.869  -4.096 6.39e-05 ***
# SA             2594.858    693.715   3.741 0.000248 ***
# Temp           2744.643    646.506   4.245 3.52e-05 ***
# I(Temp^2):SA      5.167      1.503   3.439 0.000729 ***
# SA:Temp        -234.376     64.948  -3.609 0.000400 ***
# 
# Residual standard error: 316.2 on 177 degrees of freedom
# Multiple R-squared:  0.1514,	Adjusted R-squared:  0.1275 
# F-statistic: 6.318 on 5 and 177 DF,  p-value: 2.017e-05

ASSA.q.a <- lm(AS ~ I(Temp^2) + Temp + SA, data = Start)
summary(ASSA.q.a)
AIC(ASSA.q.a) # 2651.551

LakeData2
ASSApoint <- ggplot(Start, aes(x=as.factor(Temp), y=AS, color = Pop)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)), color = "Population (Surface Area)")+
  scale_color_manual(labels = c("Finger (13.8)", "Spirit (12.2)", "Wik (8.0)", "Watson (7.5)"), values=c('darkblue','lightskyblue3','darkgrey','dodgerblue3'))+
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

MMRSApoint + RMRSApoint + ASSApoint

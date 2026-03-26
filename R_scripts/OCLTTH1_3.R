# Testing OCLTTH Predictions 1-3 at the start & end of the experiment
# Figure 3, Supplemental Figures 1 & 2

# Test the OCLTTH: Assumptions & Model Selection Criteria ---------------------------------------------------------

# OCLTTH Assumptions: 
#   1) As temp increases, MMR will plateau and/or decrease. 
#   2) As temp increases, RMR (aka SMR) will increase. 
#   3) At high temp AS will decrease.

# Model selection: 
#   Lowest AIC if delta AIC >2. 
#   Simplest model if delta AIC is <= 2


# Start of experiment ---------------------------------------


# Assumption 1 ------------------------------------------------------------

#First testing for Ecotype*Temp interaction
MMRL.i <- lm(MMR ~ Temp*Ecotype , data = Start)
plot(MMRL.i) #residuals look reasonably even
AIC(MMRL.i)# 2692.5
summary(MMRL.i)
# Residuals:
# Min      1Q  Median      3Q     Max 
# -1790.4  -236.5    -6.9   241.0   871.1 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             3.691    310.569   0.012    0.991    
# Temp                   84.650     14.132   5.990 1.13e-08 ***
#   EcotypeLimnetic       374.591    444.212   0.843    0.400    
# Temp:EcotypeLimnetic  -13.273     20.236  -0.656    0.513    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 372.9 on 179 degrees of freedom
# Multiple R-squared:  0.2586,	Adjusted R-squared:  0.2462 
# F-statistic: 20.82 on 3 and 179 DF,  p-value: 1.291e-11

confint(MMRL.i, level = 0.95)
#                     2.5 %    97.5 %
# (Intercept)          -609.15540  616.53814
# Temp                   56.76301  112.53662
# EcotypeLimnetic      -501.97441 1251.15701
# Temp:EcotypeLimnetic  -53.20477   26.65875

#Test quadratic fit
MMRQ.quad <- lm(MMR ~ I(Temp^2)*Ecotype + Temp * Ecotype, data = Start)
plot(MMRQ.quad)
AIC(MMRQ.quad) #2696.322
summary(MMRQ.quad)


#n.s.--ecotype only
MMRQ.quad.a <- lm(MMR ~ I(Temp^2) + Temp + Ecotype, data = Start)
plot(MMRQ.quad.a)
AIC(MMRQ.quad.a) #2692.939
summary(MMRQ.quad.a)



### No evidence of quadratic effect (interaction or additive), but Limnetic fish clearly appear to decline in MMR at highest temperature (26C). 
MMRpoint <- ggplot(Start, aes(x=as.factor(Temp), y=MMR, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = .4), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = expression("Temperature"*degree*C), y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = .4)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = .4))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), 
              method = "lm", 
              formula = y ~ x, 
              size = .5,
              se = TRUE, 
              alpha = .2)   +
  annotate("text", x =0.5, y = Inf, label = "A", hjust = 0, vjust = 1, size = 5, fontface = "bold")
#annotate("text", x =3.5, y = 200, label = "Temperature: p<0.0001", hjust = 0, vjust = -1, size = 2, fontface = "bold")
MMRpoint


#Testing whether adding high temp to model changes the fit.

#Test linear + high temp*Ecotype. Marginally significant Limnetic*high temp - limnetic fish showed decreased MMR at 26C at the start of the experiment
Start.high <- Start %>%
  mutate(high_temp = ifelse(Temp == max(Temp), 1, 0))

test_high <- lm(MMR ~ Temp + Ecotype * high_temp, data = Start.high)
plot(test_high)
AIC(test_high) # 2692.565
summary(test_high)

# Call:
#   lm(formula = MMR ~ Temp + Ecotype * high_temp, data = Start.high)
# 
# Residuals:
# Min       1Q   Median       3Q      Max 
# -1812.26  -216.21   -15.76   244.27   867.17 
# 
# Coefficients:
#                              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  84.44     277.03   0.305    0.761    
# Temp                         80.32      12.99   6.183 4.19e-09 ***
# EcotypeLimnetic             118.90      59.44   2.001    0.047 *  
# high_temp                    92.03     128.83   0.714    0.476    
# EcotypeLimnetic:high_temp  -236.47     157.55  -1.501    0.135    
# 
# Residual standard error: 372 on 178 degrees of freedom
# Multiple R-squared:  0.2664,	Adjusted R-squared:   0.25 
# F-statistic: 16.16 on 4 and 178 DF,  p-value: 2.611e-11



# Assumption 2 ------------------------------------------------------------
#As temp increases, RMR will increase.

Start$RMR <- as.numeric(Start$RMR)

#Significant interaction, but not as good a fit as just temp quadratic
RMRL <- lm(RMR~Temp*Ecotype, data = Start)
plot(RMRL)
AIC(RMRL) # 2450.772
summary(RMRL)


#no ecotype*quadratic interaction
RMRQ.i <- lm(RMR~ I(Temp^2)*Ecotype + Temp*Ecotype, data = Start)
plot(RMRQ.i)
AIC(RMRQ.i)#2443.451
summary(RMRQ.i)

### Best model ####
#quadratic additive--no evidence for ecotype, but strong quadratic effect of temp
RMRQ.a <- lm(RMR~ I(Temp^2) + Temp+ Ecotype, data = Start)
plot(RMRQ.a) 
AIC(RMRQ.a)#2445.096 - within 2 AIC of last model
summary(RMRQ.a)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -594.39 -117.74    1.53  123.70  531.02 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     2567.493   1007.488   2.548  0.01166 *  
# I(Temp^2)          7.494      2.145   3.494  0.00060 ***
# Temp            -246.710     93.596  -2.636  0.00913 ** 
# EcotypeLimnetic   23.910     28.137   0.850  0.39658    
# 
# Residual standard error: 189.7 on 179 degrees of freedom
# Multiple R-squared:  0.5858,	Adjusted R-squared:  0.5788 
# F-statistic: 84.37 on 3 and 179 DF,  p-value: < 2.2e-16



RMRpoint <- ggplot(Start, aes(x=as.factor(Temp), y=RMR, color=Ecotype)) +
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
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
# annotate("text", x = 5.5, y = 70, 
#          label = expression(bold(Temperature^{2}~": p<0.001")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 5.5, y = 5, 
#          label = expression(bold("Temperature: p<0.05")), 
#          hjust = 1, vjust = 0, size = 2)
RMRpoint



# Assumption 3 ------------------------------------------------------------
#At high temp AS will decrease.

hist(Start$AS)

#significant interaction
ASL <- lm(AS ~ Temp*Ecotype, data = Start)
plot(ASL)
AIC(ASL) # 2653.389
summary(ASL)
# Call:
#   lm(formula = AS ~ Temp * Ecotype, data = Start)
# 
# Residuals:
#    Min       1Q   Median       3Q      Max 
# -1140.94  -217.79    -4.47   205.58   865.35 
# 
# Coefficients:
#                       Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            653.27     279.09   2.341   0.0203 *
# Temp                    17.70      12.70   1.394   0.1650  
# EcotypeLimnetic        932.27     399.19   2.335   0.0206 *
# Temp:EcotypeLimnetic   -39.65      18.19  -2.180   0.0305 *
# 
# Residual standard error: 335.1 on 179 degrees of freedom
# Multiple R-squared:  0.03609,	Adjusted R-squared:  0.01994 
# F-statistic: 2.234 on 3 and 179 DF,  p-value: 0.08587
emmeans(ASL, pairwise ~ Ecotype, var = "Temp")
# Ecotype  emmean   SE  df lower.CL upper.CL
# Benthic    1039 35.7 179      968     1109
# Limnetic   1107 34.4 179     1040     1175
# $contrasts
# contrast           estimate   SE  df t.ratio p.value
# Benthic - Limnetic    -68.6 49.6 179  -1.383  0.1683

#quadratic interaction
ASQ <- lm(AS ~ I(Temp^2)*Ecotype +Temp*Ecotype, data = Start)
plot(ASQ)
AIC(ASQ)#2653.723
summary(ASQ)


#quadratic additive
ASQ.a <- lm(AS ~ I(Temp^2) + Temp + Ecotype, data = Start)
plot(ASQ.a)
AIC(ASQ.a)#2654.31
summary(ASQ.a)


ASpoint <- ggplot(Start, aes(x=as.factor(Temp), y=AS, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
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
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
# annotate("text", x = 5.5, y = 70, 
#          label = expression(bold(Temperature^{2}~": p<0.05")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 5.5, y = 5, 
#          label = expression(bold("Temperature: p<0.05")), 
#          hjust = 1, vjust = 0, size = 2)
ASpoint


ASpointL <- ggplot(Start, aes(x=as.factor(Temp), y=AS, color=Ecotype)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = .4), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))+
  scale_color_manual( values=c('black','blue'))+
  stat_summary(fun = mean, na.rm = TRUE, 
               geom = "point", shape = "circle",
               size = 1, 
               position = position_dodge(width = .4)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE, 
               geom = "errorbar", width = .3,
               position = position_dodge(width = .4))+
  stat_smooth(aes(x = as.numeric(factor(Temp))), 
              method = "lm", 
              formula = y ~ x, 
              size = .5,
              se = TRUE, 
              alpha = .2)   +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1, size = 5, fontface = "bold")
# annotate("text", x = 5.5, y = 70, 
#          label = expression(bold(Ecotype~": p<0.05")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 5.5, y = 5, 
#          label = expression(bold("Temperature*Ecotype: p<0.05")), 
#          hjust = 1, vjust = 0, size = 2)
ASpointL



# End of Experiment,secondary metabolic rates -------------------------------------------------------


# Assumption 1 ------------------------------------------------------------

# test linear interaction
MMREndL <- lm(MMR ~ Temp * Ecotype + Sex, data = End)
plot(MMREndL)
AIC(MMREndL) #1049.426
summary(MMREndL)

# test quadratic interaction 
MMREndQ.i <- lm(MMR ~ I(Temp^2)* Ecotype + Temp * Ecotype + Sex, data = End)
plot(MMREndQ.i)
AIC(MMREndQ.i) # 1042.542
summary(MMREndQ.i)

# test quadratic additive 
MMREndQ.a<- lm(MMR ~ I(Temp^2) + Temp + Ecotype + Sex, data = End)
plot(MMREndQ.a)
AIC(MMREndQ.a) #1041.022
summary(MMREndQ.a)

###### Best model ######
#No sex effect, test quadratic additive 
MMREndQ.a.noSex<- lm(MMR ~ I(Temp^2) + Temp + Ecotype, data = End)
plot(MMREndQ.a.noSex)
AIC(MMREndQ.a.noSex) #1039.209
summary(MMREndQ.a.noSex)
# Call:
#   lm(formula = MMR ~ I(Temp^2) + Temp + Ecotype, data = End)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1018.30  -284.38     1.26   296.23  1472.45 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     14998.928   3252.320   4.612 1.92e-05 ***
# I(Temp^2)          29.854      7.145   4.179 8.94e-05 ***
# Temp            -1283.462    307.558  -4.173 9.11e-05 ***
# EcotypeLimnetic   288.780    113.649   2.541   0.0135 *  
# 
# Residual standard error: 432.2 on 65 degrees of freedom
# Multiple R-squared:  0.241,	Adjusted R-squared:  0.206 
# F-statistic:  6.88 on 3 and 65 DF,  p-value: 0.0004276

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
# annotate("text", x = 3.5, y = 250, 
#          label = expression(bold(Temperature^{2}~": p<0.0001")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 3.5, y = 130, 
#          label = expression(bold("Temperature: p<0.0001")), 
#          hjust = 1, vjust = 0, size = 2)+
# annotate("text", x = 3.5, y = 5, 
#          label = expression(bold("Ecotype: p<0.05")), 
#          hjust = 1, vjust = 0, size = 2)
MMREnd

MMRSex <- ggplot(End, aes(x=factor(Temp), y=MMR, color=Sex)) +
  geom_point(position =  position_jitterdodge(jitter.width = 0.0001, jitter.height = 0.0001, dodge.width = 1), alpha= 0.1) +
  theme_classic() +
  theme(legend.position = "right") +
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



# Assess the effects of fibrosis on MMR ----------------------------

###### Best model ######
# Add fibrosis to the best model for MMR
MMREndQ.fib <- lm(MMR ~ I(Temp^2) +Temp + Ecotype + Sex + Fibrosis, data = End)
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
#No effect of temp

# Did fibrosis rates vary by ecotype?
Fib.ecotype <- glm(Fibrosis ~ Ecotype,family = "binomial", data = Fib.data)
summary(Fib.ecotype)
#                 Estimate Std. Error z value Pr(>|z|)   
# (Intercept)      -0.5878     0.2108  -2.788   0.0053 **
# EcotypeLimnetic   0.5351     0.3116   1.717   0.0859 . 
# 
# Null deviance: 236.02  on 173  degrees of freedom
# Residual deviance: 233.05  on 172  degrees of freedom
# AIC: 237.05


emmeans(Fib.ecotype,~Ecotype)
# Ecotype   emmean    SE  df asymp.LCL asymp.UCL
# Benthic  -0.5878 0.211 Inf    -1.001    -0.175
# Limnetic -0.0526 0.229 Inf    -0.502     0.397
# 
# Results are given on the logit (not the response) scale. 
# Confidence level used: 0.95

exp(-0.5878) / (1 + exp(-0.5878)) #0.3571398 = Benthic fibrosis rate
exp(-0.0526) / (1 + exp(-0.0526)) #0.486853 = Limnetic fibrosis rate

effectsize::standardize_parameters(Fib.ecotype, exp = TRUE)
# # Standardization method: refit
# 
# Parameter          | Std_Odds_Ratio |       95% CI
# --------------------------------------------------
# (Intercept)        |           0.56 | [0.36, 0.83]
# Ecotype [Limnetic] |           1.71 | [0.93, 3.16]
#
# - Response is unstandardized.



# Assumption 2 ------------------------------------------------------------

End$RMR <- as.numeric(End$RMR)
hist(End$RMR)

# test linear interaction
RMREndL <- lm(RMR~Temp*Ecotype + Sex, data = End)
plot(RMREndL)
AIC(RMREndL) #910.7636
summary(RMREndL)


# test quadratic interaction
RMREndQ.i <- lm(RMR ~ I(Temp^2)*Ecotype +Temp*Ecotype + Sex, data = End)
plot(RMREndQ.i)
AIC(RMREndQ.i) # 902.283
summary(RMREndQ.i)

###### Best model ######
# test quadratic additive - sig effect of quadratic and linear temp but not ecotype or sex
RMREndQ.a <- lm(RMR ~ I(Temp^2) + Temp + Ecotype + Sex, data = End)
plot(RMREndQ.a)
AIC(RMREndQ.a) #901.9981
summary(RMREndQ.a)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -397.03  -91.15    5.62   83.35  473.54 
# 
# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     5797.192   1195.327   4.850 8.25e-06 ***
# I(Temp^2)         12.859      2.626   4.896 6.95e-06 ***
# Temp            -520.241    113.039  -4.602 2.03e-05 ***
# EcotypeLimnetic   51.954     41.846   1.242    0.219    
# SexM             -62.240     38.723  -1.607    0.113    
# 
# Residual standard error: 158.8 on 64 degrees of freedom
# Multiple R-squared:  0.4464,	Adjusted R-squared:  0.4117 
# F-statistic:  12.9 on 4 and 64 DF,  p-value: 9.283e-08

# test quadratic additive without sex since it was n.s.
RMREndQ.nosex <- lm(RMR ~ I(Temp^2) + Temp + Ecotype, data = End)
plot(RMREndQ.nosex)
AIC(RMREndQ.nosex) #902.7287
summary(RMREndQ.nosex)


RMREnd <- ggplot(End, aes(x=as.factor(Temp), y=RMR, color=Ecotype)) +
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
# annotate("text", x = 3.5, y = 50, 
#          label = expression(bold(Temperature^{2}~": p<0.0001")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 3.5, y = 5, 
#          label = expression(bold("Temperature: p<0.0001")), 
#          hjust = 1, vjust = 0, size = 2)
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

# test linear interaction
ASEndL <- lm(AS ~ Temp*Ecotype + Sex, data = End)
plot(ASEndL)
AIC(ASEndL) #1007.751
summary(ASEndL)


# test quadratic interaction
ASEndQ.i <- lm(AS ~ I(Temp^2) * Ecotype +Temp * Ecotype + Sex, data = End)
plot(ASEndQ.i)
AIC(ASEndQ.i) #1003.939
summary(ASEndQ.i)


# test quadratic additive--model with sex, but sex is n.s.
ASEndQ.a <- lm(AS ~ I(Temp^2) + Temp+  Ecotype + Sex, data = End)
plot(ASEndQ.a)
AIC(ASEndQ.a) #1002.022
summary(ASEndQ.a)


###### Best model ######
# test quadratic additive, no sex.
ASEndQ.a.noSex <- lm(AS ~ I(Temp^2) + Temp+  Ecotype, data = End)
plot(ASEndQ.a.noSex)
AIC(ASEndQ.a.noSex) #1001.906
summary(ASEndQ.a.noSex)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -597.00 -209.61   38.07  180.63 1151.97 
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     9224.817   2481.992   3.717 0.000422 ***
# I(Temp^2)         16.907      5.452   3.101 0.002854 ** 
# Temp            -760.723    234.711  -3.241 0.001879 ** 
# EcotypeLimnetic  232.677     86.731   2.683 0.009249 ** 
# 
# Residual standard error: 329.9 on 65 degrees of freedom
# Multiple R-squared:  0.2676,	Adjusted R-squared:  0.2338 
# F-statistic: 7.915 on 3 and 65 DF,  p-value: 0.0001409

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
# annotate("text", x = 3.5, y = 2400, 
#          label = expression(bold(Temperature^{2}~": p<0.01")), 
#          hjust = 1, vjust = 0, size = 2) +
# annotate("text", x = 3.5, y = 2300, 
#          label = expression(bold("Temperature: p<0.01")), 
#          hjust = 1, vjust = 0, size = 2)+
# annotate("text", x = 3.5, y = 2200, 
#          label = expression(bold("Ecotype: p<0.01")), 
#          hjust = 1, vjust = 0, size = 2)
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



# Assess the effects of fibrosis on AS ----------------------------

# Add fibrosis to the best model
ASEndQ.fib <- lm(AS ~ I(Temp^2) + Temp+  Ecotype + Fibrosis, data = End)
plot(ASEndQ.fib)
AIC(ASEndQ.fib) # 998.106
summary(ASEndQ.fib)
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -629.39 -218.70   14.45  181.09 1066.24 
# 
# Coefficients:
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

End$Fibrosis <- as.factor(End$Fibrosis)
End$Fibrosis
ASfib <- ggplot(End, aes(x=factor(Temp), y=AS, color=Fibrosis)) +
  geom_point(aes(shape=Ecotype, group=interaction(Fibrosis, Ecotype)),
             position = position_jitterdodge(jitter.width = 0.5, jitter.height = 0.0001, dodge.width = .5), 
             size = 3, alpha = 0.2) +
  theme_classic() +
  scale_color_manual(values = c('black', "blue"), labels = c("Absent", "Present")) +
  labs(x = expression("Temperature"*degree*C), y = expression(Aerobic~Scope~(mgO[2]/kg/hr))) +
  stat_summary(aes(shape=Ecotype, group=interaction(Fibrosis, Ecotype)),
               fun = mean, na.rm = TRUE, 
               geom = "point",
               size = 3, 
               position = position_dodge(width = 0.5)) +
  stat_summary(aes(group=interaction(Fibrosis, Ecotype)),
               fun.data = mean_se, na.rm = TRUE,
               geom = "errorbar", width = .3,
               position = position_dodge(width = 0.5)) +
  stat_smooth(aes(x = as.numeric(factor(Temp)), group=Fibrosis), 
              method = "lm", 
              formula = y ~ poly(x, 2, raw = TRUE), 
              linewidth = 1, se = TRUE, alpha = .1)  
ASfib

ASfib <- ggplot(End, aes(x=Fibrosis, y=AS, color=Fibrosis)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) +
  theme_classic() +
  scale_color_manual(values = c('black', "blue")) +
  labs(x = "Fibrosis", y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))
ASfib

FibLabels <- c("Absent", "Present")

ASfib <- ggplot(End, aes(x=Fibrosis, y=AS, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
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
  labs(x = "Fibrosis", y = expression(Aerobic~Scope~(mgO[2]/kg/hr)))
ASfib

SMRfib <- ggplot(End, aes(x=Fibrosis, y=RMR, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(legend.position = "none") +
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
  labs(x = "Fibrosis", y = expression(Standard~Metabolic~Rate~(mgO[2]/kg/hr)))
SMRfib

MMRfib <- ggplot(End, aes(x=Fibrosis, y=MMR, color=Ecotype)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitterdodge(jitter.width=0.2, dodge.width=0.75)) +
  theme_classic() +
  theme(legend.position = "none")+
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
  labs(x = "Fibrosis", y = expression(Maximum~Metabolic~Rate~(mgO[2]/kg/hr)))
MMRfib


# Make figure -------------------------------------------------------------



pdf(file = "/home/ekerns/ThermTol/Figures/Figure3Q.pdf",
    width = 7,
    height = 6.25)

MMRpoint + RMRpoint + ASpoint + MMREnd + RMREnd + ASEnd

dev.off()

pdf(file = "/home/ekerns/ThermTol/Figures/Figure3L.pdf",
    width = 7,
    height = 6.25)

MMRpoint + RMRpoint + ASpointL + MMREnd + RMREnd + ASEnd

dev.off()

pdf(file = "/home/ekerns/ThermTol/Figures/SuppFigure2.pdf",
    width = 6.6, 
    height = 4)

MMRSex + RMRSex + ASSex

dev.off()

pdf(file = "/home/ekerns/ThermTol/Figures/SuppFigure1.pdf",
    width = 6.5, 
    height = 4)

MMRSApoint + RMRSApoint + ASSApoint

dev.off()

pdf(file = "/home/ekerns/ThermTol/Figures/SuppFigure1b.pdf",
    width = 6.5, 
    height = 4)

MMRfib+ SMRfib + ASfib 

dev.off()




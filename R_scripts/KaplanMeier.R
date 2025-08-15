# Mortality
# Figure 2

######### Load Packages #########
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
library(lme4)
library(emmeans)
library(lmerTest)
library(lmtest)
library(ranger)
library(ggfortify)
library(survival)
library(survminer)
library(ggsurvfit)

# import data
Survival <- read_xlsx("R:/ThermTol_MultStress/Sasser_Ch3/Ch3/ThermalToleranceProject_SuspectDataRemoved.xlsx", sheet = "Survival")
head(Survival)

# covert to factors
Survival$Ecotype <- factor(Survival$Ecotype, 
                          levels = c("B", "L"), 
                          labels = c("Benthic", "Limnetic"))

SurvBen <- subset(Survival, Ecotype == "Benthic")
SurvLim <- subset(Survival, Ecotype == "Limnetic")



# actual code
surv_object <- Surv(time = Survival$Time, event = Survival$Status)
surv_object 

surv_Ben <- Surv(time = SurvBen$Time, event = SurvBen$Status)
surv_Lim <- Surv(time = SurvLim$Time, event = SurvLim$Status)


# ~ temp alone
fit1 <- survfit(surv_object ~ Temp, data = Survival)
summary(fit1)

ggsurvplot(fit1, data = Survival, pval = TRUE)

# ~ temp (ben)
fitB <- survfit(surv_Ben ~ Temp, data = SurvBen)
summary(fitB)

ggsurvplot(fitB, data = SurvBen, pval = TRUE)

# ~ temp (lim)
fitL <- survfit(surv_Lim ~ Temp, data = SurvLim)
summary(fitL)

ggsurvplot(fitL, data = SurvLim, pval = TRUE)



# Fit a Cox proportional hazards model
Temp.cox <- coxph(Surv(Time, Status) ~ Temp, data = Survival)
Temp.cox
# coef exp(coef) se(coef)     z        p
# Temp 0.44803   1.56522  0.08491 5.277 1.32e-07
# 
# Likelihood ratio test=38.48  on 1 df, p=5.544e-10
# n= 211, number of events= 36 

Eco.cox <- coxph(Surv(Time, Status) ~ Ecotype, data = Survival)
Eco.cox
# coef exp(coef) se(coef)     z       p
# EcotypeLimnetic 1.0899    2.9740   0.3722 2.928 0.00341
# 
# Likelihood ratio test=9.67  on 1 df, p=0.001872
# n= 211, number of events= 36 

ET.cox <- coxph(Surv(Time, Status) ~ Ecotype*Temp, data = Survival)
ET.cox

# coef  exp(coef)   se(coef)      z      p
# EcotypeLimnetic      -7.1193092  0.0008093  4.3218498 -1.647 0.0995
# Temp                  0.2757683  1.3175425  0.1318690  2.091 0.0365
# EcotypeLimnetic:Temp  0.3476727  1.4157688  0.1783450  1.949 0.0512
# 
# Likelihood ratio test=56.54  on 3 df, p=3.215e-12
# n= 211, number of events= 36 

## New plot attempt
# all together
p <- survfit2(Surv(Time, Status) ~ Temp, data = Survival) |>
  ggsurvfit(linewidth = 1) +
  scale_color_manual(values = c('dodgerblue4','deepskyblue','forestgreen', 'darkorange','#F00505')) +
  scale_ggsurvfit() +
  theme_classic(base_size = 20) +
  labs(
    y = "Percentage survival",
    x = "Time (days)"
  )+
  add_legend_title("Temperature") +
  theme(legend.position = "bottom") 
p

#Print
pdf(file = "TTKaplanMeierCombined.pdf", useDingbats = FALSE)
p
dev.off()

# benthic
p1 <- survfit2(Surv(Time, Status) ~ Temp, data = SurvBen) |>
  ggsurvfit(linewidth = 1) +
  scale_color_manual(values = c('dodgerblue4','deepskyblue','forestgreen', 'darkorange','#F00505')) +
  scale_ggsurvfit() +
  theme_classic(base_size = 20) +
  labs(title = "A",
    y = "Percentage survival",
    x = "Time (days)"
  )+
  add_legend_title("Temperature") +
  theme(legend.position = "bottom") 
p1

#Limnetic
p2 <- survfit2(Surv(Time, Status) ~ Temp, data = SurvLim) |>
  ggsurvfit(linewidth = 1) +
  scale_color_manual(values = c('dodgerblue4','deepskyblue','forestgreen', 'darkorange','#F00505')) +
  scale_ggsurvfit() +
  theme_classic(base_size = 20) +
  labs(title = "B",
       y = "Percentage survival",
       x = "Time (days)"
  )+
  add_legend_title("Temperature") +
  theme(legend.position = "bottom") 
p2

#Print
pdf(file = "TTKaplanMeier.pdf", useDingbats = FALSE)
KapMei <- ggarrange(p1,p2 ,ncol=2, nrow=1, common.legend = TRUE, legend = "bottom")
KapMei
dev.off()

pdf(file = "KMBenthic.pdf", useDingbats = FALSE)
p1
dev.off()

pdf(file = "KMLimnetic.pdf", useDingbats = FALSE)
p2
dev.off()


###############################################
Surv18 <- subset(Survival, Temp == "18")
Surv20 <- subset(Survival, Temp == "20")
Surv22 <- subset(Survival, Temp == "22")
Surv24 <- subset(Survival, Temp == "24")
Surv26 <- subset(Survival, Temp == "26")



surv_18 <- Surv(time = Surv18$Time, event = Surv18$Status)
surv_18

surv_20 <- Surv(time = Surv20$Time, event = Surv20$Status)
surv_20

surv_22 <- Surv(time = Surv22$Time, event = Surv22$Status)
surv_22

surv_24 <- Surv(time = Surv24$Time, event = Surv24$Status)
surv_24

surv_26 <- Surv(time = Surv26$Time, event = Surv26$Status)
surv_26

# by temp
# ecotype at 18
E18 <- survfit(surv_18 ~ Ecotype, data = Surv18)
summary(E18)

ggsurvplot(E18, data = Surv18, pval = TRUE)

# ecotype at 20
E20 <- survfit(surv_20 ~ Ecotype, data = Surv20)
summary(E20)

ggsurvplot(E20, data = Surv20, pval = TRUE)

# ecotype at 22
E22 <- survfit(surv_22 ~ Ecotype, data = Surv22)
summary(E22)

ggsurvplot(E22, data = Surv22, pval = TRUE)

# ecotype at 24
E24 <- survfit(surv_24 ~ Ecotype, data = Surv24)
summary(E24)

ggsurvplot(E24, data = Surv24, pval = TRUE)

# ecotype at 26
E26 <- survfit(surv_26 ~ Ecotype, data = Surv26)
summary(E26)

ggsurvplot(E26, data = Surv26, pval = TRUE)


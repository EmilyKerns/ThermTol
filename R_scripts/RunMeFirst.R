#Run me first

#This script loads packages and organizes the data needed for all analyses and graphs in the publication.

# Load Packages -----------------------------------------------------------
{
library(dplyr)
library(tidyverse)
library(readxl)
library(writexl)
library(patchwork)
library(emmeans)
library(nlme)
library(pgirmess)
library(lme4)
library(ordinal)
}

# Load data ---------------------------------------------------------------

{
  ThermTol <- read_xlsx("ThermalToleranceProject_OG.xlsx")
  head(ThermTol)
  
  ThermTol$Ecotype <- as.factor(ThermTol$Ecotype)
  ThermTol$Fibrosis <- as.factor(ThermTol$Fibrosis)
  ThermTol$Trial <- as.factor(ThermTol$Trial)
}


# Remove fish that died in the rig, don't have AS value, or with negative MMR/RMR --------

{
which(ThermTol$RMR<0) # none
which(ThermTol$MMR<0) # 11 173 231
ThermTol[11,"FishID"] #WK18S_4
ThermTol[173,"FishID"] #WK24S_7
ThermTol[231,"FishID"] #WT24S_15

todrop.1 <- which(ThermTol$FishID %in% c("WK24S_7","SL26S_6","SL26S_9","SL26S_10","SL26S_11","FG26S_2","WT_24S_15","WT22S_2","WT24S_15","WK18S_4"))

ThermTol <- ThermTol[-todrop.1,]
ThermTol <- ThermTol[!is.na(ThermTol$AS),]
ThermTol$RMR <- as.numeric(ThermTol$RMR)
ThermTol$MMR <- as.numeric(ThermTol$MMR)
ThermTol$AS <- as.numeric(ThermTol$AS)
}

# Remove outliers from RMR (aka SMR) and MMR -----------------------------------


# Remove RMR outliers based on >2SD above the mean ------------------------

# For AS, outliers are the mean +/- >2SD. For RMR, outliers are the mean + 2SD. For MMR, outliers are mean - 2SD. Because population/temperature/timepoint are small sample sizes, outliers are assesed by grouping the whole population together regardless of treatment. Therefore, if 2 points are identified as outliers but they are from the same temperature treatment, they will not be excluded to avoid throwing away temperature effects on metabolism. 

# one tailed outlier test - removing excessively high RMR values within each population

{
fg <- ThermTol[ThermTol$Pop == "FG",]
wt <- ThermTol[ThermTol$Pop == "WT",]
sl <- ThermTol[ThermTol$Pop == "SL",]
wk <- ThermTol[ThermTol$Pop == "WK",]
}

{
# Finger
fg$RMR <- as.numeric(fg$RMR)
hist(fg$RMR)
mean(fg$RMR) #800.0645
sd(fg$RMR) #252.6455

# Outliers = >2 SD above the mean
800.0645 + (2*252.6455) #1305.355

which(fg$RMR>1305.355) #19
fg[19,"FishID"] # FG22S_6

# Watson
wt$RMR <- as.numeric(wt$RMR)
hist(wt$RMR)
mean(wt$RMR) #762.8983
sd(wt$RMR) #269.4358

# Outliers = >2 SD above the mean
762.8983 + (2*269.4358) #1301.77

which(wt$RMR>1301.77) #2 53 54 62 63
wt[c(2, 53, 54, 62, 63),"FishID"] # WT20S_10, WT26S_3, WT26S_4, WT26S_13, WT26S_6 - retaining WT26S fish since they are all the same temperature treatment, it is likely an effect of temp

# Spirit
sl$RMR <- as.numeric(sl$RMR)
hist(sl$RMR)
mean(sl$RMR) #745.4838
sd(sl$RMR) #293.9025

# Outliers = >2 SD above the mean
745.4838 + (2*293.9025) #1333.289

which(sl$RMR>1333.289) #51 53 54 58 60
sl[c(51, 53, 54, 58, 60),"FishID"] # SL24S_4, SL24S_6, SL24S_7 ,SL24S_11, SL24S_13 - as with WT, retaining SL24S because it is likely an effect of temp 

# Wik
wk$RMR <- as.numeric(wk$RMR)
hist(wk$RMR)
mean(wk$RMR) #865.8102
sd(wk$RMR) #354.9285

# Outliers = >2 SD above the mean
865.8102 + (2*354.9285) #1575.667

which(wk$RMR>1575.667) #12 15
wk[c(12, 15),"FishID"] # WK26S_3, WK26S_6 - again, likely an effect of temp, keeping for further anlaysis

todrop.2 <- which(ThermTol$FishID %in% c("WT20S_10", "FG22S_6"))
ThermTol <- ThermTol[-todrop.2,]
hist(ThermTol$RMR)
}

# Remove MMR outliers based on >2SD below the mean ------------------------

# one tailed outlier test - removing excessively low MMR values within each population

{
# Finger
fg$MMR <- as.numeric(fg$MMR)
hist(fg$MMR)
mean(fg$MMR) #1765.955
sd(fg$MMR) #273.8971

# Outliers = >2 SD below the mean
1765.955 - (2*273.8971) #1218.161

which(fg$MMR<1218.161) #14
fg[14,"FishID"] # FG22S_1

# Watson
wt$MMR <- as.numeric(wt$MMR)
hist(wt$MMR)
mean(wt$MMR) #1732.038
sd(wt$MMR) #453.289

# Outliers = >2 SD below the mean
1732.038 - (2*453.289) #825.46

which(wt$MMR<825.46) #26 33 34
wt[c(26, 33, 34),"FishID"] # WT22F_7 , WT22F_14, WT22F_15 - retaining WT22 fish since they are all the same temperature treatment, it is likely an effect of temp

# Spirit
sl$MMR <- as.numeric(sl$MMR)
hist(sl$MMR)
mean(sl$MMR) #1807.167
sd(sl$MMR) #491.9084

# Outliers = >2 SD below the mean
1807.167 - (2*491.9084) #823.3502

which(sl$MMR<823.3502) #6 27 28
sl[c(6, 27, 28),"FishID"] # SL22S_6, SL22F_7, SL22F_8 - as with WT, retaining SL22 because it is likely an effect of temp 

# Wik
wk$MMR <- as.numeric(wk$MMR)
hist(wk$MMR)
mean(wk$MMR) #2095.462
sd(wk$MMR) #480.8518

# Outliers = >2 SD below the mean
2095.462 - (2*480.8518) #1133.758

which(wk$MMR<1133.758) #1 3
wk[c(1, 3),"FishID"] # WK18S_1, WK18S_3 - again, both in same temperature treatment so retaining for further analysis

todrop.3 <- which(ThermTol$FishID %in% c("FG22S_1"))
ThermTol <- ThermTol[-todrop.3,]
hist(ThermTol$MMR)
}


# final sample sizes
ftable(addmargins(table(ThermTol$Trial, ThermTol$Pop, ThermTol$Temp)))
#           18  20  22  24  26 Sum
# 
# F   FG     0   0   0   0   0   0
# SL        15   0  13   0   0  28
# WK        0   0   0   0   0   0
# WT        15   0  15   0  11  41
# Sum       30   0  28   0  11  69
# S   FG     7   7   6   6   0  26
# SL        15  15  14  15   6  65
# WK         3   7   7   6   7  30
# WT        16   3  15  15  13  62
# Sum       41  32  42  42  26 183
# Sum FG     7   7   6   6   0  26
# SL        30  15  27  15   6  93
# WK         3   7   7   6   7  30
# WT        31   3  30  15  24 103
# Sum       71  32  70  42  37 252


# Save cleaned data ----------------------------------

write_xlsx(ThermTol, "ThermTol_OutliersRemoved.xlsx")


# Organize data for downstream analysis -----------------------------------

EcotypeMortality <- read_xlsx("EcotypeMortality.xlsx")
head(EcotypeMortality)

Density <- read_xlsx("ThermalToleranceProject_SuspectDataRemoved.xlsx", sheet = "ThermTolData")

{
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
  colnames(LakeData)[10] <- "SA_Old"
  colnames(LakeData)[7] <- "SA"
  colnames(LakeData)[1] <- "Pop"
  head(LakeData)
  LakeSummary <- LakeData %>% 
    dplyr::group_by(Pop) %>% 
    dplyr::summarise(MinDepth = max(DepthM, na.rm = TRUE), 
              MaxDepth = min(DepthM, na.rm = TRUE),
              MedDepth = median(DepthM, na.rm = TRUE),
              MeanDepth = mean(DepthM, na.rm = TRUE),
              SA = first(SA),
              .groups = 'drop')
  LakeSummary
}

{
  ThermTol<- left_join(ThermTol, Density[ , c("FishID", "CCD")], by = c("FishID"))
  ThermTol <- merge(ThermTol, LakeSummary, by = "Pop")
  head(ThermTol)
   ThermTol$Trial <- factor(ThermTol$Trial, levels = c("S", "F"), 
                           labels = c("Start", "End"))
  
  #If considering temp as factor
  ThermTolTAF <- ThermTol 
  ThermTolTAF$Temp <- as.factor(ThermTolTAF$Temp)
  
  EcotypeMortalityTAF <- EcotypeMortality 
  EcotypeMortalityTAF$Temp <- as.factor(EcotypeMortalityTAF$Temp)
  
}

# Clean up health metrics -------------------------------------------------
#original and log transforms were used to check residual fits in the models below. The non-transformed values are used in the publication.
{
  hist(ThermTol$GonadMass)
  #ThermTol$logGonadMass <- log(ThermTol$GonadMass)
  #hist(ThermTol$logGonadMass)
  
  
  hist(ThermTol$Mass)
  #ThermTol$LogMass <- log(ThermTol$Mass)
  #hist(ThermTol$Mass)
  
  hist(ThermTol$SL)
  #ThermTol$LogSL <- log(ThermTol$SL)
  #hist(ThermTol$LogSL)
  
  #all indexes are represented as %s, body condition is Fulton's index K for fish of stickleback size
  
  ThermTol$BodyCond <- ((ThermTol$Mass/(ThermTol$SL)^3)*100000)
  hist(ThermTol$BodyCond)
  
  hist(ThermTol$HSI)
  ThermTol$HSI <- (ThermTol$LiverMass/ThermTol$Mass)*100
  hist(ThermTol$HSI)
  
  hist(ThermTol$SSI)
  ThermTol$SSI <- (ThermTol$SpleenMass/ThermTol$Mass)*100
  hist(ThermTol$SSI)
  
  hist(ThermTol$GSI)
  ThermTol$GSI <- (ThermTol$GonadMass/ThermTol$Mass)*100
  hist(ThermTol$GSI)
}



## Subset data ------------------------------------------------------------
{
  ThermTol$Trial <- sub(".*([SF])_.*", "\\1", ThermTol$FishID)
  ThermTol <- ThermTol %>%
    mutate(
      Trial = case_when(
        Trial == "S" ~ "Start",
        Trial == "F" ~ "End",
        TRUE ~ "Yer missing something"  # A default value if none of the conditions match
      )
    ) 
  head(ThermTol)
  ThermTol$Trial
  
  Start <-subset(ThermTol, Trial== "Start" )
  End <-subset(ThermTol, Trial == "End" )
  StartTAF <-subset(ThermTolTAF, Trial== "Start" )
  EndTAF <-subset(ThermTolTAF, Trial == "End" )
  Limno <-subset(ThermTol, Ecotype== "Limnetic" )
  Benno <-subset(ThermTol, Ecotype== "Benthic" )
}

{
  # a second dataframe to retain all fish for analysis that doesn't involve the respirometry data
  ThermTol2 <- read_xlsx("ThermalToleranceProject_OG.xlsx")
  head(ThermTol)
  hist(ThermTol2$MMR)
  
  ThermTol2$Ecotype <- as.factor(ThermTol2$Ecotype)
  ThermTol2$Fibrosis <- as.factor(ThermTol2$Fibrosis)
  ThermTol2$Trial <- as.factor(ThermTol2$Trial)
  
  ThermTol2<- left_join(ThermTol2, Density[ , c("FishID", "CCD")], by = c("FishID"))
  ThermTol2$Trial <- factor(ThermTol2$Trial, levels = c("S", "F"), 
                           labels = c("Start", "End"))
  
  #If considering temp as factor
  ThermTolTAF <- ThermTol 
  ThermTolTAF$Temp <- as.factor(ThermTolTAF$Temp)
  
  EcotypeMortalityTAF <- EcotypeMortality 
  EcotypeMortalityTAF$Temp <- as.factor(EcotypeMortalityTAF$Temp)
  
    ThermTol2$Trial <- sub(".*([SF])_.*", "\\1", ThermTol2$FishID)
    ThermTol2 <- ThermTol2 %>%
      mutate(
        Trial = case_when(
          Trial == "S" ~ "Start",
          Trial == "F" ~ "End",
          TRUE ~ "Yer missing something"  # A default value if none of the conditions match
        )
      ) 
    ThermTol2$BodyCond <- ((ThermTol2$Mass/(ThermTol2$SL)^3)*100000)
    hist(ThermTol2$BodyCond)
    
    hist(ThermTol2$HSI)
    ThermTol2$HSI <- (ThermTol2$LiverMass/ThermTol2$Mass)*100
    hist(ThermTol2$HSI)
    
    hist(ThermTol2$SSI)
    ThermTol2$SSI <- (ThermTol2$SpleenMass/ThermTol2$Mass)*100
    hist(ThermTol2$SSI)
    
    hist(ThermTol2$GSI)
    ThermTol2$GSI <- (ThermTol2$GonadMass/ThermTol2$Mass)*100
    hist(ThermTol2$GSI)
    head(ThermTol2)
    ThermTol2$Trial
    
    Start2 <-subset(ThermTol2, Trial== "Start" )
    End2 <-subset(ThermTol2, Trial == "End" )
    Limno2 <-subset(ThermTol2, Ecotype== "Limnetic" )
    Benno2 <-subset(ThermTol2, Ecotype== "Benthic" )
  }


dataE <- read.csv2("dataExploratory.csv")
dataC <- read.csv2("dataConfirmatory.csv")

library(lavaan)
library(dplyr)
library(psych)
library(lmtest)
library(kableExtra)
library(corrplot)
library(tidyverse)


# Descriptives ------------------------------------------------------------

describe(select(dataE, PCL1:PCL20))
lowerCor(select(dataE, PCL1:PCL20))
describe(select(dataC, PCL1:PCL20))
lowerCor(select(dataC, PCL1:PCL20))


##Descriptives on a combined dataset
dataFull <- rbind(dataE, dataC)
describe(select(dataFull, PCL1:PCL20))
lowerCor(select(dataFull, PCL1:PCL20))
a <- select(dataFull, PCL1:PCL20)
a <- a %>% rename(B1 = PCL1,
                  B2 = PCL2,
                  B3 = PCL3,
                  B4 = PCL4,
                  B5 = PCL5,
                  C1 = PCL6,
                  C2 = PCL7,
                  D1 = PCL8,
                  D2 = PCL9,
                  D3 = PCL10,
                  D4 = PCL11,
                  D5 = PCL12,
                  D6 = PCL13,
                  D7 = PCL14,
                  E1 = PCL15,
                  E2 = PCL16,
                  E3 = PCL17,
                  E4 = PCL18,
                  E5 = PCL19,
                  E6 = PCL20)
a <- lowerCor(a)
corrplot::corrplot.mixed(a)

dataFull$PTSD <-  rowSums(select(dataFull, PCL1:PCL20), na.rm = TRUE)
dataFull <- dataFull %>% 
  rowwise() %>% 
  mutate(QIDS = max(QIDS1, QIDS2, QIDS3, QIDS4) + QIDS5 + max(QIDS6, QIDS7, QIDS8, QIDS9) + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + max(QIDS15, QIDS16))
dataFull$RISKY <-  rowSums(select(dataFull, risky1:risky14), na.rm = TRUE)
dataFull$GAD <-  rowSums(select(dataFull, GAD1:GAD7), na.rm = TRUE)
dataFull$AUDIT <-  rowSums(select(dataFull, AUDIT1:AUDIT10), na.rm = TRUE)
describe(select(dataFull, PCL1:PCL20))
lowerCor(select(dataFull, PCL1:PCL20))
describe(select(dataFull, PTSD, QIDS, RISKY, GAD, AUDIT))
lowerCor(select(dataFull, PTSD, QIDS, RISKY, GAD, AUDIT))
print(corr.test(select(dataFull, PTSD, QIDS, RISKY, GAD, AUDIT)), short = FALSE)
#reliability - omega
omega(select(dataFull, PCL1:PCL20), nfactors = 1, poly = TRUE)
omega(select(dataFull, QIDS1:QIDS16), nfactors = 1, poly = TRUE)
omega(select(dataFull, risky1:risky14), nfactors = 1, poly = FALSE)
omega(select(dataFull, GAD1:GAD7), nfactors = 1, poly = TRUE)
omega(select(dataFull, AUDIT1:AUDIT10), nfactors = 1, poly = TRUE)



# One-factor model --------------------------------------------------------

oneFactor <- 'ptsd =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5 + PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + 
                      PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedOneFactor <- cfa(oneFactor, dataE, meanstructure = TRUE,
                       std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                       orthogonal = FALSE, missing = "ML",
                       bootstrap = 5000)
fitmeasures(fittedOneFactor, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedOneFactor, "all")
summary(fittedOneFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedOneFactor)$cov
residuals(fittedOneFactor)$cov > .10
modindices(fittedOneFactor)[order(-modindices(fittedOneFactor)$mi), ]

# DSM-5 model -------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]; 
# 2. Avoidance [C1, C2]; 
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]; 
# 4. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])


dsm5 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
         Avoidance =~ PCL6 + PCL7
         NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
         AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedDsm5 <- cfa(dsm5, dataE, meanstructure = TRUE,
                  std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                  orthogonal = FALSE, missing = "ML",
                  bootstrap = 5000)
fitmeasures(fittedDsm5, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedDsm5, "all")
summary(fittedDsm5, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5)$cov
residuals(fittedDsm5)$cov > .10
modindices(fittedDsm5)[order(-modindices(fittedDsm5)$mi), ]


# Three-factor model ------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance/Negative alternation in cognition and mood [C1, C2, D1, D2, D3, D4, D5, D6, D7]
# 3. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])

threeFactor <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                AvoidNegAlternation =~ PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedThreeFactor <- cfa(threeFactor, dataE, meanstructure = TRUE,
                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                         orthogonal = FALSE, missing = "ML",
                         bootstrap = 5000)
fitmeasures(fittedThreeFactor, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedThreeFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedThreeFactor)$cov
residuals(fittedThreeFactor)$cov > .10
modindices(fittedThreeFactor)[order(-modindices(fittedThreeFactor)$mi), ]


# DSM-5 Dysphoria model ---------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Dyshoria [D1, D2, D3, D4, D5, D6, D7, E1, E2, E5, E6]
# 4. Alternation in arousal and reactivity [E3, E4])

dsm5dysphoria <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL19 + PCL20
                  AltArousal =~ PCL17 + PCL18'
fittedDsm5dysphoria <- cfa(dsm5dysphoria, dataE, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedDsm5dysphoria, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphoria, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphoria)$cov
residuals(fittedDsm5dysphoria)$cov > .10
modindices(fittedDsm5dysphoria)[order(-modindices(fittedDsm5dysphoria)$mi), ]


# DSM-5 Dysphoric arousal model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]
# 4. Dysphoric arousal [E1, E2, E5, E6]
# 5. Anxious arousal [E3, E4])

dsm5dysphAruosal <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18'
fittedDsm5dysphArousal <- cfa(dsm5dysphAruosal, dataE, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedDsm5dysphArousal, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphArousal, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphArousal)$cov
residuals(fittedDsm5dysphArousal)$cov > .10
modindices(fittedDsm5dysphArousal)[order(-modindices(fittedDsm5dysphArousal)$mi), ]


# Anhedonia model ---------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Dysphoric arousal [E1, E2, E5, E6]
# 6. Anxious arousal [E3, E4])

anhedonia <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
              Avoidance =~ PCL6 + PCL7
              NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
              Anhedonia =~ PCL12 + PCL13 + PCL14
              DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
              AnxArousal =~ PCL17 + PCL18'
fittedAnhedonia <- cfa(anhedonia, dataE, meanstructure = TRUE,
                       std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                       orthogonal = FALSE, missing = "ML",
                       bootstrap = 5000)
fitmeasures(fittedAnhedonia, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedAnhedonia, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedAnhedonia)$cov
residuals(fittedAnhedonia)$cov > .10
modindices(fittedAnhedonia)[order(-modindices(fittedAnhedonia)$mi), ]


# Externalizing behaviors model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Numbing [D1, D2, D3, D4, D5, D6, D7]
# 4. Externalizing behaviors [E1, E2]
# 5. Anxious arousal [E3, E4]
# 6. Dysphoric arousal [E5, E6])

extBehaviors <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                 Avoidance =~ PCL6 + PCL7
                 Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                 ExtBehaviors =~ PCL15 + PCL16
                 AnxArousal =~ PCL17 + PCL18 
                 DysphArousal =~ PCL19 + PCL20'
fittedExtBehaviors <- cfa(extBehaviors, dataE, meanstructure = TRUE,
                          std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                          orthogonal = FALSE, missing = "ML",
                          bootstrap = 5000)
fitmeasures(fittedExtBehaviors, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedExtBehaviors, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedExtBehaviors)$cov
residuals(fittedExtBehaviors)$cov > .10
modindices(fittedExtBehaviors)[order(-modindices(fittedExtBehaviors)$mi), ]


# Hybrid model ------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Externalizing behaviours [E1, E2]
# 6. Anxious arousal [E3, E4]
# 7. Dysphoric arousal [E5, E6])

hybridModel <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                Avoidance =~ PCL6 + PCL7
                NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
                Anhedonia =~ PCL12 + PCL13 + PCL14
                ExtBehavior =~ PCL15 + PCL16 
                AnxArousal =~ + PCL17 + PCL18
                DysphArousal =~ PCL19 + PCL20'
fittedHybridModel <- cfa(hybridModel, dataE, meanstructure = TRUE,
                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                         orthogonal = FALSE, missing = "ML",
                         bootstrap = 5000)
fitmeasures(fittedHybridModel, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedHybridModel, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedHybridModel)$cov
residuals(fittedHybridModel)$cov > .10
modindices(fittedHybridModel)[order(-modindices(fittedHybridModel)$mi), ]


# MIs for CFAs ------------------------------------------------------------

##Apart from the one-factor model
# PCL12 ~~ PCL14
# PCL15 ~~ PCL16
# PCL2 ~~ PCL15


# Fit measures summary ----------------------------------------------------

model <- c("One-factor model", "DSM-5 model", "Three-factor model", "DSM-5 Dysphoria model", "DSM-5 Dysphoric arousal model",
           "Anhedonia model", "Externalizing behaviors model", "Hybrid model")
chisq <- c(fitmeasures(fittedOneFactor, "chisq"), fitmeasures(fittedDsm5, "chisq"), 
           fitmeasures(fittedThreeFactor, "chisq"), fitmeasures(fittedDsm5dysphoria, "chisq"), 
           fitmeasures(fittedDsm5dysphArousal, "chisq"), fitmeasures(fittedAnhedonia, "chisq"), 
           fitmeasures(fittedExtBehaviors, "chisq"), fitmeasures(fittedHybridModel, "chisq")) 
chisq <- round(chisq, 3)
df <- c(fitmeasures(fittedOneFactor, "df"), fitmeasures(fittedDsm5, "df"), 
        fitmeasures(fittedThreeFactor, "df"), fitmeasures(fittedDsm5dysphoria, "df"), 
        fitmeasures(fittedDsm5dysphArousal, "df"), fitmeasures(fittedAnhedonia, "df"), 
        fitmeasures(fittedExtBehaviors, "df"), fitmeasures(fittedHybridModel, "df")) 
df <- round(df, 3)
pvalue <- c(fitmeasures(fittedOneFactor, "pvalue"), fitmeasures(fittedDsm5, "pvalue"), 
            fitmeasures(fittedThreeFactor, "pvalue"), fitmeasures(fittedDsm5dysphoria, "pvalue"),
            fitmeasures(fittedDsm5dysphArousal, "pvalue"), fitmeasures(fittedAnhedonia, "pvalue"), 
            fitmeasures(fittedExtBehaviors, "pvalue"), fitmeasures(fittedHybridModel, "pvalue")) 
pvalue <- round(pvalue, 3)
cfi <- c(fitmeasures(fittedOneFactor, "cfi"), fitmeasures(fittedDsm5, "cfi"), 
         fitmeasures(fittedThreeFactor, "cfi"), fitmeasures(fittedDsm5dysphoria, "cfi"), 
         fitmeasures(fittedDsm5dysphArousal, "cfi"), fitmeasures(fittedAnhedonia, "cfi"), 
         fitmeasures(fittedExtBehaviors, "cfi"), fitmeasures(fittedHybridModel, "cfi")) 
cfi <- round(cfi, 3)
tli <- c(fitmeasures(fittedOneFactor, "tli"), fitmeasures(fittedDsm5, "tli"), 
         fitmeasures(fittedThreeFactor, "tli"), fitmeasures(fittedDsm5dysphoria, "tli"), 
         fitmeasures(fittedDsm5dysphArousal, "tli"), fitmeasures(fittedAnhedonia, "tli"), 
         fitmeasures(fittedExtBehaviors, "tli"), fitmeasures(fittedHybridModel, "tli")) 
tli <- round(tli, 3)
rmsea <- c(fitmeasures(fittedOneFactor, "rmsea"), fitmeasures(fittedDsm5, "rmsea"), 
           fitmeasures(fittedThreeFactor, "rmsea"), fitmeasures(fittedDsm5dysphoria, "rmsea"), 
           fitmeasures(fittedDsm5dysphArousal, "rmsea"), fitmeasures(fittedAnhedonia, "rmsea"), 
           fitmeasures(fittedExtBehaviors, "rmsea"), fitmeasures(fittedHybridModel, "rmsea"))
rmsea <- round(rmsea, 3)
rmseaLI <- c(fitmeasures(fittedOneFactor, "rmsea.ci.lower"), fitmeasures(fittedDsm5, "rmsea.ci.lower"), 
             fitmeasures(fittedThreeFactor, "rmsea.ci.lower"), fitmeasures(fittedDsm5dysphoria, "rmsea.ci.lower"), 
             fitmeasures(fittedDsm5dysphArousal, "rmsea.ci.lower"), fitmeasures(fittedAnhedonia, "rmsea.ci.lower"),
             fitmeasures(fittedExtBehaviors, "rmsea.ci.lower"), fitmeasures(fittedHybridModel, "rmsea.ci.lower")) 
rmseaLI <- round(rmseaLI, 3)
rmseaUI<- c(fitmeasures(fittedOneFactor, "rmsea.ci.upper"), fitmeasures(fittedDsm5, "rmsea.ci.upper"), 
            fitmeasures(fittedThreeFactor, "rmsea.ci.upper"), fitmeasures(fittedDsm5dysphoria, "rmsea.ci.upper"), 
            fitmeasures(fittedDsm5dysphArousal, "rmsea.ci.upper"), fitmeasures(fittedAnhedonia, "rmsea.ci.upper"), 
            fitmeasures(fittedExtBehaviors, "rmsea.ci.upper"), fitmeasures(fittedHybridModel, "rmsea.ci.upper")) 
rmseaUI <- round(rmseaUI, 3)
srmr <- c(fitmeasures(fittedOneFactor, "srmr"), fitmeasures(fittedDsm5, "srmr"), 
          fitmeasures(fittedThreeFactor, "srmr"), fitmeasures(fittedDsm5dysphoria, "srmr"), 
          fitmeasures(fittedDsm5dysphArousal, "srmr"), fitmeasures(fittedAnhedonia, "srmr"), 
          fitmeasures(fittedExtBehaviors, "srmr"), fitmeasures(fittedHybridModel, "srmr")) 
srmr <- round(srmr, 3)
aic <- c(fitmeasures(fittedOneFactor, "aic"), fitmeasures(fittedDsm5, "aic"), 
         fitmeasures(fittedThreeFactor, "aic"), fitmeasures(fittedDsm5dysphoria, "aic"), 
         fitmeasures(fittedDsm5dysphArousal, "aic"), fitmeasures(fittedAnhedonia, "aic"), 
         fitmeasures(fittedExtBehaviors, "aic"), fitmeasures(fittedHybridModel, "aic")) 
aic <- round(aic, 3)
bic <- c(fitmeasures(fittedOneFactor, "bic"), fitmeasures(fittedDsm5, "bic"), 
         fitmeasures(fittedThreeFactor, "bic"), fitmeasures(fittedDsm5dysphoria, "bic"), 
         fitmeasures(fittedDsm5dysphArousal, "bic"), fitmeasures(fittedAnhedonia, "bic"), 
         fitmeasures(fittedExtBehaviors, "bic"), fitmeasures(fittedHybridModel, "bic")) 
bic <- round(bic, 3)

modelFits <- data.frame(model, chisq, df, pvalue, cfi, tli, rmsea, rmseaLI, rmseaUI, srmr, aic, bic)


kbl(modelFits)
modelFits %>% kbl() %>% kable_material()


lowerCor(select(dataE, PCL1:PCL20))
mean(lowerCor(select(dataE, PCL1:PCL20)))
sd(lowerCor(select(dataE, PCL1:PCL20)))



# QIDS --------------------------------------------------------------------

qids <- 'qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
                 QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16'
fittedQids <- cfa(qids, dataE, meanstructure = TRUE,
                  std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                  orthogonal = FALSE, missing = "ML",
                  bootstrap = 5000)
fitmeasures(fittedQids, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedQids, "all")
summary(fittedQids, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedQids)$cov
residuals(fittedQids)$cov > .10
modindices(fittedQids)[order(-modindices(fittedQids)$mi), ]


qids2 <- 'qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
                  QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                  QIDS6 ~~ QIDS8
                  QIDS7 ~~ QIDS9
                  QIDS2 ~~ QIDS3'
fittedQids2 <- cfa(qids2, dataE, meanstructure = TRUE,
                   std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                   orthogonal = FALSE, missing = "ML",
                   bootstrap = 5000)
fitmeasures(fittedQids2, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedQids2, "all")
summary(fittedQids2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedQids2)$cov
residuals(fittedQids2)$cov > .10
modindices(fittedQids2)[order(-modindices(fittedQids2)$mi), ]


# GAD ---------------------------------------------------------------------

gad <- 'gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7'
fittedGad <- cfa(gad, dataE, meanstructure = TRUE,
                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                 orthogonal = FALSE, missing = "ML",
                 bootstrap = 5000)
fitmeasures(fittedGad, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGad, "all")
summary(fittedGad, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGad)$cov
residuals(fittedGad)$cov > .10
modindices(fittedGad)[order(-modindices(fittedGad)$mi), ]


gad2 <- 'gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                GAD1 ~~ GAD6
                GAD5 ~~ GAD7'
fittedGad2 <- cfa(gad2, dataE, meanstructure = TRUE,
                  std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                  orthogonal = FALSE, missing = "ML",
                  bootstrap = 5000)
fitmeasures(fittedGad2, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGad2, "all")
summary(fittedGad2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGad2)$cov
residuals(fittedGad2)$cov > .10
modindices(fittedGad2)[order(-modindices(fittedGad2)$mi), ]


# Risky behavior ----------------------------------------------------------

risky <- 'risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                   risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14'
fittedRisky <- cfa(risky, dataE, meanstructure = TRUE,
                   std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                   orthogonal = FALSE, missing = "ML",
                   bootstrap = 5000)
fitmeasures(fittedRisky, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedRisky, "all")
summary(fittedRisky, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedRisky)$cov
residuals(fittedRisky)$cov > .10
modindices(fittedRisky)[order(-modindices(fittedRisky)$mi), ]

risky2 <- 'risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                    risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                    risky13 ~~ risky14
                    risky2 ~~ risky7'
fittedRisky2 <- cfa(risky2, dataE, meanstructure = TRUE,
                    std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                    orthogonal = FALSE, missing = "ML",
                    bootstrap = 5000)
fitmeasures(fittedRisky2, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedRisky2, "all")
summary(fittedRisky2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedRisky2)$cov
residuals(fittedRisky2)$cov > .10
modindices(fittedRisky2)[order(-modindices(fittedRisky2)$mi), ]

# AUDIT -------------------------------------------------------------------

audit <- 'audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                   AUDIT8 + AUDIT9 + AUDIT10'
fittedAudit <- cfa(audit, dataE, meanstructure = TRUE,
                   std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                   orthogonal = FALSE, missing = "ML",
                   bootstrap = 5000)
fitmeasures(fittedAudit, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedAudit, "all")
summary(fittedAudit, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedAudit)$cov
residuals(fittedAudit)$cov > .10
modindices(fittedAudit)[order(-modindices(fittedAudit)$mi), ]

audit2 <- 'audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                    AUDIT8 + AUDIT9 + AUDIT10
                    AUDIT2 ~~ AUDIT3
                    AUDIT1 ~~ AUDIT2'
fittedAudit2 <- cfa(audit2, dataE, meanstructure = TRUE,
                    std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                    orthogonal = FALSE, missing = "ML",
                    bootstrap = 5000)
fitmeasures(fittedAudit2, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedAudit2, "all")
summary(fittedAudit2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedAudit2)$cov
residuals(fittedAudit2)$cov > .10
modindices(fittedAudit2)[order(-modindices(fittedAudit2)$mi), ]


# Global One-Factor Model -------------------------------------------------

# globalOneFactor <- 'ptsd =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5 + PCL6 + PCL7 + PCL8 + PCL9 + PCL10 +
#                             PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
#                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
#                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                              AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalOneFactor <- sem(globalOneFactor, dataE, meanstructure = TRUE,
#                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                              orthogonal = FALSE, missing = "ML",
#                              bootstrap = 5000)
# fitmeasures(fittedGlobalOneFactor,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalOneFactor, "all")
# summary(fittedGlobalOneFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalOneFactor)$cov
# residuals(fittedGlobalOneFactor)$cov > .10
# globalOneFactorMI <- modindices(fittedGlobalOneFactor)[order(-modindices(fittedGlobalOneFactor)$mi), ]
# globalOneFactorMI

globalOneFactor2 <- 'ptsd =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5 + PCL6 + PCL7 + PCL8 + PCL9 + PCL10 +
                             PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2'
fittedGlobalOneFactor2 <- sem(globalOneFactor2, dataE, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalOneFactor2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalOneFactor2, "all")
summary(fittedGlobalOneFactor2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalOneFactor2)$cov
residuals(fittedGlobalOneFactor2)$cov > .10
globalOneFactor2MI <- modindices(fittedGlobalOneFactor2)[order(-modindices(fittedGlobalOneFactor2)$mi), ]
globalOneFactor2MI


####Uncomment the lines no. 390-663 to run all general models 

# # Global DSM-5 model ------------------------------------------------------
# 
# globalDsm5 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                Avoidance =~ PCL6 + PCL7
#                NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
#                AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
#                qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                        QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                         risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                         AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalDsm5 <- sem(globalDsm5, dataE, meanstructure = TRUE,
#                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                         orthogonal = FALSE, missing = "ML",
#                         bootstrap = 5000)
# fitmeasures(fittedGlobalDsm5,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalDsm5, "all")
# summary(fittedGlobalDsm5, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalDsm5)$cov
# residuals(fittedGlobalDsm5)$cov > .10
# globalDsm5MI <- modindices(fittedGlobalDsm5)[order(-modindices(fittedGlobalDsm5)$mi), ]
# globalDsm5MI
# #                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 701           audit =~  risky1 187.701  0.438   0.438    0.597    0.597
# # 2890         AUDIT2 ~~  AUDIT3  94.882  0.245   0.245    0.452    0.452
# # 1798          PCL20 ~~   QIDS1  83.773  0.265   0.265    0.389    0.389
# # 2860        risky13 ~~ risky14  80.414  0.047   0.047    0.405    0.405
# # 2320         QIDS12 ~~ risky14  65.638  0.050   0.050    0.352    0.352
# # 2066          QIDS6 ~~   QIDS8  59.140  0.087   0.087    0.318    0.318
# # 2107          QIDS7 ~~   QIDS9  47.363  0.122   0.122    0.282    0.282
# # 2319         QIDS12 ~~ risky13  39.966  0.041   0.041    0.272    0.272
# # 709           audit =~  risky9  37.959 -0.117  -0.117   -0.283   -0.283
# # 2677         risky2 ~~  risky7  37.230  0.023   0.023    0.268    0.268
# # 475      AltArousal =~  risky6  37.075  0.193   0.193    0.263    0.263
# # 527            qids =~  risky7  37.046 -0.080  -0.080   -0.269   -0.269

globalDsm52 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                Avoidance =~ PCL6 + PCL7
                NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
                qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                        QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                        QIDS6 ~~ QIDS8
                        QIDS7 ~~ QIDS9
                        QIDS2 ~~ QIDS3
                gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                       GAD1 ~~ GAD6
                       GAD5 ~~ GAD7
                risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                         risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                         risky13 ~~ risky14
                         risky2 ~~ risky7
                audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                         AUDIT8 + AUDIT9 + AUDIT10
                         AUDIT2 ~~ AUDIT3
                         AUDIT1 ~~ AUDIT2'
fittedGlobalDsm52 <- sem(globalDsm52, dataE, meanstructure = TRUE,
                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                         orthogonal = FALSE, missing = "ML",
                         bootstrap = 5000)
fitmeasures(fittedGlobalDsm52,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalDsm52, "all")
summary(fittedGlobalDsm52, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalDsm52)$cov
residuals(fittedGlobalDsm52)$cov > .10
globalDsm52MI <- modindices(fittedGlobalDsm52)[order(-modindices(fittedGlobalDsm52)$mi), ]
globalDsm52MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 710           audit =~  risky1 185.228  0.445   0.445    0.606    0.606
# 1807          PCL20 ~~   QIDS1  84.288  0.266   0.266    0.390    0.390
# 718           audit =~  risky9  43.272 -0.128  -0.128   -0.310   -0.310
# 2326         QIDS12 ~~ risky14  37.466  0.035   0.035    0.241    0.241
# 484      AltArousal =~  risky6  34.581  0.190   0.190    0.259    0.259
# 2810         risky9 ~~ risky11  34.165  0.027   0.027    0.258    0.258
# 2779         risky7 ~~ risky11  32.517  0.019   0.019    0.236    0.236
# 595             gad =~  risky6  31.303  0.174   0.174    0.238    0.238
# 535            qids =~  risky6  31.010  0.174   0.174    0.238    0.238
# 2615           GAD6 ~~ risky10  30.627  0.087   0.087    0.227    0.227
# 2661         risky1 ~~  risky9  30.210 -0.048  -0.048   -0.249   -0.249
# 423  NegAlternation =~  risky6  29.838  0.173   0.173    0.237    0.237


# # Global Three-factor model -----------------------------------------------
# 
# globalThreeFactor <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                       AvoidNegAlternation =~ PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
#                       AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
#                       qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                               QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                       gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                       risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                                risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                       audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                                AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalThreeFactor <- sem(globalThreeFactor, dataE, meanstructure = TRUE,
#                                std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                                orthogonal = FALSE, missing = "ML",
#                                bootstrap = 5000)
# fitmeasures(fittedGlobalThreeFactor,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalThreeFactor, "all")
# summary(fittedGlobalThreeFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalThreeFactor)$cov
# residuals(fittedGlobalThreeFactor)$cov > .10
# globalThreeFactorMI <- modindices(fittedGlobalThreeFactor)[order(-modindices(fittedGlobalThreeFactor)$mi), ]
# globalThreeFactorMI
# #                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 701           audit =~  risky1 187.701  0.438   0.438    0.597    0.597
# # 2890         AUDIT2 ~~  AUDIT3  94.882  0.245   0.245    0.452    0.452
# # 1798          PCL20 ~~   QIDS1  83.773  0.265   0.265    0.389    0.389
# # 2860        risky13 ~~ risky14  80.414  0.047   0.047    0.405    0.405
# # 2320         QIDS12 ~~ risky14  65.638  0.050   0.050    0.352    0.352
# # 2066          QIDS6 ~~   QIDS8  59.140  0.087   0.087    0.318    0.318
# # 2107          QIDS7 ~~   QIDS9  47.363  0.122   0.122    0.282    0.282
# # 2319         QIDS12 ~~ risky13  39.966  0.041   0.041    0.272    0.272
# # 709           audit =~  risky9  37.959 -0.117  -0.117   -0.283   -0.283
# # 2677         risky2 ~~  risky7  37.230  0.023   0.023    0.268    0.268
# # 475      AltArousal =~  risky6  37.075  0.193   0.193    0.263    0.263
# # 527            qids =~  risky7  37.046 -0.080  -0.080   -0.269   -0.269


globalThreeFactor2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                       AvoidNegAlternation =~ PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                       AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
                       qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                               QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                               QIDS6 ~~ QIDS8
                               QIDS7 ~~ QIDS9
                               QIDS2 ~~ QIDS3
                       gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                              GAD1 ~~ GAD6
                              GAD5 ~~ GAD7
                       risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                risky13 ~~ risky14
                                risky2 ~~ risky7
                       audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                AUDIT8 + AUDIT9 + AUDIT10
                                AUDIT2 ~~ AUDIT3
                                AUDIT1 ~~ AUDIT2'
fittedGlobalThreeFactor2 <- sem(globalThreeFactor2, dataE, meanstructure = TRUE,
                                std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                orthogonal = FALSE, missing = "ML",
                                bootstrap = 5000)
fitmeasures(fittedGlobalThreeFactor2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalThreeFactor2, "all")
summary(fittedGlobalThreeFactor2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalThreeFactor2)$cov
residuals(fittedGlobalThreeFactor2)$cov > .10
globalThreeFactor2MI <- modindices(fittedGlobalThreeFactor2)[order(-modindices(fittedGlobalThreeFactor2)$mi), ]
globalThreeFactor2MI
#                      lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 634                audit =~  risky1 185.033  0.445   0.445    0.606    0.606
# 1731               PCL20 ~~   QIDS1  84.155  0.266   0.266    0.390    0.390
# 642                audit =~  risky9  43.069 -0.128  -0.128   -0.309   -0.309
# 968                 PCL6 ~~    PCL7  39.031  0.201   0.201    0.267    0.267
# 2250              QIDS12 ~~ risky14  37.552  0.035   0.035    0.241    0.241
# 408           AltArousal =~  risky6  34.616  0.190   0.190    0.259    0.259
# 2734              risky9 ~~ risky11  34.268  0.027   0.027    0.259    0.259
# 2703              risky7 ~~ risky11  32.529  0.019   0.019    0.236    0.236
# 519                  gad =~  risky6  31.307  0.174   0.174    0.238    0.238
# 459                 qids =~  risky6  31.093  0.174   0.174    0.238    0.238
# 2539                GAD6 ~~ risky10  30.681  0.088   0.088    0.228    0.228
# 2585              risky1 ~~  risky9  30.242 -0.048  -0.048   -0.249   -0.249


# # Global DSM-5 Dysphoria model --------------------------------------------
# 
# globalDysphoria <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                     Avoidance =~ PCL6 + PCL7
#                     NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
#                     AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
#                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                              AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalDysphoria <- sem(globalDysphoria, dataE, meanstructure = TRUE,
#                             std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                             orthogonal = FALSE, missing = "ML",
#                             bootstrap = 5000)
# fitmeasures(fittedGlobalDysphoria,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalDysphoria, "all")
# summary(fittedGlobalDysphoria, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalDysphoria)$cov
# residuals(fittedGlobalDysphoria)$cov > .10
# globalDysphoriaMI <- modindices(fittedGlobalDysphoria)[order(-modindices(fittedGlobalDysphoria)$mi), ]
# globalDysphoriaMI
# #                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 701           audit =~  risky1 187.701  0.438   0.438    0.597    0.597
# # 2890         AUDIT2 ~~  AUDIT3  94.882  0.245   0.245    0.452    0.452
# # 1798          PCL20 ~~   QIDS1  83.773  0.265   0.265    0.389    0.389
# # 2860        risky13 ~~ risky14  80.414  0.047   0.047    0.405    0.405
# # 2320         QIDS12 ~~ risky14  65.638  0.050   0.050    0.352    0.352
# # 2066          QIDS6 ~~   QIDS8  59.140  0.087   0.087    0.318    0.318
# # 2107          QIDS7 ~~   QIDS9  47.363  0.122   0.122    0.282    0.282
# # 2319         QIDS12 ~~ risky13  39.966  0.041   0.041    0.272    0.272
# # 709           audit =~  risky9  37.959 -0.117  -0.117   -0.283   -0.283
# # 2677         risky2 ~~  risky7  37.230  0.023   0.023    0.268    0.268
# # 475      AltArousal =~  risky6  37.075  0.193   0.193    0.263    0.263
# # 527            qids =~  risky7  37.046 -0.080  -0.080   -0.269   -0.269


globalDysphoria2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                     AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2'
fittedGlobalDysphoria2 <- sem(globalDysphoria2, dataE, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalDysphoria2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalDysphoria2, "all")
summary(fittedGlobalDysphoria2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalDysphoria2)$cov
residuals(fittedGlobalDysphoria2)$cov > .10
globalDysphoria2MI <- modindices(fittedGlobalDysphoria2)[order(-modindices(fittedGlobalDysphoria2)$mi), ]
globalDysphoria2MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 710           audit =~  risky1 185.228  0.445   0.445    0.606    0.606
# 1807          PCL20 ~~   QIDS1  84.288  0.266   0.266    0.390    0.390
# 718           audit =~  risky9  43.272 -0.128  -0.128   -0.310   -0.310
# 2326         QIDS12 ~~ risky14  37.466  0.035   0.035    0.241    0.241
# 484      AltArousal =~  risky6  34.581  0.190   0.190    0.259    0.259
# 2810         risky9 ~~ risky11  34.165  0.027   0.027    0.258    0.258
# 2779         risky7 ~~ risky11  32.517  0.019   0.019    0.236    0.236
# 595             gad =~  risky6  31.303  0.174   0.174    0.238    0.238
# 535            qids =~  risky6  31.010  0.174   0.174    0.238    0.238
# 2615           GAD6 ~~ risky10  30.627  0.087   0.087    0.227    0.227
# 2661         risky1 ~~  risky9  30.210 -0.048  -0.048   -0.249   -0.249
# 423  NegAlternation =~  risky6  29.838  0.173   0.173    0.237    0.237


# # Global DSM-5 Dysphoric Arousal model ------------------------------------
# 
# globalDysphoricArousal <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                            Avoidance =~ PCL6 + PCL7
#                            NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
#                            DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
#                            AnxArousal =~ PCL17 + PCL18
#                            qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                                    QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                            gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                            risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                                     risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                            audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                                     AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalDysphoricArousal <- sem(globalDysphoricArousal, dataE, meanstructure = TRUE,
#                                     std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                                     orthogonal = FALSE, missing = "ML",
#                                     bootstrap = 5000)
# fitmeasures(fittedGlobalDysphoricArousal,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalDysphoricArousal, "all")
# summary(fittedGlobalDysphoricArousal, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalDysphoricArousal)$cov
# residuals(fittedGlobalDysphoricArousal)$cov > .10
# globalDysphoricArousalMI <- modindices(fittedGlobalDysphoricArousal)[order(-modindices(fittedGlobalDysphoricArousal)$mi), ]
# globalDysphoricArousalMI
# # lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 778           audit =~  risky1 187.389  0.438   0.438    0.596    0.596
# # 2967         AUDIT2 ~~  AUDIT3  94.915  0.245   0.245    0.453    0.453
# # 1780          PCL20 ~~   QIDS1  84.216  0.268   0.268    0.389    0.389
# # 2937        risky13 ~~ risky14  81.412  0.048   0.048    0.406    0.406
# # 2397         QIDS12 ~~ risky14  66.115  0.050   0.050    0.353    0.353
# # 2143          QIDS6 ~~   QIDS8  59.159  0.087   0.087    0.318    0.318
# # 2184          QIDS7 ~~   QIDS9  47.358  0.122   0.122    0.282    0.282
# # 2396         QIDS12 ~~ risky13  40.413  0.041   0.041    0.274    0.274
# # 552      AnxArousal =~  risky6  39.716  0.193   0.193    0.264    0.264
# # 786           audit =~  risky9  38.958 -0.118  -0.118   -0.286   -0.286
# # 2754         risky2 ~~  risky7  37.701  0.023   0.023    0.270    0.270
# # 604            qids =~  risky7  36.756 -0.080  -0.080   -0.268   -0.268


globalDysphoricArousal2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                           Avoidance =~ PCL6 + PCL7
                           NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                           DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                           AnxArousal =~ PCL17 + PCL18
                           qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                   QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                   QIDS6 ~~ QIDS8
                                   QIDS7 ~~ QIDS9
                                   QIDS2 ~~ QIDS3
                           gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                                  GAD1 ~~ GAD6
                                  GAD5 ~~ GAD7
                           risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                    risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                    risky13 ~~ risky14
                                    risky2 ~~ risky7
                           audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                    AUDIT8 + AUDIT9 + AUDIT10
                                    AUDIT2 ~~ AUDIT3
                                    AUDIT1 ~~ AUDIT2'
fittedGlobalDysphoricArousal2 <- sem(globalDysphoricArousal2, dataE, meanstructure = TRUE,
                                     std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                     orthogonal = FALSE, missing = "ML",
                                     bootstrap = 5000)
fitmeasures(fittedGlobalDysphoricArousal2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalDysphoricArousal2, "all")
summary(fittedGlobalDysphoricArousal2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalDysphoricArousal2)$cov
residuals(fittedGlobalDysphoricArousal2)$cov > .10
globalDysphoricArousal2MI <- modindices(fittedGlobalDysphoricArousal2)[order(-modindices(fittedGlobalDysphoricArousal2)$mi), ]
globalDysphoricArousal2MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 787           audit =~  risky1 184.585  0.444   0.444    0.604    0.604
# 1789          PCL20 ~~   QIDS1  84.675  0.269   0.269    0.390    0.390
# 795           audit =~  risky9  44.603 -0.130  -0.130   -0.314   -0.314
# 2403         QIDS12 ~~ risky14  37.580  0.035   0.035    0.241    0.241
# 561      AnxArousal =~  risky6  37.163  0.190   0.190    0.259    0.259
# 270  Reexperiencing =~   PCL11  34.827  0.757   0.757    0.642    0.642
# 496    DysphArousal =~  risky6  33.789  0.194   0.194    0.266    0.266
# 2887         risky9 ~~ risky11  33.295  0.027   0.027    0.255    0.255
# 2856         risky7 ~~ risky11  32.654  0.019   0.019    0.236    0.236
# 672             gad =~  risky6  31.751  0.175   0.175    0.240    0.240
# 2738         risky1 ~~  risky9  31.504 -0.048  -0.048   -0.254   -0.254
# 612            qids =~  risky6  31.210  0.175   0.175    0.239    0.239


# # Global Anhedonia model --------------------------------------------------
# 
# globalAnhedonia <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                      Avoidance =~ PCL6 + PCL7
#                      NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
#                      Anhedonia =~ PCL12 + PCL13 + PCL14
#                      DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
#                      AnxArousal =~ PCL17 + PCL18
#                      qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                              QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                      gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                      risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                               risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                      audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                               AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalAnhedonia <- sem(globalAnhedonia, dataE, meanstructure = TRUE,
#                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                              orthogonal = FALSE, missing = "ML",
#                              bootstrap = 5000)
# fitmeasures(fittedGlobalAnhedonia,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalAnhedonia, "all")
# summary(fittedGlobalAnhedonia, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalAnhedonia)$cov
# residuals(fittedGlobalAnhedonia)$cov > .10
# globalAnhedoniaMI <- modindices(fittedGlobalAnhedonia)[order(-modindices(fittedGlobalAnhedonia)$mi), ]
# globalAnhedoniaMI
# #                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 856           audit =~  risky1 187.889  0.439   0.439    0.597    0.597
# # 3045         AUDIT2 ~~  AUDIT3  94.791  0.245   0.245    0.452    0.452
# # 1858          PCL20 ~~   QIDS1  84.467  0.268   0.268    0.389    0.389
# # 3015        risky13 ~~ risky14  81.391  0.048   0.048    0.406    0.406
# # 2475         QIDS12 ~~ risky14  66.018  0.050   0.050    0.353    0.353
# # 2221          QIDS6 ~~   QIDS8  59.056  0.087   0.087    0.318    0.318
# # 2262          QIDS7 ~~   QIDS9  47.329  0.122   0.122    0.282    0.282
# # 2474         QIDS12 ~~ risky13  40.535  0.041   0.041    0.274    0.274
# # 630      AnxArousal =~  risky6  39.786  0.193   0.193    0.264    0.264
# # 864           audit =~  risky9  39.551 -0.119  -0.119   -0.288   -0.288
# # 2832         risky2 ~~  risky7  37.530  0.023   0.023    0.269    0.269
# # 682            qids =~  risky7  37.449 -0.080  -0.080   -0.270   -0.270


globalAnhedonia2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                     Anhedonia =~ PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2'
fittedGlobalAnhedonia2 <- sem(globalAnhedonia2, dataE, meanstructure = TRUE,

                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)

fitmeasures(fittedGlobalAnhedonia2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalAnhedonia2, "all")
summary(fittedGlobalAnhedonia2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalAnhedonia2)$cov
residuals(fittedGlobalAnhedonia2)$cov > .10
globalAnhedonia2MI <- modindices(fittedGlobalAnhedonia2)[order(-modindices(fittedGlobalAnhedonia2)$mi), ]
globalAnhedonia2MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 185.054  0.444   0.444    0.605    0.605
# 1867          PCL20 ~~   QIDS1  84.890  0.269   0.269    0.390    0.390
# 873           audit =~  risky9  45.168 -0.131  -0.131   -0.316   -0.316
# 2481         QIDS12 ~~ risky14  37.470  0.035   0.035    0.240    0.240
# 639      AnxArousal =~  risky6  37.227  0.190   0.190    0.259    0.259
# 574    DysphArousal =~  risky6  34.758  0.196   0.196    0.268    0.268
# 2965         risky9 ~~ risky11  33.248  0.027   0.027    0.255    0.255
# 2934         risky7 ~~ risky11  32.567  0.019   0.019    0.236    0.236
# 750             gad =~  risky6  31.561  0.175   0.175    0.239    0.239
# 2816         risky1 ~~  risky9  31.541 -0.048  -0.048   -0.254   -0.254
# 690            qids =~  risky6  31.468  0.175   0.175    0.240    0.240
# 447       NegAffect =~  risky6  30.746  0.176   0.176    0.240    0.240


# # Global Externalizing Behaviors model ------------------------------------
# 
# globalExtBehaviors <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                        Avoidance =~ PCL6 + PCL7
#                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
#                        ExtBehaviors =~ PCL15 + PCL16
#                        AnxArousal =~ PCL17 + PCL18 
#                        DysphArousal =~ PCL19 + PCL20
#                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                                 AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalExtBehaviors <- sem(globalExtBehaviors, dataE, meanstructure = TRUE,
#                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                                 orthogonal = FALSE, missing = "ML",
#                                 bootstrap = 5000)
# fitmeasures(fittedGlobalExtBehaviors,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalExtBehaviors, "all")
# summary(fittedGlobalExtBehaviors, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalExtBehaviors)$cov
# residuals(fittedGlobalExtBehaviors)$cov > .10
# globalExtBehaviorsMI <- modindices(fittedGlobalExtBehaviors)[order(-modindices(fittedGlobalExtBehaviors)$mi), ]
# globalExtBehaviorsMI
# #                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# # 856           audit =~  risky1 187.592  0.438   0.438    0.597    0.597
# # 3045         AUDIT2 ~~  AUDIT3  94.791  0.245   0.245    0.452    0.452
# # 1953          PCL20 ~~   QIDS1  83.318  0.263   0.263    0.398    0.398
# # 3015        risky13 ~~ risky14  82.785  0.048   0.048    0.409    0.409
# # 2475         QIDS12 ~~ risky14  64.482  0.050   0.050    0.348    0.348
# # 2221          QIDS6 ~~   QIDS8  59.110  0.087   0.087    0.318    0.318
# # 2262          QIDS7 ~~   QIDS9  47.554  0.122   0.122    0.283    0.283
# # 455    ExtBehaviors =~    PCL2  41.489  0.487   0.487    0.523    0.523
# # 864           audit =~  risky9  40.304 -0.120  -0.120   -0.290   -0.290
# # 2474         QIDS12 ~~ risky13  40.239  0.041   0.041    0.273    0.273
# # 565      AnxArousal =~  risky6  39.642  0.193   0.193    0.264    0.264
# # 2832         risky2 ~~  risky7  38.233  0.023   0.023    0.271    0.271


globalExtBehaviors2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                        Avoidance =~ PCL6 + PCL7
                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                        ExtBehaviors =~ PCL15 + PCL16
                        AnxArousal =~ PCL17 + PCL18
                        DysphArousal =~ PCL19 + PCL20
                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                QIDS6 ~~ QIDS8
                                QIDS7 ~~ QIDS9
                                QIDS2 ~~ QIDS3
                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                               GAD1 ~~ GAD6
                               GAD5 ~~ GAD7
                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                 risky13 ~~ risky14
                                 risky2 ~~ risky7
                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                 AUDIT8 + AUDIT9 + AUDIT10
                                 AUDIT2 ~~ AUDIT3
                                 AUDIT1 ~~ AUDIT2'
fittedGlobalExtBehaviors2 <- sem(globalExtBehaviors2, dataE, meanstructure = TRUE,
                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                 orthogonal = FALSE, missing = "ML",
                                 bootstrap = 5000)
fitmeasures(fittedGlobalExtBehaviors2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalExtBehaviors2, "all")
summary(fittedGlobalExtBehaviors2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalExtBehaviors2)$cov
residuals(fittedGlobalExtBehaviors2)$cov > .10
globalExtBehaviors2MI <- modindices(fittedGlobalExtBehaviors2)[order(-modindices(fittedGlobalExtBehaviors2)$mi), ]
globalExtBehaviors2MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 184.825  0.444   0.444    0.605    0.605
# 1962          PCL20 ~~   QIDS1  84.094  0.265   0.265    0.400    0.400
# 873           audit =~  risky9  46.301 -0.132  -0.132   -0.319   -0.319
# 464    ExtBehaviors =~    PCL2  40.172  0.475   0.475    0.510    0.510
# 574      AnxArousal =~  risky6  37.115  0.190   0.190    0.259    0.259
# 2481         QIDS12 ~~ risky14  36.233  0.035   0.035    0.236    0.236
# 281  Reexperiencing =~   PCL11  35.512  0.764   0.764    0.647    0.647
# 639    DysphArousal =~  risky6  33.716  0.183   0.183    0.250    0.250
# 2816         risky1 ~~  risky9  32.994 -0.049  -0.049   -0.260   -0.260
# 2934         risky7 ~~ risky11  32.002  0.019   0.019    0.234    0.234
# 750             gad =~  risky6  31.817  0.176   0.176    0.240    0.240
# 690            qids =~  risky6  31.300  0.174   0.174    0.238    0.238


# # Global Hybrid model -----------------------------------------------------
# 
# globalHybrid <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
#                  Avoidance =~ PCL6 + PCL7
#                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
#                  Anhedonia =~ PCL12 + PCL13 + PCL14
#                  ExtBehavior =~ PCL15 + PCL16 
#                  AnxArousal =~ + PCL17 + PCL18
#                  DysphArousal =~ PCL19 + PCL20
#                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 + 
#                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
#                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
#                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
#                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
#                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
#                           AUDIT8 + AUDIT9 + AUDIT10'
# fittedGlobalHybrid <- sem(globalHybrid, dataE, meanstructure = TRUE,
#                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
#                           orthogonal = FALSE, missing = "ML",
#                           bootstrap = 5000)
# fitmeasures(fittedGlobalHybrid,
#             c("chisq", "df", "pvalue", "cfi", "tli",
#               "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
# fitmeasures(fittedGlobalHybrid, "all")
# summary(fittedGlobalHybrid, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
# residuals(fittedGlobalHybrid)$cov
# residuals(fittedGlobalHybrid)$cov > .10
# globalHybridMI <- modindices(fittedGlobalHybrid)[order(-modindices(fittedGlobalHybrid)$mi), ]
# globalHybridMI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 935           audit =~  risky1 188.103  0.439   0.439    0.598    0.598
# 3124         AUDIT2 ~~  AUDIT3  94.646  0.245   0.245    0.452    0.452
# 2032          PCL20 ~~   QIDS1  83.564  0.263   0.263    0.397    0.397
# 3094        risky13 ~~ risky14  82.772  0.048   0.048    0.408    0.408
# 2554         QIDS12 ~~ risky14  64.417  0.050   0.050    0.348    0.348
# 2300          QIDS6 ~~   QIDS8  59.015  0.087   0.087    0.318    0.318
# 2341          QIDS7 ~~   QIDS9  47.519  0.122   0.122    0.283    0.283
# 943           audit =~  risky9  40.862 -0.121  -0.121   -0.292   -0.292
# 2553         QIDS12 ~~ risky13  40.330  0.041   0.041    0.273    0.273
# 644      AnxArousal =~  risky6  39.734  0.193   0.193    0.264    0.264
# 534     ExtBehavior =~    PCL2  39.693  0.469   0.469    0.504    0.504
# 2911         risky2 ~~  risky7  38.146  0.023   0.023    0.271    0.271


globalHybrid2 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2'
fittedGlobalHybrid2 <- sem(globalHybrid2, dataE, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybrid2,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybrid2, "all")
summary(fittedGlobalHybrid2, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalHybrid2)$cov
residuals(fittedGlobalHybrid2)$cov > .10
globalHybrid2MI <- modindices(fittedGlobalHybrid2)[order(-modindices(fittedGlobalHybrid2)$mi), ]
globalHybrid2MI
#                 lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 944           audit =~  risky1 185.282  0.445   0.445    0.606    0.606
# 2041          PCL20 ~~   QIDS1  84.283  0.265   0.265    0.399    0.399
# 952           audit =~  risky9  46.808 -0.132  -0.132   -0.320   -0.320
# 543     ExtBehavior =~    PCL2  38.472  0.458   0.458    0.492    0.492
# 653      AnxArousal =~  risky6  37.197  0.190   0.190    0.259    0.259
# 2560         QIDS12 ~~ risky14  36.154  0.035   0.035    0.235    0.235
# 718    DysphArousal =~  risky6  34.711  0.185   0.185    0.253    0.253
# 2895         risky1 ~~  risky9  33.003 -0.049  -0.049   -0.260   -0.260
# 3013         risky7 ~~ risky11  31.964  0.019   0.019    0.233    0.233
# 829             gad =~  risky6  31.617  0.175   0.175    0.239    0.239
# 769            qids =~  risky6  31.525  0.175   0.175    0.239    0.239
# 3044         risky9 ~~ risky11  31.241  0.026   0.026    0.247    0.247

# MIs for global models ---------------------------------------------------

##Same order in all models
# 1. audit =~ risky1; risky1 = problematic alcohol use; makes sense
# 2. AUDIT2 ~~ AUDIT3; AUDIT2 = how many drinks, AUDIT3 = how often more than 6 drinks; makes sense
# 3. PCL20 ~~ QIDS1; PCL20 = trouble falling/staying asleep, QIDS1 = problems with falling asleep; makes sense
# 4. risky13 ~~ risky14; risky13 = deliberately injuring, risky14 = suicidal behaviors; makes sense
# 5. QIDS12 ~~ risky14; QIDS12 = thoughts on death/suicide; risky14 = suicidal behaviors; makes sense 
# 6. QIDS6 ~~ QIDS8; QIDS6 = decreased appetite, QIDS8 = decreased weight; makes sense
# 7. QIDS7 ~~ QIDS9; QIDS7 = increased appetite, QIDS9 = increased weight; makes sense
##Slightly different order amongst models
# 8. QIDS12 ~~ risky13; QIDS12 = thoughts on death/suicide, risky13 = deliberately injuring; makes sense
# 9. audit =~ risky9; risky9 = physically aggressive behaviors; not so much sense
# 10. ExtBehavior =~ PCL2; PCL2 = dreams about the stressful experience; does not make sense
# 11. AnxArousal =~ risky6; risky6 = problematic eating; does not make sense
# 12. risky2 ~~ risky7; risky2 = drug use, risky7 = illegal behavior; makes sense
# 13. qids =~ risky7; risky7 = illegal behavior; does not make sense
# 14. AltArousal =~ risky6; risky6 = problematic eating; does not make sense



# Models comparison -------------------------------------------------------

##List of nested and non-nested CFA models

# NESTED
# DSM-5 vs. Three-factor
# DSM-5 vs. Dysphoric arousal
# DSM-5 vs. Anhedonia
# DSM-5 vs. Externalizing behaviours
# Three-factor vs. Dysphoric arousal
# Three-factor vs. Anhedonia
# Three-factor vs. Externalizing behaviours
# Dysphoria vs. Dysphoric arousal
# Dysphoria vs. Anhedonia
# Dysphoria vs. Externalizing behaviours
# Dysphoric arousal vs. Anhedonia
# Dysphoric arousal vs. Externalizing Behaviors
# 1F vs ALL
# Hybrid vs ALL

# NON-NESTED
# DSM-5 vs. Dysphoria
# Three-factor vs. Dysphoria
# Anhedonia vs. Externalizing behaviours

##LRT for testing nested models
##AIC/BIC for testing non-nested models


# One-factor - 7 comparisons
# DSM-5 - 7 comparisons
# Three-factor - 7 comparisons
# Dysphoria - 7 comparisons
# Dysphoric Arousal - 7 comparisons
# Externalizing behaviors - 7 comparisons
# Hybrid - 7 comparisons


# Nested models comparisons -----------------------------------------------

##One-factor comparisons 

anova(fittedOneFactor, fittedDsm5, fittedThreeFactor, fittedDsm5dysphoria, 
      fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors, fittedHybridModel)
anova(fittedOneFactor, fittedDsm5)
anova(fittedOneFactor, fittedThreeFactor)
anova(fittedOneFactor, fittedDsm5dysphoria)
anova(fittedOneFactor, fittedDsm5dysphArousal)
anova(fittedOneFactor, fittedAnhedonia)
anova(fittedOneFactor, fittedExtBehaviors)
anova(fittedOneFactor, fittedHybridModel)


##DSM-5 comparisons 

anova(fittedDsm5, fittedThreeFactor, fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors)
anova(fittedDsm5, fittedThreeFactor)
anova(fittedDsm5, fittedDsm5dysphArousal)
anova(fittedDsm5, fittedAnhedonia)
anova(fittedDsm5, fittedExtBehaviors)

##Three-factor comparisons 

anova(fittedThreeFactor, fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors)
anova(fittedThreeFactor, fittedDsm5dysphArousal)
anova(fittedThreeFactor, fittedAnhedonia)
anova(fittedThreeFactor, fittedExtBehaviors)

##Dysphoria comparisons 

anova(fittedDsm5dysphoria, fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors)
anova(fittedDsm5dysphoria, fittedDsm5dysphArousal)
anova(fittedDsm5dysphoria, fittedAnhedonia)
anova(fittedDsm5dysphoria, fittedExtBehaviors)

##Dysphoric arousal comparisons 

anova(fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors)
anova(fittedDsm5dysphArousal, fittedAnhedonia)
anova(fittedDsm5dysphArousal, fittedExtBehaviors)

##Hybrid model comparisons

anova(fittedHybridModel, fittedOneFactor, fittedDsm5, fittedThreeFactor, 
      fittedDsm5dysphoria, fittedDsm5dysphArousal, fittedAnhedonia, fittedExtBehaviors)
anova(fittedHybridModel, fittedOneFactor)
anova(fittedHybridModel, fittedDsm5)
anova(fittedHybridModel, fittedThreeFactor)
anova(fittedHybridModel, fittedDsm5dysphoria)
anova(fittedHybridModel, fittedDsm5dysphArousal)
anova(fittedHybridModel, fittedAnhedonia)
anova(fittedHybridModel, fittedExtBehaviors)

# Non-nested models comparisons -------------------------------------------

##DSM-5 vs. Dysphoria
anova(fittedDsm5, fittedDsm5dysphoria)
lrtest(fittedDsm5, fittedDsm5dysphoria)

##Three-factor vs. Dysphoria

anova(fittedThreeFactor, fittedDsm5dysphoria)
lrtest(fittedThreeFactor, fittedDsm5dysphoria)

##Anhedonia vs. Externalizing behaviors

anova(fittedAnhedonia, fittedExtBehaviors)
lrtest(fittedAnhedonia, fittedExtBehaviors)


# A summary of the comparison ---------------------------------------------

# One-factor: W (-) L (DSM-5, Three-factor, Dysphoria, Dysphoric Arousal, Anhedonia, Externalizing Behaviors, Hybrid)
# DSM-5: W (One-factor, Three-factor) L (Anhedonia, Externalizing Behaviors)
# Three-factor: W (One-factor) L (DSM-5, DSM-5 Dysphoricl Arousal, Externalizing Behaviors)
# Dysphoria: W (One-factor) L (Dysphoric Arousal, Externalizing Behaviors, Anhedonia)
# Dysphoric Arousal: W (One-factor, Three-factor) L (Anhedonia, Externalizing behaviors)
# Hybrid: W (One-factor, DSM-5, Three-factor, Dysphoria, Dysphoric Arousal, Anhedonia, Externalizing Behaviors)


# Final ranking -----------------------------------------------------------

# 1. Hybrid model
# 2. Anhedonia
# 3. Externalizing behavior


# Final MIs ---------------------------------------------------------------

##Hybrid

# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 944           audit =~  risky1 185.282  0.445   0.445    0.606    0.606
# 2041          PCL20 ~~   QIDS1  84.283  0.265   0.265    0.399    0.399
# 952           audit =~  risky9  46.808 -0.132  -0.132   -0.320   -0.320
# 543     ExtBehavior =~    PCL2  38.472  0.458   0.458    0.492    0.492
# 653      AnxArousal =~  risky6  37.197  0.190   0.190    0.259    0.259
# 2560         QIDS12 ~~ risky14  36.154  0.035   0.035    0.235    0.235
# 718    DysphArousal =~  risky6  34.711  0.185   0.185    0.253    0.253
# 2895         risky1 ~~  risky9  33.003 -0.049  -0.049   -0.260   -0.260
# 3013         risky7 ~~ risky11  31.964  0.019   0.019    0.233    0.233
# 829             gad =~  risky6  31.617  0.175   0.175    0.239    0.239
# 769            qids =~  risky6  31.525  0.175   0.175    0.239    0.239
# 3044         risky9 ~~ risky11  31.241  0.026   0.026    0.247    0.247
# 459       NegAffect =~  risky6  30.738  0.176   0.176    0.240    0.240
# 2849           GAD6 ~~ risky10  29.563  0.085   0.085    0.224    0.224

##Externalizing behavior
# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 184.825  0.444   0.444    0.605    0.605
# 1962          PCL20 ~~   QIDS1  84.094  0.265   0.265    0.400    0.400
# 873           audit =~  risky9  46.301 -0.132  -0.132   -0.319   -0.319
# 464    ExtBehaviors =~    PCL2  40.172  0.475   0.475    0.510    0.510
# 574      AnxArousal =~  risky6  37.115  0.190   0.190    0.259    0.259
# 2481         QIDS12 ~~ risky14  36.233  0.035   0.035    0.236    0.236
# 281  Reexperiencing =~   PCL11  35.512  0.764   0.764    0.647    0.647
# 639    DysphArousal =~  risky6  33.716  0.183   0.183    0.250    0.250
# 2816         risky1 ~~  risky9  32.994 -0.049  -0.049   -0.260   -0.260
# 2934         risky7 ~~ risky11  32.002  0.019   0.019    0.234    0.234
# 750             gad =~  risky6  31.817  0.176   0.176    0.240    0.240
# 690            qids =~  risky6  31.300  0.174   0.174    0.238    0.238
# 2965         risky9 ~~ risky11  31.259  0.026   0.026    0.247    0.247
# 513    ExtBehaviors =~ risky10  30.172  0.232   0.232    0.298    0.298

##Anhedonia

# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 185.054  0.444   0.444    0.605    0.605
# 1867          PCL20 ~~   QIDS1  84.890  0.269   0.269    0.390    0.390
# 873           audit =~  risky9  45.168 -0.131  -0.131   -0.316   -0.316
# 2481         QIDS12 ~~ risky14  37.470  0.035   0.035    0.240    0.240
# 639      AnxArousal =~  risky6  37.227  0.190   0.190    0.259    0.259
# 574    DysphArousal =~  risky6  34.758  0.196   0.196    0.268    0.268
# 2965         risky9 ~~ risky11  33.248  0.027   0.027    0.255    0.255
# 2934         risky7 ~~ risky11  32.567  0.019   0.019    0.236    0.236
# 750             gad =~  risky6  31.561  0.175   0.175    0.239    0.239
# 2816         risky1 ~~  risky9  31.541 -0.048  -0.048   -0.254   -0.254
# 690            qids =~  risky6  31.468  0.175   0.175    0.240    0.240
# 447       NegAffect =~  risky6  30.746  0.176   0.176    0.240    0.240
# 2770           GAD6 ~~ risky10  30.296  0.087   0.087    0.226    0.226
# 853           audit =~  QIDS12  29.325  0.112   0.112    0.223    0.223




# Final models - exploratory dataset --------------------------------------

globalHybrid3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1'
fittedGlobalHybrid3 <- sem(globalHybrid3, dataE, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybrid3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybrid3, "all")
summary(fittedGlobalHybrid3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalHybrid3)$cov
residuals(fittedGlobalHybrid3)$cov > .10
globalHybrid3MI <- modindices(fittedGlobalHybrid3)[order(-modindices(fittedGlobalHybrid3)$mi), ]
globalHybrid3MI


globalAnhedonia3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                     Anhedonia =~ PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2
                     audit =~ risky1
                     audit =~ risky9
                     PCL20 ~~ QIDS1'
fittedGlobalAnhedonia3 <- sem(globalAnhedonia3, dataE, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalAnhedonia3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalAnhedonia3, "all")
summary(fittedGlobalAnhedonia3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalAnhedonia3)$cov
residuals(fittedGlobalAnhedonia3)$cov > .10
globalAnhedonia3MI <- modindices(fittedGlobalAnhedonia3)[order(-modindices(fittedGlobalAnhedonia3)$mi), ]
globalAnhedonia3MI


globalExtBehaviors3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                        Avoidance =~ PCL6 + PCL7
                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                        ExtBehaviors =~ PCL15 + PCL16
                        AnxArousal =~ PCL17 + PCL18
                        DysphArousal =~ PCL19 + PCL20
                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                QIDS6 ~~ QIDS8
                                QIDS7 ~~ QIDS9
                                QIDS2 ~~ QIDS3
                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                               GAD1 ~~ GAD6
                               GAD5 ~~ GAD7
                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                 risky13 ~~ risky14
                                 risky2 ~~ risky7
                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                 AUDIT8 + AUDIT9 + AUDIT10
                                 AUDIT2 ~~ AUDIT3
                                 AUDIT1 ~~ AUDIT2
                        audit =~ risky1
                        audit =~ risky9
                        PCL20 ~~ QIDS1'
fittedGlobalExtBehaviors3 <- sem(globalExtBehaviors3, dataE, meanstructure = TRUE,
                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                 orthogonal = FALSE, missing = "ML",
                                 bootstrap = 5000)
fitmeasures(fittedGlobalExtBehaviors3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalExtBehaviors3, "all")
summary(fittedGlobalExtBehaviors3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalExtBehaviors3)$cov
residuals(fittedGlobalExtBehaviors3)$cov > .10
globalExtBehaviors3MI <- modindices(fittedGlobalExtBehaviors3)[order(-modindices(fittedGlobalExtBehaviors3)$mi), ]
globalExtBehaviors3MI



# Confirmatory dataset ----------------------------------------------------

oneFactorC <- 'ptsd =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5 + PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + 
                      PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedOneFactorC <- cfa(oneFactorC, dataC, meanstructure = TRUE,
                        std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                        orthogonal = FALSE, missing = "ML",
                        bootstrap = 5000)
fitmeasures(fittedOneFactorC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedOneFactorC, "all")
summary(fittedOneFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedOneFactorC)$cov
residuals(fittedOneFactorC)$cov > .10
modindices(fittedOneFactorC)[order(-modindices(fittedOneFactorC)$mi), ]

# DSM-5 model -------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]; 
# 2. Avoidance [C1, C2]; 
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]; 
# 4. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])


dsm5C <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
         Avoidance =~ PCL6 + PCL7
         NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
         AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedDsm5C <- cfa(dsm5C, dataC, meanstructure = TRUE,
                   std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                   orthogonal = FALSE, missing = "ML",
                   bootstrap = 5000)
fitmeasures(fittedDsm5C, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedDsm5C, "all")
summary(fittedDsm5C, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5C)$cov
residuals(fittedDsm5C)$cov > .10
modindices(fittedDsm5C)[order(-modindices(fittedDsm5C)$mi), ]


# Three-factor model ------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance/Negative alternation in cognition and mood [C1, C2, D1, D2, D3, D4, D5, D6, D7]
# 3. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])

threeFactorC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                AvoidNegAlternation =~ PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedThreeFactorC <- cfa(threeFactorC, dataC, meanstructure = TRUE,
                          std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                          orthogonal = FALSE, missing = "ML",
                          bootstrap = 5000)
fitmeasures(fittedThreeFactorC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedThreeFactorC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedThreeFactorC)$cov
residuals(fittedThreeFactorC)$cov > .10
modindices(fittedThreeFactorC)[order(-modindices(fittedThreeFactorC)$mi), ]


# DSM-5 Dysphoria model ---------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Dyshoria [D1, D2, D3, D4, D5, D6, D7, E1, E2, E5, E6]
# 4. Alternation in arousal and reactivity [E3, E4])

dsm5dysphoriaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL19 + PCL20
                  AltArousal =~ PCL17 + PCL18'
fittedDsm5dysphoriaC <- cfa(dsm5dysphoria, dataC, meanstructure = TRUE,
                            std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                            orthogonal = FALSE, missing = "ML",
                            bootstrap = 5000)
fitmeasures(fittedDsm5dysphoriaC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphoriaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphoriaC)$cov
residuals(fittedDsm5dysphoriaC)$cov > .10
modindices(fittedDsm5dysphoriaC)[order(-modindices(fittedDsm5dysphoriaC)$mi), ]


# DSM-5 Dysphoric arousal model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]
# 4. Dysphoric arousal [E1, E2, E5, E6]
# 5. Anxious arousal [E3, E4])

dsm5dysphAruosalC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18'
fittedDsm5dysphArousalC <- cfa(dsm5dysphAruosalC, dataC, meanstructure = TRUE,
                               std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                               orthogonal = FALSE, missing = "ML",
                               bootstrap = 5000)
fitmeasures(fittedDsm5dysphArousalC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphArousalC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphArousalC)$cov
residuals(fittedDsm5dysphArousalC)$cov > .10
modindices(fittedDsm5dysphArousalC)[order(-modindices(fittedDsm5dysphArousalC)$mi), ]

# Anhedonia model ---------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Dysphoric arousal [E1, E2, E5, E6]
# 6. Anxious arousal [E3, E4])

anhedoniaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
              Avoidance =~ PCL6 + PCL7
              NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
              Anhedonia =~ PCL12 + PCL13 + PCL14
              DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
              AnxArousal =~ PCL17 + PCL18'
fittedAnhedoniaC <- cfa(anhedoniaC, dataC, meanstructure = TRUE,
                        std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                        orthogonal = FALSE, missing = "ML",
                        bootstrap = 5000)
fitmeasures(fittedAnhedoniaC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedAnhedoniaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedAnhedoniaC)$cov
residuals(fittedAnhedoniaC)$cov > .10
modindices(fittedAnhedoniaC)[order(-modindices(fittedAnhedoniaC)$mi), ]


# Externalizing behaviors model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Numbing [D1, D2, D3, D4, D5, D6, D7]
# 4. Externalizing behaviors [E1, E2]
# 5. Anxious arousal [E3, E4]
# 6. Dysphoric arousal [E5, E6])

extBehaviorsC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                 Avoidance =~ PCL6 + PCL7
                 Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                 ExtBehaviors =~ PCL15 + PCL16
                 AnxArousal =~ PCL17 + PCL18 
                 DysphArousal =~ PCL19 + PCL20'
fittedExtBehaviorsC <- cfa(extBehaviorsC, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedExtBehaviorsC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedExtBehaviorsC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedExtBehaviorsC)$cov
residuals(fittedExtBehaviorsC)$cov > .10
modindices(fittedExtBehaviorsC)[order(-modindices(fittedExtBehaviorsC)$mi), ]


# Hybrid model ------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Externalizing behaviours [E1, E2]
# 6. Anxious arousal [E3, E4]
# 7. Dysphoric arousal [E5, E6])

hybridModelC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                Avoidance =~ PCL6 + PCL7
                NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
                Anhedonia =~ PCL12 + PCL13 + PCL14
                ExtBehavior =~ PCL15 + PCL16 
                AnxArousal =~ + PCL17 + PCL18
                DysphArousal =~ PCL19 + PCL20'
fittedHybridModelC <- cfa(hybridModelC, dataC, meanstructure = TRUE,
                          std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                          orthogonal = FALSE, missing = "ML",
                          bootstrap = 5000)
fitmeasures(fittedHybridModelC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedHybridModelC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedHybridModelC)$cov
residuals(fittedHybridModelC)$cov > .10
modindices(fittedHybridModelC)[order(-modindices(fittedHybridModelC)$mi), ]


# MIs for CFAs ------------------------------------------------------------

##Apart from the one-factor model
# PCL12 ~~ PCL14
# PCL15 ~~ PCL16
# PCL2 ~~ PCL15


# Fit measures summary ----------------------------------------------------

model <- c("One-FactorC model", "DSM-5 model", "Three-FactorC model", "DSM-5 Dysphoria model", "DSM-5 Dysphoric arousal model",
           "Anhedonia model", "Externalizing behaviors model", "Hybrid model")
chisq <- c(fitmeasures(fittedOneFactorC, "chisq"), fitmeasures(fittedDsm5C, "chisq"), 
           fitmeasures(fittedThreeFactorC, "chisq"), fitmeasures(fittedDsm5dysphoriaC, "chisq"), 
           fitmeasures(fittedDsm5dysphArousalC, "chisq"), fitmeasures(fittedAnhedoniaC, "chisq"), 
           fitmeasures(fittedExtBehaviorsC, "chisq"), fitmeasures(fittedHybridModelC, "chisq")) 
chisq <- round(chisq, 3)
df <- c(fitmeasures(fittedOneFactorC, "df"), fitmeasures(fittedDsm5C, "df"), 
        fitmeasures(fittedThreeFactorC, "df"), fitmeasures(fittedDsm5dysphoriaC, "df"), 
        fitmeasures(fittedDsm5dysphArousalC, "df"), fitmeasures(fittedAnhedoniaC, "df"), 
        fitmeasures(fittedExtBehaviorsC, "df"), fitmeasures(fittedHybridModelC, "df")) 
df <- round(df, 3)
pvalue <- c(fitmeasures(fittedOneFactorC, "pvalue"), fitmeasures(fittedDsm5C, "pvalue"), 
            fitmeasures(fittedThreeFactorC, "pvalue"), fitmeasures(fittedDsm5dysphoriaC, "pvalue"),
            fitmeasures(fittedDsm5dysphArousal, "pvalue"), fitmeasures(fittedAnhedoniaC, "pvalue"), 
            fitmeasures(fittedExtBehaviorsC, "pvalue"), fitmeasures(fittedHybridModelC, "pvalue")) 
pvalue <- round(pvalue, 3)
cfi <- c(fitmeasures(fittedOneFactorC, "cfi"), fitmeasures(fittedDsm5, "cfi"), 
         fitmeasures(fittedThreeFactorC, "cfi"), fitmeasures(fittedDsm5dysphoriaC, "cfi"), 
         fitmeasures(fittedDsm5dysphArousal, "cfi"), fitmeasures(fittedAnhedoniaC, "cfi"), 
         fitmeasures(fittedExtBehaviorsC, "cfi"), fitmeasures(fittedHybridModelC, "cfi")) 
cfi <- round(cfi, 3)
tli <- c(fitmeasures(fittedOneFactorC, "tli"), fitmeasures(fittedDsm5, "tli"), 
         fitmeasures(fittedThreeFactorC, "tli"), fitmeasures(fittedDsm5dysphoriaC, "tli"), 
         fitmeasures(fittedDsm5dysphArousalC, "tli"), fitmeasures(fittedAnhedoniaC, "tli"), 
         fitmeasures(fittedExtBehaviorsC, "tli"), fitmeasures(fittedHybridModelC, "tli")) 
tli <- round(tli, 3)
rmsea <- c(fitmeasures(fittedOneFactorC, "rmsea"), fitmeasures(fittedDsm5C, "rmsea"), 
           fitmeasures(fittedThreeFactorC, "rmsea"), fitmeasures(fittedDsm5dysphoriaC, "rmsea"), 
           fitmeasures(fittedDsm5dysphArousalC, "rmsea"), fitmeasures(fittedAnhedoniaC, "rmsea"), 
           fitmeasures(fittedExtBehaviorsC, "rmsea"), fitmeasures(fittedHybridModelC, "rmsea"))
rmsea <- round(rmsea, 3)
rmseaLI <- c(fitmeasures(fittedOneFactorC, "rmsea.ci.lower"), fitmeasures(fittedDsm5, "rmsea.ci.lower"), 
             fitmeasures(fittedThreeFactorC, "rmsea.ci.lower"), fitmeasures(fittedDsm5dysphoriaC, "rmsea.ci.lower"), 
             fitmeasures(fittedDsm5dysphArousalC, "rmsea.ci.lower"), fitmeasures(fittedAnhedoniaC, "rmsea.ci.lower"),
             fitmeasures(fittedExtBehaviorsC, "rmsea.ci.lower"), fitmeasures(fittedHybridModelC, "rmsea.ci.lower")) 
rmseaLI <- round(rmseaLI, 3)
rmseaUI<- c(fitmeasures(fittedOneFactorC, "rmsea.ci.upper"), fitmeasures(fittedDsm5, "rmsea.ci.upper"), 
            fitmeasures(fittedThreeFactorC, "rmsea.ci.upper"), fitmeasures(fittedDsm5dysphoriaC, "rmsea.ci.upper"), 
            fitmeasures(fittedDsm5dysphArousalC, "rmsea.ci.upper"), fitmeasures(fittedAnhedoniaC, "rmsea.ci.upper"), 
            fitmeasures(fittedExtBehaviorsC, "rmsea.ci.upper"), fitmeasures(fittedHybridModelC, "rmsea.ci.upper")) 
rmseaUI <- round(rmseaUI, 3)
srmr <- c(fitmeasures(fittedOneFactorC, "srmr"), fitmeasures(fittedDsm5, "srmr"), 
          fitmeasures(fittedThreeFactorC, "srmr"), fitmeasures(fittedDsm5dysphoriaC, "srmr"), 
          fitmeasures(fittedDsm5dysphArousalC, "srmr"), fitmeasures(fittedAnhedoniaC, "srmr"), 
          fitmeasures(fittedExtBehaviorsC, "srmr"), fitmeasures(fittedHybridModelC, "srmr")) 
srmr <- round(srmr, 3)
aic <- c(fitmeasures(fittedOneFactorC, "aic"), fitmeasures(fittedDsm5C, "aic"), 
         fitmeasures(fittedThreeFactorC, "aic"), fitmeasures(fittedDsm5dysphoriaC, "aic"), 
         fitmeasures(fittedDsm5dysphArousalC, "aic"), fitmeasures(fittedAnhedoniaC, "aic"), 
         fitmeasures(fittedExtBehaviorsC, "aic"), fitmeasures(fittedHybridModelC, "aic")) 
aic <- round(aic, 3)
bic <- c(fitmeasures(fittedOneFactorC, "bic"), fitmeasures(fittedDsm5C, "bic"), 
         fitmeasures(fittedThreeFactorC, "bic"), fitmeasures(fittedDsm5dysphoriaC, "bic"), 
         fitmeasures(fittedDsm5dysphArousalC, "bic"), fitmeasures(fittedAnhedoniaC, "bic"), 
         fitmeasures(fittedExtBehaviorsC, "bic"), fitmeasures(fittedHybridModelC, "bic")) 
bic <- round(bic, 3)

modelFits <- data.frame(model, chisq, df, pvalue, cfi, tli, rmsea, rmseaLI, rmseaUI, srmr, aic, bic)

kbl(modelFits)
modelFits %>% kbl() %>% kable_material()


# Nested models comparisons -----------------------------------------------

anova(fittedOneFactorC, fittedDsm5C, fittedThreeFactorC, fittedDsm5dysphoriaC, 
      fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC, fittedHybridModelC)
anova(fittedOneFactorC, fittedDsm5C)
anova(fittedOneFactorC, fittedThreeFactorC)
anova(fittedOneFactorC, fittedDsm5dysphoriaC)
anova(fittedOneFactorC, fittedDsm5dysphArousalC)
anova(fittedOneFactorC, fittedAnhedoniaC)
anova(fittedOneFactorC, fittedExtBehaviorsC)
anova(fittedOneFactorC, fittedHybridModelC)


##DSM-5 comparisons 

anova(fittedDsm5C, fittedThreeFactorC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5C, fittedThreeFactorC)
anova(fittedDsm5C, fittedDsm5dysphArousalC)
anova(fittedDsm5C, fittedAnhedoniaC)
anova(fittedDsm5C, fittedExtBehaviorsC)

##Three-factor comparisons 

anova(fittedThreeFactorC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedThreeFactorC, fittedDsm5dysphArousalC)
anova(fittedThreeFactorC, fittedAnhedoniaC)
anova(fittedThreeFactorC, fittedExtBehaviorsC)

##Dysphoria comparisons 

anova(fittedDsm5dysphoriaC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5dysphoriaC, fittedDsm5dysphArousalC)
anova(fittedDsm5dysphoriaC, fittedAnhedoniaC)
anova(fittedDsm5dysphoriaC, fittedExtBehaviorsC)

##Dysphoric arousal comparisons 

anova(fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5dysphArousalC, fittedAnhedoniaC)
anova(fittedDsm5dysphArousalC, fittedExtBehaviorsC)

##Hybrid model comparisons

anova(fittedHybridModelC, fittedOneFactorC, fittedDsm5C, fittedThreeFactorC, 
      fittedDsm5dysphoriaC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedHybridModelC, fittedOneFactorC)
anova(fittedHybridModelC, fittedDsm5C)
anova(fittedHybridModelC, fittedThreeFactorC)
anova(fittedHybridModelC, fittedDsm5dysphoriaC)
anova(fittedHybridModelC, fittedDsm5dysphArousalC)
anova(fittedHybridModelC, fittedAnhedoniaC)
anova(fittedHybridModelC, fittedExtBehaviorsC)

# Non-nested models comparisons -------------------------------------------

##DSM-5 vs. Dysphoria
anova(fittedDsm5C, fittedDsm5dysphoriaC)
lrtest(fittedDsm5C, fittedDsm5dysphoriaC)

##Three-factor vs. Dysphoria

anova(fittedThreeFactorC, fittedDsm5dysphoriaC)
lrtest(fittedThreeFactorC, fittedDsm5dysphoriaC)

##Anhedonia vs. Externalizing behaviors

anova(fittedAnhedoniaC, fittedExtBehaviorsC)
lrtest(fittedAnhedoniaC, fittedExtBehaviorsC)





# Final models - confirmatory dataset -------------------------------------

globalHybridC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1'
fittedGlobalHybridC <- sem(globalHybridC, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybridC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybridC, "all")
summary(fittedGlobalHybridC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalHybridC)$cov
residuals(fittedGlobalHybridC)$cov > .10
globalHybridCMI <- modindices(fittedGlobalHybridC)[order(-modindices(fittedGlobalHybridC)$mi), ]
globalHybridCMI


globalAnhedoniaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                     Anhedonia =~ PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2
                     audit =~ risky1
                     audit =~ risky9
                     PCL20 ~~ QIDS1'
fittedGlobalAnhedoniaC <- sem(globalAnhedoniaC, dataC, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalAnhedoniaC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalAnhedoniaC, "all")
summary(fittedGlobalAnhedoniaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalAnhedoniaC)$cov
residuals(fittedGlobalAnhedoniaC)$cov > .10
globalAnhedoniaCMI <- modindices(fittedGlobalAnhedoniaC)[order(-modindices(fittedGlobalAnhedoniaC)$mi), ]
globalAnhedoniaCMI


globalExtBehaviorsC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                        Avoidance =~ PCL6 + PCL7
                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                        ExtBehaviors =~ PCL15 + PCL16
                        AnxArousal =~ PCL17 + PCL18
                        DysphArousal =~ PCL19 + PCL20
                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                QIDS6 ~~ QIDS8
                                QIDS7 ~~ QIDS9
                                QIDS2 ~~ QIDS3
                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                               GAD1 ~~ GAD6
                               GAD5 ~~ GAD7
                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                 risky13 ~~ risky14
                                 risky2 ~~ risky7
                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                 AUDIT8 + AUDIT9 + AUDIT10
                                 AUDIT2 ~~ AUDIT3
                                 AUDIT1 ~~ AUDIT2
                        audit =~ risky1
                        audit =~ risky9
                        PCL20 ~~ QIDS1'
fittedGlobalExtBehaviorsC <- sem(globalExtBehaviorsC, dataC, meanstructure = TRUE,
                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                 orthogonal = FALSE, missing = "ML",
                                 bootstrap = 5000)
fitmeasures(fittedGlobalExtBehaviorsC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalExtBehaviorsC, "all")
summary(fittedGlobalExtBehaviorsC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalExtBehaviorsC)$cov
residuals(fittedGlobalExtBehaviorsC)$cov > .10
globalExtBehaviorsCMI <- modindices(fittedGlobalExtBehaviorsC)[order(-modindices(fittedGlobalExtBehaviorsC)$mi), ]
globalExtBehaviorsCMI









# Equality of regression coefficients -------------------------------------

##Test equality of regression coefficients using Wald test with lavaan (Klopp, 2020)
lwt <- function(object, constraints = NULL, verbose = FALSE, std = FALSE) {
  
  if(object@optim$npar > 0L && !object@optim$converged)
    stop("lavaan ERROR: model did not converge")
  
  if(is.null(constraints) || nchar(constraints) == 0L) {
    stop("lavaan ERROR: constraints are empty")
  }
  
  # remove == constraints from parTable
  PT <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
  eq.idx <- which(PT$op == "==")
  if(length(eq.idx) > 0L) {
    PT <- PT[-eq.idx,]
  }
  partable <- as.list(PT)
  
  # parse constraints
  FLAT <- lavParseModelString( constraints ); CON <- attr(FLAT, "constraints")
  LIST <- list()
  if(length(CON) > 0L) {
    lhs = unlist(lapply(CON, "[[", "lhs"))
    op = unlist(lapply(CON, "[[",  "op"))
    rhs = unlist(lapply(CON, "[[", "rhs"))
    LIST$lhs        <- c(LIST$lhs,        lhs)
    LIST$op         <- c(LIST$op,         op)
    LIST$rhs        <- c(LIST$rhs,        rhs)
  } else {
    stop("lavaan ERROR: no equality constraints found in constraints argument")
  }
  
  # theta = free parameters only
  # changed: to include the std argument
  if (!std) {
    theta <- object@optim$x
  }
  
  if(std) {
    theta <- lavaan:::lav_standardize_all(object)
    fixed <- which(unclass(partable(object))$free %in% c(0))
    if (length(fixed) !=0 ) theta <- theta[-fixed]
  }
  
  # build constraint function
  ceq.function <- lav_partable_constraints_ceq(partable = partable,
                                               con = LIST, debug = FALSE)
  
  # compute jacobian restrictions
  JAC <- try(lav_func_jacobian_complex(func = ceq.function, x = theta), silent=TRUE)
  if(inherits(JAC, "try-error")) { # eg. pnorm()
    JAC <- lav_func_jacobian_simple(func = ceq.function, x = theta)
  }
  
  if(verbose) {
    cat("Restriction matrix (jacobian):\n"); print(JAC); cat("\n")
  }
  
  # linear restriction
  theta.r <- ceq.function( theta )
  
  if(verbose) {
    cat("Restricted theta values:\n"); print(theta.r); cat("\n")
  }
  
  # get VCOV
  # VCOV <- vcov(object, labels = FALSE)
  # avoid S4 dispatch
  # changed: to include the std argument
  if (!std) {
    # changed: addad lavaan::: because function is not exported
    VCOV <- lavaan:::lav_object_inspect_vcov(object, standardized = FALSE,
                                             free.only = TRUE,
                                             add.labels = FALSE,
                                             add.class = FALSE,
                                             remove.duplicated = FALSE)
  }
  
  if (std) {
    # changed: addad lavaan::: because function is not exported
    VCOV <- lavaan:::lav_object_inspect_vcov(object, standardized = TRUE,
                                             free.only = TRUE,
                                             add.labels = FALSE,
                                             add.class = FALSE,
                                             remove.duplicated = FALSE)
  }
  
  
  # restricted vcov
  info.r  <- JAC %*% VCOV %*% t(JAC)
  
  # Wald test statistic
  Wald <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)
  
  # df
  Wald.df <- nrow(JAC)
  
  # p-value based on chisq
  Wald.pvalue <- 1 - pchisq(Wald, df=Wald.df)
  
  list(stat=Wald, df=Wald.df, p.value=Wald.pvalue, se=object@Options$se)
}

##Hypotheses: 
##1. Risky behavior will have the strongest correlation with the externalizing behavior cluster (in comparison to all the other clusters)
##2. The relationship between anhedonia and depression will be stronger than the relationship between anhedonia and anxiety
##3. The correlations between negative affect and anxiety/depression will not differ significantly
##4. The relationship between anxious arousal and anxiety will be stronger than the relationship between anxious arousal and depression
##5. The correlations between dysphoric arousal and anxiety/depression will not differ significantly

globalHybridC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1
                  
                  risky ~~ a*ExtBehavior
                  risky ~~ b*Reexperiencing
                  risky ~~ c*Avoidance
                  risky ~~ d*NegAffect
                  risky ~~ e*Anhedonia
                  risky ~~ f*AnxArousal
                  risky ~~ g*DysphArousal
                  qids ~~ h*Anhedonia
                  gad ~~ i*Anhedonia
                  qids ~~ j*NegAffect
                  gad ~~ k*NegAffect
                  qids ~~ l*AnxArousal
                  gad ~~ m*AnxArousal
                  qids ~~ n*DysphArousal
                  gad ~~ o*DysphArousal

'
fittedGlobalHybridC <- sem(globalHybridC, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybridC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybridC, "all")
summary(fittedGlobalHybridC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)

unlist(lwt(fittedGlobalHybridC, 'a == b', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'a == c', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'a == d', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'a == e', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'a == f', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'a == g', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'h == i', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'j == k', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'l == m', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'n == o', std = TRUE)[c(1,3)])








# Descriptives on a combined dataset --------------------------------------






# Final ranking -----------------------------------------------------------

# 1. Hybrid model
# 2. Anhedonia
# 3. Externalizing behavior


# Final MIs ---------------------------------------------------------------

##Hybrid

# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 944           audit =~  risky1 185.282  0.445   0.445    0.606    0.606
# 2041          PCL20 ~~   QIDS1  84.283  0.265   0.265    0.399    0.399
# 952           audit =~  risky9  46.808 -0.132  -0.132   -0.320   -0.320
# 543     ExtBehavior =~    PCL2  38.472  0.458   0.458    0.492    0.492
# 653      AnxArousal =~  risky6  37.197  0.190   0.190    0.259    0.259
# 2560         QIDS12 ~~ risky14  36.154  0.035   0.035    0.235    0.235
# 718    DysphArousal =~  risky6  34.711  0.185   0.185    0.253    0.253
# 2895         risky1 ~~  risky9  33.003 -0.049  -0.049   -0.260   -0.260
# 3013         risky7 ~~ risky11  31.964  0.019   0.019    0.233    0.233
# 829             gad =~  risky6  31.617  0.175   0.175    0.239    0.239
# 769            qids =~  risky6  31.525  0.175   0.175    0.239    0.239
# 3044         risky9 ~~ risky11  31.241  0.026   0.026    0.247    0.247
# 459       NegAffect =~  risky6  30.738  0.176   0.176    0.240    0.240
# 2849           GAD6 ~~ risky10  29.563  0.085   0.085    0.224    0.224

##Externalizing behavior
# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 184.825  0.444   0.444    0.605    0.605
# 1962          PCL20 ~~   QIDS1  84.094  0.265   0.265    0.400    0.400
# 873           audit =~  risky9  46.301 -0.132  -0.132   -0.319   -0.319
# 464    ExtBehaviors =~    PCL2  40.172  0.475   0.475    0.510    0.510
# 574      AnxArousal =~  risky6  37.115  0.190   0.190    0.259    0.259
# 2481         QIDS12 ~~ risky14  36.233  0.035   0.035    0.236    0.236
# 281  Reexperiencing =~   PCL11  35.512  0.764   0.764    0.647    0.647
# 639    DysphArousal =~  risky6  33.716  0.183   0.183    0.250    0.250
# 2816         risky1 ~~  risky9  32.994 -0.049  -0.049   -0.260   -0.260
# 2934         risky7 ~~ risky11  32.002  0.019   0.019    0.234    0.234
# 750             gad =~  risky6  31.817  0.176   0.176    0.240    0.240
# 690            qids =~  risky6  31.300  0.174   0.174    0.238    0.238
# 2965         risky9 ~~ risky11  31.259  0.026   0.026    0.247    0.247
# 513    ExtBehaviors =~ risky10  30.172  0.232   0.232    0.298    0.298

##Anhedonia

# lhs op     rhs      mi    epc sepc.lv sepc.all sepc.nox
# 865           audit =~  risky1 185.054  0.444   0.444    0.605    0.605
# 1867          PCL20 ~~   QIDS1  84.890  0.269   0.269    0.390    0.390
# 873           audit =~  risky9  45.168 -0.131  -0.131   -0.316   -0.316
# 2481         QIDS12 ~~ risky14  37.470  0.035   0.035    0.240    0.240
# 639      AnxArousal =~  risky6  37.227  0.190   0.190    0.259    0.259
# 574    DysphArousal =~  risky6  34.758  0.196   0.196    0.268    0.268
# 2965         risky9 ~~ risky11  33.248  0.027   0.027    0.255    0.255
# 2934         risky7 ~~ risky11  32.567  0.019   0.019    0.236    0.236
# 750             gad =~  risky6  31.561  0.175   0.175    0.239    0.239
# 2816         risky1 ~~  risky9  31.541 -0.048  -0.048   -0.254   -0.254
# 690            qids =~  risky6  31.468  0.175   0.175    0.240    0.240
# 447       NegAffect =~  risky6  30.746  0.176   0.176    0.240    0.240
# 2770           GAD6 ~~ risky10  30.296  0.087   0.087    0.226    0.226
# 853           audit =~  QIDS12  29.325  0.112   0.112    0.223    0.223




# Final models - exploratory dataset --------------------------------------

globalHybrid3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1'
fittedGlobalHybrid3 <- sem(globalHybrid3, dataE, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybrid3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybrid3, "all")
summary(fittedGlobalHybrid3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalHybrid3)$cov
residuals(fittedGlobalHybrid3)$cov > .10
globalHybrid3MI <- modindices(fittedGlobalHybrid3)[order(-modindices(fittedGlobalHybrid3)$mi), ]
globalHybrid3MI


globalAnhedonia3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                     Anhedonia =~ PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2
                     audit =~ risky1
                     audit =~ risky9
                     PCL20 ~~ QIDS1'
fittedGlobalAnhedonia3 <- sem(globalAnhedonia3, dataE, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalAnhedonia3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalAnhedonia3, "all")
summary(fittedGlobalAnhedonia3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalAnhedonia3)$cov
residuals(fittedGlobalAnhedonia3)$cov > .10
globalAnhedonia3MI <- modindices(fittedGlobalAnhedonia3)[order(-modindices(fittedGlobalAnhedonia3)$mi), ]
globalAnhedonia3MI


globalExtBehaviors3 <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                        Avoidance =~ PCL6 + PCL7
                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                        ExtBehaviors =~ PCL15 + PCL16
                        AnxArousal =~ PCL17 + PCL18
                        DysphArousal =~ PCL19 + PCL20
                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                QIDS6 ~~ QIDS8
                                QIDS7 ~~ QIDS9
                                QIDS2 ~~ QIDS3
                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                               GAD1 ~~ GAD6
                               GAD5 ~~ GAD7
                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                 risky13 ~~ risky14
                                 risky2 ~~ risky7
                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                 AUDIT8 + AUDIT9 + AUDIT10
                                 AUDIT2 ~~ AUDIT3
                                 AUDIT1 ~~ AUDIT2
                        audit =~ risky1
                        audit =~ risky9
                        PCL20 ~~ QIDS1'
fittedGlobalExtBehaviors3 <- sem(globalExtBehaviors3, dataE, meanstructure = TRUE,
                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                 orthogonal = FALSE, missing = "ML",
                                 bootstrap = 5000)
fitmeasures(fittedGlobalExtBehaviors3,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalExtBehaviors3, "all")
summary(fittedGlobalExtBehaviors3, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalExtBehaviors3)$cov
residuals(fittedGlobalExtBehaviors3)$cov > .10
globalExtBehaviors3MI <- modindices(fittedGlobalExtBehaviors3)[order(-modindices(fittedGlobalExtBehaviors3)$mi), ]
globalExtBehaviors3MI



# Confirmatory dataset ----------------------------------------------------

oneFactorC <- 'ptsd =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5 + PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + 
                      PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedOneFactorC <- cfa(oneFactorC, dataC, meanstructure = TRUE,
                       std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                       orthogonal = FALSE, missing = "ML",
                       bootstrap = 5000)
fitmeasures(fittedOneFactorC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedOneFactorC, "all")
summary(fittedOneFactor, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedOneFactorC)$cov
residuals(fittedOneFactorC)$cov > .10
modindices(fittedOneFactorC)[order(-modindices(fittedOneFactorC)$mi), ]

# DSM-5 model -------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]; 
# 2. Avoidance [C1, C2]; 
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]; 
# 4. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])


dsm5C <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
         Avoidance =~ PCL6 + PCL7
         NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
         AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedDsm5C <- cfa(dsm5C, dataC, meanstructure = TRUE,
                  std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                  orthogonal = FALSE, missing = "ML",
                  bootstrap = 5000)
fitmeasures(fittedDsm5C, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedDsm5C, "all")
summary(fittedDsm5C, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5C)$cov
residuals(fittedDsm5C)$cov > .10
modindices(fittedDsm5C)[order(-modindices(fittedDsm5C)$mi), ]


# Three-factor model ------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance/Negative alternation in cognition and mood [C1, C2, D1, D2, D3, D4, D5, D6, D7]
# 3. Alterations in arousal and reactivity [E1, E2, E3, E4, E5, E6])

threeFactorC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                AvoidNegAlternation =~ PCL6 + PCL7 + PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                AltArousal =~ PCL15 + PCL16 + PCL17 + PCL18 + PCL19 + PCL20'
fittedThreeFactorC <- cfa(threeFactorC, dataC, meanstructure = TRUE,
                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                         orthogonal = FALSE, missing = "ML",
                         bootstrap = 5000)
fitmeasures(fittedThreeFactorC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedThreeFactorC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedThreeFactorC)$cov
residuals(fittedThreeFactorC)$cov > .10
modindices(fittedThreeFactorC)[order(-modindices(fittedThreeFactorC)$mi), ]


# DSM-5 Dysphoria model ---------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Dyshoria [D1, D2, D3, D4, D5, D6, D7, E1, E2, E5, E6]
# 4. Alternation in arousal and reactivity [E3, E4])

dsm5dysphoriaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14 + PCL15 + PCL16 + PCL19 + PCL20
                  AltArousal =~ PCL17 + PCL18'
fittedDsm5dysphoriaC <- cfa(dsm5dysphoria, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedDsm5dysphoriaC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphoriaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphoriaC)$cov
residuals(fittedDsm5dysphoriaC)$cov > .10
modindices(fittedDsm5dysphoriaC)[order(-modindices(fittedDsm5dysphoriaC)$mi), ]


# DSM-5 Dysphoric arousal model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative alternation in cognition and mood [D1, D2, D3, D4, D5, D6, D7]
# 4. Dysphoric arousal [E1, E2, E5, E6]
# 5. Anxious arousal [E3, E4])

dsm5dysphAruosalC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAlternation =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18'
fittedDsm5dysphArousalC <- cfa(dsm5dysphAruosalC, dataC, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedDsm5dysphArousalC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedDsm5dysphArousalC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedDsm5dysphArousalC)$cov
residuals(fittedDsm5dysphArousalC)$cov > .10
modindices(fittedDsm5dysphArousalC)[order(-modindices(fittedDsm5dysphArousalC)$mi), ]

# Anhedonia model ---------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Dysphoric arousal [E1, E2, E5, E6]
# 6. Anxious arousal [E3, E4])

anhedoniaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
              Avoidance =~ PCL6 + PCL7
              NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
              Anhedonia =~ PCL12 + PCL13 + PCL14
              DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
              AnxArousal =~ PCL17 + PCL18'
fittedAnhedoniaC <- cfa(anhedoniaC, dataC, meanstructure = TRUE,
                       std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                       orthogonal = FALSE, missing = "ML",
                       bootstrap = 5000)
fitmeasures(fittedAnhedoniaC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedAnhedoniaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedAnhedoniaC)$cov
residuals(fittedAnhedoniaC)$cov > .10
modindices(fittedAnhedoniaC)[order(-modindices(fittedAnhedoniaC)$mi), ]


# Externalizing behaviors model -------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Numbing [D1, D2, D3, D4, D5, D6, D7]
# 4. Externalizing behaviors [E1, E2]
# 5. Anxious arousal [E3, E4]
# 6. Dysphoric arousal [E5, E6])

extBehaviorsC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                 Avoidance =~ PCL6 + PCL7
                 Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                 ExtBehaviors =~ PCL15 + PCL16
                 AnxArousal =~ PCL17 + PCL18 
                 DysphArousal =~ PCL19 + PCL20'
fittedExtBehaviorsC <- cfa(extBehaviorsC, dataC, meanstructure = TRUE,
                          std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                          orthogonal = FALSE, missing = "ML",
                          bootstrap = 5000)
fitmeasures(fittedExtBehaviorsC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedExtBehaviorsC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedExtBehaviorsC)$cov
residuals(fittedExtBehaviorsC)$cov > .10
modindices(fittedExtBehaviorsC)[order(-modindices(fittedExtBehaviorsC)$mi), ]


# Hybrid model ------------------------------------------------------------

# 1. Reexperiencing [B1, B2, B3, B4, B5]
# 2. Avoidance [C1, C2]
# 3. Negative affect [D1, D2, D3, D4]
# 4. Anhedonia [D5, D6, D7]
# 5. Externalizing behaviours [E1, E2]
# 6. Anxious arousal [E3, E4]
# 7. Dysphoric arousal [E5, E6])

hybridModelC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                Avoidance =~ PCL6 + PCL7
                NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11 
                Anhedonia =~ PCL12 + PCL13 + PCL14
                ExtBehavior =~ PCL15 + PCL16 
                AnxArousal =~ + PCL17 + PCL18
                DysphArousal =~ PCL19 + PCL20'
fittedHybridModelC <- cfa(hybridModelC, dataC, meanstructure = TRUE,
                         std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                         orthogonal = FALSE, missing = "ML",
                         bootstrap = 5000)
fitmeasures(fittedHybridModelC, 
            c("chisq", "df", "pvalue", "cfi", "tli", 
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
summary(fittedHybridModelC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedHybridModelC)$cov
residuals(fittedHybridModelC)$cov > .10
modindices(fittedHybridModelC)[order(-modindices(fittedHybridModelC)$mi), ]


# MIs for CFAs ------------------------------------------------------------

##Apart from the one-factor model
# PCL12 ~~ PCL14
# PCL15 ~~ PCL16
# PCL2 ~~ PCL15


# Fit measures summary ----------------------------------------------------

model <- c("One-FactorC model", "DSM-5 model", "Three-FactorC model", "DSM-5 Dysphoria model", "DSM-5 Dysphoric arousal model",
           "Anhedonia model", "Externalizing behaviors model", "Hybrid model")
chisq <- c(fitmeasures(fittedOneFactorC, "chisq"), fitmeasures(fittedDsm5C, "chisq"), 
           fitmeasures(fittedThreeFactorC, "chisq"), fitmeasures(fittedDsm5dysphoriaC, "chisq"), 
           fitmeasures(fittedDsm5dysphArousalC, "chisq"), fitmeasures(fittedAnhedoniaC, "chisq"), 
           fitmeasures(fittedExtBehaviorsC, "chisq"), fitmeasures(fittedHybridModelC, "chisq")) 
chisq <- round(chisq, 3)
df <- c(fitmeasures(fittedOneFactorC, "df"), fitmeasures(fittedDsm5C, "df"), 
        fitmeasures(fittedThreeFactorC, "df"), fitmeasures(fittedDsm5dysphoriaC, "df"), 
        fitmeasures(fittedDsm5dysphArousalC, "df"), fitmeasures(fittedAnhedoniaC, "df"), 
        fitmeasures(fittedExtBehaviorsC, "df"), fitmeasures(fittedHybridModelC, "df")) 
df <- round(df, 3)
pvalue <- c(fitmeasures(fittedOneFactorC, "pvalue"), fitmeasures(fittedDsm5C, "pvalue"), 
            fitmeasures(fittedThreeFactorC, "pvalue"), fitmeasures(fittedDsm5dysphoriaC, "pvalue"),
            fitmeasures(fittedDsm5dysphArousal, "pvalue"), fitmeasures(fittedAnhedoniaC, "pvalue"), 
            fitmeasures(fittedExtBehaviorsC, "pvalue"), fitmeasures(fittedHybridModelC, "pvalue")) 
pvalue <- round(pvalue, 3)
cfi <- c(fitmeasures(fittedOneFactorC, "cfi"), fitmeasures(fittedDsm5, "cfi"), 
         fitmeasures(fittedThreeFactorC, "cfi"), fitmeasures(fittedDsm5dysphoriaC, "cfi"), 
         fitmeasures(fittedDsm5dysphArousal, "cfi"), fitmeasures(fittedAnhedoniaC, "cfi"), 
         fitmeasures(fittedExtBehaviorsC, "cfi"), fitmeasures(fittedHybridModelC, "cfi")) 
cfi <- round(cfi, 3)
tli <- c(fitmeasures(fittedOneFactorC, "tli"), fitmeasures(fittedDsm5, "tli"), 
         fitmeasures(fittedThreeFactorC, "tli"), fitmeasures(fittedDsm5dysphoriaC, "tli"), 
         fitmeasures(fittedDsm5dysphArousalC, "tli"), fitmeasures(fittedAnhedoniaC, "tli"), 
         fitmeasures(fittedExtBehaviorsC, "tli"), fitmeasures(fittedHybridModelC, "tli")) 
tli <- round(tli, 3)
rmsea <- c(fitmeasures(fittedOneFactorC, "rmsea"), fitmeasures(fittedDsm5C, "rmsea"), 
           fitmeasures(fittedThreeFactorC, "rmsea"), fitmeasures(fittedDsm5dysphoriaC, "rmsea"), 
           fitmeasures(fittedDsm5dysphArousalC, "rmsea"), fitmeasures(fittedAnhedoniaC, "rmsea"), 
           fitmeasures(fittedExtBehaviorsC, "rmsea"), fitmeasures(fittedHybridModelC, "rmsea"))
rmsea <- round(rmsea, 3)
rmseaLI <- c(fitmeasures(fittedOneFactorC, "rmsea.ci.lower"), fitmeasures(fittedDsm5, "rmsea.ci.lower"), 
             fitmeasures(fittedThreeFactorC, "rmsea.ci.lower"), fitmeasures(fittedDsm5dysphoriaC, "rmsea.ci.lower"), 
             fitmeasures(fittedDsm5dysphArousalC, "rmsea.ci.lower"), fitmeasures(fittedAnhedoniaC, "rmsea.ci.lower"),
             fitmeasures(fittedExtBehaviorsC, "rmsea.ci.lower"), fitmeasures(fittedHybridModelC, "rmsea.ci.lower")) 
rmseaLI <- round(rmseaLI, 3)
rmseaUI<- c(fitmeasures(fittedOneFactorC, "rmsea.ci.upper"), fitmeasures(fittedDsm5, "rmsea.ci.upper"), 
            fitmeasures(fittedThreeFactorC, "rmsea.ci.upper"), fitmeasures(fittedDsm5dysphoriaC, "rmsea.ci.upper"), 
            fitmeasures(fittedDsm5dysphArousalC, "rmsea.ci.upper"), fitmeasures(fittedAnhedoniaC, "rmsea.ci.upper"), 
            fitmeasures(fittedExtBehaviorsC, "rmsea.ci.upper"), fitmeasures(fittedHybridModelC, "rmsea.ci.upper")) 
rmseaUI <- round(rmseaUI, 3)
srmr <- c(fitmeasures(fittedOneFactorC, "srmr"), fitmeasures(fittedDsm5, "srmr"), 
          fitmeasures(fittedThreeFactorC, "srmr"), fitmeasures(fittedDsm5dysphoriaC, "srmr"), 
          fitmeasures(fittedDsm5dysphArousalC, "srmr"), fitmeasures(fittedAnhedoniaC, "srmr"), 
          fitmeasures(fittedExtBehaviorsC, "srmr"), fitmeasures(fittedHybridModelC, "srmr")) 
srmr <- round(srmr, 3)
aic <- c(fitmeasures(fittedOneFactorC, "aic"), fitmeasures(fittedDsm5C, "aic"), 
         fitmeasures(fittedThreeFactorC, "aic"), fitmeasures(fittedDsm5dysphoriaC, "aic"), 
         fitmeasures(fittedDsm5dysphArousalC, "aic"), fitmeasures(fittedAnhedoniaC, "aic"), 
         fitmeasures(fittedExtBehaviorsC, "aic"), fitmeasures(fittedHybridModelC, "aic")) 
aic <- round(aic, 3)
bic <- c(fitmeasures(fittedOneFactorC, "bic"), fitmeasures(fittedDsm5C, "bic"), 
         fitmeasures(fittedThreeFactorC, "bic"), fitmeasures(fittedDsm5dysphoriaC, "bic"), 
         fitmeasures(fittedDsm5dysphArousalC, "bic"), fitmeasures(fittedAnhedoniaC, "bic"), 
         fitmeasures(fittedExtBehaviorsC, "bic"), fitmeasures(fittedHybridModelC, "bic")) 
bic <- round(bic, 3)

modelFits <- data.frame(model, chisq, df, pvalue, cfi, tli, rmsea, rmseaLI, rmseaUI, srmr, aic, bic)

kbl(modelFits)
modelFits %>% kbl() %>% kable_material()


# Nested models comparisons -----------------------------------------------

anova(fittedOneFactorC, fittedDsm5C, fittedThreeFactorC, fittedDsm5dysphoriaC, 
      fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC, fittedHybridModelC)
anova(fittedOneFactorC, fittedDsm5C)
anova(fittedOneFactorC, fittedThreeFactorC)
anova(fittedOneFactorC, fittedDsm5dysphoriaC)
anova(fittedOneFactorC, fittedDsm5dysphArousalC)
anova(fittedOneFactorC, fittedAnhedoniaC)
anova(fittedOneFactorC, fittedExtBehaviorsC)
anova(fittedOneFactorC, fittedHybridModelC)


##DSM-5 comparisons 

anova(fittedDsm5C, fittedThreeFactorC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5C, fittedThreeFactorC)
anova(fittedDsm5C, fittedDsm5dysphArousalC)
anova(fittedDsm5C, fittedAnhedoniaC)
anova(fittedDsm5C, fittedExtBehaviorsC)

##Three-factor comparisons 

anova(fittedThreeFactorC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedThreeFactorC, fittedDsm5dysphArousalC)
anova(fittedThreeFactorC, fittedAnhedoniaC)
anova(fittedThreeFactorC, fittedExtBehaviorsC)

##Dysphoria comparisons 

anova(fittedDsm5dysphoriaC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5dysphoriaC, fittedDsm5dysphArousalC)
anova(fittedDsm5dysphoriaC, fittedAnhedoniaC)
anova(fittedDsm5dysphoriaC, fittedExtBehaviorsC)

##Dysphoric arousal comparisons 

anova(fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedDsm5dysphArousalC, fittedAnhedoniaC)
anova(fittedDsm5dysphArousalC, fittedExtBehaviorsC)

##Hybrid model comparisons

anova(fittedHybridModelC, fittedOneFactorC, fittedDsm5C, fittedThreeFactorC, 
      fittedDsm5dysphoriaC, fittedDsm5dysphArousalC, fittedAnhedoniaC, fittedExtBehaviorsC)
anova(fittedHybridModelC, fittedOneFactorC)
anova(fittedHybridModelC, fittedDsm5C)
anova(fittedHybridModelC, fittedThreeFactorC)
anova(fittedHybridModelC, fittedDsm5dysphoriaC)
anova(fittedHybridModelC, fittedDsm5dysphArousalC)
anova(fittedHybridModelC, fittedAnhedoniaC)
anova(fittedHybridModelC, fittedExtBehaviorsC)

# Non-nested models comparisons -------------------------------------------

##DSM-5 vs. Dysphoria
anova(fittedDsm5C, fittedDsm5dysphoriaC)
lrtest(fittedDsm5C, fittedDsm5dysphoriaC)

##Three-factor vs. Dysphoria

anova(fittedThreeFactorC, fittedDsm5dysphoriaC)
lrtest(fittedThreeFactorC, fittedDsm5dysphoriaC)

##Anhedonia vs. Externalizing behaviors

anova(fittedAnhedoniaC, fittedExtBehaviorsC)
lrtest(fittedAnhedoniaC, fittedExtBehaviorsC)





# Final models - confirmatory dataset -------------------------------------

globalHybridC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1'
fittedGlobalHybridC <- sem(globalHybridC, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybridC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybridC, "all")
summary(fittedGlobalHybridC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalHybridC)$cov
residuals(fittedGlobalHybridC)$cov > .10
globalHybridCMI <- modindices(fittedGlobalHybridC)[order(-modindices(fittedGlobalHybridC)$mi), ]
globalHybridCMI


globalAnhedoniaC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                     Avoidance =~ PCL6 + PCL7
                     NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                     Anhedonia =~ PCL12 + PCL13 + PCL14
                     DysphArousal =~ PCL15 + PCL16 + PCL19 + PCL20
                     AnxArousal =~ PCL17 + PCL18
                     qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                             QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                             QIDS6 ~~ QIDS8
                             QIDS7 ~~ QIDS9
                             QIDS2 ~~ QIDS3
                     gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                            GAD1 ~~ GAD6
                            GAD5 ~~ GAD7
                     risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                              risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                              risky13 ~~ risky14
                              risky2 ~~ risky7
                     audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                              AUDIT8 + AUDIT9 + AUDIT10
                              AUDIT2 ~~ AUDIT3
                              AUDIT1 ~~ AUDIT2
                     audit =~ risky1
                     audit =~ risky9
                     PCL20 ~~ QIDS1'
fittedGlobalAnhedoniaC <- sem(globalAnhedoniaC, dataC, meanstructure = TRUE,
                              std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                              orthogonal = FALSE, missing = "ML",
                              bootstrap = 5000)
fitmeasures(fittedGlobalAnhedoniaC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalAnhedoniaC, "all")
summary(fittedGlobalAnhedoniaC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalAnhedoniaC)$cov
residuals(fittedGlobalAnhedoniaC)$cov > .10
globalAnhedoniaCMI <- modindices(fittedGlobalAnhedoniaC)[order(-modindices(fittedGlobalAnhedoniaC)$mi), ]
globalAnhedoniaCMI


globalExtBehaviorsC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                        Avoidance =~ PCL6 + PCL7
                        Numbing =~ PCL8 + PCL9 + PCL10 + PCL11 + PCL12 + PCL13 + PCL14
                        ExtBehaviors =~ PCL15 + PCL16
                        AnxArousal =~ PCL17 + PCL18
                        DysphArousal =~ PCL19 + PCL20
                        qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                                QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                                QIDS6 ~~ QIDS8
                                QIDS7 ~~ QIDS9
                                QIDS2 ~~ QIDS3
                        gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                               GAD1 ~~ GAD6
                               GAD5 ~~ GAD7
                        risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                                 risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                                 risky13 ~~ risky14
                                 risky2 ~~ risky7
                        audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                                 AUDIT8 + AUDIT9 + AUDIT10
                                 AUDIT2 ~~ AUDIT3
                                 AUDIT1 ~~ AUDIT2
                        audit =~ risky1
                        audit =~ risky9
                        PCL20 ~~ QIDS1'
fittedGlobalExtBehaviorsC <- sem(globalExtBehaviorsC, dataC, meanstructure = TRUE,
                                 std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                                 orthogonal = FALSE, missing = "ML",
                                 bootstrap = 5000)
fitmeasures(fittedGlobalExtBehaviorsC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalExtBehaviorsC, "all")
summary(fittedGlobalExtBehaviorsC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)
residuals(fittedGlobalExtBehaviorsC)$cov
residuals(fittedGlobalExtBehaviorsC)$cov > .10
globalExtBehaviorsCMI <- modindices(fittedGlobalExtBehaviorsC)[order(-modindices(fittedGlobalExtBehaviorsC)$mi), ]
globalExtBehaviorsCMI



# Descriptives on a combined dataset --------------------------------------

dataFull <- rbind(dataE, dataC)
describe(select(dataFull, PCL1:PCL20))[, 3:4]
lowerCor(select(dataFull, PCL1:PCL20))






fittedDsm5@test$standard$stat - fittedThreeFactor@test$standard$stat
pchisq(abs(fittedDsm5@test$standard$stat - fittedThreeFactor@test$standard$stat),
       abs(fittedDsm5@test$standard$df - fittedThreeFactor@test$standard$df),
       lower.tail = FALSE)

anova(fittedOneFactor, fittedDsm5, fittedThreeFactor)
lrtest(fittedGlobalOneFactor, fittedGlobalDsm5, fittedGlobalThreeFactor,
       fittedGlobalDysphoria, fittedDsm5dysphArousal, fittedGlobalAnhedonia,
       fittedGlobalExtBehaviors, fittedGlobalHybrid)
lrtest(fittedAnhedonia, fittedHybridModel)


# Equality of regression coefficients -------------------------------------

##Test equality of regression coefficients using Wald test with lavaan (Klopp, 2020)
lwt <- function(object, constraints = NULL, verbose = FALSE, std = FALSE) {
  
  if(object@optim$npar > 0L && !object@optim$converged)
    stop("lavaan ERROR: model did not converge")
  
  if(is.null(constraints) || nchar(constraints) == 0L) {
    stop("lavaan ERROR: constraints are empty")
  }
  
  # remove == constraints from parTable
  PT <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
  eq.idx <- which(PT$op == "==")
  if(length(eq.idx) > 0L) {
    PT <- PT[-eq.idx,]
  }
  partable <- as.list(PT)
  
  # parse constraints
  FLAT <- lavParseModelString( constraints ); CON <- attr(FLAT, "constraints")
  LIST <- list()
  if(length(CON) > 0L) {
    lhs = unlist(lapply(CON, "[[", "lhs"))
    op = unlist(lapply(CON, "[[",  "op"))
    rhs = unlist(lapply(CON, "[[", "rhs"))
    LIST$lhs        <- c(LIST$lhs,        lhs)
    LIST$op         <- c(LIST$op,         op)
    LIST$rhs        <- c(LIST$rhs,        rhs)
  } else {
    stop("lavaan ERROR: no equality constraints found in constraints argument")
  }
  
  # theta = free parameters only
  # changed: to include the std argument
  if (!std) {
    theta <- object@optim$x
  }
  
  if(std) {
    theta <- lavaan:::lav_standardize_all(object)
    fixed <- which(unclass(partable(object))$free %in% c(0))
    if (length(fixed) !=0 ) theta <- theta[-fixed]
  }
  
  # build constraint function
  ceq.function <- lav_partable_constraints_ceq(partable = partable,
                                               con = LIST, debug = FALSE)
  
  # compute jacobian restrictions
  JAC <- try(lav_func_jacobian_complex(func = ceq.function, x = theta), silent=TRUE)
  if(inherits(JAC, "try-error")) { # eg. pnorm()
    JAC <- lav_func_jacobian_simple(func = ceq.function, x = theta)
  }
  
  if(verbose) {
    cat("Restriction matrix (jacobian):\n"); print(JAC); cat("\n")
  }
  
  # linear restriction
  theta.r <- ceq.function( theta )
  
  if(verbose) {
    cat("Restricted theta values:\n"); print(theta.r); cat("\n")
  }
  
  # get VCOV
  # VCOV <- vcov(object, labels = FALSE)
  # avoid S4 dispatch
  # changed: to include the std argument
  if (!std) {
    # changed: addad lavaan::: because function is not exported
    VCOV <- lavaan:::lav_object_inspect_vcov(object, standardized = FALSE,
                                             free.only = TRUE,
                                             add.labels = FALSE,
                                             add.class = FALSE,
                                             remove.duplicated = FALSE)
  }
  
  if (std) {
    # changed: addad lavaan::: because function is not exported
    VCOV <- lavaan:::lav_object_inspect_vcov(object, standardized = TRUE,
                                             free.only = TRUE,
                                             add.labels = FALSE,
                                             add.class = FALSE,
                                             remove.duplicated = FALSE)
  }
  
  
  # restricted vcov
  info.r  <- JAC %*% VCOV %*% t(JAC)
  
  # Wald test statistic
  Wald <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)
  
  # df
  Wald.df <- nrow(JAC)
  
  # p-value based on chisq
  Wald.pvalue <- 1 - pchisq(Wald, df=Wald.df)
  
  list(stat=Wald, df=Wald.df, p.value=Wald.pvalue, se=object@Options$se)
}


globalHybridC <- 'Reexperiencing =~ PCL1 + PCL2 + PCL3 + PCL4 + PCL5
                  Avoidance =~ PCL6 + PCL7
                  NegAffect =~ PCL8 + PCL9 + PCL10 + PCL11
                  Anhedonia =~ PCL12 + PCL13 + PCL14
                  ExtBehavior =~ PCL15 + PCL16
                  AnxArousal =~ + PCL17 + PCL18
                  DysphArousal =~ PCL19 + PCL20
                  qids =~ QIDS1 + QIDS2 + QIDS3 + QIDS4 + QIDS5 + QIDS6 + QIDS7 + QIDS8 +
                          QIDS9 + QIDS10 + QIDS11 + QIDS12 + QIDS13 + QIDS14 + QIDS15 + QIDS16
                          QIDS6 ~~ QIDS8
                          QIDS7 ~~ QIDS9
                          QIDS2 ~~ QIDS3
                  gad =~ GAD1 + GAD2 + GAD3 + GAD4 + GAD5 + GAD6 + GAD7
                         GAD1 ~~ GAD6
                         GAD5 ~~ GAD7
                  risky =~ risky1 + risky2 + risky3 + risky4 + risky5 + risky6 + risky7 +
                           risky8 + risky9 + risky10 + risky11 + risky12 + risky13 + risky14
                           risky13 ~~ risky14
                           risky2 ~~ risky7
                  audit =~ AUDIT1 + AUDIT2 + AUDIT3 + AUDIT4 + AUDIT5 + AUDIT6 + AUDIT7 +
                           AUDIT8 + AUDIT9 + AUDIT10
                           AUDIT2 ~~ AUDIT3
                           AUDIT1 ~~ AUDIT2
                  audit =~ risky1
                  audit =~ risky9
                  PCL20 ~~ QIDS1
                  
                  risky ~~ ebRisky*ExtBehavior
                  risky ~~ reRisky*Reexperiencing
                  risky ~~ avRisky*Avoidance
                  risky ~~ naRisky*NegAffect
                  risky ~~ anRisky*Anhedonia
                  risky ~~ aaRisky*AnxArousal
                  risky ~~ daRisky*DysphArousal
                  qids ~~ anQids*Anhedonia
                  qids ~~ reQids*Reexperiencing
                  qids ~~ avQids*Avoidance
                  qids ~~ naQids*NegAffect
                  qids ~~ ebQids*ExtBehavior
                  qids ~~ aaQids*AnxArousal
                  qids ~~ daQids*DysphArousal
                  gad ~~ aaGad*AnxArousal
                  gad ~~ reGad*Reexperiencing
                  gad ~~ avGad*Avoidance
                  gad ~~ naGad*NegAffect
                  gad ~~ anGad*Anhedonia
                  gad ~~ ebGad*ExtBehavior
                  gad ~~ daGad*DysphArousal
'

fittedGlobalHybridC <- sem(globalHybridC, dataC, meanstructure = TRUE,
                           std.lv = TRUE, mimic = "Mplus", estimator = "ML",
                           orthogonal = FALSE, missing = "ML",
                           bootstrap = 5000)
fitmeasures(fittedGlobalHybridC,
            c("chisq", "df", "pvalue", "cfi", "tli",
              "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic" ,"bic"))
fitmeasures(fittedGlobalHybridC, "all")
summary(fittedGlobalHybridC, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)

unlist(lwt(fittedGlobalHybridC, 'ebRisky == reRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'ebRisky == avRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'ebRisky == naRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'ebRisky == anRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'ebRisky == aaRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'ebRisky == daRisky', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == reQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == avQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == naQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == ebQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == aaQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == daQids', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == reGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == avGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == naGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == anGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == ebGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaGad == daGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'anQids == anGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'naQids == naGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'aaQids == aaGad', std = TRUE)[c(1,3)])
unlist(lwt(fittedGlobalHybridC, 'daQids == daGad', std = TRUE)[c(1,3)])

#Holm p-value adjustment (p-values are sorted from the smallest to largerst)
##externalizing behaviours
pvalues <- c(1.072233e-07,2.623170e-07,7.826395e-07,6.545618e-06,1.447224e-06,0.005855748)
p.adjust(pvalues, method = "holm")
##anhedonia
pvalues <- c(4.145189e-09,1.207538e-04,0.006308463,0.0153679,0.1712953,0.3227154)
p.adjust(pvalues, method = "holm")
##anxious arousal
pvalues <- c(0.000440757,0.001537357,0.4996941,0.5514683,0.81585572,0.8762490)
p.adjust(pvalues, method = "holm")
##anhedonia vs negative affect
pvalues <- c(0.02197239,0.05336481)
p.adjust(pvalues, method = "holm")
##anxious arousal vs dysphoric arousal
pvalues <- c(0.01016448,0.87701635)
p.adjust(pvalues, method = "holm")

# Proportion of people endorsed in each cluster in hybrid model---------------------------on full dataset

sum(is.na(select(dataFull, PCL1:PCL20)))

nrow(subset(dataFull, PCL1 > 1 | PCL2 > 1 | PCL3 > 1 | PCL4 > 1 | PCL5 > 1)) / nrow(dataFull) ##41.5%
nrow(subset(dataFull, PCL6 > 1 | PCL7 > 1)) / nrow(dataFull) ##38.2%
nrow(subset(dataFull, PCL8 > 1 | PCL9 > 1 | PCL10 > 1 | PCL11 > 1)) / nrow(dataFull) ##41.1%
nrow(subset(dataFull, PCL12 > 1 | PCL13 > 1 | PCL14 > 1)) / nrow(dataFull) ##33.6%
nrow(subset(dataFull, PCL15 > 1 | PCL16 > 1)) / nrow(dataFull) ##27.6%
nrow(subset(dataFull, PCL17 > 1 | PCL18 > 1)) / nrow(dataFull) ##35.8%
nrow(subset(dataFull, PCL19 > 1 | PCL20 > 1)) / nrow(dataFull) ##29.7%

#prevalence of PTSD (hybrid model)

dataFull$endRE1 <- ifelse((dataFull$PCL1) > 1, 1, 0)
dataFull$endRE2 <- ifelse((dataFull$PCL2) > 1, 1, 0)
dataFull$endRE3 <- ifelse((dataFull$PCL3) > 1, 1, 0)
dataFull$endRE4 <- ifelse((dataFull$PCL4) > 1, 1, 0)
dataFull$endRE5 <- ifelse((dataFull$PCL5) > 1, 1, 0)
dataFull$endRE <- ifelse((dataFull$endRE1 + dataFull$endRE2 + dataFull$endRE3 + dataFull$endRE4 + dataFull$endRE5) > 0, 1, 0)

dataFull$endAV1 <- ifelse((dataFull$PCL6) > 1, 1, 0)
dataFull$endAV2 <- ifelse((dataFull$PCL7) > 1, 1, 0)
dataFull$endAV <- ifelse((dataFull$endAV1 + dataFull$endAV2) > 0, 1, 0)

dataFull$endNA1 <- ifelse((dataFull$PCL8) > 1, 1, 0)
dataFull$endNA2 <- ifelse((dataFull$PCL9) > 1, 1, 0)
dataFull$endNA3 <- ifelse((dataFull$PCL10) > 1, 1, 0)
dataFull$endNA4 <- ifelse((dataFull$PCL11) > 1, 1, 0)
dataFull$endNA <- ifelse((dataFull$endNA1 + dataFull$endNA2 + dataFull$endNA3 + dataFull$endNA4) > 0, 1, 0)

dataFull$endAN1 <- ifelse((dataFull$PCL12) > 1, 1, 0)
dataFull$endAN2 <- ifelse((dataFull$PCL13) > 1, 1, 0)
dataFull$endAN3 <- ifelse((dataFull$PCL14) > 1, 1, 0)
dataFull$endAN <- ifelse((dataFull$endAN1 + dataFull$endAN2 + dataFull$endAN3) > 0, 1, 0)

dataFull$endEB1 <- ifelse((dataFull$PCL15) > 1, 1, 0)
dataFull$endEB2 <- ifelse((dataFull$PCL16) > 1, 1, 0)
dataFull$endEB <- ifelse((dataFull$endEB1 + dataFull$endEB2) > 0, 1, 0)

dataFull$endAA1 <- ifelse((dataFull$PCL17) > 1, 1, 0)
dataFull$endAA2 <- ifelse((dataFull$PCL18) > 1, 1, 0)
dataFull$endAA <- ifelse((dataFull$endAA1 + dataFull$endAA2) > 0, 1, 0)

dataFull$endDA1 <- ifelse((dataFull$PCL19) > 1, 1, 0)
dataFull$endDA2 <- ifelse((dataFull$PCL20) > 1, 1, 0)
dataFull$endDA <- ifelse((dataFull$endDA1 + dataFull$endDA2) > 0, 1, 0)

dataFull$prevalence <- ifelse((dataFull$endRE + dataFull$endAV + dataFull$endNA + dataFull$endAN + dataFull$endEB + dataFull$endAA + dataFull$endDA) == 7, 1,0)

table(dataFull$prevalence)[2] / nrow(dataFull) ##10.6%

#prevalence PTSD (DSM-5 model)
##Re-experiencing and avoidance same structure like hybrid model
#dataFull$endRE1 <- ifelse((dataFull$PCL1) > 1, 1, 0)
#dataFull$endRE2 <- ifelse((dataFull$PCL2) > 1, 1, 0)
#dataFull$endRE3 <- ifelse((dataFull$PCL3) > 1, 1, 0)
#dataFull$endRE4 <- ifelse((dataFull$PCL4) > 1, 1, 0)
#dataFull$endRE5 <- ifelse((dataFull$PCL5) > 1, 1, 0)
#dataFull$endRE <- ifelse((dataFull$endRE1 + dataFull$endRE2 + dataFull$endRE3 + dataFull$endRE4 + dataFull$endRE5) > 0, 1, 0)
#dataFull$endAV1 <- ifelse((dataFull$PCL6) > 1, 1, 0)
#dataFull$endAV2 <- ifelse((dataFull$PCL7) > 1, 1, 0)
#dataFull$endAV <- ifelse((dataFull$endAV1 + dataFull$endAV2) > 0, 1, 0)

dataFull$endNACM1 <- ifelse((dataFull$PCL8) > 1, 1, 0)
dataFull$endNACM2 <- ifelse((dataFull$PCL9) > 1, 1, 0)
dataFull$endNACM3 <- ifelse((dataFull$PCL10) > 1, 1, 0)
dataFull$endNACM4 <- ifelse((dataFull$PCL11) > 1, 1, 0)
dataFull$endNACM5 <- ifelse((dataFull$PCL12) > 1, 1, 0)
dataFull$endNACM6 <- ifelse((dataFull$PCL13) > 1, 1, 0)
dataFull$endNACM7 <- ifelse((dataFull$PCL14) > 1, 1, 0)
dataFull$endNACM <- ifelse((dataFull$endNACM1 + dataFull$endNACM2 + dataFull$endNACM3 + dataFull$endNACM4 + dataFull$endNACM5 + dataFull$endNACM6 + dataFull$endNACM7) > 1, 1, 0)

dataFull$endAAR1 <- ifelse((dataFull$PCL15) > 1, 1, 0)
dataFull$endAAR2 <- ifelse((dataFull$PCL16) > 1, 1, 0)
dataFull$endAAR3 <- ifelse((dataFull$PCL17) > 1, 1, 0)
dataFull$endAAR4 <- ifelse((dataFull$PCL18) > 1, 1, 0)
dataFull$endAAR5 <- ifelse((dataFull$PCL19) > 1, 1, 0)
dataFull$endAAR6 <- ifelse((dataFull$PCL20) > 1, 1, 0)
dataFull$endAAR <- ifelse((dataFull$endAAR1 + dataFull$endAAR2 + dataFull$endAAR3 + dataFull$endAAR4 + dataFull$endAAR5 + dataFull$endAAR6) > 1, 1, 0)

dataFull$prevalenceDSM5 <- ifelse((dataFull$endRE + dataFull$endAV + dataFull$endNACM + dataFull$endAAR) == 4, 1,0)

table(dataFull$endNACM)[2] / nrow(dataFull) ##32.8%
table(dataFull$endAAR)[2] / nrow(dataFull) ##30.9%
table(dataFull$prevalenceDSM5)[2] / nrow(dataFull) ##18.7%

# Endorsement of PTSD symptoms
dataFull %>% 
  select(endRE1, endRE2, endRE3, endRE4, endRE5, endAV1, endAV2, endNA1, endNA2, endNA3, endNA4, endAN1, endAN2, endAN3, endEB1, endEB2, endAA1, endAA2, endDA1, endDA2) %>% 
  map(~round(prop.table(table(.x)), 3))

#B1: Recurrent thoughts of trauma 24.4%
#B2: Recurrent dreams of trauma 15.3%
#B3:Flashbacks 18.2%
#B4: Psychological cue reactivity 28.9%
#B5: Physiological cue reactivity 25.1%
#C1: Avoidance of thoughts of trauma 31.4% The most endorsed symptom
#C2: Avoidance of reminders of trauma 24.4%
#D1: Memory impairment 18.5%
#D2: Negative beliefs 23.9%
#D3: Distorted blame 23%
#D4: Persistent negative emotional state 26.1%
#D5: Diminished interest in activities 20%
#D6: Feelings of detachment from others 20.9%
#D7: Inability to experience positive emotions 21.5%
#E1: Irritability or anger 22.4%
#E2: Reckless or self-destructive behaviour 14.4% The least endorsed symptom
#E3: Hypervigilance 28.4%
#E4: Exaggerated startle response 22.9%
#E5: Difficulty concentrating 19.1%
#E6: Sleeping difficulties 23.3%

























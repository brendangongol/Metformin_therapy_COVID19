
 
#############################################################################################
#### This is an analysis of the effect of metformin on COVID-19 morbindity and mortality ####
#############################################################################################
homedir <- "C:/Users/breng/Dropbox/COVID19 metformin"
setwd(homedir)
library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ClinicalDataMineR)
dt <- fread("./data/2020.5.13 English DM COVID19 Spreadsheet.csv")#[,c(1:41)]
dt <- dt[,c(1:12,14,16,18:32,13,15,17,33:52)]

#####################
#### format data ####
#####################
#### Numeric columns, blanks are replaced with NA
#### N/A values are replaced with NA
#### All treatments are transformed into 'Y' or 'N' binary categorical variables
dt$Outcome <- gsub("D", "0", dt$Outcome)
dt$Outcome <- gsub("S", "1", dt$Outcome)
dt$Outcome <- as.numeric(dt$Outcome)
dt$Weight <- as.numeric(dt$Weight)
dt$O2_Saturation <- gsub("%", "", dt$O2_Saturation)
dt$O2_Saturation <- as.numeric(dt$O2_Saturation)
dt$HbA1C <- gsub("%", "", dt$HbA1C)
dt$HbA1C <- as.numeric(dt$HbA1C)
colnames(dt)[2] <- "ID"
colnames(dt)
dt[dt==''|dt==' ']<-"N"
#### Replace N/A with NA
dt[dt=='N/A']<-NA
#### convert columns into a binary classification
dt[grep("[A-Za-z][A-Za-z]", dt$`Secretagogues A`),]$`Secretagogues A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Secretagogues B`),]$`Secretagogues B` <- "Y"
dt[dt$Meds_ACEI_ARB == "No meds",]$Meds_ACEI_ARB <- "N"
dt[grep("[A-Za-z][A-Za-z]", dt$Meds_ACEI_ARB),]$Meds_ACEI_ARB <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Glycosidase_inhibitors A`),]$`Glycosidase_inhibitors A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Glycosidase_inhibitors B`),]$`Glycosidase_inhibitors B` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`DPP4_inhibitor B`),]$`DPP4_inhibitor B` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`DPP4_inhibitor A`),]$`DPP4_inhibitor A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`TZD A`)]$`TZD A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`TZD B`)]$`TZD B` <- "Y"
dt[grep("[0-9]", dt$Smoking_history),]$Smoking_history <- "Y"
dt[grep("[0-9]", dt$Hypertension),]$Hypertension <- "Y"
dt[grep("[0-9]", dt$CAD_years),]$CAD_years <- "Y"
str(dt)

#### record counts ####
#######################
frequenciesdyn("dt", "MRN")
frequenciesdyn("dt", "ID")
frequenciesdyn("dt", "ID, MRN")
frequenciesdyn("dt", "Outcome")
frequenciesdyn("dt", "Intervention")
frequenciesdyn("dt", "Sex")
frequenciesdyn("dt", "Age")
frequenciesdyn("dt", "Diagnosis")
frequenciesdyn("dt", "Smoking_history")
frequenciesdyn("dt", "OSA")

mean(dt$Age)

mean(dt[!is.na(dt$BMI),]$BMI)
ave(dt[!is.na(dt$BMI),]$BMI)


##############################################################################################################################
##############################################################################################################################
#### Analysis with medications taken Before (B) admission ####################################################################
##############################################################################################################################
##############################################################################################################################

######################################################################
#### Logistic regression outcome (death or survival) as dependent ####
######################################################################
#### stepwise increase
######################
glm.fit <- glm(Outcome ~ `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose, data = dt, family = binomial);summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP + D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)

#### stepwise decrease
######################
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP + D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose + CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG + Glucose, data = dt, family = binomial);summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B` +
                 + CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + 
                 Hypertension + Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B` + Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin B`, data = dt, family = binomial); summary(glm.fit)

#### selected variations of independent variabbles ####
#######################################################
glm.fit <- glm(Outcome ~ CAD_years + `Metformin B` + Age + Hypertension + Hyperlipidemia,
               data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years + `Metformin B` + Age + Hypertension + Hyperlipidemia + 
                 O2_Saturation + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years + `Metformin B` + Age + Hypertension + Hyperlipidemia + 
                 O2_Saturation + CHOL + TG + LDL_C + Glucose, data = dt, family = binomial); summary(glm.fit)
#### with lab values ####
glm.fit <- glm(Outcome ~ CHOL + TG + HDL_C + LDL_C + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin B` + Glucose, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Glucose, data = dt, family = binomial); summary(glm.fit)

#### Each variable alone ####
#############################
glm.fit <- glm(Outcome ~ Age, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Weight, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ BMI, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ O2_Saturation, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Secretagogues B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Glucose, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ HbA1C, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Insulin B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `TZD B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `DPP4_inhibitor B`, data = dt, family = binomial); summary(glm.fit)

#### Isolate people taking metformin, insulin, all other medications, metformin + insulin, insulin + others, metformin + others ####
####################################################################################################################################
#### people taking only metformin and no other medication ####
met <- dt[#dt$`Metformin B` == "Y" &
         (!dt$`Insulin B` == "Y") &
         (!dt$`Secretagogues B` == "Y") &
         (!dt$`TZD B` == "Y") &
         (!dt$`Glycosidase_inhibitors B` == "Y") &
             # (!dt$`GLP_1 B` == "Y") & # column consists of only NA
             # (!dt$`SGLT_2_inhibitor B` == "Y") # column consists of only NA
         (!dt$`DPP4_inhibitor B` == "Y"),]
# met[met==' '| met=='']<-NA
glm.fit <- glm(Outcome ~ `Metformin B`, data = met, family = binomial); summary(glm.fit)

met$Outcome <- gsub(0, "Deceased", met$Outcome)
met$Outcome <- gsub(1, "Survived", met$Outcome)
ggplot(met, aes(`Metformin B`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking only insulin and no other medication ####
ins <- dt[(!dt$`Metformin B` == "Y") &
          # (!dt$`Insulin B` == "Y") &
          (!dt$`Secretagogues B` == "Y") &
          (!dt$`TZD B` == "Y") &
          (!dt$`Glycosidase_inhibitors B` == "Y") &
          (!dt$`DPP4_inhibitor B` == "Y"),]
# ins[ins==' '| ins=='']<-NA
glm.fit <- glm(Outcome ~ `Insulin B`, data = ins, family = binomial); summary(glm.fit)

ins$Outcome <- gsub(0, "Deceased", ins$Outcome)
ins$Outcome <- gsub(1, "Survived", ins$Outcome)
ggplot(ins, aes(`Insulin B`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking no metformin and other medication ####
nomet <- dt[(!dt$`Metformin B` == "Y"),]
nomet$other <- "no_other"
nomet[(nomet$`Insulin B` == "Y")|(nomet$`Secretagogues B` == "Y")|(nomet$`TZD B` == "Y")|(nomet$`Glycosidase_inhibitors B` == "Y")|(nomet$`DPP4_inhibitor B` == "Y"),]$other <- "other"
# nomet[nomet==' '| nomet=='']<-NA
glm.fit <- glm(Outcome ~ other, data = nomet, family = binomial); summary(glm.fit)

nomet$Outcome <- gsub(0, "Deceased", nomet$Outcome)
nomet$Outcome <- gsub(1, "Survived", nomet$Outcome)
ggplot(nomet, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()


#### people taking no insulin and other medication ####
noins <- dt[(!dt$`Insulin B` == "Y"),]
noins$other <- "no_other"
noins[(noins$`Metformin B` == "Y")|(noins$`Secretagogues B` == "Y")|(noins$`TZD B` == "Y")|(noins$`Glycosidase_inhibitors B` == "Y")|(noins$`DPP4_inhibitor B` == "Y"),]$other <- "other"
# noins[nomet==' '| noins=='']<-NA
glm.fit <- glm(Outcome ~ other, data = noins, family = binomial); summary(glm.fit)

noins$Outcome <- gsub(0, "Deceased", noins$Outcome)
noins$Outcome <- gsub(1, "Survived", noins$Outcome)
ggplot(noins, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking no insulin, no metformin, and other medication ####
noins <- dt[(!dt$`Insulin B` == "Y"),]
noins$other <- "no_other"
noins[(!noins$`Metformin B` == "Y") & ((noins$`Secretagogues B` == "Y")|(noins$`TZD B` == "Y")|(noins$`Glycosidase_inhibitors B` == "Y")|(noins$`DPP4_inhibitor B` == "Y")),]$other <- "other"
# noins[nomet==' '| noins=='']<-NA
glm.fit <- glm(Outcome ~ other, data = noins, family = binomial); summary(glm.fit)

noins$Outcome <- gsub(0, "Deceased", noins$Outcome)
noins$Outcome <- gsub(1, "Survived", noins$Outcome)
ggplot(noins, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()


#### people taking insulin and metformin, and no other medication ####
metins <- dt
metins$met_ins <- "N"
metins[((metins$`Insulin B` == "Y")&(metins$`Metformin B` == "Y")),]$met_ins <- "Y"
metins <- metins[((!metins$`Secretagogues B` == "Y")&(!metins$`TZD B` == "Y")&(!metins$`Glycosidase_inhibitors B` == "Y")&(!metins$`DPP4_inhibitor B` == "Y")),]
metins$met_ins
length(metins[metins$met_ins == "Y",]$met_ins)
glm.fit <- glm(Outcome ~ met_ins, data = metins, family = binomial); summary(glm.fit)

metins$Outcome <- gsub(0, "Deceased", metins$Outcome)
metins$Outcome <- gsub(1, "Survived", metins$Outcome)
ggplot(metins, aes(`met_ins`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#######################################################
#### Logistic regression intervention as dependent ####
#######################################################
dt$vent <- 0
dt[dt$Intervention == "NIV" | dt$Intervention == "Invasive Ventilation",]$vent <- 1

#### stepwise increase
######################
glm.fit <- glm(vent ~ `Metformin B`, data = dt, family = binomial); summary(glm.fit)

#### Each variable alone ####
#############################
glm.fit <- glm(vent ~ Age, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Weight, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ BMI, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ O2_Saturation, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ CAD_years, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ `Secretagogues B`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ Glucose, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ HbA1C, data = dt, family = binomial); summary(glm.fit)

############################################################################################################################################
############################################################################################################################################
#### Repeat with medications after admission ###############################################################################################
############################################################################################################################################
############################################################################################################################################

######################################################################
#### Logistic regression outcome (death or survival) as dependent ####
######################################################################
#### stepwise increase
######################
glm.fit <- glm(Outcome ~ `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose, data = dt, family = binomial);summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP + D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)

#### stepwise decrease
######################
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP + D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose + CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG + Glucose, data = dt, family = binomial);summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A` +
                 + CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history + `Secretagogues A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB + Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia + Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + 
                 Hypertension + Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A` + Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + CAD_years + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + O2_Saturation + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + BMI + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + Weight + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Age + `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin A`, data = dt, family = binomial); summary(glm.fit)

#### selected variations of independent variabbles ####
#######################################################
glm.fit <- glm(Outcome ~ CAD_years + `Metformin A` + Age + Hypertension + Hyperlipidemia,
               data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years + `Metformin A` + Age + Hypertension + Hyperlipidemia + 
                 O2_Saturation + CHOL + TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years + `Metformin A` + Age + Hypertension + Hyperlipidemia + 
                 O2_Saturation + CHOL + TG + LDL_C + Glucose, data = dt, family = binomial); summary(glm.fit)
#### with lab values ####
glm.fit <- glm(Outcome ~ CHOL + TG + HDL_C + LDL_C + Glucose + CRP + D_dimmer + HbA1C, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin A` + Glucose, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Glucose, data = dt, family = binomial); summary(glm.fit)

#### Each variable alone ####
#############################
glm.fit <- glm(Outcome ~ Age, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Weight, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ BMI, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ O2_Saturation, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CAD_years, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Hypertension, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Hyperlipidemia, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Meds_ACEI_ARB, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Smoking_history, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Secretagogues A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CHOL, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ TG, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ Glucose, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ CRP, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ D_dimmer, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ HbA1C, data = dt, family = binomial); summary(glm.fit)

###################################################################################################
#### Additional information on the number of people takin metformin before and after admission ####
###################################################################################################
dt$`Metformin A`
length(dt[dt$`Metformin A` == "Y",]$`Metformin A`)

#### outcomes of people taking metformin before admission
dt[dt$`Metformin B` == "Y",]$Outcome
#### outcomes of people taking metformin after admission
dt[dt$`Metformin A` == "Y",]$Outcome
#### overlap between people taking metformin before and after admission
dt[dt$`Metformin B` == "Y" & dt$`Metformin A` == "Y",]$Outcome
intersect(dt[dt$`Metformin B` == "Y",]$MRN, dt[dt$`Metformin A` == "Y",]$MRN)

#### people not taking metformin prior to admission but taking it after admission
dt[(!dt$`Metformin B` == "Y") & dt$`Metformin A` == "Y",]$Outcome
dt[(!dt$`Metformin B` == "Y") & dt$`Metformin A` == "Y",]$MRN

#### people taking metformin prior to admission but not after admission
dt[(dt$`Metformin B` == "Y") & (!dt$`Metformin A` == "Y"),]$Outcome
dt[(dt$`Metformin B` == "Y") & (!dt$`Metformin A` == "Y"),]$MRN
length(dt[dt$`Metformin A` == "Y",]$Outcome)

glm.fit <- glm(Outcome ~ `Metformin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin B`, data = dt, family = binomial); summary(glm.fit)

dt[dt$`Metformin B` == "Y",]$Outcome
dt[dt$`Metformin A` == "Y",]$Outcome

####################################################################################################################################
#### Isolate people taking metformin, insulin, all other medications, metformin + insulin, insulin + others, metformin + others ####
####################################################################################################################################
#### people taking only metformin and no other medication ####
met <- dt[#dt$`Metformin A` == "Y" &
    (!dt$`Insulin A` == "Y") &
    (!dt$`Secretagogues A` == "Y") &
    (!dt$`TZD A` == "Y") &
    (!dt$`Glycosidase_inhibitors A` == "Y") &
    # (!dt$`GLP_1 A` == "Y") & # column consists of only NA
    # (!dt$`SGLT_2_inhibitor A` == "Y") # column consists of only NA
    (!dt$`DPP4_inhibitor A` == "Y"),]
# met[met==' '| met=='']<-NA
glm.fit <- glm(Outcome ~ `Metformin A`, data = met, family = binomial); summary(glm.fit)

met$Outcome <- gsub(0, "Deceased", met$Outcome)
met$Outcome <- gsub(1, "Survived", met$Outcome)
ggplot(met, aes(`Metformin A`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()


#### people taking only insulin and no other medication ####
ins <- dt[(!dt$`Metformin A` == "Y") &
            # (!dt$`Insulin A` == "Y") &
            (!dt$`Secretagogues A` == "Y") &
            (!dt$`TZD A` == "Y") &
            (!dt$`Glycosidase_inhibitors A` == "Y") &
            (!dt$`DPP4_inhibitor A` == "Y"),]
# ins[ins==' '| ins=='']<-NA
glm.fit <- glm(Outcome ~ `Insulin A`, data = ins, family = binomial); summary(glm.fit)

ins$Outcome <- gsub(0, "Deceased", ins$Outcome)
ins$Outcome <- gsub(1, "Survived", ins$Outcome)
ggplot(ins, aes(`Insulin A`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking no metformin and other medication ####
nomet <- dt[(!dt$`Metformin A` == "Y"),]
nomet$other <- "no_other"
nomet[(nomet$`Insulin A` == "Y")|(nomet$`Secretagogues A` == "Y")|(nomet$`TZD A` == "Y")|(nomet$`Glycosidase_inhibitors A` == "Y")|(nomet$`DPP4_inhibitor A` == "Y"),]$other <- "other"
# nomet[nomet==' '| nomet=='']<-NA
glm.fit <- glm(Outcome ~ other, data = nomet, family = binomial); summary(glm.fit)

nomet$Outcome <- gsub(0, "Deceased", nomet$Outcome)
nomet$Outcome <- gsub(1, "Survived", nomet$Outcome)
ggplot(nomet, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking no insulin and other medication ####
noins <- dt[(!dt$`Insulin A` == "Y"),]
noins$other <- "no_other"
noins[(noins$`Metformin A` == "Y")|(noins$`Secretagogues A` == "Y")|(noins$`TZD A` == "Y")|(noins$`Glycosidase_inhibitors A` == "Y")|(noins$`DPP4_inhibitor A` == "Y"),]$other <- "other"
# noins[nomet==' '| noins=='']<-NA
glm.fit <- glm(Outcome ~ other, data = noins, family = binomial); summary(glm.fit)

noins$Outcome <- gsub(0, "Deceased", noins$Outcome)
noins$Outcome <- gsub(1, "Survived", noins$Outcome)
ggplot(noins, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking no insulin, no metformin, and other medication ####
noins <- dt[(!dt$`Insulin A` == "Y"),]
noins$other <- "no_other"
noins[(!noins$`Metformin A` == "Y")&((noins$`Secretagogues A` == "Y")|(noins$`TZD A` == "Y")|(noins$`Glycosidase_inhibitors A` == "Y")|(noins$`DPP4_inhibitor A` == "Y")),]$other <- "other"
# noins[nomet==' '| noins=='']<-NA
glm.fit <- glm(Outcome ~ other, data = noins, family = binomial); summary(glm.fit)

noins$Outcome <- gsub(0, "Deceased", noins$Outcome)
noins$Outcome <- gsub(1, "Survived", noins$Outcome)
ggplot(noins, aes(`other`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#### people taking insulin and metformin, and no other medication ####
metins <- dt
metins$met_ins <- "N"
metins[((metins$`Insulin A` == "Y")&(metins$`Metformin A` == "Y")),]$met_ins <- "Y"
metins <- metins[((!metins$`Secretagogues A` == "Y")&(!metins$`TZD A` == "Y")&(!metins$`Glycosidase_inhibitors A` == "Y")&(!metins$`DPP4_inhibitor A` == "Y")),]
metins$met_ins
length(metins[metins$met_ins == "Y",]$met_ins)
glm.fit <- glm(Outcome ~ met_ins, data = metins, family = binomial); summary(glm.fit)

metins$Outcome <- gsub(0, "Deceased", metins$Outcome)
metins$Outcome <- gsub(1, "Survived", metins$Outcome)
ggplot(metins, aes(`met_ins`, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#######################################################
#### Logistic regression intervention as dependent ####
#######################################################
dt$vent <- 0
dt[dt$Intervention == "NIV" | dt$Intervention == "Invasive Ventilation",]$vent <- 1

#### stepwise increase
######################
glm.fit <- glm(vent ~ `Metformin A`, data = dt, family = binomial); summary(glm.fit)

#### Each variable alone ####
#############################
glm.fit <- glm(vent ~ `Secretagogues A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ `Insulin A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ `TZD A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ `Glycosidase_inhibitors A`, data = dt, family = binomial); summary(glm.fit)
glm.fit <- glm(vent ~ `DPP4_inhibitor A`, data = dt, family = binomial); summary(glm.fit)

############################################
#### Make plots of diseases and outcome ####
############################################
gdt <- dt[,c("CAD_years", "Hypertension", "Hyperlipidemia", "Intervention", "vent", "Outcome")]
gdt$Outcome <- gsub(0, "Deceased", gdt$Outcome)
gdt$Outcome <- gsub(1, "Survived", gdt$Outcome)
gdt$vent <- gsub(0, "not_ventilated", gdt$vent)
gdt$vent <- gsub(1, "ventilated", gdt$vent)


### CAD
ggplot(gdt, aes(CAD_years, ..count..)) + 
geom_bar(aes(fill = Outcome), position = "dodge") +
theme_pubr()

### Hypertension
ggplot(gdt, aes(Hypertension, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

### Hyperlipidemia
ggplot(gdt, aes(Hyperlipidemia, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

### intervention 
ggplot(gdt, aes(Intervention, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

### vent 
ggplot(gdt, aes(vent, ..count..)) + 
  geom_bar(aes(fill = Outcome), position = "dodge") +
  theme_pubr()

#



##########################################################################################
#### Loop through all combinations of data and return the values that are significant ####
##########################################################################################

#### Subset data to look at before or after medications ####
############################################################"SGLT_2_inhibitor B", "SGLT_2_inhibitor A", "GLP_1 B", "GLP_1 A", 
dt_sub <- dt[,c("Insulin B", "Metformin B", "Secretagogues B",
                "TZD B", "Glycosidase_inhibitors B",  "DPP4_inhibitor B",
                 "Insulin A", "Metformin A",  "Secretagogues A",
                "TZD A",  "Glycosidase_inhibitors A", "DPP4_inhibitor A",
                 "Steroid use", "Intervention", "Outcome")]
colnames(dt_sub) <- gsub(" ", "_", colnames(dt_sub))
write.table(dt_sub, "./results/1 subsetted_data_table.xls", row.names=FALSE, quote=FALSE, sep="\t")

frequenciesdyn("dt_sub", "Metformin_B")


#### Generate all combinations of columns ####
comb <- do.call(CJ, replicate(13, 0:1, FALSE))
comb[comb$V1 == 1, ]$V1 <- 1
comb[comb$V2 == 1, ]$V2 <- 2
comb[comb$V3 == 1, ]$V3 <- 3
comb[comb$V4 == 1, ]$V4 <- 4
comb[comb$V5 == 1, ]$V5 <- 5
comb[comb$V6 == 1, ]$V6 <- 6
comb[comb$V7 == 1, ]$V7 <- 7
comb[comb$V8 == 1, ]$V8 <- 8
comb[comb$V9 == 1, ]$V9 <- 9
comb[comb$V10 == 1, ]$V10 <- 10
comb[comb$V11 == 1, ]$V11 <- 11
comb[comb$V12 == 1, ]$V12 <- 12
comb[comb$V13 == 1, ]$V13 <- 13
#### create key ####
comb[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb)]
k <- comb$key

#### Loop through all column combinations and return only the values that are significant for each model ####
#############################################################################################################
glmcompileR <- function(DT, key, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    #### subset data table to contain only required columns
    spl <- strsplit(key[m], split = "_")[[1]]
    spl <- spl[!spl == "0"]
    newdt <- data.table()
    for(i in 1:length(spl)){
      temp <- DT[,c(as.numeric(spl[i])), with = FALSE]
      newdt <- cbind(newdt, temp)
      }
    #### add outcome variable to last column 
    newdt$Outcome <- DT$Outcome
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", p2))
    }else{
      fmla <- as.formula(paste0(colnames(newdt)[ncol(newdt)], " ~ ", nam))
      }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
        #### subset data table to retain only significant indepenent variables 
        dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
      }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:3)]),]
  return(compiled_dt)
}

compiled_dt <- glmcompileR(DT=dt_sub, key = k, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/2 complete_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

##########################################
#### subset out metformin monotherapy ####
##########################################
met <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
              (!dt_sub$`Secretagogues_A` == "Y") &
              (!dt_sub$`TZD_A` == "Y") &
              (!dt_sub$`Glycosidase_inhibitors_A` == "Y") &
              (!dt_sub$`DPP4_inhibitor_A` == "Y") &
              (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/3 metformin_monotherapy_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_A`, data = met, family = binomial); summary(glm.fit)

met <- dt_sub[(!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_B` == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/4 metformin_monotherapy_B_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_B`, data = met, family = binomial); summary(glm.fit)

#### metformin monotherapy before and after
met <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_A` == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$Steroid_use == "Y") &
                (!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_B` == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/5 metformin_monotherapy_BA_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_B`, data = met, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin_A`, data = met, family = binomial); summary(glm.fit)


################################################################
#### people taking no metformin but taking other medication ####
################################################################
#### no metformin after
nomet <- dt_sub[(!dt_sub$`Metformin_A` == "Y"),]
nomet$other <- "no_other"
nomet[(nomet$`Insulin_A` == "Y")|(nomet$`Secretagogues_A` == "Y")|(nomet$`TZD_A` == "Y")|(nomet$`Glycosidase_inhibitors_A` == "Y")|(nomet$`DPP4_inhibitor_A` == "Y")|(nomet$Steroid_use == "Y"),]$other <- "other"
nomet
write.table(nomet, "./results/6 no_metformin_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = nomet, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
nomet <- nomet[nomet$other == "other",]
nomet2 <- nomet[,c(1,3, 5:7, 9:16), with = FALSE]
nomet2
write.table(nomet2, "./results/7 no_metformin_multivariant_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=nomet2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/8 no_met_complete_compiled_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

##########################################
#### subset out Acarbose monotherapy #####
##########################################
#### acarbose monotherapy after
aca <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$Metformin_A == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/9 Acarbose_monotherapy_dt_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_A`, data = aca, family = binomial); summary(glm.fit)

#### acarbose monotherapy before
aca <- dt_sub[(!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$Metformin_B == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/10 Acarbose_monotherapy_dt_B.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_A`, data = aca, family = binomial); summary(glm.fit)

#### Acarbose monotherapy before and after
aca <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$Metformin_A == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$Metformin_B == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/11 Acarbose_monotherapy_dt_BA.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_A`, data = aca, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_B`, data = aca, family = binomial); summary(glm.fit)

#############################################################
#### subset out no Acarbose but taking other medications ####
#############################################################
#### no Acarbose taking others after
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_A == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_A` == "Y")|
        (noaca$`Secretagogues_A` == "Y")|
        (noaca$`TZD_A` == "Y")|
        (noaca$`Glycosidase_inhibitors_A` == "Y")|
        (noaca$Metformin_A == "Y")|
        (noaca$`DPP4_inhibitor_A` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
write.table(noaca, "./results/12 no_Acarbose_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 5:9, 12:16), with = FALSE]
noaca2
write.table(noaca2, "./results/13 no_Acarbose_multivariant_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/14 no_Acarbose_complete_compiled_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#### no Acarbose but taking other medications before
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_B == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_B` == "Y")|
        (noaca$`Secretagogues_B` == "Y")|
        (noaca$`TZD_B` == "Y")|
        (noaca$`Glycosidase_inhibitors_B` == "Y")|
        (noaca$Metformin_B == "Y")|
        (noaca$`DPP4_inhibitor_B` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
write.table(noaca, "./results/15 no_Acarbose_before_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 6:9, 11:16), with = FALSE]
noaca2
write.table(noaca2, "./results/16 no_Acarbose_before_multivariant_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/17 no_Acarbosee_before_complete_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#### no Acarbose taking others before and after
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_B == "Y"),]
noaca <- noaca[(!noaca$Glycosidase_inhibitors_A == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_B` == "Y")|
        (noaca$`Secretagogues_B` == "Y")|
        (noaca$`TZD_B` == "Y")|
        (noaca$Metformin_B == "Y")|
        (noaca$`Glycosidase_inhibitors_B` == "Y")|
        (noaca$`DPP4_inhibitor_B` == "Y")|
        (noaca$`Insulin_A` == "Y")|
        (noaca$`Secretagogues_A` == "Y")|
        (noaca$`TZD_A` == "Y")|
        (noaca$Metformin_A == "Y")|
        (noaca$`Glycosidase_inhibitors_A` == "Y")|
        (noaca$`DPP4_inhibitor_A` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
# noaca
write.table(noaca, "./results/18 no_Acarbosee_before_after_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(9, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 6:9, 12:16), with = FALSE]
noaca2
write.table(noaca2, "./results/19 no_Acarbose_before_after_multivariant_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/20 no_Acarbose_before_after_complete_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#############################################################################
#### subset out only records taking glycosidase before and after therapy ####
#############################################################################
dt_sub$glyco_diff <- "NA"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$glyco_diff <- "Before_After"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "N"),]$glyco_diff <- "Before"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$glyco_diff <- "After"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "N"),]$glyco_diff <- "Never"
dt_sub[dt_sub$glyco_diff == "NA",]$glyco_diff <- "After"
write.table(dt_sub, "./results/21 glyco_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_sub$Before_After <- "N"
dt_sub$Before <- "N"
dt_sub$After <- "N"
dt_sub$Never <- "N"
dt_sub[dt_sub$glyco_diff == "Before_After",]$Before_After <- "Y"
dt_sub[dt_sub$glyco_diff == "Before",]$Before <- "Y"
dt_sub[dt_sub$glyco_diff == "After",]$After <- "Y"
dt_sub[dt_sub$glyco_diff == "Never",]$Never <- "Y"
write.table(dt_sub, "./results/22 glyco_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_glyco <- dt_sub[,c(17:20,14:15)]
dt_glyco
#### Generate all combinations of columns ####
comb <- do.call(CJ, replicate(4, 0:1, FALSE))
comb[comb$V1 == 1, ]$V1 <- 1
comb[comb$V2 == 1, ]$V2 <- 2
comb[comb$V3 == 1, ]$V3 <- 3
comb[comb$V4 == 1, ]$V4 <- 4
#### create key ####
comb[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb)]
k3 <- comb$key

compiled <- glmcompileR(DT = dt_glyco, key = k3, significant = "F")
compiled
write.table(compiled, "./results/23 glyco_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled)

#### adjust columns to include people taking glycosidase inhibitors before and after in the before or after calculations ####
#############################################################################################################################
dt_sub$pre <- "N"
dt_sub$post <- "N"
#### pre/post and pre
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y")|
         (dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "N"),]$pre <- "Y"
#### pre/post and post
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y")|
         (dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$post <- "Y"
write.table(dt_sub, "./results/24 glyco_grouped.xls", row.names=FALSE, quote=FALSE, sep="\t")

glm.fit <- glm(Outcome ~ pre, data = dt_sub, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ post, data = dt_sub, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ pre + post, data = dt_sub, family = binomial); summary(glm.fit)


###########################################################################
#### subset out only records taking metformin before and after therapy ####
###########################################################################
dt_sub$met_diff <- "NA"
dt_sub[(dt_sub$Metformin_B == "Y" & dt_sub$Metformin_A == "Y"),]$met_diff <- "Before_After"
dt_sub[(dt_sub$Metformin_B == "Y" & dt_sub$Metformin_A == "N"),]$met_diff <- "Before"
dt_sub[(dt_sub$Metformin_B == "N" & dt_sub$Metformin_A == "Y"),]$met_diff <- "After"
dt_sub[(dt_sub$Metformin_B == "N" & dt_sub$Metformin_A == "N"),]$met_diff <- "Never"
# dt_sub[dt_sub$glyco_met == "NA",]$met_diff <- "After"
write.table(dt_sub, "./results/25 met_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_sub$met_Before_After <- "N"
dt_sub$met_Before <- "N"
dt_sub$met_After <- "N"
dt_sub$met_Never <- "N"
dt_sub[dt_sub$met_diff == "Before_After",]$met_Before_After <- "Y"
dt_sub[dt_sub$met_diff == "Before",]$met_Before <- "Y"
dt_sub[dt_sub$met_diff == "After",]$met_After <- "Y"
dt_sub[dt_sub$met_diff == "Never",]$met_Never <- "Y"
write.table(dt_sub, "./results/26 met_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_met <- dt_sub[,c(24,27,14:15)]
dt_met
#### Generate all combinations of columns ####
comb <- do.call(CJ, replicate(2, 0:1, FALSE))
comb[comb$V1 == 1, ]$V1 <- 1
comb[comb$V2 == 1, ]$V2 <- 2
#### create key ####
comb[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb)]
k3 <- comb$key

compiled <- glmcompileR(DT = dt_met, key = k3, significant = "F")
compiled
write.table(compiled, "./results/27 met_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled)


###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
#### repeat with ventilation as the outcome ###################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
dt <- fread("./data/2020.5.13 English DM COVID19 Spreadsheet.csv")#[,c(1:41)]
dt <- dt[,c(1:12,14,16,18:32,13,15,17,33:52)]
# dt$Outcome <- gsub("D", "0", dt$Outcome)
# dt$Outcome <- gsub("S", "1", dt$Outcome)
# dt$Outcome <- as.numeric(dt$Outcome)
dt$Weight <- as.numeric(dt$Weight)
dt$O2_Saturation <- gsub("%", "", dt$O2_Saturation)
dt$O2_Saturation <- as.numeric(dt$O2_Saturation)
dt$HbA1C <- gsub("%", "", dt$HbA1C)
dt$HbA1C <- as.numeric(dt$HbA1C)
colnames(dt)[2] <- "ID"
colnames(dt)
dt[dt==''|dt==' ']<-"N"
#### Replace N/A with NA
dt[dt=='N/A']<-NA
#### convert columns into a binary classification
dt[grep("[A-Za-z][A-Za-z]", dt$`Secretagogues A`),]$`Secretagogues A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Secretagogues B`),]$`Secretagogues B` <- "Y"
dt[dt$Meds_ACEI_ARB == "No meds",]$Meds_ACEI_ARB <- "N"
dt[grep("[A-Za-z][A-Za-z]", dt$Meds_ACEI_ARB),]$Meds_ACEI_ARB <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Glycosidase_inhibitors A`),]$`Glycosidase_inhibitors A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`Glycosidase_inhibitors B`),]$`Glycosidase_inhibitors B` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`DPP4_inhibitor B`),]$`DPP4_inhibitor B` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`DPP4_inhibitor A`),]$`DPP4_inhibitor A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`TZD A`)]$`TZD A` <- "Y"
dt[grep("[A-Za-z][A-Za-z]", dt$`TZD B`)]$`TZD B` <- "Y"
dt[grep("[0-9]", dt$Smoking_history),]$Smoking_history <- "Y"
dt[grep("[0-9]", dt$Hypertension),]$Hypertension <- "Y"
dt[grep("[0-9]", dt$CAD_years),]$CAD_years <- "Y"
str(dt)
dt_sub <- dt[,c("Insulin B", "Metformin B", "Secretagogues B",
                "TZD B", "Glycosidase_inhibitors B",  "DPP4_inhibitor B",
                "Insulin A", "Metformin A",  "Secretagogues A",
                "TZD A",  "Glycosidase_inhibitors A", "DPP4_inhibitor A",
                "Steroid use", "Intervention", "Outcome")]
colnames(dt_sub) <- gsub(" ", "_", colnames(dt_sub))

setnames(dt_sub, colnames(dt_sub), c("Insulin_B",  "Metformin_B",  "Secretagogues_B",  "TZD_B", "Glycosidase_inhibitors_B", "DPP4_inhibitor_B", "Insulin_A", "Metformin_A",             
                             "Secretagogues_A","TZD_A", "Glycosidase_inhibitors_A", "DPP4_inhibitor_A",        
                             "Steroid_use",  "Outcome", "Survival"))
# dt_sub$Survival <- as.character(dt_sub$Survival)
# dt_sub[dt_sub$Survival == 1,]$Survival <- "survived"
# dt_sub[dt_sub$Survival == 0,]$Survival <- "deceased"

dt_sub[!duplicated(dt_sub$Outcome),]$Outcome
dt_sub[!(dt_sub$Outcome == "Invasive Ventilation" | dt_sub$Outcome == "NIV"),]$Outcome <- "0"
dt_sub[dt_sub$Outcome == "Invasive Ventilation" | dt_sub$Outcome == "NIV",]$Outcome <- "1"
dt_sub[!duplicated(dt_sub$Outcome),]$Outcome
dt_sub$Outcome <- as.numeric(dt_sub$Outcome) 

##########################################
#### subset out metformin monotherapy ####
##########################################
met <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_A` == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/28 vent_metformin_monotherapy_dt_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_A`, data = met, family = binomial); summary(glm.fit)

#### metformin monotherapy before
met <- dt_sub[(!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_B` == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/29 vent_metformin_monotherapy_dt_B.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_B`, data = met, family = binomial); summary(glm.fit)

#### metformin monotherapy before and after
met <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_A` == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$Steroid_use == "Y") &
                (!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$`Glycosidase_inhibitors_B` == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
met
write.table(met, "./results/30 vent_metformin_monotherapy_BA_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Metformin_B`, data = met, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Metformin_A`, data = met, family = binomial); summary(glm.fit)

#####################################################################
#### people taking no metformin but are taking other medications ####
#####################################################################
#### no metformin after
nomet <- dt_sub[(!dt_sub$`Metformin_A` == "Y"),]
nomet$other <- "no_other"
nomet[(nomet$`Insulin_A` == "Y")|(nomet$`Secretagogues_A` == "Y")|(nomet$`TZD_A` == "Y")|(nomet$`Glycosidase_inhibitors_A` == "Y")|(nomet$`DPP4_inhibitor_A` == "Y")|(nomet$Steroid_use == "Y"),]$other <- "other"
nomet
write.table(nomet, "./results/31 vent_no_metformin_dt_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = nomet, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
nomet <- nomet[nomet$other == "other",]
nomet2 <- nomet[,c(1,3, 5:7, 9:16), with = FALSE]
nomet2
write.table(nomet2, "./results/32 vent_no_metformin_multivariant_dt_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=nomet2, key = k2, significant = "T")
compiled_dt
# write.table(compiled_dt, "./results/33 vent_no_met_complete_compiled_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

##########################################
#### subset out Acarbose monotherapy #####
##########################################
#### acarbose monotherapy after
aca <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$Metformin_A == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/34 vent_Acarbose_monotherapy_dt_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_A`, data = aca, family = binomial); summary(glm.fit)

#### acarbose monotherapy before
aca <- dt_sub[(!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$Metformin_B == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/35 vent_Acarbose_monotherapy_dt_B.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_B`, data = aca, family = binomial); summary(glm.fit)

#### Acarbose monotherapy before and after
aca <- dt_sub[(!dt_sub$`Insulin_A` == "Y") &
                (!dt_sub$`Secretagogues_A` == "Y") &
                (!dt_sub$`TZD_A` == "Y") &
                (!dt_sub$Metformin_A == "Y") &
                (!dt_sub$`DPP4_inhibitor_A` == "Y") &
                (!dt_sub$`Insulin_B` == "Y") &
                (!dt_sub$`Secretagogues_B` == "Y") &
                (!dt_sub$`TZD_B` == "Y") &
                (!dt_sub$Metformin_B == "Y") &
                (!dt_sub$`DPP4_inhibitor_B` == "Y") &
                (!dt_sub$Steroid_use == "Y"),]
aca
write.table(aca, "./results/36 vent_Acarbose_monotherapy_AB_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_A`, data = aca, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ `Glycosidase_inhibitors_B`, data = aca, family = binomial); summary(glm.fit)

#############################################################
#### subset out no Acarbose but taking other medications ####
#############################################################
#### no Acarbose taking others after
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_A == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_A` == "Y")|
        (noaca$`Secretagogues_A` == "Y")|
        (noaca$`TZD_A` == "Y")|
        (noaca$`Glycosidase_inhibitors_A` == "Y")|
        (noaca$Metformin_A == "Y")|
        (noaca$`DPP4_inhibitor_A` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
write.table(noaca, "./results/37 vent_no_Acarbose_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 5:9, 12:16), with = FALSE]
noaca2
write.table(noaca2, "./results/38 vent_no_Acarbose_multivariant_A_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/39 vent_no_Acarbose_complete_compiled_A.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#### no Acarbose taking other medications before
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_B == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_B` == "Y")|
        (noaca$`Secretagogues_B` == "Y")|
        (noaca$`TZD_B` == "Y")|
        (noaca$`Glycosidase_inhibitors_B` == "Y")|
        (noaca$Metformin_B == "Y")|
        (noaca$`DPP4_inhibitor_B` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
write.table(noaca, "./results/40 vent_no_Acarbose_before_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(10, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
comb2[comb2$V10 == 1, ]$V10 <- 10
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 6:9, 11:16), with = FALSE]
noaca2
write.table(noaca2, "./results/41 vent_no_Acarbose_before_multivariant_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/42 vent_no_Acarbosee_before_complete_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#### no Acarbose taking other medications before and after
noaca <- dt_sub[(!dt_sub$Glycosidase_inhibitors_B == "Y"),]
noaca <- noaca[(!noaca$Glycosidase_inhibitors_A == "Y"),]
noaca$other <- "no_other"
noaca[(noaca$`Insulin_B` == "Y")|
        (noaca$`Secretagogues_B` == "Y")|
        (noaca$`TZD_B` == "Y")|
        (noaca$Metformin_B == "Y")|
        (noaca$`Glycosidase_inhibitors_B` == "Y")|
        (noaca$`DPP4_inhibitor_B` == "Y")|
        (noaca$`Insulin_A` == "Y")|
        (noaca$`Secretagogues_A` == "Y")|
        (noaca$`TZD_A` == "Y")|
        (noaca$Metformin_A == "Y")|
        (noaca$`Glycosidase_inhibitors_A` == "Y")|
        (noaca$`DPP4_inhibitor_A` == "Y")|
        (noaca$Steroid_use == "Y"),]$other <- "other"
# noaca
write.table(noaca, "./results/43 vent_no_Acarbosee_before_after_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
glm.fit <- glm(Outcome ~ other, data = noaca, family = binomial); summary(glm.fit)

#### Generate all combinations of columns ####
comb2 <- do.call(CJ, replicate(9, 0:1, FALSE))
comb2[comb2$V1 == 1, ]$V1 <- 1
comb2[comb2$V2 == 1, ]$V2 <- 2
comb2[comb2$V3 == 1, ]$V3 <- 3
comb2[comb2$V4 == 1, ]$V4 <- 4
comb2[comb2$V5 == 1, ]$V5 <- 5
comb2[comb2$V6 == 1, ]$V6 <- 6
comb2[comb2$V7 == 1, ]$V7 <- 7
comb2[comb2$V8 == 1, ]$V8 <- 8
comb2[comb2$V9 == 1, ]$V9 <- 9
#### create key ####
comb2[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb2)]
k2 <- comb2$key
noaca <- noaca[noaca$other == "other",]
noaca2 <- noaca[,c(1:3, 6:9, 12:16), with = FALSE]
noaca2
write.table(noaca2, "./results/44 vent_no_Acarbose_before_after_multivariant_dt.xls", row.names=FALSE, quote=FALSE, sep="\t")
compiled_dt <- glmcompileR(DT=noaca2, key = k2, significant = "T")
compiled_dt
write.table(compiled_dt, "./results/45 vent_no_Acarbose_before_after_complete_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled_dt)

#############################################################################
#### subset out only records taking glycosidase before and after therapy ####
#############################################################################
dt_sub$glyco_diff <- "NA"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$glyco_diff <- "Before_After"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "N"),]$glyco_diff <- "Before"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$glyco_diff <- "After"
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "N"),]$glyco_diff <- "Never"
dt_sub[dt_sub$glyco_diff == "NA",]$glyco_diff <- "After"
write.table(dt_sub, "./results/46 vent_glyco_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_sub$Before_After <- "N"
dt_sub$Before <- "N"
dt_sub$After <- "N"
dt_sub$Never <- "N"
dt_sub[dt_sub$glyco_diff == "Before_After",]$Before_After <- "Y"
dt_sub[dt_sub$glyco_diff == "Before",]$Before <- "Y"
dt_sub[dt_sub$glyco_diff == "After",]$After <- "Y"
dt_sub[dt_sub$glyco_diff == "Never",]$Never <- "Y"
write.table(dt_sub, "./results/47 vent_glyco_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_glyco <- dt_sub[,c(17:20,14:15)]
dt_glyco
#### Generate all combinations of columns ####
comb <- do.call(CJ, replicate(4, 0:1, FALSE))
comb[comb$V1 == 1, ]$V1 <- 1
comb[comb$V2 == 1, ]$V2 <- 2
comb[comb$V3 == 1, ]$V3 <- 3
comb[comb$V4 == 1, ]$V4 <- 4
#### create key ####
comb[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb)]
k3 <- comb$key

compiled <- glmcompileR(DT = dt_glyco, key = k3, significant = "F")
compiled
write.table(compiled, "./results/48 vent_glyco_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled)

#### adjust columns to include people taking glycosidase inhibitors before and after in the before or after calculations ####
#############################################################################################################################
dt_sub$pre <- "N"
dt_sub$post <- "N"
#### pre/post and pre
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y")|
         (dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "N"),]$pre <- "Y"
#### pre/post and post
dt_sub[(dt_sub$Glycosidase_inhibitors_B == "Y" & dt_sub$Glycosidase_inhibitors_A == "Y")|
         (dt_sub$Glycosidase_inhibitors_B == "N" & dt_sub$Glycosidase_inhibitors_A == "Y"),]$post <- "Y"
write.table(dt_sub, "./results/49 vent_glyco_grouped.xls", row.names=FALSE, quote=FALSE, sep="\t")

glm.fit <- glm(Outcome ~ pre, data = dt_sub, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ post, data = dt_sub, family = binomial); summary(glm.fit)
glm.fit <- glm(Outcome ~ pre + post, data = dt_sub, family = binomial); summary(glm.fit)


###########################################################################
#### subset out only records taking metformin before and after therapy ####
###########################################################################
dt_sub$met_diff <- "NA"
dt_sub[(dt_sub$Metformin_B == "Y" & dt_sub$Metformin_A == "Y"),]$met_diff <- "Before_After"
dt_sub[(dt_sub$Metformin_B == "Y" & dt_sub$Metformin_A == "N"),]$met_diff <- "Before"
dt_sub[(dt_sub$Metformin_B == "N" & dt_sub$Metformin_A == "Y"),]$met_diff <- "After"
dt_sub[(dt_sub$Metformin_B == "N" & dt_sub$Metformin_A == "N"),]$met_diff <- "Never"
# dt_sub[dt_sub$glyco_met == "NA",]$met_diff <- "After"
write.table(dt_sub, "./results/50 vent_met_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_sub$met_Before_After <- "N"
dt_sub$met_Before <- "N"
dt_sub$met_After <- "N"
dt_sub$met_Never <- "N"
dt_sub[dt_sub$met_diff == "Before_After",]$met_Before_After <- "Y"
dt_sub[dt_sub$met_diff == "Before",]$met_Before <- "Y"
dt_sub[dt_sub$met_diff == "After",]$met_After <- "Y"
dt_sub[dt_sub$met_diff == "Never",]$met_Never <- "Y"
write.table(dt_sub, "./results/51 vent_met_annotated.xls", row.names=FALSE, quote=FALSE, sep="\t")

dt_met <- dt_sub[,c(24,27,14:15)]
dt_met
#### Generate all combinations of columns ####
comb <- do.call(CJ, replicate(2, 0:1, FALSE))
comb[comb$V1 == 1, ]$V1 <- 1
comb[comb$V2 == 1, ]$V2 <- 2
#### create key ####
comb[, key := do.call(paste, c(.SD, sep = "_")), .SDcols = names(comb)]
k3 <- comb$key

compiled <- glmcompileR(DT = dt_met, key = k3, significant = "F")
compiled
write.table(compiled, "./results/52 vent_met_compiled.xls", row.names=FALSE, quote=FALSE, sep="\t")
rm(compiled)




#########################################################################################################################
#########################################################################################################################
#### look at length of hospital stay between people taking metformin or not as an outcome ###############################
#########################################################################################################################
#########################################################################################################################
dt[dt$`Metformin A` == "Y",]$`Lenth of hospital stay`

dt2 <- dt[dt$Outcome == 1,]

ggplot(dt, aes(x=`Metformin A`, y=`Lenth of hospital stay`)) + 
  geom_violin(trim=FALSE)

p<-ggplot(dt, aes(x=`Lenth of hospital stay`, color=`Metformin A`)) +
  geom_density()#+
  # geom_vline(data=mu, aes(xintercept=grp.mean, color=`Metformin A`),
  #            linetype="dashed")
p

tiff(file = "length of hospital stay and metformin.tiff", width = 800, height = 500, units = "px", res = 100)
p
dev.off()

p<-ggplot(dt2, aes(x=`Lenth of hospital stay`, color=`Metformin A`)) +
  geom_density()#+
# geom_vline(data=mu, aes(xintercept=grp.mean, color=`Metformin A`),
#            linetype="dashed")
p

tiff(file = "length of hospital stay and metformin.tiff", width = 800, height = 500, units = "px", res = 100)
p
dev.off()

test <- suppressWarnings(pairwise.wilcox.test(dt$`Lenth of hospital stay`, as.factor(dt$`Metformin A`), conf.int = TRUE))
test$p.value

test <- suppressWarnings(pairwise.wilcox.test(dt2$`Lenth of hospital stay`, as.factor(dt2$`Metformin A`), conf.int = TRUE))
test$p.value


dt[,`Metformin A`:=ifelse(`Metformin A`=="N",0,1)]
glm.fit <- glm(`Metformin A` ~ `Lenth of hospital stay`, data = dt, family = binomial); summary(glm.fit)


dt2[,`Metformin A`:=ifelse(`Metformin A`=="N",0,1)]
glm.fit <- glm(`Metformin A` ~ `Lenth of hospital stay`, data = dt2, family = binomial); summary(glm.fit)







######################################################################################################################
######################################################################################################################
#### example code ####################################################################################################
######################################################################################################################
######################################################################################################################
rep(glm.fit$formula, nrow(dtf))
summary(glm.fit)$coefficients[,"Pr(>|z|)"]
rownames(summary(glm.fit)$coefficients)


# create dummy data set
dat <- data.frame(y = rnorm(10), x1 = rnorm(10)*2, x2 = rnorm(10)*3, x3 = rnorm(10)*4, scale = c(rep(1,5), rep(2,5)))

# create data frame to store results
results <- data.frame()

# loop through the scales and each variable
for(scale in unique(dat$scale)){
  for(var in names(dat)[c(-1,-length(dat))]){
    # dynamically generate formula
    fmla <- as.formula(paste0("y ~ ", var))
    
    # fit glm model
    fit <- glm(fmla, data=dat[dat$scale = scale,])
    
    ## capture summary stats
    intercept <- coef(summary(fit))[1]
    slope <- coef(summary(fit))[2]
    p.value <- coef(summary(fit))[8]
    AIC <- AIC(fit)
    Deviance <- deviance(fit)
    
    # get coefficents of fit
    cfit <- coef(summary(fit))
    
    # create temporary data frame
    df <- data.frame(var = var, scale = scale, intercept = cfit[1],
                     slope = cfit[2], p.value = cfit[8],
                     AIC = AIC(fit), Deviance = deviance(fit), stringsAsFactors = F)
    
    # bind rows of temporary data frame to the results data frame
    results <- rbind(results, df)
  }
}




cols <- c(6, 9, 10:12, 18, 20:25, 27)
# finDT <- data.table()
# for(i in 1:length(cols)){
# print(i)
# test <- shapiro.test(as.numeric(dt[[cols[i]]]))
# 
# if(test$p.value < 0.05){
#   distribution <- "Not normally distributed"
# }else{
#   distribution <- "normally distributed"
# }
# 
# DT <- data.table(method = test$method,
#                  pvalue = test$p.value,
#                  variable = colnames(dt[,cols[i], with = FALSE]),
#                  normality = distribution
#                  )
# finDT <- rbind(finDT,DT)
# 
# }
# finDT

normalityfindR <- function(DT, cols){
  finDT <- data.table()
  for(i in 1:length(cols)){
    
    test <- shapiro.test(as.numeric(DT[[cols[i]]]))
    
    if(test$p.value < 0.05){
      distribution <- "Not normally distributed"
    }else{
      distribution <- "normally distributed"
    }
    
    DT2 <- data.table(method = test$method,
                      pvalue = test$p.value,
                      variable = colnames(dt[,cols[i], with = FALSE]),
                      normality = distribution
    )
    finDT <- rbind(finDT,DT2)
  }
  return(finDT)
}

cols <- c(6, 9, 10:12, 18, 20:25, 27)
normalityfindR(DT = dt, cols = cols)





















fname <- "./data/2020.5.13 English DM COVID19 Spreadsheet.csv"
dt <- fread(fname)[,c(1:12,14,16,18:32,13,15,17,33:52)]
colnames(dt)[3] <- "ID"
names(dt) <- gsub(" ", "_", tolower(names(dt)))
#### convert outcome to a numeric column
dt[,outcome:=ifelse(outcome=="D",0,1)]
#### convert columns to a binary numeric column
dt[,c("weight","o2_saturation","hba1c") := .(as.numeric(weight),
                                             as.numeric(gsub("%","",o2_saturation)),
                                             as.numeric(gsub("%","",hba1c)))]
#### Convert columns into binary classifiers
dt[,c("secretagogues_b", "secretagogues_a", "glycosidase_inhibitors_b",
      "glycosidase_inhibitors_a", "dpp4_inhibitor_b", "dpp4_inhibitor_a",
      "tzd_b", "tzd_a", "meds_acei_arb", "statins", "life_style_modification",
      "cad_meds") := .(ifelse(secretagogues_b!='',"Y","N"),
                       ifelse(secretagogues_a!='',"Y","N"),
                       ifelse(glycosidase_inhibitors_b!='',"Y","N"),
                       ifelse(glycosidase_inhibitors_a!='',"Y","N"),
                       ifelse(dpp4_inhibitor_b!='',"Y","N"),
                       ifelse(dpp4_inhibitor_a!='',"Y","N"),
                       ifelse(tzd_b!='',"Y","N"),
                       ifelse(tzd_a!='',"Y","N"),
                       ifelse(meds_acei_arb!='',"Y","N"),
                       ifelse(statins!='',"Y","N"),
                       ifelse(life_style_modification!='',"Y","N"),
                       ifelse(cad_meds!='',"Y","N"))]
dt[,c("smoking_history") := .(ifelse(grepl("[0-9]",smoking_history) | smoking_history=="Y","Y","N"))]
dt[,c("hypertension") := .(ifelse(grepl("[0-9]",hypertension) | hypertension=="Y","Y","N"))]
dt[,c("cad_years") := .(ifelse(grepl("[0-9]",cad_years) | cad_years=="Y","Y","N"))]
#### clean missing values based on the information provided by the individuals who procured the dataset.
dt[dt==''|dt==' ']<-"N"
dt[dt=='N/A']<-NA
#### remove these columns that are empty
remove <- c("glp_1_a", "glp_1_b", 'osa', 'sglt_2_inhibitor_a', 'sglt_2_inhibitor_b') 
dt <- dt[,! ..remove]
#colnames(dt)
#str(dt)



dt_sub2 <- dt[,c("Secretagogues_b", "Secretagogues_a", "glycosidase_inhibitors_b",
                 "glycosidase_inhibitors_a", "dpp4_inhibitor_b", "dpp4_inhibitor_a",
                 "tzd_b", "tzd_a", "meds_acei_arb", "statins", "life_style_modification",
                 "cad_meds", "smoking_history", "hypertension", "cad_years","procalcitonin", 
                 "hyperlipidemia", "insulin_b", "insulin_a", "steroid_use") := .(ifelse(dt$secretagogues_b=='Y',1,0),
                                                                                 ifelse(dt$secretagogues_a=='Y',1,0),
                                                                                 ifelse(dt$glycosidase_inhibitors_b=='Y',1,0),
                                                                                 ifelse(dt$glycosidase_inhibitors_a=='Y',1,0),
                                                                                 ifelse(dt$dpp4_inhibitor_b=='Y',1,0),
                                                                                 ifelse(dt$dpp4_inhibitor_a=='Y',1,0),
                                                                                 ifelse(dt$tzd_b=='Y',1,0),
                                                                                 ifelse(dt$tzd_a=='Y',1,0),
                                                                                 ifelse(dt$meds_acei_arb=='Y',1,0),
                                                                                 ifelse(dt$statins=='Y',1,0),
                                                                                 ifelse(dt$life_style_modification=='Y',1,0),
                                                                                 ifelse(dt$cad_meds=='Y',1,0),
                                                                                 ifelse(dt$smoking_history=='Y',1,0),
                                                                                 ifelse(dt$hypertension=='Y',1,0),
                                                                                 ifelse(dt$cad_years=='Y',1,0),
                                                                                 ifelse(dt$procalcitonin=='Y',1,0),
                                                                                 ifelse(dt$hyperlipidemia=='Y',1,0),
                                                                                 ifelse(dt$insulin_b=='Y',1,0),
                                                                                 ifelse(dt$insulin_a=='Y',1,0),
                                                                                 ifelse(dt$steroid_use=='Y',1,0))]


dt_sub2 <- dt_sub2[,c(13:16, 29:44)][,c(1:9, 11:15, 17:20, 10, 16)]
dt_sub2 <- dt_sub2[, c(1:9, 11:14, 16:20)]
str(dt_sub2)
#dt_sub2[,outcome:=ifelse(outcome=="Y",1,0)] # transform metformin_b into a binary numeric column
k <- 1:(length(colnames(dt_sub2))-2)
k

DT <- dt_sub2
key <- k
independent_column <- 19
m <- 3

dependentglmR <- function(DT, key, independent_column, significant){
  compiled_dt <- data.table()
  pb <- txtProgressBar(min = 0, max = length(key), style = 3)
  for(m in 2:length(key)){
    print(m)
    #### subset data table to contain only required columns
    newdt <- DT[,c(as.numeric(key[m])), with = FALSE]
    #### add outcome variable to last column 
    newdt <- cbind(newdt, DT[,independent_column, with = FALSE])
    
    #### generate formula 
    nam <- colnames(newdt[, -(ncol(newdt)), with = FALSE])
    if(length(nam) > 1){
      plus <- "+"
      p1 <- paste(nam[-length(nam)], plus, sep = " ")
      p2 <- paste(p1, nam[length(nam)])
      fmla <- as.formula(paste0(p2, " ~ ", colnames(newdt)[ncol(newdt)]))
    }else{
      fmla <- as.formula(paste0(nam, " ~ ", colnames(newdt)[ncol(newdt)]))
    }
    
    #### fit glm model
    glm.fit <- glm(fmla, data=newdt, family = binomial)
    
    #### create data table 
    dtf <- data.table(coef(summary(glm.fit)))#[,"Pr(>|z|)"]
    dtf$names <- rownames(summary(glm.fit)$coefficients)
    chr <- as.character(glm.fit$formula)
    dtf$formula <- paste(chr[2], chr[1], chr[3:length(chr)])
    if(significant == "T"){
      #### subset data table to retain only significant indepenent variables 
      dtf <- dtf[dtf$`Pr(>|z|)` < 0.05 & (!dtf$names == "(Intercept)"),]
    }
    if(nrow(dtf) > 0){
      compiled_dt <- rbind(compiled_dt, dtf)
    }
    setTxtProgressBar(pb, m)
  }
  close(pb)
  compiled_dt <- compiled_dt[!duplicated(compiled_dt[,c(1:3)]),]
  compiled_dt <- compiled_dt[(!compiled_dt$names == "(Intercept)"),]
  
  compiled_dt$Significance <- "NA"
  compiled_dt[`Pr(>|z|)` < 0.05,]$Significance <- "Significant"
  compiled_dt[`Pr(>|z|)` > 0.05,]$Significance <- "Not-Significant"
  return(compiled_dt)
}

dependentglmR(DT = dt_sub2, key = k, independent_col = 18, significant = "F")

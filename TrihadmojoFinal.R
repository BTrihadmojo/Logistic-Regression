#SOC401 FINAL Assignment
#Bambang Trihadmojo

library(MASS)
library(psych)
library(moments)
library(tidyverse)
library(sjPlot)
library(arsenal)
library(car)
library(ggplot2)
library(FSA)
library(gvlma)
library(caret)

#Setting working firectory
setwd('')

#importing dataset
gcsvaccine <- read.csv('trihadmojofinal.csv', header = T)

#checking missing value
desc_table <- describe(gcsvaccine,na.rm = T) %>% mutate(nmiss=195-n) %>% 
  select(n,nmiss,mean,sd)
write.csv(desc_table, 'descriptive table.csv')

#descriptive statistic
tabbdesc <- CreateTableOne(vars = c("M1", "VA", "PS", "GE", "RQ", "RL",
                                    "CC", "GDPcapita", "effectiveness", 
                                    "importance", "safety", "GAVI"), 
                           factorVars = c("GAVI"), data = gcsvaccine)
tabdescprint <- print(tabbdesc2, showAllLevels = T)

write.csv(tabdescprint, "descriptive statistics.csv")

tab1 <- describe(gcsvaccine[,c("M1", "VA", "PS", "GE", "RQ", "RL",
                               "CC", "GDPcapita", "effectiveness", 
                               "importance", "safety", "GAVI")], na.rm = T,
                 interp = T, skew = T, ranges = F, trim = .1, type = 1)
write.csv(tab1, "tab1.csv")
hist(gcsvaccine$safety)

tab2 <- describe(gcsvaccine[,c("M1revlog", "VA", "PS", "GE", "RQ", "RL",
                               "CC", "GDPcapitalog", "effectivenesslog", 
                               "importancelog", "safetylog", "GAVI")], na.rm = T,
                 interp = T, skew = T, ranges = F, trim = .1, type = 1)
write.csv(tab2, "tab2.csv")

#generating histogram for skewed variables before transformation
par(mfrow = c(3,2))
hist(gcsvaccine$M1, main = NULL, xlab ="Measless vaccine coverage")
hist(gcsvaccine$GDPcapita, main = NULL, xlab = "GDP per capita")
hist(gcsvaccine$importance, main = NULL, xlab = "Perceptions of the importance of vaccines for children")
hist(gcsvaccine$safety, main = NULL, xlab = "Perceptions of the safety of vaccines")
hist(gcsvaccine$effectiveness, main = NULL, xlab ="Perceptions of the effectiveness of vaccines")

#transforming skewed variables
kurang <- function(x){(100-x)}
gcsvaccine$M1rev <- sapply(gcsvaccine$M1, kurang)
gcsvaccine$M1revlog <- log(gcsvaccine$M1rev)

gcsvaccine$GDPcapitalog <- log(gcsvaccine$GDPcapita)
gcsvaccine$importancelog <- log(gcsvaccine$importance+1)
gcsvaccine$safetylog <- log(gcsvaccine$safety+1)
gcsvaccine$effectivenesslog <- log(gcsvaccine$effectiveness+1)

#generating histogram for skewed variables after transformation
par(mfrow = c(3,2))
hist(gcsvaccine$M1revlog, main = NULL, xlab ="Measless vaccine coverage")
hist(gcsvaccine$GDPcapitalog, main = NULL, xlab = "GDP per capita")
hist(gcsvaccine$importancelog, main = NULL, xlab = "Perceptions of the importance of vaccines for children")
hist(gcsvaccine$safetylog, main = NULL, xlab = "Perceptions of the safety of vaccines")
hist(gcsvaccine$effectivenesslog, main = NULL, xlab ="Perceptions of the effectiveness of vaccines")

#Regression model for the relationship between governance and measles vaccine coverage
fit1 <- lm(M1revlog~ + VA + PS + GE + RQ + RL + CC, data = gcsvaccine)

#Regression model that include economic factor
fit2 <- lm(M1revlog~ + VA + PS + GE + RQ + RL + CC + GDPcapitalog, data = gcsvaccine)

#Regression model that includes vaccine hesitancy variables
fit3 <- lm(M1revlog~ + VA + PS + GE + RQ + RL + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog, data = gcsvaccine)

#Regression model that includes international aid for health
fit4 <- lm(M1revlog~ + VA + PS + GE + RQ + RL + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
tab_model(fit1, fit2, fit3, fit4, show.ci = .95, show.se = T, p.style = "stars")

#nonlinearity check
nl1 <- lm(M1revlog~ + VA + I(VA^2), data = gcsvaccine)
summary(nl1)

nl2 <- lm(M1revlog~ + PS + I(PS^2), data = gcsvaccine)
summary(nl2)

nl3 <- lmnl1 <- lm(M1revlog~ + GE + I(GE^2), data = gcsvaccine)
summary(nl3)

nl4 <-lm(M1revlog~ + RQ + I(RQ^2), data = gcsvaccine)
summary(nl4)
 
nl5 <- lm(M1revlog~ + RL + I(RL^2), data = gcsvaccine)
summary(nl5)

nl6 <- lm(M1revlog~ + CC + I(CC^2), data = gcsvaccine)
summary(nl6)

nl7 <- lm(M1revlog~ + GDPcapita + I(GDPcapita^2), data = gcsvaccine)
summary(nl7)

nl8 <- lm(M1revlog~ + effectiveness + I(effectiveness^2), data = gcsvaccine)
summary(nl8)

nl9 <- lm(M1revlog~ + safety + I(safety^2), data = gcsvaccine)
summary(nl9)

nl10 <- lm(M1revlog~ + importance + I(importance^2), data = gcsvaccine)
summary(nl10)

nl11 <- lm(M1revlog~ + as.factor(GAVI) + I(as.factor(GAVI^2)), data = gcsvaccine)
summary(nl11)

tab_model(nl1, nl3,  nl5, nl6, nl7, show.ci = .95, show.se = T, p.style = "stars")

fit5 <- lm(M1revlog~ + log(VA) + PS + log(GE) + RQ + log(RL) + log(CC) + 
             GDPcapitalog + effectiveness + safety + importance + 
             as.factor(GAVI), data = gcsvaccine)
tab_model(fit5, fit4, show.ci = .95, show.se = T, p.style = "stars")

#Heteroskedasicity check
gvlma(fit4)

#Multicolonierity check
car::vif(fit4)

fita <- lm(M1revlog~ + VA + PS + GE+ + RQ+ RL + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)

fitb <- lm(M1revlog~ + VA + PS + RQ + RL + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
car::vif(fitb)

fitc <- lm(M1revlog~ + VA + PS + GE + RL + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
car::vif(fitc)

fitd <- lm(M1revlog~ + VA + PS + GE + RQ + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
car::vif(fitd)

fite <- lm(M1revlog~ + VA + PS + GE + RQ+ RL + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
car::vif(fite)

fitf <- lm(M1revlog~ + VA + PS + CC + GDPcapitalog +
             effectivenesslog + safetylog + importancelog + as.factor(GAVI),
           data = gcsvaccine)
car::vif(fitf)


100-exp(coef(fit1))
100-exp(coef(fit2))
(exp(-0.32)-1) * 100
(exp(-1.16)-1) * 100
(exp(-0.38)-1) * 100
(exp(-1.10)-1) * 100
M1revlog
#VAR model

library(glmnet)
library(covidcast)
library(readr)
library(dplyr)
library(ggplot2)
library(car)
library(vars)
library(sandwich)
library(urca)
# Collect data ####


rm(list = ls())
options(error = function() { 
  library(RPushbullet)
  pbPost("note", "Error", geterrmessage())
  if(!interactive()) stop(geterrmessage())
})
#collect data
anxiety_1 <- covidcast_signal("fb-survey", "smoothed_wanxious_5d",
                              geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02")
depression_1 <- covidcast_signal("fb-survey", "smoothed_wdepressed_5d",
                                 geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02")

anxiety_2 <- covidcast_signal("fb-survey", "smoothed_wanxious_7d",
                              geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10")
depression_2 <- covidcast_signal("fb-survey", "smoothed_wdepressed_7d",
                                 geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10")


#Filter for specific state
adep_var_1 <- anxiety_1 %>% filter(geo_value == 'co')
ddep_var_1 <- depression_1 %>% filter(geo_value == 'co')
adep_var_2 <- anxiety_2 %>% filter(geo_value == 'co')
ddep_var_2 <- depression_2 %>% filter(geo_value == 'co')

###### Lockdown IV for first time frame ###

IV_1_1 <- covidcast_signal("fb-survey", "smoothed_wspent_time_1d",
                           geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02") %>% 
  filter(geo_value == 'co')
int <- intersect(adep_var_1$time_value, IV_1_1$time_value)
adep_var_1 <- subset(adep_var_1, time_value %in% int)

IV_1_2 <- covidcast_signal("fb-survey", "smoothed_wwork_outside_home_1d",
                           geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_1$time_value, IV_1_2$time_value)
adep_var_1 <- subset(adep_var_1, time_value %in% int)

IV_1_3 <- covidcast_signal("jhu-csse", "confirmed_7dav_incidence_prop",
                           geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_1$time_value, IV_1_3$time_value)
adep_var_1 <- subset(adep_var_1, time_value %in% int)

IV_1_4 <- covidcast_signal("jhu-csse", "deaths_7dav_incidence_prop",
                           geo_type="state", start_day = "2020-09-08", end_day = "2021-03-02") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_1$time_value, IV_1_4$time_value)
adep_var_1 <- subset(adep_var_1, time_value %in% int)

#compare with depression
int<- intersect(adep_var_1$time_value, ddep_var_1$time_value)
adep_var_1 <- subset(adep_var_1, time_value %in% int)

#subset once again to ensure homogeneity 
int <- intersect(adep_var_1$time_value, IV_1_1$time_value)
ddep_var_1 <- subset(ddep_var_1, time_value %in% int)
IV_1_1 <- subset(IV_1_1, time_value %in% int)
IV_1_2 <- subset(IV_1_2, time_value %in% int)
IV_1_3 <- subset(IV_1_3, time_value %in% int)
IV_1_4 <- subset(IV_1_4, time_value %in% int)

#remove duplicate dates
adep_var_1 <- adep_var_1[!duplicated(adep_var_1$time_value),]
ddep_var_1 <- ddep_var_1[!duplicated(ddep_var_1$time_value),]
IV_1_1 <- IV_1_1[!duplicated(IV_1_1$time_value),]
IV_1_2 <- IV_1_2[!duplicated(IV_1_2$time_value),]
IV_1_3 <- IV_1_3[!duplicated(IV_1_3$time_value),]
IV_1_4 <- IV_1_4[!duplicated(IV_1_4$time_value),]


#create variables for analysis
spentTime <- IV_1_1 %>% dplyr::select(value) %>% rename(spentTime = value)
workOutsideHome <- IV_1_2 %>% dplyr::select(value) %>% rename(workOutsideHome = value)
incidenceProp <- IV_1_3 %>% dplyr::select(value) %>% rename(incidenceProp = value)
deathCumProp <- IV_1_4 %>% dplyr::select(value) %>% rename(deathCumProp = value)
anxiety1 <- adep_var_1 %>% dplyr::select(value) %>% rename(anxiety = value)
depression1 <- ddep_var_1 %>% dplyr::select(value) %>% rename(depression = value)
dates1 <- IV_1_1 %>% dplyr::select(time_value) %>% rename(dates1 = time_value)

#first order differencing to ensure stationary time series
spentTime <- data.matrix(spentTime)
workOutsideHome <- data.matrix(workOutsideHome)
incidenceProp <- data.matrix(incidenceProp)
deathCumProp <- data.matrix(deathCumProp)
anxiety1 <- data.matrix(anxiety1)
depression1 <- data.matrix(depression1)


#test stationarity of data using auto correlation function,
#most of the values should lie within the blue dashed lines
# acf(spentTime)
# acf(workOutsideHome)
# acf(incidenceProp)
# acf(deathCumProp)
# acf(anxiety1)
# acf(depression1)

###### Lockdown IV for second time frame #
IV_2_1 <- covidcast_signal("fb-survey", "smoothed_wcovid_vaccinated",
                           geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_2$time_value, IV_2_1$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

IV_2_2 <- covidcast_signal("fb-survey", "smoothed_wspent_time_indoors_1d",
                           geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_2$time_value, IV_2_2$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

IV_2_3 <- covidcast_signal("fb-survey", "smoothed_wwork_outside_home_indoors_1d",
                           geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_2$time_value, IV_2_3$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

IV_2_4 <- covidcast_signal("jhu-csse", "confirmed_7dav_incidence_prop",
                           geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_2$time_value, IV_2_4$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

IV_2_5 <- covidcast_signal("jhu-csse", "deaths_7dav_incidence_prop",
                           geo_type="state", start_day = "2021-03-02", end_day = "2022-01-10") %>%
  filter(geo_value == 'co')
int <- intersect(adep_var_2$time_value, IV_2_5$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

#calc intersect for depression and anxiety
int <- intersect(adep_var_2$time_value, ddep_var_2$time_value)
adep_var_2 <- subset(adep_var_2, time_value %in% int)

#subset everything once again to ensure homogeneity 
int <- intersect(adep_var_2$time_value, IV_2_5$time_value)
ddep_var_2 <- subset(ddep_var_2, time_value %in% int)
IV_2_1 <- subset(IV_2_1, time_value %in% int)
IV_2_2 <- subset(IV_2_2, time_value %in% int)
IV_2_3 <- subset(IV_2_3, time_value %in% int)
IV_2_4 <- subset(IV_2_4, time_value %in% int)
IV_2_5 <- subset(IV_2_5, time_value %in% int)

#remove duplicate rows
adep_var_2 <- adep_var_2[!duplicated(adep_var_2$time_value),]
ddep_var_2 <- ddep_var_2[!duplicated(ddep_var_2$time_value),]
IV_2_1 <- IV_2_1[!duplicated(IV_2_1$time_value),]
IV_2_2 <- IV_2_2[!duplicated(IV_2_2$time_value),]
IV_2_3 <- IV_2_3[!duplicated(IV_2_3$time_value),]
IV_2_4 <- IV_2_4[!duplicated(IV_2_4$time_value),]
IV_2_5 <- IV_2_5[!duplicated(IV_2_5$time_value),]


covidVaccinated <- IV_2_1 %>% dplyr::select(value) %>% rename(covidVaccinated = value)
spentTimeIndoors <- IV_2_2 %>% dplyr::select(value) %>% rename(spentTimeIndoors = value)
workOutsideHomeIndoors <- IV_2_3 %>% dplyr::select(value) %>% rename(workOutsideHomeIndoors = value)
incidenceProp_2 <- IV_2_4 %>% dplyr::select(value) %>% rename(incidenceProp_2 = value)
deathCumProp_2 <- IV_2_5 %>% dplyr::select(value) %>% rename(deathCumProp_2 = value)
anxiety2 <- adep_var_2 %>% dplyr::select(value) %>% rename(anxiety = value)
depression2 <- ddep_var_2 %>% dplyr::select(value) %>% rename(depression = value)
dates2 <- IV_2_1 %>% dplyr::select(time_value) %>% rename(dates2 = time_value)

covidVaccinated <- data.matrix(covidVaccinated)
spentTimeIndoors <- data.matrix(spentTimeIndoors)
workOutsideHomeIndoors <- data.matrix(workOutsideHomeIndoors)
incidenceProp_2 <- data.matrix(incidenceProp_2)
deathCumProp_2 <- data.matrix(deathCumProp_2)
anxiety2 <- data.matrix(anxiety2)
depression2 <- data.matrix(depression2)

#test stationarity of data using auto correlation function,
#most of the values should lie within the blue dashed lines
# acf(covidVaccinated)
# acf(spentTimeIndoors)
# acf(workOutsideHomeIndoors)
# acf(incidenceProp_2)
# acf(deathCumProp_2)
# acf(anxiety2)
# acf(depression2)

######  Define IVs ###

anxiety1 <- data.frame(anxiety1)
depression1 <- data.frame(depression1)
anxiety2 <- data.frame(anxiety2)
depression2 <- data.frame(depression2)

dep_var1_death <- data.frame(anxiety1, depression1, spentTime, workOutsideHome, 
                       incidenceProp, deathCumProp)
dep_var2_death <- data.frame(anxiety2, depression2, covidVaccinated, spentTimeIndoors,
                       workOutsideHomeIndoors, incidenceProp_2, deathCumProp_2)
dep_var1_incidence <- data.frame(anxiety1, depression1, spentTime, workOutsideHome, 
                                 deathCumProp, incidenceProp)
dep_var2_incidence <- data.frame(anxiety2, depression2, covidVaccinated, spentTimeIndoors,
                             workOutsideHomeIndoors, deathCumProp_2, incidenceProp_2)
dep_var1_work <- data.frame(anxiety1, depression1, spentTime,
                            deathCumProp, incidenceProp, workOutsideHome)
dep_var2_work <- data.frame(anxiety2, depression2, covidVaccinated, spentTimeIndoors,
                            deathCumProp_2, incidenceProp_2, workOutsideHomeIndoors)
dep_var1_spent <- data.frame(anxiety1, depression1, workOutsideHome, 
                             deathCumProp, incidenceProp, spentTime)
dep_var2_spent <- data.frame(anxiety2, depression2, covidVaccinated,
                             workOutsideHomeIndoors, deathCumProp_2, incidenceProp_2, spentTimeIndoors)
dep_var2_vac <- data.frame(anxiety2, depression2,
                           workOutsideHomeIndoors, deathCumProp_2, incidenceProp_2, spentTimeIndoors, covidVaccinated)

#FOR SAVING THE DATA
tp1data <- data.frame(dates1, anxiety1, depression1, spentTime, workOutsideHome, 
                      incidenceProp, deathCumProp) 
tp2data <- data.frame(dates2, anxiety2, depression2, covidVaccinated, spentTimeIndoors,
                      workOutsideHomeIndoors, incidenceProp_2, deathCumProp_2)

write.csv(tp1data, "C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/data/Colorado1.csv",row.names = FALSE)
write.csv(tp2data, "C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/data/Colorado2.csv",row.names = FALSE)
###### Run models ##
####FIRST TIME PERIOD###
#save J test results
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/Johansen tests")
sink('ColoradoJtest1.txt')
summary(ca.jo(dep_var1_death, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun'))
sink()
#Specify the model and its order
vecmodel1_death <- ca.jo(dep_var1_death, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl1_death = vars::vec2var(vecmodel1_death, r = 1)
vecmodel1_incidence <- ca.jo(dep_var1_incidence, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl1_incidence = vars::vec2var(vecmodel1_incidence, r = 1)
vecmodel1_work <- ca.jo(dep_var1_work, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl1_work = vars::vec2var(vecmodel1_work, r = 1)
vecmodel1_spent <- ca.jo(dep_var1_spent, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl1_spent = vars::vec2var(vecmodel1_spent, r = 1)
#collect irf results for plotting later, if you get a singularity error, add 'boot = FALSE'
irfresult1_death = irf(mdl1_death, n.ahead = 20, impulse = 'deathCumProp')
irfresult1_incidence = irf(mdl1_incidence, n.ahead = 20, impulse = 'incidenceProp')
irfresult1_work = irf(mdl1_work, n.ahead = 20, impulse = 'workOutsideHome')
irfresult1_spent = irf(mdl1_spent, n.ahead = 20, impulse = 'spentTime')
#irfresult1 = irf(mdl1, boot = FALSE)
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/IRF coeffs")
sink('irfColorado1_death.txt')
print(irfresult1_death)
sink()
#
sink('irfColorado1_incidence.txt')
print(irfresult1_incidence)
sink()
#
sink('irfColorado1_work.txt')
print(irfresult1_work)
sink()
#
sink('irfColorado1_spent.txt')
print(irfresult1_spent)
sink()
#save the irf plots. You'll need to hit enter in the console a few times to get all of the plots to save
#after running line: plot(irfresult1)
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/IRF")
pdf('irfColorado1_death.pdf')
plot(irfresult1_death)
dev.off()
# tiff('irfColorado1_death.tiff', res = 300)
# plot(irfresult1_death)
# dev.off()
#
pdf('irfColorado1_incidence.pdf')
plot(irfresult1_incidence)
dev.off()
# tiff('irfColorado1_incidence.tiff', res = 300)
# plot(irfresult1_incidence)
# dev.off()
#
pdf('irfColorado1_work.pdf')
plot(irfresult1_work)
dev.off()
# tiff('irfColorado1_work.tiff', res = 300)
# plot(irfresult1_work)
# dev.off()
#
pdf('irfColorado1_spent.pdf')
plot(irfresult1_spent)
dev.off()
# tiff('irfColorado1_spent.tiff', res = 300)
# plot(irfresult1_spent)
# dev.off()
#save model coefficients
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/Coefficients")
sink('CoeffColorado1.txt')
print(c("Number of days",length(data.matrix(anxiety1)),"VECM EQN ESTIMATES",summary(cajools(vecmodel1_death)),"VAR EQN ESTIMATES",mdl1_death$A))
sink()

####SECOND TIME PERIOD###
#save J test results
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/Johansen tests")
sink('ColoradoJtest2.txt')
summary(ca.jo(dep_var2_death, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun'))
sink()
#Specify the model and its order
vecmodel2_death <- ca.jo(dep_var2_death, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl2_death = vars::vec2var(vecmodel2_death, r = 1)
vecmodel2_incidence <- ca.jo(dep_var2_incidence, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl2_incidence = vars::vec2var(vecmodel2_incidence, r = 1)
vecmodel2_work <- ca.jo(dep_var2_work, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl2_work = vars::vec2var(vecmodel2_work, r = 1)
vecmodel2_spent <- ca.jo(dep_var2_spent, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl2_spent = vars::vec2var(vecmodel2_spent, r = 1)
vecmodel2_vac <- ca.jo(dep_var2_vac, type = 'eigen', ecdet = 'const', K = 2, spec = 'longrun')
mdl2_vac = vars::vec2var(vecmodel2_vac, r = 1)

#collect irf results, if you get a singularity error, add 'boot = FALSE'
irfresult2_death = irf(mdl2_death, n.ahead = 20, impulse = 'deathCumProp_2')
irfresult2_incidence = irf(mdl2_incidence, n.ahead = 20, impulse = 'incidenceProp_2')
irfresult2_work = irf(mdl2_work, n.ahead = 20, impulse = 'workOutsideHomeIndoors')
irfresult2_spent = irf(mdl2_spent, n.ahead = 20, impulse = 'spentTimeIndoors')
irfresult2_vac = irf(mdl2_vac, n.ahead = 20, impulse = 'covidVaccinated')
#irfresult1 = irf(mdl1, boot = FALSE)
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/IRF coeffs")
sink('irfColorado2_death.txt')
print(irfresult2_death)
sink()
#
sink('irfColorado2_incidence.txt')
print(irfresult2_incidence)
sink()
#
sink('irfColorado2_work.txt')
print(irfresult2_work)
sink()
#
sink('irfColorado2_spent.txt')
print(irfresult2_spent)
sink()
#
sink('irfColorado2_vac.txt')
print(irfresult2_vac)
sink()

#save the irf plots. You'll need to hit enter in the console a few times to get all of the plots to save
#after running line: plot(irfresult1)
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/IRF")
pdf('irfColorado2_death.pdf')
plot(irfresult2_death)
dev.off()
# tiff('irfColorado2_death.tiff', res = 300)
# plot(irfresult2_death)
# dev.off()
#
pdf('irfColorado2_incidence.pdf')
plot(irfresult2_incidence)
dev.off()
# tiff('irfColorado2_incidence.tiff', res = 300)
# plot(irfresult2_incidence)
# dev.off()
#
pdf('irfColorado2_work.pdf')
plot(irfresult2_work)
dev.off()
# tiff('irfColorado2_work.tiff', res = 300)
# plot(irfresult2_work)
# dev.off()
#
pdf('irfColorado2_spent.pdf')
plot(irfresult2_spent)
dev.off()
# tiff('irfColorado2_spent.tiff', res = 300)
# plot(irfresult2_spent)
# dev.off()
#
pdf('irfColorado2_vac.pdf')
plot(irfresult2_vac)
dev.off()
# tiff('irfColorado2_vac.tiff', res = 300)
# plot(irfresult2_vac)
# dev.off()

#save model coefficients
setwd("C:/Users/a426f262/OneDrive - University of Kansas/Desktop/School/Research/COVID-19 and Mental Health/code + results/results/Coefficients")
sink('CoeffColorado2.txt')
print(c("Number of days",length(data.matrix(anxiety2)),"VECM EQN ESTIMATES",summary(cajools(vecmodel2_death)),"VAR EQN ESTIMATES",mdl2_death$A))
sink()
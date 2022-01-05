options(digits = 22)
library(lmerTest)
library(merTools)
library(nortest)
library(MuMIn)
library(performance)
library(car)
library(plyr)

warnings()

# R statistical analysis
#The majority of this code comes from the r_calculation_10.R with #'s referencing models in that file
#For simplicity this code shows the derivation of;
# 1) A linear model of the PWGD
# 2) the two reference constant models J21 and J23.
# 3) model J22 with ionic strength as an additional fixed effect


all_FULL_HEADER<- read.csv("/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/bin/all_FULL_HEADER.txt", header = TRUE,sep="\t") #all_FULL_HEADER_300_AGE

#linear model of PWGD
m_J1 <- lm(Log_Ca_Mg ~ 1 + kelvin + I(kelvin^-1),data=all_FULL_HEADER)  #m_lm

m_J21 <- lmer(Log_Ca_Mg ~ 1 + kelvin + I(kelvin^-1) + (1|FIELD) + (1|DepthID) + (1|FORMATION) + (1|BASIN) + (1|LITHOLOGY) + (1|AGE_VAL) + (1|AGE_VAL_2) ,data=all_FULL_HEADER)  #m_6_AGE_3
m_J23 <- lmer(dsp_val ~ 1 + kelvin + I(kelvin^-1) + (1|FIELD) + (1|DepthID) + (1|FORMATION) + (1|BASIN) + (1|LITHOLOGY) + (1|AGE_VAL) + (1|AGE_VAL_2),data=all_FULL_HEADER) #m_6_dsp

m_J22 <- lmer(Log_Ca_Mg ~ 1 + kelvin + I(kelvin^-1) + mu + (1|FIELD) + (1|DepthID) + (1|FORMATION) + (1|BASIN) + (1|LITHOLOGY) + (1|AGE_VAL) + (1|AGE_VAL_2) ,data=all_FULL_HEADER) #m_6_mu


compare_AIC <- AIC(m_J1,m_J21,m_J23,m_J22)
print(compare_AIC)

###### MODEL J1 OUTPUT #####

deg_25 <- data.frame(kelvin = 298.15)
deg_200 <- data.frame(kelvin = 473.15)

sink("/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J1.txt")
print("m_J1")
print(icc(m_J1))
print(summary(m_J1))
print(r.squaredGLMM(m_J1))
predict_25 <- predict(m_J1,newdata=deg_25,interval='confidence')
predict_200 <- predict(m_J1,newdata=deg_200,interval='confidence')
print('25 degc')
print(predict_25)
val_a <- (2 * -8.48) - predict_25[1]
print(val_a)
val_upr <- (2*-8.48) - predict_25[3]
val_lwr <- (2*-8.48) - predict_25[2]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
print('200 degc')
print(predict_200)
val_b <- (2 * -11.049) - predict_200[1]
print(val_b)
val_upr <- (2 * -11.049) - predict_200[3]
val_lwr <- (2 * -11.049) - predict_200[2]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
sink()


###### MODEL J21 OUTPUT #####

deg_25 <- data.frame(kelvin = 298.15,DepthID=100,BASIN=100,AGE_VAL=100,AGE_VAL_2=100,LITHOLOGY=100,FORMATION=100,FIELD=100)
deg_200 <- data.frame(kelvin = 473.15,DepthID=100,BASIN=100,AGE_VAL=100,AGE_VAL_2=100,LITHOLOGY=100,FORMATION=100,FIELD=100)


sink("/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J21.txt")
print("m_J21")
print(icc(m_J21))
print(summary(m_J21))
print(r.squaredGLMM(m_J21))
predict_25 <- predict(m_J21,newdata=deg_25,re.form=NA)
predict_200 <- predict(m_J21,newdata=deg_200,re.form=NA)
predict_25_conf <- predictInterval(m_J21,newdata=deg_25, n.sims=20000,level=0.95, which="fixed")
predict_200_conf <- predictInterval(m_J21,newdata=deg_200, n.sims=20000,level=0.95, which="fixed")
predict_25_conf_68 <- predictInterval(m_J21,newdata=deg_25, n.sims=20000,level=0.68, which="fixed")
predict_200_conf_68 <- predictInterval(m_J21,newdata=deg_200, n.sims=20000,level=0.68, which="fixed")
print('25 degc')
val_a <- (2 * -8.480) - predict_25
print(val_a)
val_upr <- (2*-8.480) - predict_25_conf[2]
val_lwr <- (2*-8.480) - predict_25_conf[3]
val_upr_68 <- (2*-8.480) - predict_25_conf_68[2]
val_lwr_68 <- (2*-8.480) - predict_25_conf_68[3]
dif <- val_upr[1] - val_lwr[1]
dif_68 <- val_upr_68[1] - val_lwr_68[1]
print(val_upr)
print(val_lwr)
print(dif)
print(val_upr_68)
print(val_lwr_68)
print(dif_68)
print('200 degc')
val_b <- (2 * -11.049) - predict_200
print(val_b)
val_upr <- (2 * -11.049) - predict_200_conf[2]
val_lwr <- (2 * -11.049) - predict_200_conf[3]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
print(coef(m_J21)$AGE_VAL)
print(coef(m_J21)$AGE_VAL_2)
print(coef(m_J21)$LITHOLOGY)
print(coef(m_J21)$BASIN)
print(coef(m_J21)$FIELD)
print(coef(m_J21)$FORMATION)
print(coef(m_J21)$DepthID)
sink()
temp_ext <- data.frame(kelvin = seq(273.15, 573.15,by=0.1),DepthID=100,BASIN=100,AGE_VAL_2=100,AGE_VAL=100,LITHOLOGY=100,FORMATION=100,FIELD=100)
m_J21_pred <- predict(m_J21,temp_ext,allow.new.levels = TRUE,  re.form=NA)  #m_6_ext_pred
m_J21_pred <- cbind(m_J21_pred,temp_ext)
write.csv(m_J21_pred,"/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J21_pred.txt") #m_6_ext_pred_AGE_3

m_J21_pred_Interval <- predictInterval(m_J21,newdata=temp_ext, n.sims=20000,level=0.95, which="fixed") #predictInterval_m_6_ext_pred
m_6_ext <- cbind(m_J21_pred_Interval,temp_ext)
write.csv(m_6_ext,"/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J21_pred_Interval.txt")


###### MODEL J23 OUTPUT #####

deg_25 <- data.frame(kelvin = 298.15,DepthID=100,BASIN=100,LITHOLOGY=100,FORMATION=100,FIELD=100,AGE_VAL_2=100,AGE_VAL=100)
deg_200 <- data.frame(kelvin = 473.15,DepthID=100,BASIN=100,LITHOLOGY=100,FORMATION=100,FIELD=100,AGE_VAL_2=100,AGE_VAL=100)

sink("/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J23.txt")
print("m_J23")
print(icc(m_J23))
print(summary(m_J23))
print(r.squaredGLMM(m_J23))
predict_25 <- predict(m_J23,newdata=deg_25,re.form=NA)
predict_200 <- predict(m_J23,newdata=deg_200,re.form=NA)
predict_25_conf <- predictInterval(m_J23,newdata=deg_25, n.sims=20000,level=0.95, which="fixed")
predict_200_conf <- predictInterval(m_J23,newdata=deg_200, n.sims=20000,level=0.95, which="fixed")
print('25 degc')
val_a <- predict_25
print(val_a)
val_upr <- predict_25_conf[2]
val_lwr <- predict_25_conf[3]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
print('200 degc')
val_b <- predict_200
print(val_b)
val_upr <- predict_200_conf[2]
val_lwr <- predict_200_conf[3]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
print(coef(m_J23)$LITHOLOGY)
print(coef(m_J23)$AGE)
print(coef(m_J23)$BASIN)
print(coef(m_J23)$FORMATION)
print(coef(m_J23)$FIELD)
print(coef(m_J23)$AGE_VAL)
print(coef(m_J23)$AGE_VAL_2)
print(coef(m_J23)$DepthID)
sink()
temp_ext <- data.frame(kelvin = seq(273.15, 573.15,by=1),DepthID=100,BASIN=100,AGE_VAL=100,AGE_VAL_2=100,LITHOLOGY=100,FORMATION=100,FIELD=100)
m_J23_dsp_pred <- predict(m_J23,temp_ext,allow.new.levels = TRUE,  re.form=NA)
m_J23_dsp_pred <- cbind(m_J23_dsp_pred,temp_ext)
write.csv(m_J23_dsp_pred,"/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J23_dsp_pred.txt")

###### Model J22 OUTPUT ######

deg_25 <- data.frame(kelvin = 298.15,mu=0.0,DepthID=100,BASIN=100,PERIOD=100,LITHOLOGY=100,FORMATION=100,FIELD=100,AGE_VAL=100,AGE_VAL_2=100)
deg_200 <- data.frame(kelvin = 473.15,mu=0.0,DepthID=100,BASIN=100,PERIOD=100,LITHOLOGY=100,FORMATION=100,FIELD=100,AGE_VAL=100,AGE_VAL_2=100)

sink("/Users/hamish/Documents/AWSprojects/data/AJS_Robertson_2022/R_CODE/model_output/m_J22.txt")
print("m_J22")
print(icc(m_J22))
print(summary(m_J22))
print(r.squaredGLMM(m_J22))
predict_25 <- predict(m_J22,newdata=deg_25,re.form=NA)
predict_200 <- predict(m_J22,newdata=deg_200,re.form=NA)
predict_25_conf <- predictInterval(m_J22,newdata=deg_25, n.sims=20000,level=0.95, which="fixed")
predict_200_conf <- predictInterval(m_J22,newdata=deg_200, n.sims=20000,level=0.95, which="fixed")
print('25 degc')
val_a <- (2 * -8.48) - predict_25
print(val_a)
val_upr <- (2*-8.48) - predict_25_conf[2]
val_lwr <- (2*-8.48) - predict_25_conf[3]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
print('200 degc')
val_b <- (2 * -11.049) - predict_200
print(val_b)
val_upr <- (2 * -11.049) - predict_200_conf[2]
val_lwr <- (2 * -11.049) - predict_200_conf[3]
dif <- val_upr[1] - val_lwr[1]
print(val_upr)
print(val_lwr)
print(dif)
sink()

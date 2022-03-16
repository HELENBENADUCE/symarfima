source("symarfimamodel.R")

# reading the data
H=read.csv(file="Osorio_2020.csv")
C=read.csv(file="SP_2020.csv")
# formating the data 
D=C[which(C$Minute==0),] 
a=which(D$Month==10 & D$Day==20& D$Hour==0) # october 20th, 0 hours 
b=which( D$Month==12 & D$Day==3& D$Hour==5) # december 3rd, 5AM
ot = H[a:b,]
st = D[a:b,]
# pt contains the desired time series
pt=ot$Wind.Speed-st$Wind.Speed
length(pt)
# Unit roots tests: p-values smaller than 0.01
PP.test(pt)
tseries::adf.test(pt)

# final model 
gh=symarfima.fit(pt,p=1, q=48,  index1 = 10, fit.alpha = T, cons=T, 
           max.fun=10000, fixed.ma=c(NA,0,0,0,rep(0,19),NA,rep(0,23),NA),h=24)
# 
# phi1     theta1    theta24    theta48          d      alpha 
# 0.77844692 0.44167273 0.10983944 0.10042120 0.42249026 0.40452866 
# varphi 
# 0.07523346 
# [1] "------"
# Estimate Std. Error z value  Pr(>|z|)    
#   phi1     0.77845    0.03029 25.7009 < 2.2e-16 ***
#   theta1   0.44167    0.02860 15.4449 < 2.2e-16 ***
#   theta24  0.10984    0.02423  4.5330 5.814e-06 ***
#   theta48  0.10042    0.02269  4.4254 9.627e-06 ***
#   d        0.42249    0.03859 10.9490 < 2.2e-16 ***
#   alpha    0.40453    0.13435  3.0111  0.002603 ** 
#   varphi   0.07523    0.00389 19.3515 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# AIC: 582.7444            
# BIC: 617.5198            
# Log-likelihood: -284.3722


# Martingale difference tests
library(vrtest)
set.seed(123)
AutoBoot.test(gh$residual, 500, wild="Mammen",prob=c(0.025,0.975))
# pval = 0.546
set.seed(123)
AutoBoot.test(gh$residual, 500, wild="Rademacher",prob=c(0.025,0.975))
# pval = 0.632
set.seed(123)
AutoBoot.test(gh$residual, 500, wild="Normal",prob=c(0.025,0.975))
# pval = 0.506
set.seed(123)
DL.test(gh$residual,500,p=50)
# Cp_pval = 0.52;  Kp_pval = 0.274
set.seed(123)
Gen.Spec.Test(gh$residual,500)
# pval = 0.434

# Ljung-Box test
Box.test(gh$residual, lag=50, fitdf=5,type="L")
# pval = 0.01293

#For forecasting purposes
af = which( D$Month==12 & D$Day==3& D$Hour==6) # december 03, 6AM
bf = which( D$Month==12 & D$Day==4& D$Hour==5) # december 04, 5AM
# contains the observed data to be compared to the forecasts:
yy = H$Wind.Speed[af:bf] - D$Wind.Speed[af:bf] 
## simple plots:
# plot.ts(yy, type="b", pch=16, ylim=c(min(yy,gh$prev),max(yy,gh$prev))) 
# lines(1:24,gh$prev, type="b", col="red")


 
library(forecast)
library(ggplot2)
library(dplyr)

# function HoltWinters requires a ts object - start is dummy
hws = HoltWinters(ts(pt, start=c(2020,10),frequency=24)) 
forecasthws<-as.vector(forecast(hws,h =  24)$mean)


# out-of-sample forecast - accuracy measures
mse1=mse2=numeric(24)
mae1=mae2=numeric(24)
for(i in 1:24){
h1=i
mse1[i] = mean((gh$prev[1:h1]-yy[1:h1])^2)
mse2[i] = mean((forecasthws[1:h1]-yy[1:h1])^2)
mae1[i] = mean(abs((gh$prev[1:h1]-yy[1:h1])/yy[1:h1]))
mae2[i] = mean(abs((forecasthws[1:h1]-yy[1:h1])/yy[1:h1]))
}
cbind(mse1,mse2,mae1,mae2)

# in sample forecast - accuracy measures.
mean((gh$mut_hat-pt[1:1062])^2) #0.118
mean((hw$fitted[1:1038]-pt[25:1062])^2) # 0.262




#crude plot
plot.ts(yy, type="b", pch=16, ylim=c(min(yy,gh$prev,forecasthws),max(yy,gh$prev,forecasthws))) 
lines(gh$prev, col="blue", pch=16, type="b")
lines(forecasthws, col="red", pch=16, type="b")
 
############ PLOTS ###############
#setting up for a fancy time series plot

J = cbind(1:1062,H[a:b,c(1,2,3,4,23)],D[a:b,23],H[a:b,23]-D[a:b,23])
colnames(J)[c(1,6,7,8)]=c("Index", "Osorio","SP","diff")
rownames(J)=NULL


p <- ggplot(J, aes(x=Index,y=diff)) +
   geom_line() + 
   #  geom_line(aes(y=Osorio), colour="blue")+
   #  geom_line(aes(y=SP), colour="red")+
   xlab("Time Index")+ylab("Wind Speed (m/s)")+theme_bw()+
   theme(axis.text.x=element_text(size=rel(2.2)),axis.text.y=element_text(size=rel(2.2)), 
         axis.title.y=element_text(size=rel(2),vjust=2),
         axis.title.x=element_text(size=rel(2),vjust=-1),
         plot.title=element_text(size=rel(2.3),hjust=0.5, vjust=2),
         axis.ticks.x.bottom=element_line()
   ) +
   ggtitle("Difference between the wind speed in Osório/RS and São Paulo/SP")+
   scale_x_continuous(breaks=c(0,265,530,795,1062))
p


wf<-15
hf<-6
dev.off()
pdf(file="F:/Dropbox/alunos/Helen/Helen/Paper/JSPI/figuras/TS_plot.pdf", width=wf,height=hf)
p
dev.off()

## ACF and PACF plots

w1<-7
h1<-5

pdf(file = "F:/Dropbox/alunos/Helen/Helen/Paper/JSPI/figuras/acf_pacf.pdf",width = w1, height = h1,family = "Times")
{
   par(mfrow=c(2,1), mai= c(0.7,1,0.5,1))
   
   acf(pt, lag=100, main="", xlab="")
   pacf(pt, lag=100, main="", xlab="")
}
dev.off()

##########################
########## Forecast plots


# #crude plot
# plot.ts(yy, type="b", pch=16, ylim=c(min(yy,gh$prev,forecasthws),max(yy,gh$prev,forecasthws))) 
# lines(gh$prev, col="blue", pch=16, type="b")
# lines(forecasthws, col="red", pch=16, type="b")


# Fancy forecast plots

#setting up the data 
obs1=c(pt[1012:1062],yy[1],rep(NA,23))
obs2=c(rep(NA, 51), yy[1:24])
insample = c(gh$mut_hat[1012:1062], rep(NA,24))
outsample = c(rep(NA,51), gh$prev[1:24])
ohw = c(rep(NA,51), forecasthws[1:24])
ind=1012:1086
bogus1=c(rep(NA,50),pt[1062],gh$prev[1],rep(NA,23))
bogus2=c(rep(NA,50),pt[1062],forecasthws[1],rep(NA,23))
M=data.frame(cbind(ind,obs1,obs2,insample, outsample, ohw, bogus1, bogus2))

names(M)=c("Index","obs.IS","obs.OS","prev.IS", "prev.OS", "prev.hw", "bogus1","bogus2")

# actual plot
tam=0.8
s <- ggplot(M, aes(x=Index,y=obs.IS),na.rm=T) + theme_bw() +
   geom_line(na.rm=T, size=tam) + 
   geom_line(aes(y=obs.OS, color="Observed"),na.rm=T, size=tam)+ geom_point(aes(y=obs.OS), colour="black",na.rm=T, size=2) +
   geom_line(aes(y=prev.IS, color="In-sample"),  linetype="dashed",na.rm=T, size=tam) +
   geom_line(aes(y=prev.OS, color="SYMMARFIMA"), na.rm=T, size=tam) + geom_point(aes(y=prev.OS), colour="blue",na.rm=T, size=2) +
   geom_line(aes(y=prev.hw, color="Holt-Winters"), na.rm=T, size=tam) + geom_point(aes(y=prev.hw), colour="red",na.rm=T, size=2) +
   geom_line(aes(y=bogus1), color="blue",na.rm=T, size=tam) +
   geom_line(aes(y=bogus2), color="red",na.rm=T, size=tam) +
   geom_vline(xintercept=1062, col="gray", size=tam) +
   xlab("Time") + ylab("Wind Speed (m/s)") + 
   theme(axis.text.x=element_text(size=rel(2.2)),axis.text.y=element_text(size=rel(2.2)), 
         axis.title.y=element_text(size=rel(2),vjust=2),
         axis.title.x=element_text(size=rel(2),vjust=-1),
         plot.title=element_text(size=rel(2.3),hjust=0.5, vjust=2),
         legend.text = element_text(color="black", size=rel(1.4)),
         legend.position=c(0.1,0.8),
         axis.ticks.x.bottom=element_line()
   ) +
   ggtitle("In-sample and out-of-sample forecasts for the wind speed data") +  
   scale_x_continuous(breaks=c(1012,1037,1062,1086))+
   scale_colour_manual("",values = c("Observed"="black", "In-sample"="green", 
                                     "SYMMARFIMA"="blue", "Holt-Winters"="red"),limits=c("Observed", "In-sample", "SYMMARFIMA","Holt-Winters")) 

s  

w2<-15
h2<-5

pdf(file = "F:/Dropbox/alunos/Helen/Helen/Paper/JSPI/figuras/forecast.pdf",width = w2, height = h2)
{
s
}
dev.off()


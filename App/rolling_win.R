source("symarfimamodel.R")

H=read.csv(file="Osorio_2020.csv")
C=read.csv(file="SP_2020.csv")
D=C[which(C$Minute==0),]

# data preparation
a=which(D$Month==10 & D$Day==20& D$Hour==0)
b=which( D$Month==12 & D$Day==10& D$Hour==0) # december 10, 0AM
ot = H[a:b,]
st = D[a:b,]
y = ot$Wind.Speed-st$Wind.Speed
length(y)

# Rolling window exercise 
# gh=symarfima.fit(y[1:1062],p=5, q=48,  index1 = 10, fit.alpha = T, cons=T, 
#                  max.fun=10000, fixed.ar = c(NA,0,0,0,0), fixed.ma=c(NA,0,0,0,rep(0,19),NA,rep(0,23),NA),h=24)
# est_coefs = matrix(rep(0,7*60), ncol=7)
# pred = matrix(rep(0,24*60), ncol=24)
# colnames(pred)=1:24
# colnames(est_coefs)=c("phi1","theta1","theta24","theta48","d","alpha","varphi")
# #pred[1,]=gh$prev
# #est_coefs[1,]=gh$coef
# 
# for(i in 0:99){
# yt = y[1:(1062+i)] 
# gh=symarfima.fit(y[1:1062],p=5, q=48,  index1 = 10, fit.alpha = T, cons=T, 
#                  max.fun=10000, fixed.ar = c(NA,0,0,0,0), fixed.ma=c(NA,0,0,0,rep(0,19),NA,rep(0,23),NA),h=24) est_coefs[i+1,]=gh$coef
# pred[i+1,] = gh$prev
# cat(i,"\n", sep="")
# }
# write.table(est_coefs, sep=",",file = "est_coef_RW.csv",row.names = F, 
#             col.names = T)
# write.table(pred, sep=",",file = "pred_RW.csv",row.names = F, 
#             col.names = T)

# 

# coefs=read.csv(file = "est_coef_RW.csv", header=T)

pred=as.matrix(read.csv(file = "pred_RW.csv", header=T))
colnames(pred)=1:24
dim(pred)
dat = stack(as.data.frame(pred))


ov=0*pred
for(i in 1:24){
  ov[,i]= y[(1062+i):(1062+99+i)]
}


err=(pred-ov)^2
round(colMeans(err),3)
err2 = (pred-ov)


library(ggplot2)
dat = stack(as.data.frame(err2))
p <- ggplot(dat) + 
  geom_boxplot(aes(x=ind,y=values), fill="lightgrey")+
  stat_boxplot(geom ='errorbar', width = 0.6, aes(x=ind,y=values)) +
  theme_bw()+ geom_hline(yintercept=0, linetype="dashed", color = "red", size=1.2)+
  theme(axis.text.x=element_text(size=rel(2.2)),axis.text.y=element_text(size=rel(2.2)), 
        axis.title.y=element_text(size=rel(2),vjust=2),
        panel.grid = element_blank(),
        axis.title.x=element_text(size=rel(2),vjust=-1),
        plot.title=element_text(size=rel(2.3),hjust=0.5, vjust=2),
        axis.ticks.x.bottom=element_line()
  ) +
  xlab("Number of steps ahead forecast")+ylab("predicted - observed")+
  ggtitle("Difference between predicted and observed values - rolling window exercise")
p


wf<-15
hf<-6
dev.off()
pdf(file="F:/Dropbox/alunos/Helen/Helen/Paper/JSPI/figuras/roll_win.pdf", width=wf,height=hf)
p
dev.off()

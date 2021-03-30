####### Result combiner #####
setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/case1_2/FInal piece of codes")
load("mymodel_rcl_pred.RData")


uncert_p<-function(pr,truev,predval,predvar)
{
  li<-qnorm(p=(1-pr)/2,mean =predval,sd=sqrt(predvar) )
  ui<-qnorm(p=(1+pr)/2,mean =predval,sd=sqrt(predvar) )
  kj_p<-as.numeric(truev<ui & truev>li)
  k_p<-mean(kj_p)
  temp<-kj_p*(ui-li)
  w_p<-mean(temp)/k_p
  return(list(kjp=kj_p,kp=k_p,wp=w_p))
}

#### Goodness statistics ##
G_stat<-function(p,kp)
{
  ap<-as.numeric(kp>p)
  s1<-3*ap-2
  s2<-kp-p
  s<-s1*s2
  dp<-p[2]-p[1]
  return(1-dp*sum(s))
}


mypredres
myrmse_int<-numeric()
myrmse_forc<-numeric()
myG_pred<-numeric()
myG_forc<-numeric()
mypredvalue<-list()
mypredvariance<-list()
myforcval<-list()
myforcvariance<-list()
testpred<-list()
testforc<-list()
my_crps_pred<-my_crps_forc<-my_logs_pred<-my_logs_forc<-numeric()
for(i in (1:nsim))
{
  myrmse_int[i]<-mypredres[[i]]$RMSE_pred
  myrmse_forc[i]<-mypredres[[i]]$RMSE_forc
  myG_pred[i]<-mypredres[[i]]$myG_pred
  myG_forc[i]<-mypredres[[i]]$myG_forc
  mypredvalue[[i]]<-mypredres[[i]]$pred.val_pred
  mypredvariance[[i]]<-mypredres[[i]]$pred.var_pred
  myforcval[[i]]<-mypredres[[i]]$pred.val_forc
  myforcvariance[[i]]<-mypredres[[i]]$pred.var_forc
  testpred[[i]]<-mypredres[[i]]$testval_pred
  testforc[[i]]<-mypredres[[i]]$testval_forc
  my_crps_pred[i]<-mypredres[[i]]$mCRPS_pred
  my_crps_forc[i]<-mypredres[[i]]$mCRPS_forc
  my_logs_pred[i]<-mypredres[[i]]$mLogS_pred
  my_logs_forc[i]<-mypredres[[i]]$mLogS_forc
}

mypredvalue<-do.call(rbind,mypredvalue)
mypredvariance<-do.call(rbind,mypredvariance)
myforcval<-do.call(rbind,myforcval)
myforcvariance<-do.call(rbind,myforcvariance)
testpred<-do.call(rbind,testpred)
testforc<-do.call(rbind,testforc)


boxplot(myG_pred,myG_forc)

prob_seq<-seq(0.01,0.99,length.out = 4*99)

mywp_pred<-mykp_pred<-gmwp_pred<-gmkp_pred<-swp_pred<-skp_pred<-numeric(length = length(prob_seq))
mywp_forc<-mykp_forc<-gmwp_forc<-gmkp_forc<-swp_forc<-skp_forc<-numeric(length = length(prob_seq))

for(i in 1:length(prob_seq))
{
  temp1<-uncert_p(p=prob_seq[i],truev = c(testpred),predval = c(mypredvalue),predvar = c(mypredvariance))
  temp2<-uncert_p(p=prob_seq[i],truev = c(testforc),predval = c(myforcval),predvar = c(myforcvariance))
  mywp_pred[i]<-temp1$wp
  mywp_forc[i]<-temp2$wp
  mykp_pred[i]<-temp1$kp
  mykp_forc[i]<-temp2$kp
  

}


plot(prob_seq,mywp_forc,cex=0.5,pch=19)
points(prob_seq,mywp_pred,cex=0.5,pch=19,col=2)

plot(prob_seq,mykp_forc,cex=0.5,pch=19)
points(prob_seq,mykp_pred,cex=0.5,pch=19,col=2)
lines(prob_seq,prob_seq)


load("gneitingmodel_rcl_pred.RData")

gmrmse_int<-numeric()
gmrmse_forc<-numeric()
gmG_pred<-numeric()
gmG_forc<-numeric()
gmpredvalue<-list()
gmpredvariance<-list()
gmforcval<-list()
gmforcvariance<-list()
testpred<-list()
testforc<-list()
gm_crps_pred<-gm_crps_forc<-gm_logs_pred<-gm_logs_forc<-numeric()

for(i in 1:nsim)
{
  gmrmse_int[i]<-mypredres[[i]]$RMSE_pred
  gmrmse_forc[i]<-mypredres[[i]]$RMSE_forc
  gmG_pred[i]<-mypredres[[i]]$myG_pred
  gmG_forc[i]<-mypredres[[i]]$myG_forc
  gmpredvalue[[i]]<-mypredres[[i]]$pred.val_pred
  gmpredvariance[[i]]<-mypredres[[i]]$pred.var_pred
  gmforcval[[i]]<-mypredres[[i]]$pred.val_forc
  gmforcvariance[[i]]<-mypredres[[i]]$pred.var_forc
  testpred[[i]]<-mypredres[[i]]$testval_pred
  testforc[[i]]<-mypredres[[i]]$testval_forc
  gm_crps_pred[i]<-mypredres[[i]]$mCRPS_pred
  gm_crps_forc[i]<-mypredres[[i]]$mCRPS_forc
  gm_logs_pred[i]<-mypredres[[i]]$mLogS_pred
  gm_logs_forc[i]<-mypredres[[i]]$mLogS_forc
  
  
  
}

gmpredvalue<-do.call(rbind,gmpredvalue)
gmpredvariance<-do.call(rbind,gmpredvariance)
gmforcval<-do.call(rbind,gmforcval)
gmforcvariance<-do.call(rbind,gmforcvariance)
testpred<-do.call(rbind,testpred)
testforc<-do.call(rbind,testforc)


boxplot(gmG_pred,gmG_forc)

prob_seq<-seq(0.01,0.99,length.out = 4*99)


for(i in 1:length(prob_seq))
{
  temp1<-uncert_p(p=prob_seq[i],truev = c(testpred),predval = c(gmpredvalue),predvar = c(gmpredvariance))
  temp2<-uncert_p(p=prob_seq[i],truev = c(testforc),predval = c(gmforcval),predvar = c(gmforcvariance))
  gmwp_pred[i]<-temp1$wp
  gmwp_forc[i]<-temp2$wp
  gmkp_pred[i]<-temp1$kp
  gmkp_forc[i]<-temp2$kp
  
  
}


plot(prob_seq,gmwp_forc,cex=0.5,pch=19)
points(prob_seq,gmwp_pred,cex=0.5,pch=19,col=2)

plot(prob_seq,gmkp_forc,cex=0.5,pch=19)
points(prob_seq,gmkp_pred,cex=0.5,pch=19,col=2)
lines(prob_seq,prob_seq)



load("separablemodel_rcl_pred.RData")

sprmse_int<-numeric()
sprmse_forc<-numeric()
spG_pred<-numeric()
spG_forc<-numeric()
sppredvalue<-list()
sppredvariance<-list()
spforcval<-list()
spforcvariance<-list()
testpred<-list()
testforc<-list()
sp_crps_pred<-sp_crps_forc<-sp_logs_pred<-sp_logs_forc<-numeric()

for(i in 1:nsim)
{
  sprmse_int[i]<-mypredres[[i]]$RMSE_pred
  sprmse_forc[i]<-mypredres[[i]]$RMSE_forc
  spG_pred[i]<-mypredres[[i]]$myG_pred
  spG_forc[i]<-mypredres[[i]]$myG_forc
  sppredvalue[[i]]<-mypredres[[i]]$pred.val_pred
  sppredvariance[[i]]<-mypredres[[i]]$pred.var_pred
  spforcval[[i]]<-mypredres[[i]]$pred.val_forc
  spforcvariance[[i]]<-mypredres[[i]]$pred.var_forc
  testpred[[i]]<-mypredres[[i]]$testval_pred
  testforc[[i]]<-mypredres[[i]]$testval_forc
  sp_crps_pred[i]<-mypredres[[i]]$mCRPS_pred
  sp_crps_forc[i]<-mypredres[[i]]$mCRPS_forc
  sp_logs_pred[i]<-mypredres[[i]]$mLogS_pred
  sp_logs_forc[i]<-mypredres[[i]]$mLogS_forc
}

sppredvalue<-do.call(rbind,sppredvalue)
sppredvariance<-do.call(rbind,sppredvariance)
spforcval<-do.call(rbind,spforcval)
spforcvariance<-do.call(rbind,spforcvariance)
testpred<-do.call(rbind,testpred)
testforc<-do.call(rbind,testforc)


boxplot(spG_pred,spG_forc)

prob_seq<-seq(0.01,0.99,length.out = 4*99)


for(i in 1:length(prob_seq))
{
  temp1<-uncert_p(p=prob_seq[i],truev = c(testpred),predval = c(sppredvalue),predvar = c(sppredvariance))
  temp2<-uncert_p(p=prob_seq[i],truev = c(testforc),predval = c(spforcval),predvar = c(spforcvariance))
  swp_pred[i]<-temp1$wp
  swp_forc[i]<-temp2$wp
  skp_pred[i]<-temp1$kp
  skp_forc[i]<-temp2$kp
  
  
}


plot(prob_seq,swp_forc,cex=0.5,pch=19)
points(prob_seq,swp_pred,cex=0.5,pch=19,col=2)

plot(prob_seq,skp_forc,cex=0.5,pch=19)
points(prob_seq,skp_pred,cex=0.5,pch=19,col=2)
lines(prob_seq,prob_seq)


##### Comparing prediction G-stat, RMSE, wp, and kp.
par(mfrow=c(3,3))
boxplot(myG_pred,gmG_pred,spG_pred,myG_forc,gmG_forc,spG_forc, main="G-stat for prediction (blue) and forecasting (grey)",names = c("Tvar","Gneit-mat","Sep","Tvar","Gneit-mat","Sep"),col = c("blue","blue","blue","grey","grey","grey"))
points(x=c(1,2,3,4,5,6),y=c(mean(myG_pred),mean(gmG_pred),mean(spG_pred),mean(myG_forc),mean(gmG_forc),mean(spG_forc)),pch=19,col="green",cex=0.7)
legend("bottomright",pch=19,col=c("green"),"Mean value")


boxplot(myrmse_int,gmrmse_int,sprmse_int, main="RMSE prediction",names = c("Tvar","Gneit-mat","Sep"),col = "blue")
points(x=c(1,2,3),y=c(mean(myrmse_int),mean(gmrmse_int),mean(sprmse_int)),pch=19,col="green",cex=0.7)
legend("bottomright",pch=19,col=c("green"),"Mean value")

boxplot(myrmse_forc,gmrmse_forc,sprmse_forc, main="RMSE forecast",names = c("Tvar","Gneit-mat","Sep"),col = "grey")
points(x=c(1,2,3),y=c(mean(myrmse_forc),mean(gmrmse_forc),mean(sprmse_forc)),pch=19,col="green",cex=0.7)
legend("bottomright",pch=19,col=c("green"),"Mean value")


#cbind(mean(myrmse_int),mean(gmrmse_int),mean(sprmse_int))

plot(prob_seq,mykp_pred,type="l",col="blue",lwd=2,xlab="p-Probability Interval", ylab = "Proportion within Probability Interval",main="Prediction")
lines(prob_seq,gmkp_pred,lty=2,col="red",lwd=2)
lines(prob_seq,skp_pred,lty=3,col="green",lwd=2)
lines(prob_seq,prob_seq,lwd=2,lty=4)
#lines(prob_seq,mykp_pred,type="l",col="blue",lwd=2)

legend("bottomright",lwd=c(2,2,2,2),lty=c(1,2,3,4),col=c("blue","red","green","black"),c("Tvar","Gneit-mat","Sep","Diagonal line"))




plot(prob_seq,mywp_pred,type="l",col="blue",lwd=2,xlab="p-Probability Interval", ylab = "Probability Interval Width",main="Prediction",ylim = c(0,max(mywp_pred,gmwp_pred,swp_pred)))
lines(prob_seq,gmwp_pred,lty=2,col="red",lwd=2)
lines(prob_seq,swp_pred,lty=3,col="green",lwd=2)
#lines(prob_seq,prob_seq,lwd=2,lty=4)
#lines(prob_seq,mykp_pred,type="l",col="blue",lwd=2)

legend("bottomright",lwd=c(2,2,2),lty=c(1,2,3),col=c("blue","red","green"),c("Tvar","Gneit-mat","Sep"))



plot(prob_seq,mykp_forc,type="l",col="blue",lwd=2,xlab="p-Probability Interval", ylab = "Proportion within Probability Interval",main="Forecast")
lines(prob_seq,gmkp_forc,lty=2,col="red",lwd=2)
lines(prob_seq,skp_forc,lty=3,col="green",lwd=2)
lines(prob_seq,prob_seq,lwd=2,lty=4)
#lines(prob_seq,mykp_pred,type="l",col="blue",lwd=2)

legend("bottomright",lwd=c(2,2,2,2),lty=c(1,2,3,4),col=c("blue","red","green","black"),c("Tvar","Gneit-mat","Sep","Diagonal line"))




plot(prob_seq,mywp_forc,type="l",col="blue",lwd=2,xlab="p-Probability Interval", ylab = "Probability Interval Width",main="Forecast",ylim = c(0,max(mywp_forc,gmwp_forc,swp_forc)))
lines(prob_seq,gmwp_forc,lty=2,col="red",lwd=2)
lines(prob_seq,swp_forc,lty=3,col="green",lwd=2)
#lines(prob_seq,prob_seq,lwd=2,lty=4)
#lines(prob_seq,mykp_pred,type="l",col="blue",lwd=2)

legend("bottomright",lwd=c(2,2,2),lty=c(1,2,3),col=c("blue","red","green"),c("Tvar","Gneit-mat","Sep"))


#plot(prob_seq,mywp_pred,pch=19,col="blue",cex=0.5,xlab="p-Probability Interval", ylab = "Probability Interval Width",main="Forecast",ylim=c(0,max(mywp_pred,gmwp_pred,swp_pred)))
#points(prob_seq,gmwp_pred,pch=19,col="red",cex=0.5)
#points(prob_seq,swp_pred,col="green",cex=0.5)
#legend("topleft",pch=c(19,19,1),col=c("blue","red","green"),c("Our model","Gneiting's-Matern","Separable Model"))






##### Comparing forecasting G-stat, RMSE, wp, and kp.

boxplot(myG_forc,gmG_forc,spG_forc, main="G-stat forecasting",names = c("Tvar","Gneit-mat","Sep"))
cbind(mean(myG_forc),mean(gmG_forc),mean(spG_forc))

boxplot(myrmse_forc,gmrmse_forc,sprmse_forc, main="RMSE forecasting",names = c("Tvar","Gneit-mat","Sep"))
cbind(mean(myrmse_forc),mean(gmrmse_forc),mean(sprmse_forc))

plot(prob_seq,mykp_forc,pch=19,col="blue",cex=0.5,xlab="p-Probability Interval", ylab = "Proportion within Probability Interval",main="forecasting")
points(prob_seq,gmkp_forc,pch=19,col="red",cex=0.5)
points(prob_seq,skp_forc,col="green",cex=0.5)
lines(prob_seq,prob_seq,lwd=2)
legend("bottomright",pch=c(19,19,1),col=c("blue","red","green"),c("Our model","Gneiting's-Matern","Separable Model"))


plot(prob_seq,mywp_forc,pch=19,col="blue",cex=0.5,xlab="p-Probability Interval", ylab = "Probability Interval Width",main="forecasting",ylim=c(0,max(mywp_forc,gmwp_forc,swp_forc)))
points(prob_seq,gmwp_forc,pch=19,col="red",cex=0.5)
points(prob_seq,swp_forc,col="green",cex=0.5)
legend("bottomright",pch=c(19,19,1),col=c("blue","red","green"),c("Our model","Gneiting's-Matern","Separable Model"))



########################################
####### Probabilistic scores ###########
########################################



par(mfrow=c(1,1))
boxplot(my_crps_pred,gm_crps_pred,sp_crps_pred,names = c("Tvar","Gneiting","Separable"),col="blue")
boxplot(my_crps_forc,gm_crps_forc,sp_crps_forc,names = c("Tvar","Gneiting","Separable"),col="blue")

boxplot(my_logs_pred,gm_logs_pred,sp_logs_pred,names = c("Tvar","Gneiting","Separable"),col="blue")
boxplot(my_logs_forc,gm_logs_forc,sp_logs_forc,names = c("Tvar","Gneiting","Separable"),col="blue")




############ Saving all the summary data with c1prefix #######

c1my_crps_pred<-my_crps_pred
c1gm_crps_pred<-gm_crps_pred
c1sp_crps_pred<-sp_crps_pred

c1my_crps_forc<-my_crps_forc
c1gm_crps_forc<-gm_crps_forc
c1sp_crps_forc<-sp_crps_forc


c1my_logs_pred<-my_logs_pred
c1gm_logs_pred<-gm_logs_pred
c1sp_logs_pred<-sp_logs_pred

c1my_logs_forc<-my_logs_forc
c1gm_logs_forc<-gm_logs_forc
c1sp_logs_forc<-sp_logs_forc




c1myG_pred<-myG_pred
c1gmG_pred<-gmG_pred
c1spG_pred<-spG_pred

c1myG_forc<-myG_forc
c1gmG_forc<-gmG_forc
c1spG_forc<-spG_forc



c1myrmse_int<-myrmse_int
c1gmrmse_int<-gmrmse_int
c1sprmse_int<-sprmse_int


c1myrmse_forc<-myrmse_forc
c1gmrmse_forc<-gmrmse_forc
c1sprmse_forc<-sprmse_forc


c1mykp_pred<-mykp_pred
c1gmkp_pred<-gmkp_pred
c1skp_pred<-skp_pred


c1mywp_pred<-mywp_pred
c1gmwp_pred<-gmwp_pred
c1swp_pred<-swp_pred

c1mykp_forc<-mykp_forc
c1gmkp_forc<-gmkp_forc
c1skp_forc<-skp_forc

c1mywp_forc<-mywp_forc
c1gmwp_forc<-gmwp_forc
c1swp_forc<-swp_forc



undelete<-c("c1my_crps_pred",
            "c1gm_crps_pred",
            "c1sp_crps_pred",
            "c1my_crps_forc",
            "c1gm_crps_forc",
            "c1sp_crps_forc",
            "c1my_logs_pred",
            "c1gm_logs_pred",
            "c1sp_logs_pred",
            "c1my_logs_forc",
            "c1gm_logs_forc",
            "c1sp_logs_forc",
            "c1myG_pred",
            "c1gmG_pred",
            "c1spG_pred",
            "c1myG_forc",
            "c1gmG_forc",
            "c1spG_forc",
            "c1myrmse_int",
            "c1gmrmse_int",
            "c1sprmse_int",
            "c1myrmse_forc",
            "c1gmrmse_forc",
            "c1sprmse_forc",
            "c1mykp_pred",
            "c1gmkp_pred",
            "c1skp_pred",
            "c1mywp_pred",
            "c1gmwp_pred",
            "c1swp_pred",
            "c1mykp_forc",
            "c1gmkp_forc",
            "c1skp_forc",
            "c1mywp_forc",
            "c1gmwp_forc",
            "c1swp_forc",
            "prob_seq"
)
rm(list = ls()[!ls()%in%undelete])
save.image("c1visdata.RData")
dir()

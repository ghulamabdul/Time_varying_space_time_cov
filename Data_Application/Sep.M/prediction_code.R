######################################
#### Setting working dirctory ########
######################################

setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Data Application/SeparableModel")
load("sprcl_ests_gmini4.RData")

######################################
######## Loading Libraries ###########
######################################

library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
library(sp)
library(gstat)
library(geoR)
library(fields)
library(MASS)
library(mvtnorm)


ncores<-detectCores()
registerDoParallel(cores = ncores-6)
getDoParWorkers()




########### Reading parameter estimates ##########

mat.a<-rep(rcl.est$par[1],times=tdays)
mat.nu<-rep(rcl.est$par[2],times=tdays)
mat.sigma<-rcl.est$par[3]
bern.a<-rcl.est$par[4]
bern.alpha<-rcl.est$par[5]
bern.beta1<-0
bern.beta2<-rcl.est$par[6]


###################################################

mytest.data1
mytest.data2

###################################################
###### Lets compute the dynamic Sigma22 ###########
###################################################
tr.l
days.window<-12
dynamic.sig22<-list()

#myfas<-f_alphas(t=tseq[1:tdays],alphas = mat.a,a=bern.a,alpha=bern.alpha,beta=bern.beta1,c.p = mean(mat.a))

#image.plot(tseq[1:tdays],tseq[1:tdays],myfas)
#myzeta<-zeta(t=tseq[1:tdays],nus =mat.nu[1:tdays],alphas = mat.a[1:tdays],a=bern.a,alpha=bern.alpha,beta=bern.beta1,beta2 =bern.beta2,c.p = mean(mat.a) )

#image.plot(tseq[1:tdays],tseq[1:tdays],myzeta)
getDoParWorkers()
pred_summary<-foreach(index.day = 1:tdays,.packages=c("fields"))%dopar%
{
dyn.seq<-(index.day-(days.window)/2):(index.day+(days.window)/2)
dyn.seq<-dyn.seq[dyn.seq>=1 & dyn.seq<=tdays]

index.day.loc<-which(dyn.seq==index.day)
###### Computing Fas #########

myfas<-f_alphas(t=tseq[dyn.seq],alphas = mat.a[dyn.seq],a=bern.a,alpha=bern.alpha,beta=bern.beta1,c.p = mean(mat.a))
myzeta<-zeta(t=tseq[dyn.seq],nus =mat.nu[dyn.seq],alphas = mat.a[dyn.seq],a=bern.a,alpha=bern.alpha,beta=bern.beta1,beta2 =bern.beta2,c.p = mean(mat.a) )
mynus<-outer(mat.nu[dyn.seq],mat.nu[dyn.seq],function(x,y) (x+y)/2)

dynamic.sig22<-matrix(NA,nrow = length(dyn.seq)*tr.l,ncol=length(dyn.seq)*tr.l)

tloc_dist<-rdist(cbind(my.res.data[,1],my.res.data[,2]))

for(i in 1:length(dyn.seq))
{
  for(j in i:length(dyn.seq))
  {
    dynamic.sig22[((i-1)*tr.l+1):((i-1)*tr.l+tr.l),((j-1)*tr.l+1):((j-1)*tr.l+tr.l)]<-dynamic.sig22[((j-1)*tr.l+1):((j-1)*tr.l+tr.l),((i-1)*tr.l+1):((i-1)*tr.l+tr.l)]<-myzeta[i,j]*my.matern(h=tloc_dist,
                                                                                                                                                                                                            a=myfas[i,j],
                                                                                                                                                                                                            nu=mynus[i,j],
                                                                                                                                                                                                            sigma=mat.sigma)
  }
}

sig22_inv<-solve(dynamic.sig22)

##### Now computing Sig11 #########

tstloc_dist<-rdist(cbind(mytest.data1$Longitude,mytest.data1$Latitude))


sig11<-myzeta[index.day.loc,index.day.loc]*my.matern(h=tstloc_dist,
                             a=myfas[index.day.loc,index.day.loc],
                             nu=mynus[index.day.loc,index.day.loc],
                             sigma=mat.sigma)


sig12<-matrix(NA,nrow =nrow(sig11) ,ncol = nrow(tloc_dist)*length(dyn.seq))

cross.dist<-rdist(cbind(mytest.data1$Longitude,mytest.data1$Latitude),cbind(my.res.data[,1],my.res.data[,2]))
for(i in 1:length(dyn.seq))
{
  sig12[1:nrow(sig11),((i-1)*nrow(tloc_dist)+1):(i*nrow(tloc_dist))]<-myzeta[index.day.loc,i]*my.matern(h=cross.dist,
                                                                                                                    a=myfas[index.day.loc,i],
                                                                                                                    nu=mynus[index.day.loc,i],
                                                                                                                    sigma=mat.sigma)
}




###########################
### Prediction value ######
###########################

mywts<-sig12%*%sig22_inv
mypred.val<-mywts%*%c(res.mat[,dyn.seq])
mypred.var<-diag(sig11-mywts%*%t(sig12))

data.frame(prediction=mypred.val,variance=mypred.var)
}


sum(do.call(rbind,pred_summary)[,2]<0)
interp.value<-do.call(rbind,pred_summary)[,1]
interp.variance<-do.call(rbind,pred_summary)[,2]

###### Now we do the forecasting #########
index.day<-320
dyn.seq<-index.day:tdays
index.day.loc<-which(dyn.seq==index.day)
###### Computing Fas #########

myfas<-f_alphas(t=tseq[dyn.seq],alphas = mat.a[dyn.seq],a=bern.a,alpha=bern.alpha,beta=bern.beta1,c.p = mean(mat.a))
myzeta<-zeta(t=tseq[dyn.seq],nus =mat.nu[dyn.seq],alphas = mat.a[dyn.seq],a=bern.a,alpha=bern.alpha,beta=bern.beta1,beta2 =bern.beta2,c.p = mean(mat.a) )
mynus<-outer(mat.nu[dyn.seq],mat.nu[dyn.seq],function(x,y) (x+y)/2)

dynamic.sig22forc<-matrix(NA,nrow = length(dyn.seq)*tr.l,ncol=length(dyn.seq)*tr.l)

tloc_dist<-rdist(cbind(my.res.data[,1],my.res.data[,2]))

for(i in 1:length(dyn.seq))
{
  for(j in i:length(dyn.seq))
  {
    dynamic.sig22forc[((i-1)*tr.l+1):((i-1)*tr.l+tr.l),((j-1)*tr.l+1):((j-1)*tr.l+tr.l)]<-dynamic.sig22forc[((j-1)*tr.l+1):((j-1)*tr.l+tr.l),((i-1)*tr.l+1):((i-1)*tr.l+tr.l)]<-myzeta[i,j]*my.matern(h=tloc_dist,
                                                                                                                                                                                              a=myfas[i,j],
                                                                                                                                                                                              nu=mynus[i,j],
                                                                                                                                                                                              sigma=mat.sigma)
  }
}

sig22forc_inv<-solve(dynamic.sig22forc)

##### Now computing Sig11 #########
tstloc_dist<-rdist(cbind(mytest.data2$Longitude,mytest.data2$Latitude))

myfas.forc<-f_alphas(t=tseq,alphas = c(mat.a,rep(mat.a[1],times=10)),a=bern.a,alpha=bern.alpha,beta=bern.beta1,c.p = mean(mat.a))
myzeta.forc<-zeta(t=tseq,nus =c(mat.nu,rep(mat.nu[1],times=10)),alphas = c(mat.a,rep(mat.a[1],times=10)),a=bern.a,alpha=bern.alpha,beta=bern.beta1,beta2 =bern.beta2,c.p = mean(mat.a) )
mynus.forc<-outer(c(mat.nu,rep(mat.nu[1],times=10)),c(mat.nu,rep(mat.nu[1],times=10)),function(x,y) (x+y)/2)

trl.forc<-length(mytest.data2$Longitude)
sig11forc<-matrix(NA,nrow=10*length(mytest.data2$Longitude),ncol = 10*length(mytest.data2$Longitude))

forc.index<-(tdays+1):365
for(i in 1:10)
{
  for(j in i:10)
  {
    sig11forc[((i-1)*trl.forc+1):((i-1)*trl.forc+trl.forc),((j-1)*trl.forc+1):((j-1)*trl.forc+trl.forc)]<-sig11forc[((j-1)*trl.forc+1):((j-1)*trl.forc+trl.forc),((i-1)*trl.forc+1):((i-1)*trl.forc+trl.forc)]<-myzeta.forc[forc.index[i],forc.index[j]]*my.matern(h=tstloc_dist,
                                                                                                                                                                                              a=myfas.forc[forc.index[i],forc.index[j]],
                                                                                                                                                                                              nu=mynus.forc[forc.index[i],forc.index[j]],
                                                                                                                                                                                              sigma=mat.sigma)
  }
}

sum(is.na(sig11forc))
cross.distforc<-rdist(cbind(mytest.data2$Longitude,mytest.data2$Latitude),cbind(my.res.data[,1],my.res.data[,2]))
dim(cross.distforc)[2]*length(dyn.seq)
sig12forc<-matrix(NA,nrow =nrow(sig11forc) ,ncol = ncol(cross.distforc)*length(dyn.seq))
dim(sig22forc_inv)
dim(sig12forc)

for(i in 1:length(forc.index))
{
  for(j in 1:length(dyn.seq))
  {
    sig12forc[((i-1)*nrow(cross.distforc)+1):(i*nrow(cross.distforc)),((j-1)*ncol(cross.distforc)+1):(j*ncol(cross.distforc))]<-myzeta.forc[forc.index[i],dyn.seq[j]]*my.matern(h=cross.distforc,
                                                                                                                                                                  a=myfas.forc[forc.index[i],dyn.seq[j]],
                                                                                                                                                                  nu=mynus.forc[forc.index[i],dyn.seq[j]],
                                                                                                                                                                  sigma=mat.sigma)
  }
}

dim(sig12forc)
sum(is.na(sig12forc))

myforc.wts<-sig12forc%*%sig22forc_inv
myforc.val<-myforc.wts%*%c(res.mat[,dyn.seq])
myforc.variance<-diag(sig11forc-myforc.wts%*%t(sig12forc))
sum(myforc.variance<0)
quilt.plot(mytest.data2$Longitude,mytest.data2$Latitude,mytest.data2$z356)
i=1
quilt.plot(mytest.data2$Longitude,mytest.data2$Latitude,myforc.val[((i-1)*805+1):(i*805)])


############### Computing prediction scores ###############

interp.true<-c(as.matrix(mytest.data1[,-c(1,2)]))
forecast.true<-c(as.matrix(mytest.data2[,-c(1,2)]))

library(scoringRules)

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


pred_stats<-function(true,predvalue,predvar)
{
  rmse<-sqrt(mean((true-predvalue)^2))
  mcrps<-mean(crps_norm(y=true,mean = predvalue,sd=sqrt(predvar)))
  mlogs<-mean(logs_norm(y=true,mean = predvalue,sd=sqrt(predvar)))
  cat("RMSE=", rmse, "\n","mCRPS=",mcrps, "\n","mLogs=", mlogs)
  
  prob_seq<-seq(0.01,0.99,length.out = 2*99)
  
  mywp<-mykp<-numeric()
  for(i in 1:length(prob_seq))
  {
    temp<-uncert_p(pr=prob_seq[i],truev = true,predval = predvalue,predvar =predvar )
    mywp[i]<-temp$wp
    mykp[i]<-temp$kp
    
  }
  
  myG<-G_stat(p=prob_seq,kp=mykp)
  cat("\nGstat=", myG,"\n")
  
  par(mfrow=c(1,2))
  {
    plot(prob_seq,mykp,type="l",lwd=2,xlab="p-Probability Interval", ylab = "Proportion within Probability Interval")
    plot(prob_seq,mywp,type="l",lwd=2,xlab="p-Probability Interval", ylab = "Probability Interval Width")
    
  }
  return(list(RMSE=rmse,MCRPS=mcrps,MLOGS=mlogs,WP=mywp,KP=mykp,G=myG))
}

forecast_summary<-pred_stats(true = forecast.true,predvalue = myforc.val,predvar = myforc.variance)
interpolation_summary<-pred_stats(true = interp.true,predvalue =interp.value,predvar = interp.variance )
rm(dynamic.sig22forc,myforc.wts,sig11forc,sig12forc,sig22forc_inv)
save.image("Separable_Prediction_summary6d320f.RData")

######### Daywise forecasting scores ###########
forc_dayw_true<-as.matrix(mytest.data2[,-c(1,2)])
forc_dayw_value<-matrix(myforc.val,ncol = ncol(forc_dayw_true),byrow = F)
forc_dayw_variance<-matrix(myforc.variance,ncol = ncol(forc_dayw_true),byrow = F)


forc_summary_daywise<-list()
forc_dw_rmse<-forc_dw_mcrps<-forc_dw_mlogs<-forc_dw_G<-numeric(10)
for(i in 1:10)
{
  forc_summary_daywise[[i]]<-pred_stats(true = forc_dayw_true[,i],predvalue = forc_dayw_value[,i],predvar = forc_dayw_variance[,i])
  forc_dw_rmse[i]<-forc_summary_daywise[[i]]$RMSE
  forc_dw_mcrps[i]<-forc_summary_daywise[[i]]$MCRPS
  forc_dw_mlogs[i]<-forc_summary_daywise[[i]]$MLOGS
  forc_dw_G[i]<-forc_summary_daywise[[i]]$G
}
save.image("Separable_Prediction_summary20fv2.RData")



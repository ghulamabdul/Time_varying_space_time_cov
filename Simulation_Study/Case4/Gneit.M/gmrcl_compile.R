setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case4_2")
load("simulated_data.RData")
# setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Case4_2/FInal piece of codes/GneitingModel")
library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
ncores<-detectCores()
registerDoParallel(cores = ncores-4)
setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case4_2/GneitingModel")

#### myrcl_compile ####
load("gmrclp1.RData")
gmrcl1<-gmrcl
load("gmrclp2.RData")
gmrcl2<-gmrcl
load("gmrclp3.RData")
gmrcl3<-gmrcl

gmrcl<-do.call(c,list(gmrcl1,gmrcl2,gmrcl3))
#getwd()


##### Extracting parameter values ####
alphas.rcl<-nus.rcl<-bern_alpha.rcl<-bern_beta1.rcl<-bern_beta2.rcl<-sigma.rcl<-numeric(length = nsim)
count_f<-numeric(length = nsim)
myest.a.rcl<-myest.nu.rcl<-matrix(NA,nrow = nsim,ncol = tl)
for(i in 1:nsim)
{
  
  alphas.rcl[i]<-gmrcl[[i]]$estimates$par[1]
  nus.rcl[i]<-gmrcl[[i]]$estimates$par[2]
  sigma.rcl[i]<-gmrcl[[i]]$estimates$par[3]
  
  bern_alpha.rcl[i]<-gmrcl[[i]]$estimates$par[4]
  bern_beta1.rcl[i]<-gmrcl[[i]]$estimates$par[5]
  bern_beta2.rcl[i]<-gmrcl[[i]]$estimates$par[6]
  
  myest.a.rcl[i,]<-alphas.rcl[i]
  myest.nu.rcl[i,]<-nus.rcl[i]
  count_f[i]<-gmrcl[[i]]$estimates$counts[1]
}


hist(count_f)

sigma.rcl

round(mean(sigma.rcl),2)
round(sd(sigma.rcl),2)


hist(bern_alpha.rcl)
round(mean(bern_alpha.rcl),2)
round(sd(bern_alpha.rcl),2)


hist(bern_beta1.rcl)

round(mean(bern_beta1.rcl),2)
round(sd(bern_beta1.rcl),2)


hist(bern_beta2.rcl)

round(mean(bern_beta2.rcl),2)
round(sd(bern_beta2.rcl),2)


hist(alphas.rcl)

round(mean(alphas.rcl),2)
round(sd(alphas.rcl),2)


hist(nus.rcl)


round(mean(nus.rcl),2)
round(sd(nus.rcl),2)



####### Now we do the spatia-temporal prediction ########
save.image("pointa.RData")

mypredres<-foreach(simkey=1:nsim,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm","scoringRules"))%dopar%{
  #bern.est.alpha
  myest.spcov<-matrix(NA,nrow = tl*sl,ncol = tl*sl)
  
  runi.fas<-f_alphas(t=tseq,alphas = myest.a.rcl[simkey,],a=mya_bern,alpha=bern_alpha.rcl[simkey],beta=bern_beta1.rcl[simkey],c.p = mean(myest.a.rcl[simkey,]))
  runi.zeta<-zeta(t=tseq,nus = myest.nu.rcl[simkey,],alphas = myest.a.rcl[simkey,],a=mya_bern,alpha=bern.est.alpha[simkey],beta=bern_beta1.rcl[simkey],c.p = mean(myest.a.rcl[simkey,]),beta2 = bern_beta2.rcl[simkey])
  #par(mfrow=c(1,2))
  #image.plot(tseq,tseq,runi.fas,xlab="t",ylab="t",main=bquote(f[as]~"("~t[i]~","~t[j]~")"))
  #image.plot(tseq,tseq,runi.zeta,xlab="t",ylab="t",main=bquote(zeta~"("~t[i]~","~t[j]~")"))
  runi.nus<-outer(myest.nu.rcl[simkey,],myest.nu.rcl[simkey,],function(a,b) (a+b)/2)
  
  
  
  ##########################################################
  ##### Creating  spatio-temporal covariance function ######
  ##########################################################
  
  temp<-matrix(NA,ncol = sl,nrow = sl)
  
  
  system.time(for(i in 1:tl)
  {
    for(j in i:tl)
    {
      tmp1<-runi.zeta[i,j]*my.matern(h=uniq.dist,
                                     a=runi.fas[i,j],
                                     nu=runi.nus[i,j],
                                     sigma=sigma.rcl[simkey])
      temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
      myest.spcov[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-myest.spcov[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
    }
  }
  )
  
  ###### Now computing sigma11, sigma12 and sigma22 ######
  pred.index1<-NULL
  tmp3<-NULL
  for(i in 1:length(train.myspat[1,,1]))
  {
    tmp3<-rmiss.matrix[,simkey]+(i-1)*sl
    pred.index1<-c(pred.index1,tmp3)
  }
  
  ### appending the indices of the last two frames ####
  pred.index1<-c(pred.index1,(sl*length(train.myspat[1,,1])+1):(tl*sl))
  
  SIGMA.11<-myest.spcov[pred.index1,pred.index1]
  SIGMA.12<-myest.spcov[pred.index1,-pred.index1]
  SIGMA.22<-myest.spcov[-pred.index1,-pred.index1]
  rm(myest.spcov)
  SIGMA.INV<-solve(SIGMA.22)
  OBSDATA<-sim1[simkey,-pred.index1]
  WTS<-SIGMA.12%*%SIGMA.INV
  
  PREDS_ST<-WTS%*%OBSDATA
  PREDS_VAR<-diag(SIGMA.11-SIGMA.12%*%SIGMA.INV%*%t(SIGMA.12))
  INTERP.PRED<-PREDS_ST[1:(test_l*length(train.myspat[1,,1]))]
  FORECAST.PRED<-PREDS_ST[-(1:(test_l*length(train.myspat[1,,1])))]
  
  INTERP.VAR<-PREDS_VAR[1:(test_l*length(train.myspat[1,,1]))]
  FORECAST.VAR<-PREDS_VAR[-(1:(test_l*length(train.myspat[1,,1])))]
  
  rm(WTS,SIGMA.11,SIGMA.22,SIGMA.12,SIGMA.INV)
  my_interp.mat<-matrix(INTERP.PRED,ncol=length(train.myspat[1,,1]))
  my_interp.var.mat<-matrix(INTERP.VAR,ncol=length(train.myspat[1,,1]))
  #image.plot(1:125,1:19,my_interp.mat)
  #image.plot(1:125,1:19,test.myspat[,1:19,simkey])
  
  myforecast.mat<-matrix(FORECAST.PRED,ncol=2)
  myforecast.var.mat<-matrix(FORECAST.VAR,ncol=2)
  
  #quilt.plot(grid,full.myspat[,20,simkey],nx=25,ny=25)
  #quilt.plot(grid,myforecast.mat[,1],nx=25,ny=25)
  
  RMSE_interp<-sqrt(mean((my_interp.mat-test.myspat[,1:19,simkey])^2))
  
  mcrps_interp<-mean(crps_norm(y=c(test.myspat[,1:19,simkey]),mean = c(my_interp.mat),sd=sqrt(c(my_interp.var.mat))))
  mlogs_interp<-mean(logs_norm(y=c(test.myspat[,1:19,simkey]),mean = c(my_interp.mat),sd=sqrt(c(my_interp.var.mat))))
  
  RMSE_forecast<-sqrt(mean((myforecast.mat-full.myspat[,20:21,simkey])^2))
  
  mcrps_forc<-mean(crps_norm(y=c(full.myspat[,20:21,simkey]),mean = c(myforecast.mat),sd=sqrt(c(myforecast.var.mat))))
  mlogs_forc<-mean(logs_norm(y=c(full.myspat[,20:21,simkey]),mean = c(myforecast.mat),sd=sqrt(c(myforecast.var.mat))))
  
  
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
  
  
  #uncert_p(pr=0.70,truev = c(test.myspat[,1:19,simkey]),predval = c(my_interp.mat),predvar =c(my_interp.var.mat) )
  
  #### Interpolation ####
  prob_seq<-seq(0.01,0.99,length.out = 2*99)
  
  wp<-kp<-numeric()
  for(i in 1:length(prob_seq))
  {
    wp[i]<-uncert_p(pr=prob_seq[i],truev = c(test.myspat[,1:19,simkey]),predval = c(my_interp.mat),predvar =c(my_interp.var.mat) )$wp
    kp[i]<-uncert_p(pr=prob_seq[i],truev = c(test.myspat[,1:19,simkey]),predval = c(my_interp.mat),predvar =c(my_interp.var.mat) )$kp
    
  }
  
  my_wp_pred<-wp
  my_kp_pred<-kp
  
  #plot(prob_seq,my_wp_pred)
  #plot(prob_seq,my_kp_pred)
  #lines(prob_seq,prob_seq)
  myG_pred<-G_stat(p=prob_seq,kp=my_kp_pred)
  
  
  
  wp2<-kp2<-numeric()
  for(i in 1:length(prob_seq))
  {
    wp2[i]<-uncert_p(pr=prob_seq[i],truev = c(full.myspat[,20:21,simkey]),predval = c(myforecast.mat),predvar =c(myforecast.var.mat) )$wp
    kp2[i]<-uncert_p(pr=prob_seq[i],truev = c(full.myspat[,20:21,simkey]),predval = c(myforecast.mat),predvar =c(myforecast.var.mat) )$kp
    
  }
  
  
  my_wp_forc<-wp2
  my_kp_forc<-kp2
  
  #plot(prob_seq,my_wp_forc,cex=0.5)
  #plot(prob_seq,my_kp_forc,cex=0.5)
  #lines(prob_seq,prob_seq)
  myG_forc<-G_stat(p=prob_seq,kp=my_kp_forc)
  
  fin.output<-list(my_wp_pred=my_wp_pred,my_kp_pred=my_kp_pred,myG_pred=myG_pred,
                   my_wp_forc=my_wp_forc,my_kp_forc=my_kp_forc,myG_forc=myG_forc,
                   RMSE_pred=RMSE_interp,RMSE_forc=RMSE_forecast,
                   mLogS_pred=mlogs_interp,mLogS_forc=mlogs_forc,
                   mCRPS_pred=mcrps_interp,mCRPS_forc=mcrps_forc,
                   pred.val_pred=my_interp.mat,pred.var_pred=my_interp.var.mat,
                   pred.val_forc=myforecast.mat,pred.var_forc=myforecast.var.mat,
                   testval_pred=test.myspat[,1:19,simkey],testval_forc=full.myspat[,20:21,simkey])
  
  fin.output
}


save.image("gneitingmodel_rcl_pred.RData")

load("gneitingmodel_rcl_pred.RData")

library(scoringRules)
gmpred_crps<-gmpred_logs<-gmforc_crps<-gmforc_logs<-numeric()

for(i in 1:100)
{
  gmpred_crps[i]<-mean(crps_norm(y=c(mypredres[[i]]$testval_pred),mean = c(mypredres[[i]]$pred.val_pred),sd=sqrt(c(mypredres[[i]]$pred.var_pred))))
  gmpred_logs[i]<-mean(logs_norm(y=c(mypredres[[i]]$testval_pred),mean = c(mypredres[[i]]$pred.val_pred),sd=sqrt(c(mypredres[[i]]$pred.var_pred))))

  gmforc_crps[i]<-mean(crps_norm(y=c(mypredres[[i]]$testval_forc),mean = c(mypredres[[i]]$pred.val_forc),sd=sqrt(c(mypredres[[i]]$pred.var_forc))))
  gmforc_logs[i]<-mean(logs_norm(y=c(mypredres[[i]]$testval_forc),mean = c(mypredres[[i]]$pred.val_forc),sd=sqrt(c(mypredres[[i]]$pred.var_forc))))

}


save.image("gmrcl_case1.RData")

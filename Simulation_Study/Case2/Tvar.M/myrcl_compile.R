setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case2_2")
load("simulated_data.RData")
setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case2_2/MyModel")
#setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Case2_2/FInal piece of codes/MyModel")
library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
ncores<-detectCores()
registerDoParallel(cores = ncores-4)

#### myrcl_compile ####
load("myrclestsp1.RData")
myrcl_estims1<-myrcl_estims
load("myrclestsp2.RData")
myrcl_estims2<-myrcl_estims
load("myrclestsp3.RData")
myrcl_estims3<-myrcl_estims

myrcl_f<-do.call(c,list(myrcl_estims1,myrcl_estims2,myrcl_estims3))
getwd()


##### Extracting parameter values ####
a.0.rcl<-a.1.rcl<-a.2.rcl<-nu.0.rcl<-nu.1.rcl<-nu.2.rcl<-bern_alpha.rcl<-bern_beta1.rcl<-bern_beta2.rcl<-sigma.rcl<-numeric(length = nsim)
count_f<-numeric(length = nsim)
myest.a.rcl<-myest.nu.rcl<-matrix(NA,nrow = nsim,ncol = tl)
for(i in 1:nsim)
{
  
  a.0.rcl[i]<-myrcl_f[[i]]$estimates$par[1]
  a.1.rcl[i]<-myrcl_f[[i]]$estimates$par[2]
  a.2.rcl[i]<-myrcl_f[[i]]$estimates$par[3]
  nu.0.rcl[i]<-myrcl_f[[i]]$estimates$par[4]
  nu.1.rcl[i]<-myrcl_f[[i]]$estimates$par[5]
  nu.2.rcl[i]<-myrcl_f[[i]]$estimates$par[6]
  sigma.rcl[i]<-myrcl_f[[i]]$estimates$par[7]

  bern_alpha.rcl[i]<-myrcl_f[[i]]$estimates$par[8]
  bern_beta1.rcl[i]<-myrcl_f[[i]]$estimates$par[9]
  bern_beta2.rcl[i]<-myrcl_f[[i]]$estimates$par[10]
  
  myest.a.rcl[i,]<-exp(a.0.rcl[i]+a.1.rcl[i]*tseq[1:tl]+a.2.rcl[i]*((tseq[1:tl])^2))
  myest.nu.rcl[i,]<-exp(nu.0.rcl[i]+nu.1.rcl[i]*tseq[1:tl]+nu.2.rcl[i]*(tseq[1:tl]^2))
  count_f[i]<-myrcl_f[[i]]$estimates$counts[1]
}

plot(tseq,alphas,type="l",col="grey",ylim=c(min(myest.a),max(myest.a)))
plot(tseq,alphas,type="l",col="grey",ylim=c(5,35),xlab = bquote(t[i]),ylab = bquote(alpha[s](t[i])))

points(tseq,alphas,pch=19)
myest.a.rcl.sd<-apply(myest.a.rcl, 2, sd)
myest.a.rcl.mean<-apply(myest.a.rcl, 2, mean)
lines(tseq,myest.a.rcl.mean,col="green",lwd=1)
points(tseq,myest.a.rcl.mean,col="blue",lwd=1,pch=19,cex=0.7)
lines(tseq,myest.a.rcl.mean+1.96*myest.a.rcl.sd,col="red",lwd=1,lty=2)
lines(tseq,myest.a.rcl.mean-1.96*myest.a.rcl.sd,col="red",lwd=1,lty=2)
legend("topright",pch = c(19,19),c("True","Average (estimate)"),col=c("black","blue"))
#lines(tseq,alphas_rcl[11,])

plot(tseq,nus,type="l",col="grey",ylim=c(min(myest.nu),max(myest.nu)),xlab = bquote(t[i]),ylab = bquote(nu[s](t[i])))
points(tseq,nus,pch=19)
myest.nu.rcl.sd<-apply(myest.nu.rcl, 2, sd)
myest.nu.rcl.mean<-apply(myest.nu.rcl, 2, mean)
lines(tseq,myest.nu.rcl.mean,col="green",lwd=1)
points(tseq,myest.nu.rcl.mean,col="blue",lwd=1,pch=19,cex=0.7)
lines(tseq,myest.nu.rcl.mean+1.96*myest.nu.rcl.sd,col="red",lwd=1,lty=2)
lines(tseq,myest.nu.rcl.mean-1.96*myest.nu.rcl.sd,col="red",lwd=1,lty=2)
legend("topright",pch = c(19,19),c("True","Average (estimate)"),col=c("black","blue"))


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
c4true_alpha<-alphas
c4myalpha.mean<-myest.a.rcl.mean
c4myalpha.sd<-myest.a.rcl.sd
c4true_nus<-nus
c4mynu.mean<-myest.nu.rcl.mean
c4mynu.sd<-myest.nu.rcl.sd

#rm(list = ls()[!ls()%in%c("tseq","c4true_alpha","c4myalpha.mean","c4myalpha.sd","c4true_nus","c4mynu.mean","c4mynu.sd")])
#save.image("c4myalphnu.RData")



#myest.a.pred<-myest.nu.pred<-matrix(NA,nrow = nsim,ncol=tl-length(train.myspat[1,,1]))
#### First interpolating estimated alphas using AR 1 process
#for(i in 1:nsim)
#{
 # a.0<-mymargins[[i]]$estimate$par[1]
 #  a.1<-mymargins[[i]]$estimate$par[2]
 #  a.2<-mymargins[[i]]$estimate$par[3]
 #  nu.0<-mymargins[[i]]$estimate$par[4]
 #  nu.1<-mymargins[[i]]$estimate$par[5]
 #  nu.2<-mymargins[[i]]$estimate$par[6]
  
#  myest.a.pred[i,]<-exp(a.0+a.1*tseq[(length(train.myspat[1,,1])+1):tl]+a.2*tseq[(length(train.myspat[1,,1])+1):tl]^2)
#  myest.nu.pred[i,]<-exp(nu.0+nu.1*tseq[(length(train.myspat[1,,1])+1):tl]+nu.2*tseq[(length(train.myspat[1,,1])+1):tl]^2)
  
#}


myest.full.a<-myest.a.rcl
myest.full.nu<-myest.nu.rcl

full_cmean_a<-colMeans(myest.full.a)
full_cmean_nu<-colMeans(myest.full.nu)

full_sd_a<-apply(myest.full.a,2,sd)
full_sd_nu<-apply(myest.full.nu,2,sd)

plot(tseq,alphas,col="grey",type="l",ylim=range(c(myest.full.a)))
points(tseq,alphas,pch=19)
points(tseq,full_cmean_a,pch=19,col="blue")
lines(tseq,full_cmean_a+1.96*full_sd_a,pch=19,col="blue",lty=2)
lines(tseq,full_cmean_a-1.96*full_sd_a,pch=19,col="blue",lty=2)







plot(tseq,nus,col="grey",type="l",ylim=range(c(nus,full_cmean_nu)))
points(tseq,nus,pch=19)
points(tseq,full_cmean_nu,pch=19,col="blue")
lines(tseq,full_cmean_nu+1.96*full_sd_nu,pch=19,col="blue",lty=2)
lines(tseq,full_cmean_nu-1.96*full_sd_nu,pch=19,col="blue",lty=2)






####### Now we do the spatia-temporal prediction ########


mypredres<-foreach(simkey=1:nsim,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm","scoringRules"))%dopar%{
  #bern.est.alpha
  myest.spcov<-matrix(NA,nrow = tl*sl,ncol = tl*sl)
  
  runi.fas<-f_alphas(t=tseq,alphas = myest.full.a[simkey,],a=mya_bern,alpha=bern_alpha.rcl[simkey],beta=bern_beta1.rcl[simkey],c.p = mean(myest.a.rcl[simkey,1:length(train.myspat[1,,1])]))
  runi.zeta<-zeta(t=tseq,nus = myest.full.nu[simkey,],alphas = myest.full.a[simkey,],a=mya_bern,alpha=bern_alpha.rcl[simkey],beta=bern_beta1.rcl[simkey],c.p = mean(myest.a.rcl[simkey,1:length(train.myspat[1,,1])]),beta2 = bern_beta2.rcl[simkey])
  #par(mfrow=c(1,2))
  #image.plot(tseq,tseq,runi.fas,xlab="t",ylab="t",main=bquote(f[as]~"("~t[i]~","~t[j]~")"))
  #image.plot(tseq,tseq,runi.zeta,xlab="t",ylab="t",main=bquote(zeta~"("~t[i]~","~t[j]~")"))
  runi.nus<-outer(myest.full.nu[simkey,],myest.full.nu[simkey,],function(a,b) (a+b)/2)
  
  
  
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
  RMSE_forecast<-sqrt(mean((myforecast.mat-full.myspat[,20:21,simkey])^2))
  mcrps_interp<-mean(crps_norm(y=c(test.myspat[,1:19,simkey]),mean = c(my_interp.mat),sd=sqrt(c(my_interp.var.mat))))
  mlogs_interp<-mean(logs_norm(y=c(test.myspat[,1:19,simkey]),mean = c(my_interp.mat),sd=sqrt(c(my_interp.var.mat))))
  
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


save.image("mymodel_rcl_pred.RData")





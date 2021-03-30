setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case4_2")
load("simulated_data.RData")

library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
ncores<-detectCores()
registerDoParallel(cores = ncores-4)
### functions ####

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*a)^nu
  num3<-besselK(x=(h*a), nu=nu)
  return(num1*num2*num3)
}



psi<-function(t,a=1,alpha=1,beta=1)
{
  temp<-(a*(t^alpha)+1)^beta
  return(temp)
}



##### Creating function for f_alphas #######
f_alphas<-function(t,alphas,a=1,alpha=1,beta=1,c.p)
{
  alpha_mean<-c.p
  len<-length(alphas)
  al.creat<-function(a1,a2)
  {
    return(((1/(a1^2))+(1/(a2^2)))/2)
  }
  d2<-outer(alphas,alphas,al.creat)
  tmat<-outer(t,t,function(t1,t2) (abs(t1-t2))^2)
  d1<-psi(t=tmat,a=a,alpha = alpha,beta=beta)/(alpha_mean^2)
  d3<-matrix(psi(t=0,a=a,alpha = alpha,beta=beta)/(alpha_mean^2),nrow=nrow(d1),ncol=ncol(d1))
  return(1/sqrt(d1+d2-d3))
}


zeta<-function(t,nus,alphas,a=0.001,alpha=0.9,beta=0.5,d=2,beta2=0.5,c.p)
{
  n2<-outer(nus,nus,function(a,b) gamma((a+b)/2))
  f<-(f_alphas(t=t,alphas = alphas,a=a,alpha=alpha,beta=beta,c.p = c.p))
  n1<-f^d
  dig<-(diag(f))^(d/2)
  #d1<-dig%*%dig
  d1<-outer(dig,dig,function(a,b) a*b)
  d2<-outer(nus,nus,function(a,b) sqrt(gamma(a)*gamma(b)))
  #n2<-outer(alphas,alphas,function(a,b) ((1/a)^(d/4))+((1/b)^(d/4)))
  #d2<-(f_alphas(t=t,alphas = alphas,a=a,alpha=alpha,beta=beta))^(d/2)
  tmat<-outer(t,t,function(t1,t2) (abs(t1-t2))^2)
  d3<-psi(t=tmat,a=a,alpha = alpha,beta=beta2)
  
  return((n1*n2)/(d1*d2*d3))
}



##### Marginal estimation #####


############## stepwise MLE ##########
mymargins<-foreach(simkey=1:nsim,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm"))%dopar%
  {
    
    start.time<-proc.time()
    ls.a<-ls.nu<-ls.sigma<-numeric(length = (tl-2))
    for(i in 1:(tl-2))
    {
      mytrial<-cbind(train.locs[,,simkey],train.myspat[,i,simkey])
      mytrial<-as.geodata(mytrial)
      myvario<-variog(mytrial)
      
      
      my.var.loss<-function(p,vardist=myvario$u,gamma_hat=myvario$v,nbins=myvario$n)
      {
        a<-p[1]
        nu<-p[2]
        sigma<-p[3]
        if(sum(p<=0,p[2]>4,p[1]>30)!=0)
        {
          return(Inf)
        }
        else
        {
          gamma_theta<-(my.matern(h=0,a=a,nu=nu,sigma = sigma)-my.matern(h=vardist,a=a,nu=nu,sigma = sigma))
          nk<-nbins
          temp<-nk*(((gamma_hat-gamma_theta)/gamma_theta)^2)
          return(sum(temp))
        }
      }
      temp<-optim(my.var.loss,par = c(runif(n=1,min=0.1,max=29),0.5,sd(train.myspat[,i,simkey])))
      ls.a[i]<-temp$par[1]
      ls.nu[i]<-temp$par[2]
      ls.sigma[i]<-temp$par[3]
      
      
      
    }
    
    
    
    
    
    
    ########## Using composite loglikelihood ###########
    
    my.comp.ini2<-c(mean(ls.a),mean(ls.nu),mean(ls.sigma))
    
    
    #### Computing mat.index for train.locs ####
    
    mytrain.dist<-rdist(train.locs[,,simkey])
    
    ######### Optimized cov computing ##############
    uniq.dist.tr<-sort(unique(c(mytrain.dist)))
    
    ##### Unique indexing #####
    mat.index.tr<-NULL
    for(i in 1:length(uniq.dist.tr))
    {
      mat.index.tr[[i]]<-cbind(which(mytrain.dist==uniq.dist.tr[i],arr.ind = T),i)
    }
    #cbind(which(dist.sp==uniq.dist[2],arr.ind = T),2)
    mat.index.tr<-do.call(rbind,mat.index.tr)
    
    
    
    
    ####### Composite MLE ######
    
    
    ##### Now creating maximum likelihood function ####
    
    mle.comp.all<-function(locs,z,p,tl=(tl-2),sl=train_l,mat.index=mat.index.tr,uniq.dist=uniq.dist.tr)
    {
      #distmat<-rdist(locs)
      a<-p[1]
      nu<-p[2]
      sigma<-p[3]
      
      if(sum(c(p<=0))!=0)
      {
        return(list(mlv=Inf))
      }
      else
      {
        nlogl<-0
        for(i in 1:tl)
        {
          C<-matrix(NA,nrow = sl,ncol = sl)
          tspec.d<-z[((i-1)*sl+1):((i-1)*sl+sl)]  
          Ctemp<-my.matern(h=uniq.dist,a=a,nu=nu,sigma = sigma)
          C[mat.index[,-3]]<-Ctemp[mat.index[,3]]
          nlogl<-nlogl-dmvnorm(x=tspec.d,mean=(rep(0,times=length(tspec.d))),sigma = C,log = T)
        }
        return(list(mlv=nlogl))
      }
    }
    
    #ini.v<-c(ini.a,ini.nu,ini.sig)
    #mle.all(locs=cbind(uk.train[[1]]$s1,uk.train[[1]]$s2),z=lmfit$residuals,p=ini.v)
    mle.comp_only_mlv<-function(par)
    {
      return(mle.comp.all(locs=train.locs[,,simkey],z=c(train.myspat[,,simkey]),p=par,tl=(tl-2),sl=train_l,mat.index=mat.index.tr,uniq.dist=uniq.dist.tr)$mlv
      )
    }
    #mle_only_mlv(c(alphas,nus,1))
    #mle.comp_only_mlv(c(alphas,nus,1))
    #system.time(mle.comp_only_mlv(my.comp.ini2))
    optim_marg_comp.loglik <- function(par){
      optim(par=par,
            fn = mle.comp_only_mlv,
            hessian=FALSE,
            control=list(trace=6,
                         pgtol=0,
                         parscale=rep(0.1,length(par)),
                         maxit=7000))
    }
    
    
    
    mysp.est<-optim_marg_comp.loglik(my.comp.ini2)
    end.time<-proc.time()
    
    list(estimate=mysp.est,time=(end.time-start.time))
    
  }

save.image("separablemarginals.RData")
setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case4_2/SeparableModel")
save.image("separablemarginals.RData")

########## Temporal correlation ###########

loc.batch<-list()
########## Now we do this 10 times ######
est.rep<-15
nbatch<-20
set.seed(123)

loc.batch.ind<-array(NA,dim = c(train_l/nbatch,est.rep*nbatch,nsim))
for(k in 1:nsim)
{
  for(i in 1:est.rep)
  {
    loc.order<-sample(1:train_l,train_l)
    loc.batch[[i]]<-matrix(loc.order,nrow =train_l/nbatch )
  }
  loc.batch.ind[,,k]<-do.call(cbind,loc.batch)
}

#### Extracting marginal estimates ######

myest.sigma<-numeric(length = nsim)
myest.a<-myest.nu<-matrix(nrow = nsim,ncol = length(train.myspat[1,,1]))

for(i in 1:nsim)
{
  myest.sigma[i]<-mymargins[[i]]$estimate$par[3]
  myest.a[i,]<-mymargins[[i]]$estimate$par[1]
  myest.nu[i,]<-mymargins[[i]]$estimate$par[2]
}

#############################################
##### Plotting marginal estimates ###########

#plot(tseq,alphas,type="l",col="grey",ylim=c(min(myest.a),max(myest.a)))
#points(tseq,alphas,pch=19)
#myest.a.sd<-apply(myest.a, 2, sd)
#myest.a.mean<-apply(myest.a, 2, mean)
#lines(tseq[1:length(train.myspat[1,,1])],myest.a.mean,col="green",lwd=1)
#points(tseq[1:length(train.myspat[1,,1])],myest.a.mean,col="blue",lwd=1,pch=19)
#lines(tseq[1:length(train.myspat[1,,1])],myest.a.mean+1.96*myest.a.sd,col="red",lwd=1,lty=2)
#lines(tseq[1:length(train.myspat[1,,1])],myest.a.mean-1.96*myest.a.sd,col="red",lwd=1,lty=2)

#plot(tseq,nus,type="l",col="grey",ylim=c(min(myest.nu),max(myest.nu)))
#points(tseq,nus,pch=19)
#myest.nu.sd<-apply(myest.nu, 2, sd)
#myest.nu.mean<-apply(myest.nu, 2, mean)
#lines(tseq[1:length(train.myspat[1,,1])],myest.nu.mean,col="green",lwd=1)
#points(tseq[1:length(train.myspat[1,,1])],myest.nu.mean,col="blue",lwd=1,pch=19)
#lines(tseq[1:length(train.myspat[1,,1])],myest.nu.mean+1.96*myest.nu.sd,col="red",lwd=1,lty=2)
#lines(tseq[1:length(train.myspat[1,,1])],myest.nu.mean-1.96*myest.nu.sd,col="red",lwd=1,lty=2)


#ag.est.sigma<-mean(est.sigma)
#mytemporal[[1]]$estimates
mytemporal<-foreach(simkey=1:nsim,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm"))%dopar%{
start.time<-proc.time()
mle.temporal.all<-function(locs,myspat1,p,t.points=tseq[1:length(train.myspat[1,,1])],a.est=myest.a[simkey,],nu.est=myest.nu[simkey,],sigma.est=myest.sigma[simkey],mysl=train_l,mytl=length(train.myspat[1,,1]),mynbatch=nbatch,myloc.batch=loc.batch.ind[,,simkey])
{
  #distmat<-rdist(locs)
  mya_bern1<-10
  myalpha_bern1<-p[1]
  mybeta_bern1<-0
  mybeta2_bern1<-p[2]
  
  if(sum(c(mya_bern1<=0,myalpha_bern1<=0,myalpha_bern1>1,mybeta_bern1<0,mybeta_bern1>1,mybeta2_bern1<0))!=0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    
    
    myfas.e<-f_alphas(t=t.points,alphas = a.est,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,c.p = mean(a.est))
    
    myzeta.e<-zeta(t=t.points,nus = nu.est,alphas = a.est,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,beta2 =mybeta2_bern1,c.p = mean(a.est) )
    
    indic<-ncol(myloc.batch)
    myloglv<-0
    for(sr.batch in 1:indic)
    {
    mynus.e<-outer(nu.est,nu.est,function(a,b) (a+b)/2)
    grid=locs
    spcov.e<-matrix(NA,ncol = mysl*mytl/mynbatch,nrow=mysl*mytl/mynbatch)
    
    dist.sp.e<-rdist(grid[myloc.batch[,sr.batch],])
    locs.e<-grid[myloc.batch[,sr.batch],]
    sl.sub<-length(locs.e[,1])
    myspat.e<-myspat1[myloc.batch[,sr.batch],]
    z<-c(myspat.e)
    ##### Creating covariance matrix #####
    for(i in 1:mytl)
    {
      for(j in i:mytl)
      {
        spcov.e[((i-1)*sl.sub+1):((i-1)*sl.sub+sl.sub),((j-1)*sl.sub+1):((j-1)*sl.sub+sl.sub)]<-spcov.e[((j-1)*sl.sub+1):((j-1)*sl.sub+sl.sub),((i-1)*sl.sub+1):((i-1)*sl.sub+sl.sub)]<-myzeta.e[i,j]*my.matern(h=dist.sp.e,
                                                                                                                                                                                                                a=myfas.e[i,j],
                                                                                                                                                                                                                nu=mynus.e[i,j],
                                                                                                                                                                                                                sigma=sigma.est)
      }
    }
    
    
    
    
    
    C<-spcov.e
    nlogl<--dmvnorm(x=z,mean=(rep(0,times=length(z))),sigma = C,log = T)
    myloglv<-myloglv+nlogl
    }
    return(list(mlv=myloglv))
  }
}
#for(i in 1:10)
#print(mle.temporal.all(locs = grid,myspat1 =myspat,p=c(mya_bern,myalpha_bern,mybeta_bern),sr.batch = i))
#system.time(mle.temporal.all(locs = grid,myspat1 =myspat,p=c(mya_bern,myalpha_bern,mybeta_bern),sr.batch = i))


#### Finding initial value for temporal correlation ####
emp.temp.cor<-cor(train.myspat[,,simkey])
temp.ini<-function(p)
{
  mya<-10
  myalpha<-p[1]
  mybeta<-0
  mybeta2<-p[2]
  if(sum(c(mya<=0,myalpha<=0,myalpha>1,mybeta<0,mybeta>1,mybeta2<0))!=0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    t1<-zeta(t=tseq[1:length(train.myspat[1,,1])],nus = myest.nu[simkey,],alphas = myest.a[simkey,],a=mya,alpha=myalpha,beta=mybeta,beta2 =mybeta2,c.p = mean(myest.a[simkey,]) )
    return(mean(sum((t1-emp.temp.cor)^2)))
  }
}
#temp.ini(c(mya_bern,myalpha_bern,mybeta_bern))
set.seed(1)
ls.ini<-runif(2,min = 0.1,max = 0.9)
my.temp.ini<-optim(par=ls.ini,temp.ini)


mle.temp_only_mlv<-function(par)
{
  return(mle.temporal.all(locs=train.locs[,,simkey],myspat1=train.myspat[,,simkey],p=par,t.points=tseq[1:length(train.myspat[1,,1])],a.est=myest.a[simkey,],nu.est=myest.nu[simkey,],sigma.est=myest.sigma[simkey],mysl=train_l,mytl=length(train.myspat[1,,1]),mynbatch=nbatch,myloc.batch=loc.batch.ind[,,simkey])
         
  )
}


optim_marg_temp.loglik <- function(par){
  optim(par=par,
        fn = mle.temp_only_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=1000))
}

#simkey

tempporal.est<-optim_marg_temp.loglik(par=my.temp.ini$par)


end.time<-proc.time()
list(estimates=tempporal.est,time=end.time-start.time)

}


save.image("separabletemporal.RData")



bern.est.alpha<-bern.est.beta1<-bern.est.beta2<-numeric(length = nsim )
for(i in 1:nsim)
{
  bern.est.alpha[i]<-mytemporal[[i]]$estimates$par[1]
  bern.est.beta1[i]<-0
  bern.est.beta2[i]<-mytemporal[[i]]$estimates$par[2]
}

mean(bern.est.alpha)
mean(bern.est.beta1)
mean(bern.est.beta2)


sd(bern.est.alpha)
sd(bern.est.beta1)
sd(bern.est.beta2)


hist(bern.est.alpha)
hist(bern.est.beta1)
hist(bern.est.beta2)
median(bern.est.beta1)
median(bern.est.beta2)




myest.full.a<-cbind(myest.a,myest.a[,1:2])
myest.full.nu<-cbind(myest.nu,myest.nu[,1:2])

#full_cmean_a<-colMeans(myest.full.a)
#full_cmean_nu<-colMeans(myest.full.nu)

#full_sd_a<-apply(myest.full.a,2,sd)
#full_sd_nu<-apply(myest.full.nu,2,sd)

#plot(tseq,alphas,col="grey",type="l",ylim=range(c(myest.full.a)))
#points(tseq,alphas,pch=19)
#points(tseq,full_cmean_a,pch=19,col="blue")
#lines(tseq,full_cmean_a+1.96*full_sd_a,pch=19,col="blue",tly=2)
#lines(tseq,full_cmean_a-1.96*full_sd_a,pch=19,col="blue",tly=2)







#plot(tseq,nus,col="grey",type="l",ylim=range(c(nus,full_cmean_nu)))
#points(tseq,nus,pch=19)
#points(tseq,full_cmean_nu,pch=19,col="blue")
#lines(tseq,full_cmean_nu+1.96*full_sd_nu,pch=19,col="blue",tly=2)
#lines(tseq,full_cmean_nu-1.96*full_sd_nu,pch=19,col="blue",tly=2)







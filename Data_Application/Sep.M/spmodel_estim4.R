setwd("/ibex/scratch/qadirga/Project_4/Data_Application/SeparableModel")
#setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Data Application/SeparableModel")
load("gmrcl_ests.RData")
rcl.est
mygm_ini<-rcl.est$par
rm(list = ls()[-which(ls()=="mygm_ini")])

load("Myprelimests_oregon.RData")
#rm(list=ls()[-c(16,19,21,22,27,33,35)])
library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
ncores<-detectCores()
registerDoParallel(cores = ncores-1)
getDoParWorkers()
### functions ####

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*a)^nu
  num3<-besselK(x=(h*a), nu=nu)
  return(num1*num2*num3)
}

ls.a<-ls.nu<-ls.sigma<-numeric()

for(i in 1:355)
{
  ls.a[i]<-my_exp[[i]]$par[1]
  ls.nu[i]<-my_exp[[i]]$par[2]
  ls.sigma[i]<-my_exp[[i]]$par[3]
}

par(mfrow=c(1,3))
plot(1:355,ls.a,type="l")
plot(1:355,ls.nu,type="l")
plot(1:355,ls.sigma,type="l")




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
library(sp)
library(gstat)
library(geoR)
library(fields)
library(MASS)
library(mvtnorm)

############## stepwise MLE ##########

ini.sd<-sd(res.mat)

alpha.func<-function(t, params)
{
  temp<-numeric(length(t))
  for(i in 1:length(t))
  {
    temp2<-0
    for(j in 1:length(params))
    {
      
      temp2<-temp2+t[i]^(j-1)*params[j]
      
    }
    temp[i]<-exp(sum(temp2))
    
  }
  return(temp)
}


nu.func<-function(t, params)
{
  temp<-numeric(length(t))
  for(i in 1:length(t))
  {
    temp2<-0
    for(j in 1:length(params))
    {
      
      temp2<-temp2+t[i]^(j-1)*params[j]
      
    }
    temp[i]<-exp(sum(temp2))
    
  }
  return(temp)
}

############ Fitting least square model ###########
##### scaling the time to 0-1 ###

t.points<-1:365
tseq<-(t.points-(min(t.points)))/(max(t.points)-min(t.points)) ###scaled

tr.days<-355
alpha.poly.order<-nu.poly.order<-1

par(mfrow=c(1,2))
plot(tseq[1:tr.days],ls.a,ylab="alpha_s(t)",xlab="t (days)",main="Spatial scale",type="l")
a.ini<-mean(ls.a)
lines(tseq[1:tr.days],rep(a.ini,tdays),col="green")  
#plot(tseq[1:tr.days],log(alpha.func(t=tseq[1:tr.days],params = alpha.lmestim$coefficients)))

#min(log(alpha.func(t=1:355,params = alpha.lmestim$coefficients)))
#max(log(alpha.func(t=1:355,params = alpha.lmestim$coefficients)))
plot(tseq[1:tr.days],ls.nu,ylab="nu_s(t)",xlab="t (days)",main="Spatial smoothness",type="l")
nu.ini<-mean(ls.nu)
lines(tseq[1:tr.days],rep(nu.ini,tdays),col="green")  
quilt.plot(my.res.data[,1],my.res.data[,2],my.res.data[,100])
map("state",regions = 'oregon',add=T)
ini.v<-c(a.ini,nu.ini,sd(res.mat))

#names(ini.v)<-NULL
train.dist<-rdist(cbind(my.res.data[,1],my.res.data[,2]))
##### Now creating maximum likelihood function ####
#p<-ini.v
#dim(res.mat)
mle.comp.all<-function(tr.dist=train.dist,z=c(res.mat),p,tl=tdays,sl=tr.l,tseq.arg=tseq)
{
  
  
  a.s<-p[1]
  nu.s<-p[2]
  sigma<-p[3]
  
  #a.s<-exp(p[1]+p[2]*(1:355))
  #nu.s<-exp(p[4]+p[5]*(1:355))
  
  loc.n<-nrow(tr.dist)
  if(sum(c(sigma<=0,a.s<=0,nu.s<=0))!=0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    #myC<-array(dim = c(nrow(tr.dist),ncol(tr.dist),tl))
    #myC<-foreach(i = tseq.arg[1:tl]) %dopar%
    #{
    # my.matern(h=tr.dist,a=a.s[i],nu=nu.s[i],sigma = sigma)
    #}
    
    myC<-my.matern(h=tr.dist,a=a.s,nu=nu.s,sigma = sigma)
    daywise.nlogl<-foreach(i=1:tl,.packages = c("mvtnorm"))%dopar%
      {
        
        tspec.d<-z[((i-1)*loc.n+1):(i*loc.n)] 
        #myC<-my.matern(h=tr.dist,a=a.s[i],nu=nu.s[i],sigma = sigma)
        temp<--dmvnorm(x=tspec.d,mean=(rep(0,times=length(tspec.d))),
                       sigma =myC ,log = T,checkSymmetry = F)
        #rm(myC,tspec.d)
        
        temp
      }
    #rm(myC)
    nlogl<-do.call(sum,daywise.nlogl)
    cat("done computing nlog",p,"\n")
    return(list(mlv=nlogl))
  }
}



system.time(mle.comp.all(p=ini.v))


mle_only_mlv<-function(par)
{
  return(mle.comp.all(p=par)$mlv
  )
}

mle_only_mlv(ini.v)

optim_marg_loglik <- function(par){
  optim(par=par,
        fn = mle_only_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=2000))
  
}


#myspace.estims<-optim_marg_loglik(ini.v)

#############################################
#save.image("Oregan_sp_p01.RData")
#load("Oregan_sp_p01.RData")

######## Now reporting preliminary estimates ######

#ml1.a<-ml1.nu<-ml1.sigma<-numeric()

#ml1.a<-rep(myspace.estims$par[1],times=tdays)
#ml1.nu<-rep(myspace.estims$par[2],times=tdays)
#ml1.sigma<-myspace.estims$par[3]


#plot(tseq[1:tdays],ls.a,type="l",xlab="time",ylab="spatial scale")
#lines(tseq[1:tdays],ml1.a,col="blue",lwd=2)


#plot(tseq[1:tdays],ls.nu,type="l",xlab="time",ylab="spatial smoothness")
#lines(tseq[1:tdays],ml1.nu,col="blue",lwd=2)




########## Temporal correlation ###########

loc.batch<-list()
########## Now we do this 10 times ######
est.rep<-4  ##### No. of repeated subsampling
nbatch<-46
set.seed(123)
loc.batch.ind<-matrix(NA,nrow=tr.l/nbatch,ncol=est.rep*nbatch)

for(i in 1:est.rep)
{
  loc.order<-sample(1:tr.l,tr.l)
  loc.batch[[i]]<-matrix(loc.order,nrow =tr.l/nbatch )
}
loc.batch.ind<-do.call(cbind,loc.batch)



mle.temporal.all<-function(locs=cbind(my.res.data[,1],my.res.data[,2]),myspat1=res.mat,p,t.points=tseq[1:tdays],a.est=ml1.a,nu.est=ml1.nu,sigma.est=ml1.sigma,mysl=length(my.res.data[,1]),mytl=tdays,mynbatch=nbatch,myloc.batch=loc.batch.ind)
{
  #distmat<-rdist(locs)
  mya_bern1<-p[1]
  myalpha_bern1<-p[2]
  mybeta_bern1<-0
  mybeta2_bern1<-p[3]
  
  if(sum(c(mya_bern1<=0,mya_bern1<=0,myalpha_bern1<=0,myalpha_bern1>1,mybeta_bern1<0,mybeta_bern1>1,mybeta2_bern1<0))!=0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    
    
    myfas.e<-f_alphas(t=t.points,alphas = a.est,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,c.p = mean(a.est))
    
    myzeta.e<-zeta(t=t.points,nus = nu.est,alphas = a.est,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,beta2 =mybeta2_bern1,c.p = mean(a.est) )
    
    indic<-ncol(myloc.batch)
    #myloglv<-0
    full.dist<-rdist(locs)
    myloglv<-foreach(sr.batch = 1:indic,.combine = '+',.packages = c("mvtnorm"))%dopar%
      {
        mynus.e<-outer(nu.est,nu.est,function(a,b) (a+b)/2)
        grid=locs
        spcov.e<-matrix(NA,ncol = mysl*mytl/mynbatch,nrow=mysl*mytl/mynbatch)
        
        dist.sp.e<-full.dist[myloc.batch[,sr.batch],myloc.batch[,sr.batch]]
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
        #myloglv<-myloglv+nlogl
        #cat("done with ",sr.batch,"\n")
        nlogl
      }
    
    cat("done with ",p,"\n")
    return(list(mlv=myloglv))
  }
}

#temp.init<-c(10,0.9,0.5,0.5)



#system.time(mle.temporal.all(p=c(10,0.9,0.5)))


#for(i in 1:10)
#print(mle.temporal.all(locs = grid,myspat1 =myspat,p=c(mya_bern,myalpha_bern,mybeta_bern),sr.batch = i))
#system.time(mle.temporal.all(locs = grid,myspat1 =myspat,p=c(mya_bern,myalpha_bern,mybeta_bern),sr.batch = i))


#### Finding initial value for temporal correlation ####
emp.temp.cor<-cor(res.mat)
image.plot(tseq[1:tdays],tseq[1:tdays],emp.temp.cor)
temp.ini<-function(p)
{
  mya<-p[1]
  myalpha<-p[2]
  mybeta<-0
  mybeta2<-p[3]
  if(sum(c(mya<=0,myalpha<=0,myalpha>1,mybeta<0,mybeta>1,mybeta2<0))!=0)
  {
    return(Inf)
  }
  else
  {
    t1<-zeta(t=tseq[1:tdays],nus = ml1.nu,alphas = ml1.a,a=mya,alpha=myalpha,beta=mybeta,beta2 =mybeta2,c.p = mean(ml1.a) )
    return(mean(sum((t1-emp.temp.cor)^2)))
  }
}
#temp.ini(c(mya_bern,myalpha_bern,mybeta_bern))
set.seed(1)
ls.ini<-runif(3,min = 0.1,max = 0.9)
#my.temp.ini<-optim(par=ls.ini,temp.ini,control = list(maxit=1000,trace=6))

#initv.temp<-my.temp.ini$par


mle.temp_only_mlv<-function(par)
{
  return(mle.temporal.all(p=par)
         
  )
}


optim_marg_temp.loglik <- function(par){
  optim(par=par,
        fn = mle.temp_only_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=2000))
}

#simkey

#tempporal.est<-optim_marg_temp.loglik(par=initv.temp)

#save.image("sptemp_mlep1.RData")

#load("sptemp_mlep1.RData")


#########################################################
######### Random composite likelihood ###################
#########################################################


#### Plotting the temporal correlation so far #####

#adhoc_temp<-zeta(t=tseq[1:tdays],nus=rep(myspace.estims$par[2],tdays),alphas = rep(myspace.estims$par[1],tdays),a=tempporal.est$par[1],alpha =tempporal.est$par[2],beta = 0,
 #                beta2 = tempporal.est$par[3],c.p = mean(rep(myspace.estims$par[1],tdays)))

#par(mfrow=c(1,2))
#image.plot(tseq[1:tdays],tseq[1:tdays],emp.temp.cor,xlab="t",ylab="t",main="Empirical")
#image.plot(tseq[1:tdays],tseq[1:tdays],adhoc_temp,xlab="t",ylab="t",main="Empirical")


##############################################################

##### Now we define the random-composite loglikelihood function ######

start.time<-proc.time()
########## Using composite loglikelihood ###########
#my.p1.ini<-c(myspace.estims$par)


#my.p2.ini<-c(tempporal.est$par)

####### Composite MLE ######
#p=my.comp.ini2

##### Now creating maximum likelihood function ####


mle_rcl<-function(locs=cbind(my.res.data[,1],my.res.data[,2]),
                  myspat1=res.mat,p,t.points=tseq[1:tdays],mysl=length(my.res.data[,1]),
                  mytl=tdays,mynbatch=nbatch,
                  myloc.batch=loc.batch.ind,z=c(res.mat),tl=tdays,sl=tr.l,tseq.arg=tseq,tr.dist=train.dist)
{
  
  
  #distmat<-rdist(locs)
  a.0<-p[1]
  nu.0<-p[2]
  sigma<-p[3]
  mya_bern1<-p[4]
  myalpha_bern1<-p[5]
  mybeta_bern1<-0
  mybeta2_bern1<-p[6]
  
  if(sum(c(a.0<=0,nu.0<=0,sigma<=0,mya_bern1<=0,mya_bern1<=0,myalpha_bern1<=0,myalpha_bern1>1,mybeta_bern1<0,mybeta_bern1>1,mybeta2_bern1<0))!=0)
  {
    return(list(mlv=Inf))
  }
  
  else{
    loc.n<-nrow(tr.dist)
    a.s<-rep(a.0,times=tl)
    nu.s<-rep(nu.0,times=tl)
    
    myfC<-my.matern(h=tr.dist,a=a.0,nu=nu.0,sigma = sigma)
    daywise.nlogl<-foreach(i=1:tl,.packages = c("mvtnorm"))%dopar%
      {
        
        tspec.d<-z[((i-1)*loc.n+1):(i*loc.n)] 
        #myC<-my.matern(h=tr.dist,a=a.s[i],nu=nu.s[i],sigma = sigma)
        temp<--dmvnorm(x=tspec.d,mean=(rep(0,times=length(tspec.d))),
                       sigma =myfC ,log = T,checkSymmetry = F)
        #rm(myC,tspec.d)
        
        temp
      }
    #rm(myC)
    nlogl1<-do.call(sum,daywise.nlogl)
    cat("done computing nlogl part1",p,"\n")
    
    myfas.e<-f_alphas(t=t.points,alphas = a.s,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,c.p = mean(a.s))
    
    myzeta.e<-zeta(t=t.points,nus = nu.s,alphas = a.s,a=mya_bern1,alpha=myalpha_bern1,beta=mybeta_bern1,beta2 =mybeta2_bern1,c.p = mean(a.s) )
    
    indic<-ncol(myloc.batch)
    #myloglv<-0
    full.dist<-rdist(locs)
    myloglv<-foreach(sr.batch = 1:indic,.combine = '+',.packages = c("mvtnorm"))%dopar%
      {
        mynus.e<-outer(nu.s,nu.s,function(a,b) (a+b)/2)
        grid=locs
        spcov.e<-matrix(NA,ncol = mysl*mytl/mynbatch,nrow=mysl*mytl/mynbatch)
        
        dist.sp.e<-full.dist[myloc.batch[,sr.batch],myloc.batch[,sr.batch]]
        locs.e<-grid[myloc.batch[,sr.batch],]
        sl.sub<-length(locs.e[,1])
        myspat.e<-myspat1[myloc.batch[,sr.batch],]
        z2<-c(myspat.e)
        ##### Creating covariance matrix #####
        for(i in 1:mytl)
        {
          for(j in i:mytl)
          {
            spcov.e[((i-1)*sl.sub+1):((i-1)*sl.sub+sl.sub),((j-1)*sl.sub+1):((j-1)*sl.sub+sl.sub)]<-spcov.e[((j-1)*sl.sub+1):((j-1)*sl.sub+sl.sub),((i-1)*sl.sub+1):((i-1)*sl.sub+sl.sub)]<-myzeta.e[i,j]*my.matern(h=dist.sp.e,
                                                                                                                                                                                                                    a=myfas.e[i,j],
                                                                                                                                                                                                                    nu=mynus.e[i,j],
                                                                                                                                                                                                                    sigma=sigma)
          }
        }
        
        
        
        
        
        C<-spcov.e
        nlogl2<--dmvnorm(x=z2,mean=(rep(0,times=length(z2))),sigma = C,log = T)
        #myloglv<-myloglv+nlogl
        #cat("done with ",sr.batch,"\n")
        nlogl2
      }
    cat("done with nlogl2 ",p,"\n")
    return(list(mlv=(nlogl1+myloglv)/2,p1=nlogl1,p2=myloglv))
  }
  
} 

#tempporal.est
mysp_ini<-mygm_ini[-6]
mysp_ini[6]<-1
system.time(mle_rcl(p=mysp_ini))

#mle_rcl(p=c(myspace.estims$par,tempporal.est$par))


mle.rcl_only_mlv<-function(par)
{
  return(mle_rcl(p=par)$mlv
  )
}


#myrcl.ini.v<-c(my.p1.ini,my.p2.ini)   

system.time(mle.rcl_only_mlv(mysp_ini))

optim_rcl.loglik <- function(par){
  optim(par=par,
        fn = mle.rcl_only_mlv,
        hessian=FALSE,
        control=list(trace=6,
                     pgtol=0,
                     parscale=rep(0.1,length(par)),
                     maxit=3000))
}

#simkey

rcl.est<-optim_rcl.loglik(par=mysp_ini)

save.image("sprcl_ests_gmini4.RData")


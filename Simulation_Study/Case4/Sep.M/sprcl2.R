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


#save.image("separabletemporal.RData")
setwd("/ibex/scratch/qadirga/Project_4/Simulation_study_rev3/Case4_2/SeparableModel")

load("separabletemporal.RData")



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




mysep.rcl<-foreach(simkey=37:72,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm"))%dopar%
  {
    
    start.time<-proc.time()
    ########## Using composite loglikelihood ###########
    
    sp.p1.ini<-c(myest.a[simkey,1],myest.nu[simkey,1],myest.sigma[simkey])
    sp.p2.ini<-c(bern.est.alpha[simkey],bern.est.beta2[simkey])
    sp.fin.ini<-c(sp.p1.ini,sp.p2.ini)
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
    
    mle_sp_rcl<-function(locs,p,sl=train_l,mat.index=mat.index.tr,uniq.dist=uniq.dist.tr,tseq=tseq,
                         myspat1,mynbatch=nbatch,myloc.batch=loc.batch.ind[,,simkey])
    {
      
      a<-p[1]
      nu<-p[2]
      sigma<-p[3]
      mya_bern1<-10
      myalpha_bern1<-p[4]
      mybeta_bern1<-0
      mybeta2_bern1<-p[5]
      
      tl<-mytl<-length(myspat1[1,])
      mysl=sl
      t.points<-tseq[1:mytl]
      z<-c(myspat1)
      sigma.est<-sigma
      
      if(sum(c(mya_bern1<=0,myalpha_bern1<=0,myalpha_bern1>1,mybeta_bern1<0,mybeta_bern1>1,mybeta2_bern1<0,a<=0,nu<=0,sigma<=0))!=0)
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
        
        nloglp1<-nlogl
        a.est<-rep(a,times=tl)
        nu.est<-rep(nu,times=tl)
        sigma.est<-sigma
        
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
        
        
        return(list(mlv=(nloglp1+myloglv)/2))
        
        
        
        }
    }
    
    
    mle.rcl_sp_only_mlv<-function(par)
    {
      return(mle_sp_rcl(locs=train.locs[,,simkey],p=par,sl=train_l,mat.index=mat.index.tr,
                        myspat1=train.myspat[,,simkey],uniq.dist=uniq.dist.tr
                        ,tseq=tseq,mynbatch=nbatch,myloc.batch=loc.batch.ind[,,simkey])
      )
    }
    
    optim_rcl.sp.loglik <- function(par){
      optim(par=par,
            fn = mle.rcl_sp_only_mlv,
            hessian=FALSE,
            control=list(trace=6,
                         pgtol=0,
                         parscale=rep(0.1,length(par)),
                         maxit=5000))
    }
    
    #simkey
    
    rcl.est<-optim_rcl.sp.loglik(par=sp.fin.ini)
    
    
    
    end.time<-proc.time()
    
    list(estimates=rcl.est,time=end.time-start.time)
    
  
}

save.image("sprclp2.RData")




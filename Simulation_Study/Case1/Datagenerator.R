
######################################################################################
#### Creating function to generate tvarying spatio-temporal covariance function ######
######################################################################################
library(MASS)
library(mvtnorm)
library(fields)
my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*a)^nu
  num3<-besselK(x=(h*a), nu=nu)
  return(num1*num2*num3)
}


######## Parameters #######
nus.init<-0.5 ### Spatial smoothness
nus.max<-1.5    ### Temporal smoothness

#nut.init<-3   #### Temporal smoothness at spatial center
#nut.corn<-0.5  #### Temporal smoothness at spatial corner



as.init<-20 ### spatial range=0.2 (1/5)
as.max<-15  ### spatial range=0.5 (1/2)


mynus.f<-function(t_i,nus.init,nus.max)
{
  temp<-nus.init+(nus.max-nus.init)*sin(t_i*pi/time_limit) 
  return(temp)
}


as<-function(t_i,as.init,as.max)
{
  temp<-as.init+(as.max-as.init)*sin(t_i*pi/time_limit) 
  return(temp)
}




###### First we create a psi function #######
psi<-function(t,a=1,alpha=1,beta=1)
{
  temp<-(a*(t^alpha)+1)^beta
  return(temp)
}


tseq<-seq(0,10,length.out = 100)

time_limit<-1
time_points<-21

tseq<-seq(0,time_limit,length.out = time_points)

#### Plotting psi #####
plot(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.5),main=bquote(Psi~"(t)"),xlab="t",ylab=bquote(Psi~"(t)"))
lines(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.5))
lines(tseq,1/psi(t=tseq^2,a=1,alpha=1,beta = 0.6))
lines(tseq,1/psi(t=tseq^2,a=1.5,alpha=1,beta = 0.5))


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



alphas<-as(t_i=tseq,as.init = as.init,as.max = as.max)

nus<-mynus.f(t_i=tseq,nus.init=nus.init,nus.max=nus.max)
#nus<-rep(1,times=time_points)
#mean(alphas)
plot(tseq,alphas)
plot(tseq,nus,xlab="t",ylab=bquote(nu[s]~"(t)"),pch=19)
lines(tseq,nus,col="grey")


plot(tseq,alphas)
plot(tseq,alphas,xlab="t",ylab=bquote(alpha[s]~"(t)"),pch=19)
lines(tseq,alphas,col="grey")



myfas<-f_alphas(t=tseq,alphas = alphas,a=0.01,alpha=0.5,beta=1,c.p=mean(alphas))
image.plot(tseq,tseq,myfas)


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


#nus<-(1+2/(tseq+1))/2

plot(tseq,nus)


par(mfrow=c(1,3))
mya_bern<-10 #0.05
myalpha_bern<-0.6 #0.5
mybeta_bern<-0.8 #0.7
mybeta2_bern<-0.1
myfas<-f_alphas(t=tseq,alphas = alphas,a=mya_bern,alpha=myalpha_bern,beta=mybeta_bern,c.p=mean(alphas))
image.plot(tseq,tseq,myfas,xlab="t",ylab="t",main=bquote(f[as]~"("~t[i]~","~t[j]~")"))

myzeta<-zeta(t=tseq,nus = nus,alphas = alphas,a=mya_bern,alpha=myalpha_bern,beta=mybeta_bern,beta2 = mybeta2_bern,c.p=mean(alphas))
image.plot(tseq,tseq,myzeta,xlab="t",ylab="t",main=bquote(zeta~"("~t[i]~","~t[j]~")"),zlim=c(0,1))
contour(tseq,tseq,myzeta,add=T)
#plot(myzeta[5,5:21])
mynus<-outer(nus,nus,function(a,b) (a+b)/2)
#diag(myzeta)
####### Creating spatial grid #########
spres<-25
x<-y<-seq(0,1,length.out = spres)

grid<-expand.grid(x,y)
par(mfrow=c(1,1))
plot(grid,xlab="x",ylab="y")


##########################################################
##### Creating  spatio-temporal covariance function ######
##########################################################
sl<-length(grid[,1])
tl<-length(tseq)
spcov<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
dist.sp<-rdist(grid)
#system.time(for(i in 1:tl)
#{
 # for(j in 1:tl)
  #{
   # spcov[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-myzeta[i,j]*my.matern(h=dist.sp,
    #                                                                                    a=myfas[i,j],
     #                                                                                   nu=mynus[i,j],
      #                                                                                  sigma=1)
#  }
#}
#)



######### Optimized cov computing ##############
uniq.dist<-sort(unique(c(dist.sp)))

##### Unique indexing #####
mat.index<-NULL
for(i in 1:length(uniq.dist))
{
  mat.index[[i]]<-cbind(which(dist.sp==uniq.dist[i],arr.ind = T),i)
}
#cbind(which(dist.sp==uniq.dist[2],arr.ind = T),2)
mat.index<-do.call(rbind,mat.index)
spcov2<-matrix(NA,ncol = sl*tl,nrow=sl*tl)
temp<-matrix(NA,ncol = sl,nrow = sl)


system.time(for(i in 1:tl)
{
  for(j in i:tl)
  {
    tmp1<-myzeta[i,j]*my.matern(h=uniq.dist,
                                a=myfas[i,j],
                                nu=mynus[i,j],
                                sigma=1)
    temp[mat.index[,-3]]<-tmp1[mat.index[,3]]
    spcov2[((i-1)*sl+1):((i-1)*sl+sl),((j-1)*sl+1):((j-1)*sl+sl)]<-spcov2[((j-1)*sl+1):((j-1)*sl+sl),((i-1)*sl+1):((i-1)*sl+sl)]<-temp
  }
}
)

#check<-chol(spcov2)
#rm(check)
library(mvtnorm)
set.seed(123)
nsim<-100
#mysim<-mvrnorm(n=nsim,mu=rep(0,times=sl*tl),Sigma = spcov2)
mysim<-rmvnorm(n=nsim,sigma=spcov2,method="chol")
sim1<-mysim
dim(sim1)

locs<-grid
rm(spcov,spcov2)
par(mfrow=c(5,5),mar=c(2,2,1,1))
key=1
#tseq
for(i in 1:(length(tseq)))
{
  if(i<=(length(tseq)-1))
  {
    quilt.plot(locs,sim1[key,(1+length(locs[,1])*(i-1)):(length(locs[,1])*i)],main=paste("t=",round(tseq[i],2)),nx=spres,ny=spres,zlim=c(min(sim1[key,]),max(sim1[key,])),add.legend = F)
  }
  else
  {
    quilt.plot(locs,sim1[key,(1+length(locs[,1])*(i-1)):(length(locs[,1])*i)],main=paste("t=",round(tseq[i],2)),nx=spres,ny=spres,zlim=c(min(sim1[key,]),max(sim1[key,])))
    
  }
}


r.locs<-sample(1:sl,size = 4)
points(locs[r.locs,],pch="x",cex=1.5)
####### Now plotting the time series ########
for(i in r.locs)
{
  k=i
  plot(tseq,sim1[key,seq(k,(tl-1)*sl+k,by=sl)],type="l",main = paste("location"=round(locs[k,1],2),round(locs[k,2],2)),xlab="time",ylab = "value")
  points(tseq,sim1[key,seq(k,(tl-1)*sl+k,by=sl)],pch=19)
}
#20*400
####### distributing columnwise ########
myspat<-matrix(NA,nrow = sl,ncol = tl)
#emp.temp.cor<-array(NA,dim=c(tl,tl,100))
#for(k in 1:100){
#for(i in 1:tl)
#{
#  myspat[,i]<-sim1[key,((i-1)*sl+1):((i-1)*sl+sl)]
#}
#emp.temp.cor[,,k]<-cor(myspat)
#}
#rm(myspat)
setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/case1_2/FInal piece of codes")
save.image("Dataset.RData")


##### Now we divide the dataset into training and testing data #####

### 20 percent testing data ####

#### Last two frames held out for forecasting #####
#### 20 percent locations from the remaning frames are held out for testing (interpolation)
test_l<-length(grid[,1])*0.2
train_l<-length(grid[,1])*0.8
set.seed(123)
rmiss.matrix<-matrix(NA,nrow=test_l,ncol=nsim)
for(i in 1:nsim)
{
  rmiss.matrix[,i]<-sample(1:sl,size = test_l)
}

#### Now dividing test and train data set #####

full.myspat<-array(NA,dim = c(sl,tl,nsim))

for(i in 1:nsim)
{
  for(j in 1:tl)
  {
  full.myspat[,j,i]<-sim1[i,((j-1)*sl+1):((j-1)*sl+sl)]
}
}
par(mfrow=c(5,5))
for(i in 1:tl)
quilt.plot(grid,full.myspat[,i,1],nx=spres,ny=spres,zlim=c(min(full.myspat[,,1]),max(full.myspat[,,1])))

train.locs<-array(NA,dim = c(train_l,2,nsim))
test.locs<-array(NA,dim=c(test_l,2,nsim))
train.myspat<-array(NA,dim = c(train_l,tl-2,nsim))
test.myspat<-array(NA,dim = c(test_l,tl,nsim))
for(i in 1:nsim)
{
  train.locs[,,i]<-as.matrix(grid[-rmiss.matrix[,i],])
  test.locs[,,i]<-as.matrix(grid[rmiss.matrix[,i],])
  train.myspat[,,i]<-full.myspat[-rmiss.matrix[,i],(1:(tl-2)),i]
  test.myspat[,,i]<-full.myspat[rmiss.matrix[,i],,i]
  
}

key<-3
par(mfrow=c(5,5))
for(i in 1:(tl-2))
{quilt.plot(train.locs[,,key],train.myspat[,i,key],zlim=c(min(full.myspat[,,key]),max(full.myspat[,,key])),nx=spres,ny=spres)  
quilt.plot(test.locs[,,key],test.myspat[,i,key],zlim=c(min(full.myspat[,,key]),max(full.myspat[,,key])),nx=spres,ny=spres,add = T)  
}

save.image("simulated_data.RData")



########## Visualizing spatio-temporal fields ######
library(viridis)
pr.time<-character(length = length(tseq))
for(i in 1:length(tseq))
{
  pr.time[i]<-as.character(paste("t=",format(round(tseq[i], 2), nsmall = 2),sep = "" ))
}
vis.data<-data.frame(x=rep(grid$Var1,times=tl),y=rep(grid$Var2,times=tl),z=c(full.myspat[,,1]),time=rep(pr.time,each=sl))
ggplot(vis.data, aes(x=x, y=y,col=z))+geom_point(data = vis.data, aes(x, y, colour = z),pch=19, size=1)+ 
  facet_wrap(~time,ncol=7)+scale_color_viridis(bquote(Z),option = "D")+labs(x=bquote(s[1]),y=bquote(s[2]))+theme(axis.text = element_text(size=7))





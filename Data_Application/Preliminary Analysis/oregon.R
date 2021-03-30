#############################
#############################
library(maps)
library(fields)
CR_ENC <- c("minnesota", "iowa", "wisconsin","michigan")
#setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Data Application")
setwd("/ibex/scratch/qadirga/Project_4/Data_Application/WholeUS")

st.choose<-"oregon"
d.length<-numeric()
mystate.data<-list()
tmean<-numeric()
dir()
load("pmdata2017.RData")
for(i in 1:365)
{
  mystate.data[[i]]<-daily.st[[i]][daily.st[[i]]$state==st.choose,]
  d.length[i]<-length(mystate.data[[i]]$Longitude)
  tmean[i]<-mean(mystate.data[[i]]$pm25_daily_average)
}




quilt.plot(mystate.data[[5]][,c(1,2,3)])
quilt.plot(mystate.data[[5]][,c(1,2)],exp(mystate.data[[5]][,c(3)]))
map('state',regions = st.choose,add = T)
plot(1:365,tmean,type = "l")

######### Choosing 15 random days of the year #######




##### Detrending the data #####

#mean.st<-list()
#interc.st<-numeric()


#for(i in 1:365)
#{  
# mean.st[[i]]<-lm(pm25_daily_average~Longitude+Latitude, data =mystate.data[[i]] )
# interc.st[i]<-mean.st[[i]]$coefficients[1]
#}

#plot(1:365,interc.st,type = "l")


################################################
############## Making Data Matrix ##############
################################################


mydata<-matrix(NA,nrow = d.length[1],ncol = 365+2)
mydata[,1]<-c(mystate.data[[1]]$Longitude)
mydata[,2]<-c(mystate.data[[1]]$Latitude)

for(i in 1:365)
{
  mydata[,i+2]<-mystate.data[[i]]$pm25_daily_average
}

mydata<-as.data.frame(mydata)

colnames(mydata)
names(mydata)[1]<-"Longitude"
names(mydata)[2]<-"Latitude"

series.names<-c()
for(i in 1:365)
{
  series.names[i]<-paste("z",i,sep="")
}

for(i in 1:365)
{
  names(mydata)[i+2]<-series.names[i]
}

mydata

par(mfrow=c(3,3))
for(i in 1:20)
{
  quilt.plot(mydata$Longitude,mydata$Latitude,mydata[,i+2],xlab="Longitude",ylab="Latitude")
  map('state',regions = st.choose,add=T)
}


######### Plotting time series ##########
par(mfrow=c(1,1))
plot(1:365,mydata[1,-c(1,2)],col=c(rep("blue",times=31+28),rep("green",times=31+30+31),
                                   rep("yellow",times=30+31+31),rep("purple",times=30+31+30),rep("blue",times=31)),
     ylab="PM2.5")
lines(1:365,mydata[1,-c(1,2)])


plot(1:365,colMeans(mydata[,-c(1,2)]),col=c(rep("blue",times=31+28),rep("green",times=31+30+31),
                                            rep("yellow",times=30+31+31),rep("purple",times=30+31+30),rep("blue",times=31)),
     ylab="PM2.5",xlab="Day")
lines(1:365,colMeans(mydata[,-c(1,2)]))

##### Season-wise means ########
win.mean<-mean(colMeans(mydata[1,-c(1,2)])[c(1:(31+28),335:365)])
spr.mean<-mean(colMeans(mydata[1,-c(1,2)])[c(60:(59+31+30+31))])
sum.mean<-mean(colMeans(mydata[1,-c(1,2)])[c(152:(151+30+31+31))])
fall.mean<-mean(colMeans(mydata[1,-c(1,2)])[c(244:(243+30+31+30))])


lines(x=c(1:(31+28)),y=rep(win.mean,times=length(c(1:(31+28)))),col="blue",lwd=2)
lines(x=c(335:365),y=rep(win.mean,times=length(c(335:365))),col="blue",lwd=2)

lines(x=c(60:(59+31+30+31)),y=rep(spr.mean,times=length(c(60:(59+31+30+31)))),col="green",lwd=2)
lines(x=c(152:(151+30+31+31)),y=rep(sum.mean,times=length(c(152:(151+30+31+31)))),col="yellow",lwd=2)
lines(x=c(244:(243+30+31+30)),y=rep(fall.mean,times=length(c(244:(243+30+31+30)))),col="purple",lwd=2)

####################################################################################################
############## Now we divide the data into test and training set ###################################
####################################################################################################

ts.l<-round(0.6*d.length[1],0)
tr.l<-d.length[1]-ts.l

set.seed(123)

r.miss<-sample(1:d.length[1],size = ts.l)

mytrain.data<-mydata[-r.miss,1:357]
mytest.data1<-mydata[r.miss,1:357]
mytest.data2<-mydata[,c(1:2,358:367)]



######### Plotting time series (training) ##########
par(mfrow=c(1,1))
plot(1:355,mytrain.data[1,-c(1,2)],col=c(rep("blue",times=31+28),rep("green",times=31+30+31),
                                         rep("yellow",times=30+31+31),rep("purple",times=30+31+30),rep("blue",times=21)),
     ylab="PM2.5")
lines(1:355,mytrain.data[1,-c(1,2)])


plot(1:355,colMeans(mytrain.data[,-c(1,2)]),col=c(rep("blue",times=31+28),rep("green",times=31+30+31),
                                                  rep("yellow",times=30+31+31),rep("purple",times=30+31+30),rep("blue",times=21)),
     ylab="PM2.5",xlab="Day")
lines(1:355,colMeans(mytrain.data[,-c(1,2)]))

##### Season-wise means ########
win.mean<-mean(colMeans(mytrain.data[,-c(1,2)])[c(1:(31+28),335:355)])
spr.mean<-mean(colMeans(mytrain.data[,-c(1,2)])[c(60:(59+31+30+31))])
sum.mean<-mean(colMeans(mytrain.data[,-c(1,2)])[c(152:(151+30+31+31))])
fall.mean<-mean(colMeans(mytrain.data[,-c(1,2)])[c(244:(243+30+31+30))])


lines(x=c(1:(31+28)),y=rep(win.mean,times=length(c(1:(31+28)))),col="blue",lwd=2)
lines(x=c(335:355),y=rep(win.mean,times=length(c(335:355))),col="blue",lwd=2)

lines(x=c(60:(59+31+30+31)),y=rep(spr.mean,times=length(c(60:(59+31+30+31)))),col="green",lwd=2)
lines(x=c(152:(151+30+31+31)),y=rep(sum.mean,times=length(c(152:(151+30+31+31)))),col="yellow",lwd=2)
lines(x=c(244:(243+30+31+30)),y=rep(fall.mean,times=length(c(244:(243+30+31+30)))),col="purple",lwd=2)

linechar_data<-data.frame(Day=1:355,logPM2_5=c(colMeans(mytrain.data[,-c(1,2)])),Season=c(rep("Winter",times=length(1:(31+28))),
                                                                                          rep("Spring",times=length(60:(59+31+30+31))),
                                                                                          rep("Summer",times=length(152:(151+30+31+31))),
                                                                                          rep("Fall",times=length(244:(243+30+31+30))),
                                                                                          rep("Winter",times=length(335:355))    
                                                                                          ))


ggplot(data=linechar_data, aes(x=Day, y=logPM2_5,color=Season))  + geom_point()+labs(y=bquote(Y["*"](t)))+
  geom_line(data = linechar_data[linechar_data$Season %in% c("Spring","Fall","Summer"),],aes(x=Day,y=logPM2_5,col=Season))+
  geom_line(data = linechar_data[1:(31+28),],aes(x=Day,y=logPM2_5,col=Season))+
  geom_line(data = linechar_data[335:355,],aes(x=Day,y=logPM2_5,col=Season))+
  geom_line(data = data.frame(Day=1:(31+28),logPM2_5=win.mean,Season="Winter"),aes(x=Day,y=logPM2_5,col=Season),linetype="solid",size=1.5)+
  geom_line(data = data.frame(Day=335:355,logPM2_5=win.mean,Season="Winter"),aes(x=Day,y=logPM2_5,col=Season),linetype="solid",size=1.5)+
  geom_line(data = data.frame(Day=60:(59+31+30+31),logPM2_5=spr.mean,Season="Spring"),aes(x=Day,y=logPM2_5,col=Season),linetype="solid",size=1.5)+
  geom_line(data = data.frame(Day=152:(151+30+31+31),logPM2_5=sum.mean,Season="Summer"),aes(x=Day,y=logPM2_5,col=Season),linetype="solid",size=1.5)+
  geom_line(data = data.frame(Day=244:(243+30+31+30),logPM2_5=fall.mean,Season="Fall"),aes(x=Day,y=logPM2_5,col=Season),linetype="solid",size=1.5)


  #geom_line(data = linechar_data[335:355,],aes(x=Day,y=logPM2_5,col=Season),linet)+
  #geom_line(data = linechar_data[335:355,],aes(x=Day,y=logPM2_5,col=Season))+
  #geom_line(data = linechar_data[335:355,],aes(x=Day,y=logPM2_5,col=Season))+
  #geom_line(data = linechar_data[335:355,],aes(x=Day,y=logPM2_5,col=Season))





########################################################################################################
####### Now we do the linear regression with harmonic temporal terms and linear spatial terms ##########
########################################################################################################

t.points<-1:355
t.points<-(t.points-1)/(365-1)
t1.h<-sin(2*pi*t.points/0.5)
plot(t.points,t1.h,type = "l")

t2.h<-cos(2*pi*t.points/0.5)
plot(t.points,t2.h,type = "l")

t3.h<-sin(2*2*pi*t.points/0.25)
plot(t.points,t3.h,type = "l")

t4.h<-cos(2*2*pi*t.points/0.25)
plot(t.points,t4.h,type = "l")

####################################
######## Preparing dataset #########
####################################

#long_pad<-lat_pad<-rep(0,times=tr.l*355)

#long_mat<-lat_mat<-matrix(NA,nrow=tr.l*355,ncol=355)
#for(i in 1:355)
#{
 # temp<-long_pad
#  temp[((i-1)*tr.l+1):(i*tr.l)]<-mytrain.data$Longitude
#  long_mat[,i]<-temp
#  temp<-lat_pad
#  temp[((i-1)*tr.l+1):(i*tr.l)]<-mytrain.data$Latitude
#  lat_mat[,i]<-temp
  
#}

#long_mat<-as.data.frame(long_mat)
#lat_mat<-as.data.frame(lat_mat)
#for(i in 1:355)
#{
 # names(long_mat)[i]<-paste("long_d",i,sep = "")
#  names(lat_mat)[i]<-paste("lat_d",i,sep = "")
  
#}
#loc_cov<-cbind(long_mat,lat_mat)
sin1t<-rep(t1.h,each=tr.l)
cos1t<-rep(t2.h,each=tr.l)
sin2t<-rep(t3.h,each=tr.l)
cos2t<-rep(t4.h,each=tr.l)
yvals<-c(as.matrix(mytrain.data[,-c(1,2)]))
mytime<-rep(t.points,each=tr.l)
############################
###### New covariates ######
############################

#st.time<-(t.points-1)/(365-1)
long_cov0<-rep(mytrain.data$Longitude,times=355)
lat_cov0<-rep(mytrain.data$Latitude,times=355)
long_cov1<-long_cov0*mytime
lat_cov1<-lat_cov0*mytime

long_cov2<-long_cov0*(mytime^2)
lat_cov2<-lat_cov0*(mytime^2)

long_cov3<-long_cov0*(mytime^3)
lat_cov3<-lat_cov0*(mytime^3)

long_cov4<-long_cov0*(mytime^4)
lat_cov4<-lat_cov0*(mytime^4)

long_cov5<-long_cov0*(mytime^5)
lat_cov5<-lat_cov0*(mytime^5)

long_cov6<-long_cov0*(mytime^6)
lat_cov6<-lat_cov0*(mytime^6)


long_cov7<-long_cov0*(mytime^7)
lat_cov7<-lat_cov0*(mytime^7)




long_cov8<-long_cov0*(mytime^8)
lat_cov8<-lat_cov0*(mytime^8)


long_cov9<-long_cov0*(mytime^8)
lat_cov9<-lat_cov0*(mytime^8)

long_cov10<-long_cov0*(mytime^10)
lat_cov10<-lat_cov0*(mytime^10)


long_cov11<-long_cov0*(mytime^11)
lat_cov11<-lat_cov0*(mytime^11)

long_cov12<-long_cov0*(mytime^11)
lat_cov12<-lat_cov0*(mytime^11)


long_cov13<-long_cov0*(mytime^13)
lat_cov13<-lat_cov0*(mytime^13)










#tim_cov<-data.frame(t=rep(t.points,each=tr.l),sin1t=sin1t,sin2t=sin2t,cos1t=cos1t,cos2t=cos2t,logpm25=yvals)


######### Check ##########
#b0<-b1<-b2<-numeric()

#for(i in 1:355)
#{
#mytry<-lm(mytrain.data[,i+2]~mytrain.data[,1]+mytrain.data[,2],data = mytrain.data)
#b0[i]<-mytry$coefficients[1]
#b1[i]<-mytry$coefficients[2]
#b2[i]<-mytry$coefficients[3]
#}
#plot(t.points,b0,type = "l")
#plot(t.points,b1,type = "l")
#plot(t.points,b2,type = "l")

lm.mat<-data.frame(logpm25=yvals,t=rep(t.points,each=tr.l),sin1t=sin1t,sin2t=sin2t,cos1t=cos1t,cos2t=cos2t,long0=long_cov0,lat0=lat_cov0,
                   long1=long_cov1,lat1=lat_cov1,long2=long_cov2,lat2=lat_cov2,long3=long_cov3,lat3=lat_cov3,
                   long4=long_cov4,lat4=lat_cov4
                   )
#lm.mat<-cbind(loc_cov,tim_cov)

my.ns.fit1<-lm(logpm25~.,data = lm.mat)
ns.summary<-summary(my.ns.fit1)
ns.summary


hist_res.data<-data.frame(Residual=ns.summary$residuals)
ggplot(hist_res.data, aes(x=Residual))+
  geom_histogram(color="darkblue", fill="lightblue")+labs(x=bquote(epsilon[log(PM[2.5])]),y="Frequency")
###### Histogram of raw data ######
hist_raw.data<-data.frame(PM25=exp(yvals))
ggplot(hist_raw.data, aes(x=PM25))+
  geom_histogram(color="darkblue", fill="lightblue")+labs(x=bquote(PM[2.5]),y="Frequency")

hist_log.data<-data.frame(PM25=(yvals))
ggplot(hist_log.data, aes(x=PM25))+
  geom_histogram(color="darkblue", fill="lightblue")+labs(x=bquote(log(PM[2.5])),y="Frequency")


#temp.r<-which(ns.summary$coefficients[,4]>0.05,arr.ind = T)
#temp.r2<-temp.r-1
#lm.mat2<-lm.mat[,-temp.r2]
#rm(daily.st,lat_mat,lm.mat3,loc_cov,long_mat,mean.st,tim_cov)
#my.ns.fit2<-lm(logpm25~.,data = lm.mat2)
sin1ts<-rep(t1.h,each=ts.l)
cos1ts<-rep(t2.h,each=ts.l)
sin2ts<-rep(t3.h,each=ts.l)
cos2ts<-rep(t4.h,each=ts.l)
mytime.s<-rep(t.points,each=ts.l)
############################
###### New covariates ######
############################

#st.time<-(t.points-1)/(365-1)
long_cov0s<-rep(mytest.data1$Longitude,times=355)
lat_cov0s<-rep(mytest.data1$Latitude,times=355)
long_cov1s<-long_cov0s*mytime.s
lat_cov1s<-lat_cov0s*mytime.s

long_cov2s<-long_cov0s*(mytime.s^2)
lat_cov2s<-lat_cov0s*(mytime.s^2)

long_cov3s<-long_cov0s*(mytime.s^3)
lat_cov3s<-lat_cov0s*(mytime.s^3)

long_cov4s<-long_cov0s*(mytime.s^4)
lat_cov4s<-lat_cov0s*(mytime.s^4)

long_cov5s<-long_cov0s*(mytime.s^5)
lat_cov5s<-lat_cov0s*(mytime.s^5)

long_cov6s<-long_cov0s*(mytime.s^6)
lat_cov6s<-lat_cov0s*(mytime.s^6)


long_cov7s<-long_cov0s*(mytime.s^7)
lat_cov7s<-lat_cov0s*(mytime.s^7)




long_cov8s<-long_cov0s*(mytime.s^8)
lat_cov8s<-lat_cov0s*(mytime.s^8)


long_cov9s<-long_cov0s*(mytime.s^8)
lat_cov9s<-lat_cov0s*(mytime.s^8)

long_cov10s<-long_cov0s*(mytime.s^10)
lat_cov10s<-lat_cov0s*(mytime.s^10)


long_cov11s<-long_cov0s*(mytime.s^11)
lat_cov11s<-lat_cov0s*(mytime.s^11)

long_cov12s<-long_cov0s*(mytime.s^11)
lat_cov12s<-lat_cov0s*(mytime.s^11)


long_cov13s<-long_cov0s*(mytime.s^13)
lat_cov13s<-lat_cov0s*(mytime.s^13)


test1cov<-data.frame(t=rep(t.points,each=ts.l),sin1t=sin1ts,sin2t=sin2ts,cos1t=cos1ts,cos2t=cos2ts,long0=long_cov0s,lat0=lat_cov0s,
                     long1=long_cov1s,lat1=lat_cov1s,long2=long_cov2s,lat2=lat_cov2s,long3=long_cov3s,lat3=lat_cov3s,
                     long4=long_cov4s,lat4=lat_cov4s
)

test1pred<-c(predict(my.ns.fit1,newdata =test1cov))
test1predmat<-matrix(test1pred,ncol=355)
quilt.plot(mytest.data1$Longitude,mytest.data1$Latitude,test1predmat[,300])

mytest1.res<-mytest.data1[,-c(1,2)]-test1predmat
quilt.plot(mytest.data1$Longitude,mytest.data1$Latitude,mytest1.res[,355])

my.orig.testdata1<-mytest.data1
mytest.data1[,-c(1,2)]<-mytest1.res ##### Now mytest.data1 contains validation residuals for interpolation
quilt.plot(mytest.data1$Longitude,mytest.data1$Latitude,mytest.data1[,355+2])

#summary(my.ns.fit2)
#################################################################
############# Rebuilding residual matrix ########################
#################################################################

my.res.data<-matrix(NA,nrow=tr.l,ncol = 355+2)

#load("clean.data.2.RData")
#my.res.data[,1]<-mystate.data[[1]]$Longitude
#my.res.data[,2]<-mystate.data[[1]]$Latitude
my.res.data[,1]<-mytrain.data$Longitude
my.res.data[,2]<-mytrain.data$Latitude

res.mat<-matrix(my.ns.fit1$residuals,nrow = tr.l,ncol=355,byrow = F)

my.res.data[,3:357]<-res.mat
par(mfrow=c(3,3))
for(i in 1:27)
{
  quilt.plot(my.res.data[,c(1,2,i+2)])
}

mean(res.mat)
sd(res.mat)
######### Now computing residuals for forecasting ###########

f.points<-356:365
f.points<-(f.points-1)/(365-1)
f1.h<-sin(2*pi*f.points/0.5)
plot(f.points,f1.h,type = "l")

f2.h<-cos(2*pi*f.points/0.5)
plot(f.points,f2.h,type = "l")

f3.h<-sin(2*2*pi*f.points/0.25)
plot(f.points,f3.h,type = "l")

f4.h<-cos(2*2*pi*f.points/0.25)
plot(f.points,f4.h,type = "l")


tf.l<-length(mytest.data2$Longitude)
sin1tf<-rep(f1.h,each=tf.l)
cos1tf<-rep(f2.h,each=tf.l)
sin2tf<-rep(f3.h,each=tf.l)
cos2tf<-rep(f4.h,each=tf.l)
mytime.f<-rep(f.points,each=tf.l)
############################
###### New covariates ######
############################

#st.time<-(t.points-1)/(365-1)
long_cov0f<-rep(mytest.data2$Longitude,times=10)
lat_cov0f<-rep(mytest.data2$Latitude,times=10)
long_cov1f<-long_cov0f*mytime.f
lat_cov1f<-lat_cov0f*mytime.f

long_cov2f<-long_cov0f*(mytime.f^2)
lat_cov2f<-lat_cov0f*(mytime.f^2)

long_cov3f<-long_cov0f*(mytime.f^3)
lat_cov3f<-lat_cov0f*(mytime.f^3)

long_cov4f<-long_cov0f*(mytime.f^4)
lat_cov4f<-lat_cov0f*(mytime.f^4)


test2cov<-data.frame(t=mytime.f,sin1t=sin1tf,sin2t=sin2tf,cos1t=cos1tf,cos2t=cos2tf,long0=long_cov0f,lat0=lat_cov0f,
                     long1=long_cov1f,lat1=lat_cov1f,long2=long_cov2f,lat2=lat_cov2f,long3=long_cov3f,lat3=lat_cov3f,
                     long4=long_cov4f,lat4=lat_cov4f
)

test2pred<-c(predict(my.ns.fit1,newdata =test2cov))
test2predmat<-matrix(test2pred,ncol=10)
dim(test2predmat)
quilt.plot(mytest.data2$Longitude,mytest.data2$Latitude,test2predmat[,5])

mytest2.res<-mytest.data2[,-c(1,2)]-test2predmat
quilt.plot(mytest.data2$Longitude,mytest.data2$Latitude,mytest2.res[,4])

my.orig.testdata2<-mytest.data2
mytest.data2[,-c(1,2)]<-mytest2.res ##### Now mytest.data2 contains validation residuals for interpolation
quilt.plot(mytest.data2$Longitude,mytest.data2$Latitude,mytest.data2[,4+2])






imp.cont<-c("tf.l","tr.l","ts.l","res.mat","my.res.data","my.ns.fit1","ns.summary","mytest.data1","mytest.data2","tseq","st.choose","my.orig.testdata2","my.orig.testdata1",
            "test1predmat","test2predmat")

#no.delete<-which(ls() %in% imp.cont)

#rmlist<-ls()[-which(ls() %in% imp.cont)]
rm(list = ls()[-which(ls() %in% imp.cont)])

hist(my.ns.fit1$residuals)


myfile.name<-paste("Mylmresid_",st.choose,".RData",sep="")
save.image(paste(myfile.name))
########################################################################################
################ Now we do the estimation ##############################################
########################################################################################

library(MASS)
library(mvtnorm)
library(fields)
library(doParallel)
library(geoR)
ncores<-detectCores()
registerDoParallel(cores = ncores-4)
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


mle.all<-function(locs,z,p)
{
  distmat<-rdist(locs)
  a<-p[1]
  nu<-p[2]
  sigma<-p[3]
  
  if(sum(c(p[1:3]<=0))!=0)
  {
    return(list(mlv=Inf))
  }
  else
  {
    C<-my.matern(h=distmat,a=a,nu=nu,sigma = sigma)
    nlogl<--dmvnorm(x=z,mean=(rep(0,times=length(z))),sigma = C,log = T)
    return(list(mlv=nlogl))
  }
}

#for(i in 1:365)
#{
# daily.st[[i]]$CR<-NA
#  State<-daily.st[[i]]$state
# daily.st[[i]]$CR[State %in% CR_NW] <- "NW"
#  daily.st[[i]]$CR[State %in% CR_W]  <- "W"
#  daily.st[[i]]$CR[State %in% CR_SW] <- "SW"
#  daily.st[[i]]$CR[State %in% CR_WNC]<- "WNC"
#  daily.st[[i]]$CR[State %in% CR_ENC]<- "ENC"
#  daily.st[[i]]$CR[State %in% CR_S]  <- "S"
#  daily.st[[i]]$CR[State %in% CR_C]  <- "C"
#  daily.st[[i]]$CR[State %in% CR_NE] <- "NE"
#  daily.st[[i]]$CR[State %in% CR_SE] <- "SE"
#  daily.st[[i]]<-na.omit(daily.st[[i]])
#}
set.seed(123)
tdays<-355
my_exp<-foreach(day=1:tdays,.packages = c("sp","gstat","geoR","fields","MASS","mvtnorm"))%dopar%
  {
    set.seed(day)
    mindata.length<-800
    
    ls.a<-ls.nu<-ls.sigma<-numeric(length = tdays)
    #tdata<-daily.st[[day]][daily.st[[day]]$CR==rchoice,c(1,2,3)]
    mylocs<-cbind(my.res.data[,1],my.res.data[,2])
    myz<-res.mat[,day]
    
    data.length<-length(myz)
    
    speci_length<-min(mindata.length,data.length)
    set.seed(123)
    r.ind<-sample(1:data.length,size = speci_length)
    
    mylocs.sub<-mylocs[r.ind,]
    myz.sub<-myz[r.ind]
    #myz.sub<-(myz.sub-mean(myz.sub))/sd(myz.sub)
    #quilt.plot(tdata)
    #quilt.plot(mylocs.sub,myz.sub)
    
    
    mytrial<-cbind(mylocs.sub,myz.sub)
    mytrial<-as.geodata(mytrial)
    myvario<-variog(mytrial)
    
    
    my.var.loss<-function(p,vardist=myvario$u,gamma_hat=myvario$v,nbins=myvario$n)
    {
      a<-p[1]
      nu<-p[2]
      sigma<-p[3]
      if(sum(p<=0,p[2]>2,p[1]>30)!=0)
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
    
    temp<-optim(my.var.loss,par = c(runif(n=1,min=1,max=29),0.5,sd(myz.sub)))
    ls.a[day]<-temp$par[1]
    ls.nu[day]<-temp$par[2]
    ls.sigma[day]<-temp$par[3]
    
    ini.v<-c(ls.a[day],ls.nu[day],ls.sigma[day])
    
    mle_only_mlv<-function(par)
    {
      return(mle.all(locs=mylocs.sub,z=myz.sub,p=par)$mlv
      )
    }
    check.value<-mle_only_mlv(ini.v)
    if(check.value==Inf)
    {
      ini.v<-c(1,1,1)
    }
    optim_marg_loglik <- function(par){
      optim(par=par,
            fn = mle_only_mlv,
            hessian=FALSE,
            control=list(trace=6,
                         pgtol=0,
                         parscale=rep(0.1,length(par)),
                         maxit=1000))
    }
    
    
    
    optim_marg_loglik(ini.v)
  }




myfile.name2<-paste("Myprelimests_",st.choose,".RData",sep="")
save.image(paste(myfile.name2))

####### Visualization ######
my_exp[[180]]
ls.a<-ls.nu<-ls.sigma<-numeric(length = tdays)

for(i in 1:tdays)
{
  ls.a[i]<-my_exp[[i]]$par[1]
  ls.nu[i]<-my_exp[[i]]$par[2]
  ls.sigma[i]<-my_exp[[i]]$par[3]
}

par(mfrow=c(1,3))
plot(1:tdays,ls.a,xlab="time",ylab="scale",main="spatial scale",type="l")
plot(1:tdays,ls.nu,xlab="time",ylab="smoothness",main="Smoothness",type="l")
plot(1:tdays,ls.sigma,xlab="time",ylab="sigma",main="Standard deviation",type="l")


alpha.data<-data.frame(Day=1:tdays,alpha=ls.a)
nu.data<-data.frame(Day=1:tdays,nu=ls.nu)
order(ls.a)
which.max(ls.a)
ggplot(data=alpha.data[-c(180,218),], aes(x=Day, y=alpha))  + geom_point(color="lightblue")+labs(y=bquote(alpha), x="Day")+
  geom_line(color="darkblue")
ls.a[218]
ggplot(data=nu.data, aes(x=Day, y=nu))  + geom_point(color="lightblue")+labs(y=bquote(nu), x="Day")+
  geom_line(color="darkblue")

empcor<-cor(res.mat)
library(fields)
image.plot(tseq[1:tdays],tseq[1:tdays],empcor)




################################################
######## Plots for the Manuscript ##############
################################################

library(ggplot2)
library(viridis)
library(mapdata)
us_state_map = map_data("state")
us_oregon_map<-us_state_map[us_state_map$region=="oregon",]
set.seed(123)
ndays<-18
vis.day<-sample(1:365,ndays)
vis.day<-sort(vis.day)
vis.day


my_vis_data<-data()
vis_charday<-character()
for(i in 1:ndays)
{
  vis_charday[i]<-format(as.Date(vis.day[i], origin = "2016-12-31"),"%d-%m-%Y")
}
myvis_matrix<-matrix(NA,nrow=length(mystate.data[[1]][,1])*ndays,ncol=3)

for(i in 1:ndays)
{
  myvis_matrix[((i-1)*length(mystate.data[[1]][,1])+1):(i*length(mystate.data[[1]][,1])),c(1,2,3)]<-cbind(c(mystate.data[[vis.day[i]]][,c(1)]),c(mystate.data[[vis.day[i]]][,c(2)]),c(mystate.data[[vis.day[i]]][,c(3)]))
}
myvis_data<-data.frame(Longitude=myvis_matrix[,1],Latitude=myvis_matrix[,2],PM2_5=(myvis_matrix[,3]),Day=rep(vis_charday,each=length(mystate.data[[1]][,1])))
myvis_data$Day<-factor(myvis_data$Day,levels = vis_charday)
ggplot(myvis_data, aes(x=Longitude, y=Latitude,col=PM2_5))+borders("state",regions = "oregon")+geom_point(data = myvis_data, aes(Longitude, Latitude, colour = PM2_5),pch=19, size=1)+ 
  facet_wrap(~Day,ncol=6)+scale_color_viridis(bquote(log(PM[2.5])),option = "D")+labs(x="Longitude",y="Latitude")+theme(axis.text = element_text(size=9))




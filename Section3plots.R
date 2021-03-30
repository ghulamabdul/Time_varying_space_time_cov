########################################################
##### Visualization of temporal nonstationarity ########
########################################################

psi_func<-function(u,a=1,alpha=1,beta=1)
{
  return((a*(u^alpha)+1)^beta)
}


u1<-seq(0,1,length.out = 1000)
plot(u1,1/psi_func(u=u1^2,a=3,alpha = 0.5),ylim=c(0,1))




alphab<-10 #### constant alpha
alphati<-20 ### varying alpha

alphab_r<-rep(alphab,times=length(u1))
alphati_r<-rep(alphati,times=length(u1))
alphab_r[1]<-alphati

my_tempcov<-function(alphab=15,alphati=15,vuti=1,vutj=1,sig=1)
{
  alphab_r<-rep(alphab,times=length(u1))
  alphati_r<-rep(alphati,times=length(u1))
  alphab_r[1]<-alphati
  mynum<-(sig^2)*gamma((vuti+vutj)/2)
  myden1<-sqrt(gamma(vuti)*gamma(vutj))*((alphab_r^2)^(1/2))*((alphati_r^2)^(1/2))
  myden2<-(psi_func(u=u1^2,a=3,alpha = 0.5)/(alphab^2)-(1/(alphab^2))+0.5*((1/(alphab_r^2))+(1/(alphati_r^2))))
  myvalue<-mynum/(myden1*myden2)
  return(myvalue)
}

mytv<-my_tempcov()
lines(u1,mytv,col="red")

p1.baseline<-my_tempcov()
p1.case1<-my_tempcov(alphati=5)
p1.case2<-my_tempcov(alphati=10)
p1.case3<-my_tempcov(alphati=20)
p1.case4<-my_tempcov(alphati=25)


plot(u1,p1.baseline)
lines(u1,p1.case1,col=1)
lines(u1,p1.case2,col=2)
lines(u1,p1.case3,col=3)
lines(u1,p1.case4,col=4)


p1.data.frame<-data.frame(timelag=rep(u1, times=5),cov=c(p1.baseline,p1.case1,p1.case2,p1.case3,p1.case4),
                          group=rep(c("alphas=15,alphati=15","alphas=15,alphati=10","alphas=15,alphati=12",
                                      "alphas=15,alphati=18","alphas=15,alphati=20"), each=length(u1)))


rep(bquote(alpha[f]),times=10)
p1.data.frame<-p1.data.frame[p1.data.frame$timelag!=0,]
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(tidyverse)
theme_set(theme_grey())

#display.brewer.all(colorblindFriendly = TRUE)

ggplot(data = p1.data.frame,aes(x=timelag,y=cov,group=group,color=group))+
  geom_line()+labs(x=bquote("|"~t[r]-t[j]~"|"),y=bquote(C(bold("0"),t[r],t[j])),color="")+theme(axis.text = element_text(size=9),legend.title = element_blank(),legend.background = element_rect(fill = "transparent", colour = "transparent"),legend.position = c(0.75, 0.75),legend.direction="vertical")+
  scale_color_discrete()


ggplot(data = p1.data.frame,aes(x=timelag,y=cov,group=group,color=group))+
geom_line()+labs(x=bquote("|"~t[r]-t[j]~"|"),y=bquote(C(bold("0"),t[r],t[j])),color="")+theme(legend.text = element_text( size = 10),axis.text = element_text(size=9),legend.title = element_blank(),legend.background = element_rect(fill = "transparent", colour = "transparent"),legend.position = c(0.75, 0.75),legend.direction="vertical")+
  scale_color_viridis(labels=c(bquote(alpha[f]~"=15,"~alpha[r]~"=10"),
                              bquote(alpha[f]~"=15,"~alpha[r]~"=12"),
                              bquote(alpha[f]~"=15,"~alpha[r]~"=15"),
                              bquote(alpha[f]~"=15,"~alpha[r]~"=18"),
                              bquote(alpha[f]~"=15,"~alpha[r]~"=20")),discrete = T)




p2.baseline<-my_tempcov()
p2.case1<-my_tempcov(vutj=1.5)
p2.case2<-my_tempcov(vutj=2.0)
p2.case3<-my_tempcov(vutj=2.5)
p2.case4<-my_tempcov(vutj=3)




plot(u1,p2.baseline)
lines(u1,p2.case1,col=1)
lines(u1,p2.case2,col=2)
lines(u1,p2.case3,col=3)
lines(u1,p2.case4,col=4)


p2.data.frame<-data.frame(timelag=rep(u1, times=5),cov=c(p2.baseline,p2.case1,p2.case2,p2.case3,p2.case4),
                          group=rep(c("1","1.5","2.0",
                                      "2.5","3"), each=length(u1)))

p2.data.frame<-p2.data.frame[p2.data.frame$timelag!=0,]


ggplot(data = p2.data.frame,aes(x=timelag,y=cov,group=group,color=group))+
  geom_line()+labs(x=bquote("|"~t[r]-t[j]~"|"),y=bquote(C(bold("0"),t[r],t[j])),color="")+theme(legend.text = element_text( size = 10),axis.text = element_text(size=9),legend.title = element_blank(),legend.background = element_rect(fill = "transparent", colour = "transparent"),legend.position = c(0.75, 0.75),legend.direction="vertical")+
  scale_color_viridis(labels=c(bquote(nu[f]~"=1.0,"~nu[r]~"=1.0"),
                               bquote(nu[f]~"=1.0,"~nu[r]~"=1.5"),
                               bquote(nu[f]~"=1.0,"~nu[r]~"=2.0"),
                               bquote(nu[f]~"=1.0,"~nu[r]~"=2.5"),
                               bquote(nu[f]~"=1.0,"~nu[r]~"=3.0")),discrete = T)





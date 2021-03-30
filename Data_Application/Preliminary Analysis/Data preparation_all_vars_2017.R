##### Preparing a trivariate data for PM2.5, WS1, WS2 #####
####################################################
############ Loading required libraries ############
####################################################
library(sp)
library(maps)
library(maptools)
library(geosphere)
library(fields)
library(MASS)
#library(scoringRules)
#library(doParallel)
library(rgdal)
#registerDoParallel(cores = 5)
#########################################################
################ Data Creation portion ##################
#########################################################
#setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Exploratory data anaysis")
#setwd("/Users/qadirga/Documents/Project 2/Rcode/Data Analysis")
#setwd("/Users/qadirga/Documents/Project 3/Data Analysis (WS and PM 2.5)/Data preparation")
#setwd("/Users/qadirga/Documents/Project 3/Manuscript/Data Analysis/New day checking")
#setwd("D:/Projects 3/Data Analysis version 2/Data Analysis-gpawspm2.5/Summer data set")
#setwd("/Users/qadirga/Documents/Project 3/Manuscript/Data Analysis Final Version")
setwd("/Users/qadirga/Documents/Project 4-Spacetime modeling/Revisit/Simulation Study/Revisit2/Revisit 3/Data Application/2017 data analysis")
dir()
#------ PM25 Data ------#

time_pm25_day_read <- system.time(pm25_day_2017 <- read.table("2017_pm25_daily_average.txt", header = TRUE, sep = ","))
time_pm25_day_read

str(pm25_day_2017); head(pm25_day_2017)
# tmp <- pm25_day_2006[pm25_day_2006$FIPS==1001020100,]

pm25_day_2017$Month <- as.numeric(substr(as.character(pm25_day_2017$Date),6,7))
pm25_day_2017$Day<- as.numeric(substr(as.character(pm25_day_2017$Date),9,10))

mns<-unique(pm25_day_2017$Month)
unique(pm25_day_2017$Day)
names(pm25_day_2017)[5]<-"pm25_daily_average"
data.jan<-pm25_day_2017[pm25_day_2017$Month==01,]
unique(data.jan$Day)
map('state')
quilt.plot(data.jan$Longitude[data.jan$Day==01],data.jan$Latitude[data.jan$Day==01],data.jan$pm25_daily_average[data.jan$Day==01],nx=round(sqrt(length(data.jan$Day==01)),0),ny=round(sqrt(length(data.jan$Day==01)),0))
map('state',add = T)
zlowlim=min(log(data.jan$pm25_daily_average))
#zhighlim=max(data.jan$pm25_daily_average)
zhighlim=max(log(data.jan$pm25_daily_average))
hist(log(pm25_day_2017$pm25_daily_average))

### Splitting the data month-wise and day-wise to get a spatio-temporal process for the year #####
myfulldata<-pm25_day_2017[,c(3,4,5,7,8)]

########### Taking log transformation ##########
myfulldata$pm25_daily_average<-log(myfulldata$pm25_daily_average)
hist(myfulldata$pm25_daily_average)
daily.st<-list()

jan<-cbind(1,1:31)
feb<-cbind(2,1:28)
mar<-cbind(3,1:31)
apr<-cbind(4,1:30)
may<-cbind(5,1:31)
jun<-cbind(6,1:30)
jul<-cbind(7,1:31)
aug<-cbind(8,1:31)
sep<-cbind(9,1:30)
oct<-cbind(10,1:31)
nov<-cbind(11,1:30)
dec<-cbind(12,1:31)


full.cal<-rbind(jan,feb,mar,apr,may,jun,jul,aug,sep,oct,nov,dec)

for(i in 1:365)
{
  
  temp<-myfulldata[myfulldata$Month==full.cal[i,1]&myfulldata$Day==full.cal[i,2],]
  temp<-na.omit(temp)
  daily.st[[i]]<-temp
}

quilt.plot(daily.st[[6]][,c(1,2,3)])

latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}



# Assign States
#coords <- as.data.frame(cbind(total_201207$Longitude,total_201207$Latitude))
#State <- latlong2state(coords)
for(i in 1:365)
{
daily.st[[i]]$state<-latlong2state(daily.st[[i]][,c(1,2)])
}
#==============================================================================================####

for(i in 1:365)
{
  daily.st[[i]]<-na.omit(daily.st[[i]])
}



### will execute this part later #######
#unique(as.factor(State))
#summary(as.factor(State))

#id_NA <- (1:length(State))[is.na(State)]
#id_ST <- (1:length(State))[!is.na(State)]
#coords_NA <- coords[id_NA,]
#coords_ST <- coords[id_ST,]
#dim(coords_NA); dim(coords_ST)

#DMatrix <- distm(coords_NA,coords_ST); dim(DMatrix)
#tmp <- apply(DMatrix, 1, min)
#min_ind <- which(DMatrix==tmp, arr.ind=T)
#NA2ST <- id_ST[min_ind[order(min_ind[,1]),2]]

#State[id_NA] <- State[NA2ST]
#unique(as.factor(State))
#summary(as.factor(State))
#==============================================================================================####
# Assing Climatic Regions (CR)
CR_NW <- c("washington", "oregon", "idaho")
CR_W <- c("california", "nevada")
CR_SW <- c("utah", "colorado", "arizona", "new mexico")
CR_WNC <- c("montana", "wyoming", "north dakota", "south dakota", "nebraska")
CR_ENC <- c("minnesota", "iowa", "wisconsin", "michigan")
CR_S <- c("kansas", "oklahoma", "texas", "arkansas", "louisiana", "mississippi")
CR_C <- c("illinois", "indiana", "ohio", "missouri", "kentucky", "west virginia", "tennessee")
CR_NE <- c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
           "pennsylvania", "new jersey", "delaware", "maryland")
CR_SE <- c("south carolina", "georgia", "alabama", "florida", "north carolina", "virginia")


for(i in 1:365)
{
daily.st[[i]]$CR<-NA
State<-daily.st[[i]]$state
daily.st[[i]]$CR[State %in% CR_NW] <- "NW"
daily.st[[i]]$CR[State %in% CR_W]  <- "W"
daily.st[[i]]$CR[State %in% CR_SW] <- "SW"
daily.st[[i]]$CR[State %in% CR_WNC]<- "WNC"
daily.st[[i]]$CR[State %in% CR_ENC]<- "ENC"
daily.st[[i]]$CR[State %in% CR_S]  <- "S"
daily.st[[i]]$CR[State %in% CR_C]  <- "C"
daily.st[[i]]$CR[State %in% CR_NE] <- "NE"
daily.st[[i]]$CR[State %in% CR_SE] <- "SE"
daily.st[[i]]<-na.omit(daily.st[[i]])
}


rchoice<-"SE"

quilt.plot(daily.st[[i]][daily.st[[i]]$CR==rchoice,c(1,2,3)])

save.image("pmdata2017.RData")

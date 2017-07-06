# R TRAINING WORKSHOP FOR EPSCoR - June 19, 2017

# by Jory S. Hecht, Postdoctoral Associate, Vermont EPSCoR

# **********************************************************************

# SET WORKING DIRECTORY ------------------------------------------------

getwd()

# Laptop
setwd("C:/Users/Jory/ownCloud/documents/R Workshop")

# UVM computer
#setwd("C:/Users/jshecht/ownCloud/documents/R Workshop")


# LESSON 1: SOME R BASICS --------------------------------------------------------

# Create your first variable
X <- c(1,2,3,4,5)
X
class(X)    # numeric

# Determine the properties of X
str(X)
length(X) # 5

# Another way to check if X is numeric
A <- is.numeric(X)    
A           # TRUE
class(A)    # logical

# Create a vector with characters
J <- c("a","b","c","d")  #Must put characters in quote
class(J)

# Create a vector with different classes of data
K <- c(1,"b",2,TRUE)
class(K)

# Any time a vector has character data, it is classified as a character vector
# even if other data types are present

# If a vector has logical and numeric data, it is classified as a logical vector

# Create dates
today <- as.Date("2017/06/19")
class(today)

# Create dates with times (two Date-Time classes)
# Important if you're downloading results from an EPSCoR datalogger

# POSIXct counts the number of seconds since the beginning of Jan 1, 1970
# POSIXlt is a named list of vectors representing different time increments

# Find the current time
now <- Sys.time()
now
class(now)

# Express the current time as local time
now_lt <- as.POSIXlt(now)

# Let's look at some components of time
now_lt$hour
now_lt$mday
now_lt$mon
now_lt$year

# Can use functions like strptime to enter new dates

# Types of objects------------------------------------------------

# Operations 

# Create a vector of five 2's
Y <- rep(2,5)
# Recall X = 1,2,3,4,5
X

#Multiply X * Y
X*Y

# Create a 2 x 2 matrix D

# 2 4
# 5 7

D <- matrix(c(2,4,5,7),nrow=2,ncol=2)
class(D)
# [Row, Column]

# Element-wise multiplication
D*D

# Inner product
D%*%D

# Outer product
D%o%D

# Transpose
t(D)

# Summarize matrix rows and columns
D
apply(D,2,sum) #sum by column
apply(D,1,sum) #sum by row

# Built-in functions for summing columns and rows
colSums(D)
rowSums(D)

# Can use apply functions to perform any function over a row or column of a matrix
# or data frame 
# lapply and sapply functions for lists, tapply functions for subsets
browseURL("https://www.r-bloggers.com/using-apply-sapply-lapply-in-r/")

# Many matrix functions also available 
browseURL("http://www.statmethods.net/advstats/matrix.html")

# Add a column of 3's to the left of D
E <- cbind(c(3,3),D)

# Check dimensions and object type
dim(E)
str(E)

# Create a multidimensional array G (matrix with > 2 dimensions)

G <- array(c(E,E),dim=c(2,nrow(E),ncol(E)))
G
H <- array(c(E,E),dim=c(10,nrow(E),ncol(E)))

# Evaluate possible combinations of climate and land-use change scenarios

# We can have annual precipitation that is 80%, 100% or 120% of current rates
# We can have increases in annual mean temp of 0 C, 2 C and 4 C
# We can have agricultural land that is 90%, 100% or 110% of current agricultural land area

# Change in precipitation
dP <- c(0.80,1.00,1.20)

# Change in temperature
dT <- c(0,2,4)

# Change in agricultural land
dAgLand <- c(0.90,1.00,1.10)

# Compute the number of scenarios (3 x 3 x 3 = 27)
Scenarios <- expand.grid(dP,dT,dAgLand)

# Creating an empty array to store results
Results <- array(NA,dim=c(length(dP),
                          length(dT),
                          length(dAgLand)))

# This is an example of how you would great an array to store results 
# Future development of this lession should involve an exercise where 
# this array is used to store results

# Next, we're going to look at categorical data 

# You might want to identify a list of acceptable categories for data you enter
# For instance, let's say you're dealing with samples of total nitrogen and toal 
# phosphorus 

# Add factors to store unordered categorical data
WQ_levels <- c("Total Nitrogen","Total Phosphorus") #Lists acceptable categories
WQ_obs <- c("Total Phosphorus","Total Nitrogen","Helium") #Lists categories of observations
WQ_obs <- factor(WQ_obs,WQ_levels,ordered=FALSE) #Categories are not any kind of scale
is.na(WQ_obs) #Checks for N/A's, which will be assigned to invalid categories
WQ_obs[!is.na(WQ_obs)] # ! = bang operator meaning "is not"

# Add factors to store ordered categorical data
# For instance, let's say we have cloud cover from citizen science observations
# containing only qualitative descriptions of cloud cover

weather_levels <- c("Sunny","Partly Cloudy","Overcast","Stormy")
weather_obs <- c("Sunny","Partly Cloudy", "Partly Cloudy","Stormy")
weather_obs_ordered <- factor(weather_obs,weather_levels,ordered=TRUE)

# Now, let's move on to lists, which can be a very useful type of object
ls_test <- list(2,"a",TRUE)
ls_test2 <- list(2,c("a","b"),matrix(rep(1,4),nrow=2,ncol=2))
names(ls_test2) <- c("Number","Characters","Matrix")

# Extract an element of the list in which you are interested
ls_test2$Matrix

# Lists are useful when dealing with records of different length 

# Remove all the objects we just created ###############################
rm(list=ls())

# LESSON 2: DAILY STREAMFLOW ANALYSIS  -----------------------------------------------------------------

# Import USGS daily mean streamflow data
Q_imp <- read.table("QDaily_Swanton_USGS_04294000.csv",skip=30,header=TRUE,sep=",",fill=TRUE)
Q_imp[1:5,] #Display first five rows
class(Q_imp)
# Can also import directly from Excel ".xlsx" files, but must indicate sheet! 

## #Determine the class of each column

# Number of columns in the data frame
ncol(Q_imp)
# Create a new vector to store the types of data the function class returns
Q_imp_data_types <- vector(mode ="character",length=ncol(Q_imp))

# Identify the class of data in each column of the data frame
for (i in 1:ncol(Q_imp)){
  Q_imp_data_types[i]=class(Q_imp[,i])
}
Q_imp_data_types

# Can also try apply functions

# Create a new data frame
Q_Daily <- data.frame(as.Date(Q_imp[,3],format="%m/%d/%Y"),Q_imp[,4:5])
Q_Daily[1:5,]
# We don't need to use POSIX formats since we are looking at the daily mean discharges

# Create new column names
colnames(Q_Daily) <- c("Date","Discharge","DataQual")
Q_Daily[1:5,]

# Display individual columns
Q_Daily$Date
Q_Daily$Discharge

## # Check for missing data--------------------

# Make sure there are entries for all days in the record
nrow(Q_Daily)
max(Q_Daily$Date)-min(Q_Daily$Date)
as.numeric(max(Q_Daily$Date)-min(Q_Daily$Date)) # Introduce Day 1 = Jan 1, 1970

# Check for missing data
sum(is.na(Q_Daily$Discharge))

# Check for negative data
nrow(Q_Daily[Q_Daily$Discharge<0 & !is.na(Q_Daily$Discharge),])
# None 

# Find rows with missing data
Q_Daily_missing <- Q_Daily[is.na(Q_Daily$Discharge),]  
# [is.na(Q_Daily$Discharge),] refers to the rows and columns to be subset

# Display missing values
Q_Daily_missing  #  (should be seven)

# We find out that it's all the provisional data, show online

# Subset for WY 1991-2016 (Oct 1 - Sep 30)
Q_Daily_WY <- Q_Daily[Q_Daily$Date >= as.Date("1990-10-01") &
                      Q_Daily$Date <= as.Date("2016-09-30"),]

# Compare the number of rows (nrow) of the subset and original data
nrow(Q_Daily_WY)
nrow(Q_Daily)

# EXplore daily flow data
plot(Q_Daily_WY$Date,Q_Daily_WY$Discharge)
plot(Q_Daily_WY$Date,Q_Daily_WY$Discharge,typ="l",col="blue",
     xlab="Year",ylab="Discharge (cfs)",
     main="Daily Mean Discharge - \n Missisquoi River at Swanton, VT")

# See color ramp
browseURL("http://bc.bojanorama.pl/wp-content/uploads/2013/04/rcolorsheet.pdf")

# Plot the y-axis on a logarithmic scale
plot(Q_Daily_WY$Date,Q_Daily_WY$Discharge,log="y",typ="l",col="gray80",
     xlab="Year",ylab="Discharge (cfs)")
abline(h=mean(Q_Daily_WY$Discharge),lty=2,col="black",lwd=2)
leg.txt <- "Mean discharge (cfs)"
legend("bottomright",lty=2,col="black",leg.txt,bty="n",cex=0.8)

# Compute summary statistics for the entire period of record
summary(Q_Daily_WY$Discharge)
mean(Q_Daily_WY$Discharge)
median(Q_Daily_WY$Discharge)
max(Q_Daily_WY$Discharge)  #(is this peak flow? No.)
min(Q_Daily_WY$Discharge)
sd(Q_Daily_WY$Discharge)

# Compute histograms
hist(Q_Daily_WY$Discharge)
hist(Q_Daily_WY$Discharge,25)  # Make 25 bins
?hist #Read more about the "hist" function
help(hist) #Read more about the "hist" function

# Compute quantiles
quantile(Q_Daily_WY$Discharge,0.05)
quantile(Q_Daily_WY$Discharge,0.95)

# Compute many quantiles (from 0.1% to 99.9%)
Q_Daily_WY_quantiles <- vector(mode="numeric", length=999)

for(q in 1:999){
  Q_Daily_WY_quantiles[q] = quantile(Q_Daily_WY$Discharge,q/1000)
}

# Plot quantiles
plot(seq(0.001,0.999,0.001),Q_Daily_WY_quantiles)
plot(1-seq(0.001,0.999,0.001),Q_Daily_WY_quantiles)
plot(1-seq(0.001,0.999,0.001),Q_Daily_WY_quantiles,log="y",
     typ="l",lty=1,
     xlab="Exceedance Probability (%)",
     ylab="Discharge (cfs)",
     main="Flow Duration Curve") 
abline(v=0.05,lty=2,col="red")
abline(v=0.95,lty=2,col="red")
text(0.10,15000,"Q5",cex=0.8)
text(0.90,1000,"Q95",cex=0.8)

## Annual flow duration curves 

# Now compute statistics for different water years

# Modify the data frame by adding a new column indicating the water year of each obs
WY <- as.numeric(substring(Q_Daily_WY$Date+92,1,4)) #Adds new date

# Let's try making a water year function
Calc_WY <- function(Z) as.numeric(substring(Z+92,1,4))
?substring # Find out more about substring

# Let's try calculating the water year again
WY <- Calc_WY(Q_Daily_WY$Date)

# Let's attach the water year column to the data frame for easy reference
Q_Daily_WY <- cbind(Q_Daily_WY,WY) #Append new column using cbind

# Create a 99 x 26 array to store annual flow duration curves
RL <- max(Q_Daily_WY$WY)-min(Q_Daily_WY$WY)+1 #Computes number of water years in record
quantiles <- seq(1,99,1)

# Create a signle loop with an apply function
# tapply is useful for subsets
AFDC_WY <- array(NA,dim=c(length(quantiles),RL))
AFDC_WY2 <- array(NA,dim=c(length(quantiles),RL))

# Compute 1-99% exceedance probabilities (99% to 1% percentiles, note probs=(1-q/100))

for (q in 1:99){
  AFDC_WY[q,]=tapply(Q_Daily_WY$Discharge,Q_Daily_WY$WY,quantile,probs=(1-q/100))
}

# Create a double loop to illustrate another tool for estimating these quantiles
for (q in 1:99){
  for (w in 1:RL){
    AFDC_WY2[q,w]=quantile(as.numeric(Q_Daily_WY$Discharge[Q_Daily_WY$WY==1990+w]),1-q/100)
  }
}

# Check to make sure they produce identical results
max(AFDC_WY - AFDC_WY2)

# Plot annual flow duration curves
plot(seq(1,99,1),AFDC_WY[,8], #col 8 corresponds to WY 1998
     xlim=c(0,100),ylim=c(10,max(Q_Daily_WY$Discharge)),
     xlab="Exceedance Probability (%)",
     ylab="Daily Discharge (cfs)",
     log="y",typ="l",col="red")
for (i in 1:RL){
  if (i == 8)
  next
  lines(AFDC_WY[,i],col="gray80")
}
lines(AFDC_WY[,8],col="red",lwd=3)
title("1998 Water Year")

# Create a boxplot to explore interannual variability of different annual flow durations
boxplot(as.numeric(AFDC_WY[10,]),as.numeric(AFDC_WY[50,]),as.numeric(AFDC_WY[99,]),
        log="y",
        names=c("Q1","Q50","Q99"),
        main="Interannual flow variability")

# Compute mean annual discharge in each water year
Q_Annual_Mean <- apply(AFDC_WY,2,mean)

# Compare mean annual discharge with Q1
WYears <- seq(1991,2016,1)
plot(WYears,Q_Annual_Mean,typ="l",ylim=c(0,1.2*max(AFDC_WY[1,])))
lines(WYears,AFDC_WY[1,],typ="l")   
lines(WYears,AFDC_WY[50,],typ="l",col="red")

# LESSON 3: PRECIPITATION AND CORRELATION WITH STREAMFLOW -------------------------------------------------------------------------

# Import daily precipitation data
P_Daily_Enosburg <- read.csv("PDaily_Enosburg.csv")
P_Daily_Enosburg[1:5,]  # Displays the first five columns

# Format dates
P_Daily_Dates <- as.Date(P_Daily_Enosburg[,1],format="%m/%d/%Y")
P_Daily_Dates[1:5]

# Modify data frame
P_Daily_Enosburg <- data.frame(P_Daily_Dates,P_Daily_Enosburg[,2])
colnames(P_Daily_Enosburg) <- c("Date","Precip")  #Name columns

# Plot all days 
plot(P_Daily_Enosburg$Date,P_Daily_Enosburg$Precip,
     xlab="Year",ylab="Daily Precipitation (in)")

#Correlate daily precip and flows
cor_P_Q_Daily <- cor(P_Daily_Enosburg$Precip,Q_Daily_WY[,2])

# Compute the R^2 goodness-of-fit metric
r2 <- round(cor_P_Q_Daily^2,2)
plot(P_Daily_Enosburg[,2],Q_Daily_WY[,2],
     xlab="Daily Precip (in)",ylab="Daily Discharge (cfs)")
text(5.5,25000,bquote(R^2 == .(r2)),cex=0.75)
title("Daily precipitation vs. discharge")
# Less than 5% of the daily variance in the discharge is explained by the
# precipitaiton on that day. Why?
#[Seasonal factors, such as spring snowmelt and summer evapotranspiration, 
# have a much greater effect on discharge. Snow does not runoff into streams 
# immediately either!]
#[Antecedent soil moisture conditions, influenced by previous days' precipitation 
#also affects streamflow.]
#[In addition, there are differences in precipitation at Enosburg and other locations
# in the Missisquoi basin. The time it takes precipitation that falls in other 
# locations to reach Swanton should also be taken into account.]

#Plot a year of streamflow and precipitation on the same plot
Q_Daily_1998 <- Q_Daily_WY[Q_Daily_WY$WY==1998,]
P_Daily_1998 <- P_Daily_Enosburg[Calc_WY(P_Daily_Enosburg[,1])==1998,]

options(scipen=999) # Eliminates scientific notation on axes
par(mar=c(5.1,4.1,4.1,5.1))
plot(Q_Daily_1998$Date,Q_Daily_1998$Discharge,typ="l",col="blue",log="y",
     ylim=c(min(Q_Daily_1998[,2]),7*max(Q_Daily_1998[,2])),
     xlab="Month",ylab="Daily Mean Discharge (cfs)",
     main="Water Year 1998")

# Adding a second vertical axis  
par(new=T)
plot(P_Daily_1998$Date,P_Daily_1998$Precip,typ="l",col="gray60",
     xaxt="n",yaxt="n",xlab=NA,ylab=NA, cex.axis=0.8,
     ylim=rev(range(P_Daily_1998$Precip)*3))
axis(side=4) 
mtext(side=4,line=3,"Daily Precpitation (in)")
leg.txt <- c("Daily precipitation - Enosburg",
             "Daily mean discharge - Swanton")
legend("bottomleft",lty=1,col=c("gray80","blue"),leg.txt,bty="n",cex=0.6)

# Reset graphical parameters to default values
par(new=F)
par(mar=c(5.1,4.1,4.1,2.1))

# Read more about plot margins
browseURL("https://www.r-bloggers.com/mastering-r-plot-part-3-outer-margins/")

# Correlate annual precipitation and discharge 
WY_Enosburg <- Calc_WY(P_Daily_Enosburg$Date)
P_Daily_WY <- cbind(P_Daily_Enosburg,WY_Enosburg)
P_Annual <- tapply(P_Daily_WY$Precip,P_Daily_WY$WY_Enosburg,sum)

# Derive linear relationship
plot(P_Annual,Q_Annual_Mean,
     xlab = "Annual precipitation (in)",
     ylab = "Mean annual discharge (cfs)",
     main = "Annual precipitation vs. discharge")

# Compute correlation 
cor(P_Annual,Q_Annual_Mean)
r2_annual <- round(cor(P_Annual,Q_Annual_Mean)^2,2)
text(33,2800,bquote(R^2 == .(r2_annual)),cex=0.75)

# Compute linear relationship
lm.Q_P_Annual <- lm(Q_Annual_Mean ~ P_Annual)

# Extract coefficients
b0 <- summary(lm.Q_P_Annual)$coefficients[1,1]
b1 <- summary(lm.Q_P_Annual)$coefficients[2,1]

# Create trend line
min(P_Annual)
max(P_Annual)

# Plot
plot(P_Annual,Q_Annual_Mean,col="blue",xlab="Annual Precipitation (in)",
     ylab="Annual Mean Discharge (cfs)",
     main="Annual Precipitation vs. Discharge \n (WY 1991-2016)") 
lines(seq(31,56,1),b0+b1*seq(31,56,1))
residuals(lm.Q_P_Annual)
hist(residuals(lm.Q_P_Annual))
hist(Q_Annual_Mean)


# ANNUAL PEAK FLOWS ---------------------------------------------------------------------------

# Import data
Q_Peak_imp <- read.csv("Qpeak_Swanton_USGS_04294000.csv",skip=65)

# Remove unnecessary columns
Q_Peak <- data.frame(as.Date(Q_Peak_imp[,3],format="%m/%d/%Y"),Q_Peak_imp[,5:8])
colnames(Q_Peak) <- c("Date","Discharge","Q_Code","Stage","S_Code")

# Remove first entry for 1990 (not a complete water year)
nrow(Q_Peak) #27
Q_Peak <- Q_Peak[2:nrow(Q_Peak_imp),]
Q_Peak

# Identify annual floods by water year
WY_Q_Peak <- Calc_WY(Q_Peak$Date)

# What is the 100-year flood? 

# Add a package for extreme value analysis
install.packages("fitdistrplus") 
library("fitdistrplus")

# Quick check for trend
plot(Q_Peak$Discharge)
# Many more formal methods for trend detection! 
# No trend apparent 

# Explore annual flood time series and histogram
plot(Q_Peak$Date,Q_Peak$Discharge,typ="b",
     xlab="Time",
     ylab="Annual Inst. Peak Discharge (cfs)",
     main="Annual Flood Series (AFS) \n Missisquoi River at Swanton (WY 1991-2016)")

hist(Q_Peak$Discharge,xlab="Annual Inst. Peak Discharge (cfs)") # Number of bins

# Fit LN2 distribution
Q_Peak_fit_ln2 <- fitdist(as.numeric(Q_Peak[,2]),"lnorm")
Q_Peak_fit_ln2
Q_Peak_fit_ln2$estimate

# Compute values of mean and sd in log space
Q_Peak_mean_ln2 <- as.numeric(Q_Peak_fit_ln2$estimate[1])
Q_Peak_sd_ln2 <- as.numeric(Q_Peak_fit_ln2$estimate[2]) 
  
# Compute z-scores for different quantiles
qnorm(0.90,0,1)
qnorm(0.99,0,1)

# Compute floods with different recurrence intervals
Q_Peak_10Y <- exp(Q_Peak_mean_ln2+qnorm(0.90,0,1)*Q_Peak_sd_ln2)
Q_Peak_100Y <- exp(Q_Peak_mean_ln2+qnorm(0.99,0,1)*Q_Peak_sd_ln2)

# Display results
Q_Peak_10Y
Q_Peak_100Y

# Compare to flood of record
max(Q_Peak$Discharge)
# Be careful extrapolating

# Compute floods of different recurrence intervals
Q_Peak_Quantiles <- vector(mode="numeric",length=99)

# Compute quantiles for p = {0.01, 0.02,....,0.99}
for (p in 1:99){
  Q_Peak_Quantiles[p] = exp(Q_Peak_mean_ln2+qnorm(p/100,0,1)*Q_Peak_sd_ln2)
}

# Create a vector for the recurrence intervals of floods with p = {0.01,0.02,0.99}
Q_Peak_RI <- 1/(1-seq(1,99,1)/100)

# Plot curve relating recurrence intervals and quantile values
plot(Q_Peak_RI,Q_Peak_Quantiles,typ="l",
     log="x",ylim=c(0,max(Q_Peak_Quantiles)),
     xlab="Recurrence Interval (Years)",ylab="Annual Inst. Peak Flow (cfs)",
     main="Flood Frequency Curve \n Missisquoi River at Swanton")

Q_Peak_sorted <- sort(Q_Peak$Discharge,decreasing=TRUE)
Q_Peak_ranks <- rank(Q_Peak_sorted)/(length(Q_Peak$Discharge)+1)
Q_Peak_RI <- 1/Q_Peak_ranks

points(Q_Peak_RI,sort(Q_Peak[,2],decreasing=FALSE))


# loading data as matrix
library(rsample)
x <-as.matrix( read.csv("/home/saroj/Documents/assignment/X.csv", header = T))
y <-as.matrix( read.csv("/home/saroj/Documents/assignment/y.csv", header = T))
time <-read.csv("/home/saroj/Documents/assignment/time.csv", header = T)
##----------------Naming Column--------------
colnames(x) <- paste0(rep("x",ncol(x)),1:ncol(x))
colnames(y) <- "y"
colnames(time) <- "time"
##--------------------Binding All the data-------------
df<-cbind(time,y,x)
df
##--------------------Reading data of x-------------
x1 <- df$x1
x1
x2 <- df$x2
x2
y <- df$y # reading data of Y 
y
##--------------------Time series plots-------------
x.ts<-ts(x,start = c(-10,10),frequency = 1)
y.ts<-ts(y,start = c(-1,1),frequency =1)
time.ts<-ts(time,start = c(-1,1),frequency =1)

plot(x.ts, main = "Time Series Analysis ")
plot(y.ts, main = "Time Series Analysis of Y(Output)")
plot(time.ts)


##--------------------Density plots-------------

densityPlot <- density(x)
plot(densityPlot,main = "X Signal Density Plot")

lines(densityPlot,lwd=4,col="Blue")

# Add the data-points with noise in the X-axis
rug(jitter(x))

density_x1=density(x1)
hist(x1,freq = FALSE,main = "X1:Histogram and Density")
lines(density_x1,lwd=4,col="blue")
# Add the data-points with noise in the X-axis
rug(jitter(x1))

densityX2=density(x2)
hist(x2,freq = FALSE,main = "X2: Histogram and Density")
lines(densityX2,lwd=4,col="Blue")
rug(jitter(x2))

densityY=density(y)
plot(densityY,main = "Density plot of Y")
lines(densityY,lwd=4,col="Blue")
rug(jitter(y))

#Scatter/correlation plot
plot(x1,y, col = c("red","blue"), xlab = "Input signal [x1]", ylab = "Output Signal [Y]", main =  "Correlation of X1 over Y ")
cor(x1,y) #correlation of X1 over Y 
plot(x2,y, col = c("red","blue"),xlab = "Input signal [x2]", ylab = "Output Signal [Y]", main =  "Correlation of X2 over Y ")
cor(x2,y) #correlation of X2 over Y 

#box plots
boxplot(y~df$x2,main="Effect of sound categories on output brain signal",ylab="Output brain signal",ylim=c(0,70),xlab = "Sound categories")

##-----------------end of task 1----------------------------##

#TASK 2
a <- matrix(x1^3) # x1 values in matrix 
b <- matrix(x1^5) # x1 values in matrix 
c <- matrix(x2) # x2 values in matrix 


ones <- matrix(1,length(x1),1)
THM1 <- cbind(ones,a,b,c)
thetaHatM1 <- solve(t(THM1) %*% THM1) %*% t(THM1) %*% y 
thetaHatM1 ## thetahat of model1 
#--------------

d <- matrix(x1) # x1 values in matrix 
e <- matrix(x2) # x2 values in matrix 

THM2 <- cbind(ones,d,e)
thetaHatM2 <- solve(t(THM2) %*% THM2) %*% t(THM2) %*% y 
thetaHatM2 ## thetahat of model2 
#-----------------------

f <- matrix(x1) # x1 values in matrix 
g <- matrix(x1^2) # x1 values in matrix 
h <- matrix(x1^4) # x1 values in matrix 
i <- matrix(x2) # x2 values in matrix 

THM3 <- cbind(ones,f,g,h,i)
thetaHatM3 <- solve(t(THM3) %*% THM3) %*% t(THM3) %*% y 
thetaHatM3 ## thetahat of model3 
#-----------------------

j <- matrix(x1) # x1 values in matrix 
k <- matrix(x1^2) # x1 values in matrix 
l <- matrix(x1^3) # x1 values in matrix 
m <- matrix(x1^5) # x1 values in matrix 
n <- matrix(x2) # x2 values in matrix 

THM4 <- cbind(ones,j,k,l,m,n)
thetaHatM4 <- solve(t(THM4) %*% THM4) %*% t(THM4) %*% y 
thetaHatM4 ## thetahat of model4 
#-----------------------

o <- matrix(x1) # x1 values in matrix 
p <- matrix(x1^3) # x1 values in matrix 
q <- matrix(x1^4) # x1 values in matrix 
r <- matrix(x2) # x2 values in matrix 

THM5 <- cbind(ones,o,p,q,r)
thetaHatM5 <- solve(t(THM5) %*% THM5) %*% t(THM5) %*% y 
thetaHatM5 ## thetahat of model5 
#-----------------------
##------------------Completion of TASK 2.1-----
#calculating Y hat 

YhatM1 <- THM1 %*% thetaHatM1 #model 1 Yhat 
YhatM1

YhatM2 <- THM2 %*% thetaHatM2 # model 2 Yhat 
YhatM2

YhatM3 <- THM3 %*% thetaHatM3 # model 3 Yhat 
YhatM3

YhatM4 <- THM4 %*% thetaHatM4 # model 4 Yhat 
YhatM4

YhatM5 <- THM5 %*%thetaHatM5 #model 5 Yhat 
YhatM5

RSSM1 <- sum((y -YhatM1)^2)  
RSSM1 # model 1 RSS 

RSSM2 <- sum((y-YhatM2)^2)
RSSM2 # model 2 RSS

RSSM3 <- sum((y-YhatM3)^2) 
RSSM3 # model 3 RSS 

RSSM4 <- sum((y- YhatM4)^2)
RSSM4  # model 4 RSS 

RSSM5 <- sum((y-YhatM5)^2)
RSSM5 # model 5 RSS 

# task2.2 completed ----------

#Task2.3 
n <- length(y) # length of y 
n

varianceM1 <- RSSM1/(n-1) # calculating variance of model1 
likelyhoodM1 <- -n/2  * log(2* pi) - n/2 * log(varianceM1)-((1/2) * varianceM1)*RSSM1  
likelyhoodM1 #likelihood of model 1
#---------------------------------------------------------
varianceM2 <- RSSM2/(n-1) # calculating variance of Model 2 
likelyhoodM2 <- -n/2  * log(2* pi) - n/2 * log(varianceM2)-(1/2* varianceM2)*RSSM2 
likelyhoodM2 #likelihood of model 2 
#----------------------------------------------------------
varianceM3 <- RSSM3/(n-1) # calculating variance of Model 3 
likelyhoodM3 <- -n/2  * log(2* pi) - n/2 * log(varianceM3)-(1/2* varianceM3)*RSSM3 
likelyhoodM3 #likelihood of model 3 
#------------------------------------------------------------
varianceM4 <- RSSM4/(n-1) #calculating variance of model4 
likelyhoodM4 <- -n/2  * log(2* pi) - n/2 * log(varianceM4)-(1/2* varianceM4)*RSSM4 
likelyhoodM4 #likelihood of model 4 
#-----------------------------------------------------------
varianceM5 <- RSSM5/(n-1) # calculating variance of mode5 
likelyhoodM5 <- -n/2  * log(2* pi) - n/2 * log(varianceM5)-(1/2* varianceM5)*RSSM5 
likelyhoodM5 #likelihood of model 5

#---TASK2.3 Completed-----
#--TASK 2.4--------
Mk1 <- length(thetaHatM1) #value of K form model 1
AICM1 <- 2 * Mk1 - 2* likelyhoodM1 #model 1 AIC 
AICM1  #model 1 AIC 

BICM1 <- Mk1 * log(n) - 2 * likelyhoodM1 #BIC of MOdel 1
BICM1 #BIC of MOdel 1
#----------------------------------------------
Mk2 <- length(thetaHatM2) #value of k for model 2
AICM2 <- 2 * Mk2 - 2* likelyhoodM2 #model 2 AIC 
AICM2 #model 2 AIC 

BICM2 <- Mk2 * log(Mk2) - 2 * likelyhoodM2 #BIC of MOdel 2
BICM2 #BIC of MOdel 2
#----------------------------------------------
Mk3 <- length(thetaHatM3)#value of k for model 3
AICM3 <- 2 * Mk3 - 2* likelyhoodM3 #model 3 AIC
AICM3 #model 3 AIC

BICM3 <- Mk3 * log(Mk3) - 2 * likelyhoodM3 #BIC of MOdel 3
BICM3  #BIC of MOdel 3
#------------------------------------------------
Mk4 <- length(thetaHatM4) # value of K for model 4
AICM4 <- 2 * Mk4 - 2* likelyhoodM4 #model 4 AIC
AICM4  #model 4 AIC

BICM4 <- Mk4 * log(Mk4) - 2 * likelyhoodM4 #BIC of MOdel 4
BICM4 #BIC of MOdel 4
#------------------------------------------------
Mk5 <- length(thetaHatM5)#valu of k for model 5
AICM5 <-  2 * Mk5 - 2* likelyhoodM5 #model 5 AIC
AICM5 #model 5 AIC

BICM5 <- Mk5 * log(Mk5) - 2 * likelyhoodM5 #BIC of MOdel 5
BICM5 #BIC of MOdel 5

#-----TASK 2.4 completed----------
#Task 2.5-------------------------
errorM1 <- y-YhatM1 # error of model 1

qqnorm(errorM1, col = "blue") # model 1 QQ plot 
qqline(errorM1, col = "red") # model 1 QQ line 

errorM2 <- y-YhatM2 # error of model 2

qqnorm(errorM2, col = "blue") # model 2 QQ plot 
qqline(errorM2, col = "red") # model 2 QQ line

errorM3 <- y- YhatM3 #error of model 3 

qqnorm(errorM3, col = "blue") # model 3 QQ plot 
qqline(errorM3, col = "red") # model 3 QQ line 

errorM4 <- y-YhatM4 # error of model 4 

qqnorm(errorM4, col = "blue") # model 4 QQ plot 
qqline(errorM4, col = "red") # model 4 QQ line 

errorM5 <- y- YhatM5

qqnorm(errorM5, col = "blue") # model 5 QQ plot 
qqline(errorM5, col = "red") # model 5 QQ line 

#Task 2.5 completed-----
#-----TASK 2.7-------

Ysplt<-initial_split(data = as.data.frame(y),prop=.7) # spliting data into training and testing 
traningY<-training(Ysplt)
testingY<-as.matrix(testing(Ysplt))
YTD<-as.matrix(traningY) #Ytesting data 
length(testingY)


#spliting X into testing and training 
Xsplit<-initial_split(data = as.data.frame(x),prop=.7)
trainingX<-training(Xsplit)
testingX<-as.matrix(testing(Xsplit))
XTD<-as.matrix(testingX)

YTD<-as.matrix(trainingX)
head(trainingX)
tn <- length(trainingX$x1)


### Estimating model parameters using Traning set

traning_ones=matrix(1 ,tn,1)
thainingThetaHatModel<-cbind(traning_ones,trainingX[,"x1"],(trainingX[,"x1"])^2,(trainingX[,"x1"])^4,(trainingX[,"x2"]))
thainingThetaHatModel


thetahatTrain=solve(t(thainingThetaHatModel) %*% thainingThetaHatModel) %*% t(thainingThetaHatModel) %*% YTD
thetahatTrain

### Model out/Prediction
#YhatTesting = XTD %*% thetahatTrain
#YhatTesting
#TestingRss=sum((traningY-YhatTesting)^2)
#TestingRss

t.test(YTD, mu=500, alternative="two.sided", conf.level=0.95)

C_I1=-0.1540987
C_I2=0.3484372

p2 <- plot(density(YTD), col="blue", lwd=2,
           main="Distribution of Traning Data")
abline(v=C_I1,col="red", lty=2)
abline(v=C_I2,col="red", lty=2)

##-----TASK 2.7 Completed----

##Task 3 
##best moduel has been selected as model 3 and from it we select largest two values and remaining we keep it constant
#values were 
thetebias <- 10.1736961 #selected values 
thetaone <- 8.5521613 # selected values 
thetatwo <- 6.2469171 #keeping constant 
thetathree <- 4.1598833 #keeping constant
thetafour <- -0.2829631 # keeping constant
E <- RSSM3 * 2 ## fixing value of eplision 
n <- 1000 #number of iteration 
t1s <- 0
t2s <- 0

##finding out  new Y hat in order to accept or reject based on threashold..
counter <- 0
for (i in 1:n) {
  t1 <- runif(1,-10.1736961,10.1736961) # finding the range 
  t1
  t2 <- runif(1,-8.5521613,8.5521613) # finding the range....
  
  NewTH <- matrix(c(t1,t2,thetatwo,thetathree,thetafour))
  
  New_Y_Hat <- THM3 %*% NewTH ## finding New Y hat... 
  
  new_RSS <- sum((y-New_Y_Hat)^2)
  new_RSS 
  if (new_RSS < E){
    t1s[i] <- t1
    t2s [i] <- t2
    counter = counter+1
    
    a <- matrix(t1s)
    b <- matrix(t2s)
  }
}
hist(a)
hist(b)
plot(a,b, col = c("red", "blue"), main = "Posterior Distribution")
##---TASK 3 completed-------





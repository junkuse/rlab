 data(cars)
 head(cars)

 plot(cars$dist ~ cars$speed)

 lines(lowess(cars$dist ~ cars$speed))

 a <- lm(dist ~ poly(speed,1,raw=TRUE),data=cars)
 b <- lm(dist ~ poly(speed,2,raw=TRUE),data=cars)

 lines(cars$speed, predict(a),col=2)
 lines(cars$speed, predict(b),col=3)

print(b)

c <- data.frame(speed = 17)

result <-  predict(b,c)

print(result)


options(prompt="R> ")

##################
### CHAPTER 17 ###
##################

################
# Section 17.1 #
################

## 17.1.1 ##

xvals <- seq(16,28,by=0.1)
fx.samp <- dnorm(xvals,22,1.5/sqrt(5))
plot(xvals,fx.samp,type="l",lty=2,lwd=2,xlab="",ylab="")
abline(h=0,col="gray")
fx <- dnorm(xvals,22,1.5)
lines(xvals,fx,lwd=2)
legend("topright",legend=c("raw obs. distbn.","sampling distbn. (mean)"),lty=1:2,lwd=c(2,2),bty="n") # this plot is continued below

#

pnorm(21.5,mean=22,sd=1.5)

#

pnorm(21.5,mean=22,sd=1.5/sqrt(5))

#

abline(v=21.5,col="gray")
xvals.sub <- xvals[xvals<=21.5]
fx.sub <- fx[xvals<=21.5]
fx.samp.sub <- fx.samp[xvals<=21.5]
polygon(cbind(c(21.5,xvals.sub),c(0,fx.sub)),density=10)
polygon(cbind(c(21.5,xvals.sub),c(0,fx.samp.sub)),density=10,angle=120,lty=2)

#

obs <- rnorm(5,mean=22,sd=1.5)  # In the text, my five values are approximately as follows: obs <- c(22.92233,23.09505,20.98653,20.10941,22.33888)
obs

#

obs.mean <- mean(obs)
obs.mean
obs.sd <- sd(obs)
obs.sd

#

obs.mean.se <- obs.sd/sqrt(5)
obs.mean.se

#

t4 <- (21.5-obs.mean)/obs.mean.se
t4

#

pt(t4,df=4)

#

xvals <- seq(-5,5,length=100)
fx.samp.t <- dt(xvals,df=4)
plot(xvals,dnorm(xvals),type="l",lty=2,lwd=2,col="gray",xlim=c(-4,4),xlab="",ylab="")
abline(h=0,col="gray")
lines(xvals,fx.samp.t,lty=3,lwd=2)
polygon(cbind(c(t4,-5,xvals[xvals<=t4]),c(0,0,fx.samp.t[xvals<=t4])),density=10,lty=3)
legend("topright",legend=c("N(0,1) standard","t (4 df)"),col=c("gray","black"),lty=2:3,lwd=c(2,2),bty="n")


## 17.1.2 ##

p.hat <- 80/118
p.hat

#

118*p.hat
118*(1-p.hat)

#

p.se <- sqrt(p.hat*(1-p.hat)/118)
p.se

#

pvals <- seq(p.hat-5*p.se,p.hat+5*p.se,length=100)
p.samp <- dnorm(pvals,mean=p.hat,sd=p.se)
plot(pvals,p.samp,type="l",xlab="",ylab="",xlim=p.hat+c(-4,4)*p.se,ylim=c(0,max(p.samp)))
abline(h=0,col="gray")

#

pvals.sub <- pvals[pvals>=0.7 & pvals<=0.75]
p.samp.sub <- p.samp[pvals>=0.7 & pvals<=0.75]
polygon(cbind(c(0.7,pvals.sub,0.75),c(0,p.samp.sub,0)),border=NA,col="gray")

#

pnorm(0.75,mean=p.hat,sd=p.se) - pnorm(0.7,mean=p.hat,sd=p.se)


## 17.1.3 ##

# no R code



################
# Section 17.2 #
################

## 17.2.1 ##

temp.sample <- rnorm(n=5,mean=22,sd=1.5) # In the text, my five values are approximately as follows: temp.sample <- c(20.46097,21.45658,21.06410,20.49367,24.92843)
temp.sample

#

temp.mean <- mean(temp.sample)
temp.mean
temp.sd <- sd(temp.sample)
temp.sd
temp.se <- temp.sd/sqrt(5)
temp.se

#

1-0.05/2
critval <- qt(0.975,df=4)
critval

#

pt(critval,4)-pt(-critval,4)

#

temp.mean-critval*temp.se
temp.mean+critval*temp.se

#

temp.mean+c(-1,1)*qt(p=0.9,df=4)*temp.se
temp.mean+c(-1,1)*qt(p=0.995,df=4)*temp.se


## 17.2.2 ##

p.hat <- 80/118
p.hat
p.se <- sqrt(p.hat*(1-p.hat)/118)
p.se

#

qnorm(0.95)

#

p.hat+c(-1,1)*qnorm(0.95)*p.se


## 17.2.3 ##

# no R code


## 17.2.4 ##

# no R code




##################
### CHAPTER 18 ###
##################

################
# Section 18.1 #
################

dbinom(18,20,0.07) + dbinom(19,20,0.07) + dbinom(20,20,0.07)


## 18.1.1 ##

# no R code


## 18.1.2 ##

# no R code


## 18.1.3 ##

# no R code


## 18.1.4 ##

# no R code


## 18.1.5 ##

# no R code



################
# Section 18.2 #
################

## 18.2.1 ##

snacks <- c(87.7,80.01,77.28,78.76,81.52,74.2,80.71,79.5,77.87,81.94,80.7,82.32,
            75.78,80.19,83.91,79.4,77.52,77.62,81.4,74.89,82.95,73.59,77.92,77.18,
            79.83,81.23,79.28,78.44,79.01,80.47,76.23,78.89,77.14,69.94,78.54,79.7,
            82.45,77.29,75.52,77.21,75.99,81.94,80.41,77.7)

#

n <- length(snacks)
snack.mean <- mean(snacks)
snack.mean
snack.sd <- sd(snacks)
snack.sd

#

snack.se <- snack.sd/sqrt(n)
snack.se

#

snack.T <- (snack.mean-80)/snack.se
snack.T

#

pt(snack.T,df=n-1)

#

snack.mean+c(-1,1)*qt(0.975,n-1)*snack.se

#

t.test(x=snacks,mu=80,alternative="less")

#

t.test(x=snacks,mu=80,alternative="two.sided")$conf.int  


## 18.2.2 ##

snacks2 <- c(80.22,79.73,81.1,78.76,82.03,81.66,80.97,81.32,80.12,78.98,79.21,
             81.48,79.86,81.06,77.96,80.73,80.34,80.01,81.82,79.3,79.08,79.47,
             78.98,80.87,82.24,77.22,80.03,79.2,80.95,79.17,81)

#

snack2.mean <- mean(snacks2)
snack2.mean
snack2.sd <- sd(snacks2)
snack2.sd

#

t.test(x=snacks2,y=snacks,alternative="greater",conf.level=0.9)

#

(snack2.mean-snack.mean) + c(-1,1)*qt(0.95,df=60)*sqrt(snack.sd^2/44+snack2.sd^2/31)

#

snack.sd/snack2.sd

#

men <- c(102,87,101,96,107,101,91,85,108,67,85,82)
women <- c(73,81,111,109,143,95,92,120,93,89,119,79,90,126,62,92,77,106,105,111)

#

mean(men)
sd(men)
mean(women)
sd(women)

#

sd(women)/sd(men)

#

t.test(x=men,y=women,alternative="two.sided",conf.level=0.95,var.equal=TRUE)

#

rate.before <- c(52,66,89,87,89,72,66,65,49,62,70,52,75,63,65,61) 
rate.after <- c(51,66,71,73,70,68,60,51,40,57,65,53,64,56,60,59) 

#

rate.d <- rate.after-rate.before
rate.d

#

rate.dbar <- mean(rate.d)
rate.dbar
rate.sd <- sd(rate.d)
rate.sd

#

rate.T <- rate.dbar/(rate.sd/sqrt(16))
rate.T
pt(rate.T,df=15)

#

t.test(x=rate.after,y=rate.before,alternative="less",conf.level=0.95,paired=TRUE)

#

rate.dbar-qt(0.975,df=15)*(rate.sd/sqrt(16))

#

rate.dbar+qt(0.975,df=15)*(rate.sd/sqrt(16))



################
# Section 18.3 #
################

## 18.3.1 ##

sick <- c(0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1)

#

sum(sick)
p.hat <- mean(sick)
p.hat

#

29*0.2
29*0.8

#

Z <- (p.hat-0.2)/sqrt(0.2*0.8/29)
Z

#

2*(1-pnorm(Z))

#

p.hat+c(-1,1)*qnorm(0.975)*sqrt(p.hat*(1-p.hat)/29)

#

prop.test(x=sum(sick),n=length(sick),p=0.2,correct=FALSE)


## 18.3.2 ##

x1 <- 180
n1 <- 233
p.hat1 <- x1/n1
p.hat1
x2 <- 175
n2 <- 197
p.hat2 <- x2/n2
p.hat2

#

p.star <- (x1+x2)/(n1+n2)
p.star

#

Z <- (p.hat2-p.hat1)/sqrt(p.star*(1-p.star)*(1/n1+1/n2))
Z

#

1-pnorm(Z)

#

prop.test(x=c(x2,x1),n=c(n2,n1),alternative="greater",correct=FALSE)

#

(p.hat2-p.hat1)+c(-1,1)*qnorm(0.975)*sqrt(p.hat1*(1-p.hat1)/n1+p.hat2*(1-p.hat2)/n2)



################
# Section 18.4 #
################

## 18.4.1 ##

hairy <- c(2,3,2,3,2,1,3,3,2,2,3,2,2,2,3,3,3,2,3,2,2,2,1,3,2,2,2,1,2,2,3,2,2,2,2,1,2,1,1,1,2,2,2,3,1,2,1,2,1,2,1,3,3)

#

x <- seq(0,20,length=100)
plot(x,dchisq(x,df=1),type="l",xlim=c(0,15),ylim=c(0,0.5),ylab="density")
lines(x,dchisq(x,df=5),lty=2)
lines(x,dchisq(x,df=10),lty=3)
abline(h=0,col="gray")
abline(v=0,col="gray")
legend("topright",legend=c("df=1","df=5","df=10"),lty=1:3)

#

n <- length(hairy)
n
hairy.tab <- table(hairy)
hairy.tab
hairy.tab/n

#

expected <- 1/3*n
expected
hairy.matrix <- cbind(1:3,hairy.tab,expected,(hairy.tab-expected)^2/expected)
dimnames(hairy.matrix) <- list(c("clean","beard OR mous.","beard AND mous."),c("i","Oi","Ei","(Oi-Ei)^2/Ei"))
hairy.matrix

#

X2 <- sum(hairy.matrix[,4])
X2

#

1-pchisq(X2,df=2)

#

chisq.test(x=hairy.tab)

#

chisq.test(x=hairy.tab,p=c(0.25,0.5,0.25))


## 18.4.2 ##

skin <- matrix(c(20,32,8,52,9,72,8,32,16,64,30,12),4,3,dimnames=list(c("Injection","Tablet","Laser","Herbal"),c("None","Partial","Full")))
skin

#

kr <- nrow(skin)
kc <- ncol(skin)

#

rowSums(skin)
colSums(skin)

#

rep(colSums(skin),each=kr)
rep(colSums(skin),each=kr)*rowSums(skin)
rep(colSums(skin),each=kr)*rowSums(skin)/sum(skin)
skin.expected <- matrix(rep(colSums(skin),each=kr)*rowSums(skin)/sum(skin),nrow=kr,ncol=kc,dimnames=dimnames(skin))
skin.expected


#

skin.array <- array(data=cbind(skin,skin.expected,(skin-skin.expected)^2/skin.expected),dim=c(kr,kc,3),dimnames=list(dimnames(skin)[[1]],dimnames(skin)[[2]],c("O[i,j]","E[i,j]","(O[i,j]-E[i,j])^2/E[i,j]")))
skin.array

#

X2 <- sum(skin.array[,,3])
X2

#

1-pchisq(X2,df=(kr-1)*(kc-1))

#

chisq.test(x=skin)



################
# Section 18.5 #
################

## 18.5.1 ##

# No R code


## 18.5.2 ##

typeI.tester <- function(mu0,sigma,n,alpha,ITERATIONS=10000){
  pvals <- rep(NA,ITERATIONS)
  for(i in 1:ITERATIONS){
    temporary.sample <- rnorm(n=n,mean=mu0,sd=sigma)
    temporary.mean <- mean(temporary.sample)
    temporary.sd <- sd(temporary.sample)
    pvals[i] <- 1-pt((temporary.mean-mu0)/(temporary.sd/sqrt(n)),df=n-1)
  }
  return(mean(pvals<alpha))
}

#

typeI.tester(mu0=0,sigma=1,n=40,alpha=0.05)

#

typeI.tester(mu0=-4,sigma=0.3,n=60,alpha=0.01)


## 18.5.3 ##

critval <- qnorm(1-0.05,mean=0,sd=1/sqrt(30))
critval

#

pnorm(critval,mean=0.5,sd=1/sqrt(30))

#

typeII.tester <- function(mu0,muA,sigma,n,alpha,ITERATIONS=10000){
  pvals <- rep(NA,ITERATIONS)
  for(i in 1:ITERATIONS){
    temporary.sample <- rnorm(n=n,mean=muA,sd=sigma)
    temporary.mean <- mean(temporary.sample)
    temporary.sd <- sd(temporary.sample)
    pvals[i] <- 1-pt((temporary.mean-mu0)/(temporary.sd/sqrt(n)),df=n-1)
  }
  return(mean(pvals>=alpha))
}

#

typeII.tester(mu0=0,muA=0.5,sigma=1,n=30,alpha=0.05)

#

typeII.tester(mu0=0,muA=0.5,sigma=1,n=30,alpha=0.01)

#

typeII.tester(mu0=0,muA=0.5,sigma=1.1,n=30,alpha=0.01)
typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=30,alpha=0.01)

#

typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=20,alpha=0.01)
typeII.tester(mu0=0,muA=0.5,sigma=1.2,n=40,alpha=0.01)

#

typeII.tester(mu0=0,muA=0.4,sigma=1.2,n=40,alpha=0.01)
typeII.tester(mu0=0,muA=0.6,sigma=1.2,n=40,alpha=0.01)


## 18.5.4 ##

power.tester <- function(nvec,...){
  nlen <- length(nvec)
  result <- rep(NA,nlen)
  pbar <- txtProgressBar(min=0,max=nlen,style=3)
  for(i in 1:nlen){
    result[i] <- 1-typeII.tester(n=nvec[i],...)
    setTxtProgressBar(pbar,i)
  }
  close(pbar)
  return(result)
}

#

sample.sizes <- 5:100

#

pow <- power.tester(nvec=sample.sizes,mu0=0,muA=0.6,sigma=1.2,alpha=0.01,ITERATIONS=5000)
pow

#

minimum.n <- sample.sizes[min(which(pow>=0.8))]
minimum.n

#

pow2 <- power.tester(nvec=sample.sizes,mu0=0,muA=0.6,sigma=1.2,alpha=0.05,ITERATIONS=5000)
minimum.n2 <- sample.sizes[min(which(pow2>0.8))]
minimum.n2

#

plot(sample.sizes,pow,xlab="sample size n",ylab="simulated power")
points(sample.sizes,pow2,col="gray")
abline(h=0.8,lty=2)
abline(v=c(minimum.n,minimum.n2),lty=3,col=c("black","gray"))
legend("bottomright",legend=c("alpha=0.01","alpha=0.05"),col=c("black","gray"),pch=1)


## 18.5.5 ##

# no R Code




##################
### CHAPTER 19 ###
##################

################
# Section 19.1 #
################

## 19.1.1 ##

table(chickwts$feed)
chick.means <- tapply(chickwts$weight,INDEX=chickwts$feed,FUN=mean)
chick.means

#

boxplot(chickwts$weight~chickwts$feed)
points(1:6,chick.means,pch=4,cex=1.5)

#

chick.sds <- tapply(chickwts$weight,INDEX=chickwts$feed,FUN=sd)
chick.sds

#

chick.meancen <- chickwts$weight-chick.means[as.numeric(chickwts$feed)]

#

qqnorm(chick.meancen,main="Normal QQ plot of residuals")
qqline(chick.meancen)


## 19.1.2 ##

# no R code


## 19.1.3 ##

chick.anova <- aov(weight~feed,data=chickwts)

#

summary(chick.anova)

#

195556/65



################
# Section 19.2 #
################

## 19.2.1 ##

tapply(warpbreaks$breaks,INDEX=list(warpbreaks$wool,warpbreaks$tension),FUN=mean)

#

wb.means <- aggregate(warpbreaks$breaks,by=list(warpbreaks$wool,warpbreaks$tension),FUN=mean)
wb.means


## 19.2.2 ##

summary(aov(breaks~wool,data=warpbreaks))
summary(aov(breaks~tension,data=warpbreaks))

#

summary(aov(breaks~wool+tension,data=warpbreaks))

#

summary(aov(breaks~wool+tension+wool:tension,data=warpbreaks))

#

interaction.plot(x.factor=wb.means[,2],trace.factor=wb.means[,1],response=wb.means$x,trace.label="wool",xlab="tension",ylab="mean warp breaks")



################
# Section 19.3 #
################

library("MASS")
boxplot(Age~Smoke,data=survey)
age.means <- tapply(survey$Age,survey$Smoke,mean)
age.meancen <- survey$Age-age.means[as.numeric(survey$Smoke)]
qqnorm(age.meancen,main="Normal QQ plot of residuals")
qqline(age.meancen)

#

tapply(survey$Age,survey$Smoke,sd)

#

kruskal.test(Age~Smoke,data=survey)




##################
### CHAPTER 20 ###
##################

################
# Section 20.1 #
################

library("MASS")
plot(survey$Height~survey$Wr.Hnd,xlab="Writing handspan (cm)",ylab="Height (cm)")

#

cor(survey$Wr.Hnd,survey$Height,use="complete.obs")

#

incomplete.obs <- which(is.na(survey$Height)|is.na(survey$Wr.Hnd))
length(incomplete.obs)



################
# Section 20.2 #
################

## 20.2.1 ##

# no R code


## 20.2.2 ##

# no R code


## 20.2.3 ##

survfit <- lm(Height~Wr.Hnd,data=survey)

#

survfit

#

plot(survey$Height~survey$Wr.Hnd,xlab="Writing handspan (cm)",ylab="Height (cm)")
abline(survfit,lwd=2)


## 20.2.4 ##

obsA <- c(survey$Wr.Hnd[197],survey$Height[197])
obsA
obsB <- c(survey$Wr.Hnd[154],survey$Height[154])
obsB

#

names(survfit)

#

mycoefs <- coef(survfit)
mycoefs
beta0.hat <- mycoefs[1]
beta1.hat <- mycoefs[2]

#

plot(survey$Height~survey$Wr.Hnd,xlab="Writing handspan (cm)",ylab="Height (cm)")
abline(survfit,lwd=2)
segments(x0=c(obsA[1],obsB[1]),y0=beta0.hat+beta1.hat*c(obsA[1],obsB[1]),x1=c(obsA[1],obsB[1]),y1=c(obsA[2],obsB[2]),lty=2)



################
# Section 20.3 #
################

## 20.3.1 ##

summary(survfit)


## 20.3.2 ##

confint(survfit,level=0.95)


## 20.3.3 ##

rho.xy <- cor(survey$Wr.Hnd,survey$Height,use="complete.obs")
rho.xy^2


## 20.3.4 ##

names(summary(survfit))

#

summary(survfit)$sigma



################
# Section 20.4 #
################

## 20.4.1 ##

# no R code


## 20.4.2 ##

as.numeric(beta0.hat+beta1.hat*14.5)
as.numeric(beta0.hat+beta1.hat*24)

#

xvals <- data.frame(Wr.Hnd=c(14.5,24))
xvals
mypred.ci <- predict(survfit,newdata=xvals,interval="confidence",level=0.95)
mypred.ci

#

mypred.pi <- predict(survfit,newdata=xvals,interval="prediction",level=0.95)
mypred.pi


## 20.4.3 ##

plot(survey$Height~survey$Wr.Hnd,xlim=c(13,24),ylim=c(140,205),xlab="Writing handspan (cm)",ylab="Height (cm)")
abline(survfit,lwd=2)

#

points(xvals[,1],mypred.ci[,1],pch=8)
segments(x0=c(14.5,24),x1=c(14.5,24),y0=c(mypred.pi[1,2],mypred.pi[2,2]),y1=c(mypred.pi[1,3],mypred.pi[2,3]),col="gray",lwd=3)
segments(x0=c(14.5,24),x1=c(14.5,24),y0=c(mypred.ci[1,2],mypred.ci[2,2]),y1=c(mypred.ci[1,3],mypred.ci[2,3]),lwd=2)

#

xseq <- data.frame(Wr.Hnd=seq(12,25,length=100))
ci.band <- predict(survfit,newdata=xseq,interval="confidence",level=0.95)
pi.band <- predict(survfit,newdata=xseq,interval="prediction",level=0.95)

#

lines(xseq[,1],ci.band[,2],lty=2)
lines(xseq[,1],ci.band[,3],lty=2)
lines(xseq[,1],pi.band[,2],lty=2,col="gray")
lines(xseq[,1],pi.band[,3],lty=2,col="gray")
legend("topleft",legend=c("Fit","95% CI","95% PI"),lty=c(1,2,2),col=c("black","black","gray"),lwd=c(2,1,1))


## 20.4.4 ##

predict(survfit,newdata=data.frame(Wr.Hnd=50),interval="confidence",level=0.95)



################
# Section 20.5 #
################

## 20.5.1 ##

class(survey$Sex)
table(survey$Sex)

#

plot(survey$Height~survey$Sex)

#

points(survey$Height~as.numeric(survey$Sex),cex=0.5)

#

means.sex <- tapply(survey$Height,INDEX=survey$Sex,FUN=mean,na.rm=TRUE)
means.sex
points(1:2,means.sex,pch=4,cex=3)

#

survfit2 <- lm(Height~Sex,data=survey)
summary(survfit2)

#

extra.obs <- factor(c("Female","Male","Male","Male","Female"))
extra.obs
predict(survfit2,newdata=data.frame(Sex=extra.obs),interval="confidence",level=0.9)


## 20.5.2 ##

is.factor(survey$Smoke)
table(survey$Smoke)
levels(survey$Smoke)

#

boxplot(Height~Smoke,data=survey)
points(1:4,tapply(survey$Height,survey$Smoke,mean,na.rm=TRUE),pch=4)

#

survfit3 <- lm(Height~Smoke,data=survey)
summary(survfit3)

#

one.of.each <- factor(levels(survey$Smoke))
one.of.each
predict(survfit3,newdata=data.frame(Smoke=one.of.each),interval="confidence",level=0.95)


## 20.5.3 ##

SmokeReordered <- relevel(survey$Smoke,ref="Never")
levels(SmokeReordered)


## 20.5.4 ##

boxplot(mtcars$mpg~mtcars$cyl,xlab="Cylinders",ylab="MPG")

#

class(mtcars$cyl)
carfit <- lm(mpg~cyl,data=mtcars)
summary(carfit)

#

plot(mtcars$mpg~mtcars$cyl,xlab="Cylinders",ylab="MPG")
abline(carfit,lwd=2)

#

carfit <- lm(mpg~factor(cyl),data=mtcars)
summary(carfit)


## 20.5.5 ##

summary(aov(Height~Smoke,data=survey))

#

sqrt(96.27)




##################
### CHAPTER 21 ###
##################

################
# Section 21.1 #
################

# no R code


################
# Section 21.2 #
################

## 21.2.1 ##

# no R code


## 21.2.2 ##

# No R code


## 21.2.3 ##

demo.data <- data.frame(y=c(1.55,0.42,1.29,0.73,0.76,-1.09,1.41,-0.32),x1=c(1.13,-0.73,0.12,0.52,-0.54,-1.15,0.20,-1.09),x2=c(1,0,1,1,0,1,0,1))
demo.data

#

Y <- matrix(demo.data$y) 
Y
n <- nrow(demo.data)
X <- matrix(c(rep(1,n),demo.data$x1,demo.data$x2),nrow=n,ncol=3)
X

#

BETA.HAT <- solve(t(X)%*%X)%*%t(X)%*%Y
BETA.HAT



################
# Section 21.3 #
################

## 21.3.1 ##

demo.fit <- lm(y~x1+x2,data=demo.data)
coef(demo.fit)

#

library("MASS")
survmult <- lm(Height~Wr.Hnd+Sex,data=survey)
summary(survmult)

#

survmult2 <- lm(Height~Wr.Hnd+Sex+Smoke,data=survey)
summary(survmult2)


## 21.3.2 ##

# no R code


## 21.3.3 ##

survcoefs <- coef(survmult)
survcoefs
as.numeric(survcoefs[1]+survcoefs[3])

#

plot(survey$Height~survey$Wr.Hnd,col=c("gray","black")[as.numeric(survey$Sex)],pch=16,xlab="Writing handspan",ylab="Height")
abline(a=survcoefs[1],b=survcoefs[2],col="gray",lwd=2)
abline(a=survcoefs[1]+survcoefs[3],b=survcoefs[2],col="black",lwd=2)
legend("topleft",legend=levels(survey$Sex),col=c("gray","black"),pch=16)


## 21.3.4 ##

confint(survmult2)


## 21.3.5 ##

R2 <- summary(survmult2)$r.squared
R2

#

n <- nrow(survey)-30
n

#

p <- length(coef(survmult2))-1
p

#

n-p-1

#

Fstat <- (R2*(n-p-1))/((1-R2)*p)
Fstat
1-pf(Fstat,df1=p,df2=n-p-1)


## 21.3.6 ##

predict(survmult,newdata=data.frame(Wr.Hnd=16.5,Sex="Male"),interval="confidence",level=0.95)

#

predict(survmult,newdata=data.frame(Wr.Hnd=13,Sex="Female"),interval="prediction",level=0.99)

#

survey[survey$Sex=="Female" & survey$Wr.Hnd==13,c("Sex","Wr.Hnd","Height")]



################
# Section 21.4 #
################

## 21.4.1 ##

x <- seq(-4,4,length=50)
y <- x
y2 <- x + x^2
y3 <- x + x^2 + x^3

#

plot(x,y,type="l")
plot(x,y2,type="l")
plot(x,y3,type="l")

#

plot(mtcars$disp,mtcars$mpg,xlab="Displacement (cu. in.)",ylab="MPG")

#

car.order1 <- lm(mpg~disp,data=mtcars)
summary(car.order1)

#

car.order2 <- lm(mpg~disp+I(disp^2),data=mtcars)
summary(car.order2)

#

car.order3 <- lm(mpg~disp+I(disp^2)+I(disp^3),data=mtcars)
summary(car.order3)

#

plot(mtcars$disp,mtcars$mpg,xlab="Displacement (cu. in.)",ylab="MPG")
abline(car.order1)

#

disp.seq <- seq(min(mtcars$disp)-50,max(mtcars$disp)+50,length=30)

#

car.order2.pred <- predict(car.order2,newdata=data.frame(disp=disp.seq))
lines(disp.seq,car.order2.pred,lty=2)

#

car.order3.pred <- predict(car.order3,newdata=data.frame(disp=disp.seq))
lines(disp.seq,car.order3.pred,lty=3)
legend("topright",lty=1:3,legend=c("order 1 (linear)","order 2 (quadratic)","order 3 (cubic)"))


## 21.4.2 ##

plot(1:1000,log(1:1000),type="l",xlab="x",ylab="",ylim=c(-8,8))
lines(1:1000,-log(1:1000),lty=2)
legend("topleft",legend=c("log(x)","-log(x)"),lty=c(1,2))

#

plot(mtcars$hp,mtcars$mpg,pch=19,col=c("black","gray")[factor(mtcars$am)],xlab="Horsepower",ylab="MPG")
legend("topright",legend=c("auto","man"),col=c("black","gray"),pch=19)

#

car.log <- lm(mpg~log(hp)+am,data=mtcars)
summary(car.log)

#

hp.seq <- seq(min(mtcars$hp)-20,max(mtcars$hp)+20,length=30)
n <- length(hp.seq)
car.log.pred <- predict(car.log,newdata=data.frame(hp=rep(hp.seq,2),am=rep(c(0,1),each=n)))

#

lines(hp.seq,car.log.pred[1:n])
lines(hp.seq,car.log.pred[(n+1):(2*n)],col="gray")


## 21.4.3 ##

# no R code


################
# Section 21.5 #
################

## 21.5.1 ##

# No R code


## 21.5.2 ##

library("faraway")
dia.fit <- lm(chol~age+frame+age:frame,data=diabetes)
summary(dia.fit)

#

dia.coef <- coef(dia.fit)
dia.coef

#

dia.small <- c(dia.coef[1],dia.coef[2])
dia.small
dia.medium <- c(dia.coef[1]+dia.coef[3],dia.coef[2]+dia.coef[5])
dia.medium
dia.large <- c(dia.coef[1]+dia.coef[4],dia.coef[2]+dia.coef[6])
dia.large

#

cols <- c("black","darkgray","lightgray")
plot(diabetes$chol~diabetes$age,col=cols[diabetes$frame],cex=0.5,xlab="age",ylab="cholesterol")
abline(coef=dia.small,lwd=2)
abline(coef=dia.medium,lwd=2,col="darkgray")
abline(coef=dia.large,lwd=2,col="lightgray")
legend("topright",legend=c("small frame","medium frame","large frame"),lty=1,lwd=2,col=cols)


## 21.5.3 ##

warp.fit <- lm(breaks~wool*tension,data=warpbreaks)
summary(warp.fit)

#

nd <- data.frame(wool=c("A","A","B","B"),tension=c("L","H","L","M"))
predict(warp.fit,newdata=nd,interval="confidence",level=0.9)


## 21.5.4 ##

car.fit <- lm(mpg~hp*wt,data=mtcars)
summary(car.fit)


## 21.5.5 ##

library("boot")
nuc.fit <- lm(cost~cap+cum.n*ne*ct,data=nuclear)
summary(nuc.fit)




##################
### CHAPTER 22 ###
##################

################
# Section 22.1 #
################

## 22.1.1 ##

# no R code


## 22.1.2 ##

# no R code



################
# Section 22.2 #
################

## 22.2.1 ##

library("MASS")
survmult <- lm(Height~Wr.Hnd+Sex,data=survey)
survmult2 <- lm(Height~Wr.Hnd+Sex+Smoke,data=survey)
survmult
survmult2

#

anova(survmult,survmult2)

#

library("faraway")
diab <- diabetes[-which(is.na(diabetes$age) | is.na(diabetes$frame)),]

#

dia.model1 <- lm(chol~1,data=diab)
dia.model2 <- lm(chol~age,data=diab)
dia.model3 <- lm(chol~age+frame,data=diab)
dia.model4 <- lm(chol~age*frame,data=diab)

#

anova(dia.model1,dia.model2,dia.model3,dia.model4)


## 22.2.2 ##

library(boot)
nuc.0 <- lm(cost~1,data=nuclear)
summary(nuc.0)

#

add1(nuc.0,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")

#

nuc.1 <- update(nuc.0,formula=.~.+date)
summary(nuc.1)

#

add1(nuc.1,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")

#

nuc.2 <- update(nuc.1,formula=.~.+cap)
add1(nuc.2,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")

#

nuc.3 <- update(nuc.2,formula=.~.+pt)
add1(nuc.3,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")

#

nuc.4 <- update(nuc.3,formula=.~.+ne)

#

add1(nuc.4,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,test="F")

#

summary(nuc.4)


## 22.2.3 ##

nuc.0 <- lm(cost~date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,data=nuclear)
summary(nuc.0)

#

drop1(nuc.0,test="F")

#

nuc.1 <- update(nuc.0,.~.-bw)

#

drop1(nuc.1,test="F")

#

nuc.2 <- update(nuc.1,.~.-pt)
drop1(nuc.2,test="F")

#

nuc.3 <- update(nuc.2,.~.-t1)
drop1(nuc.3,test="F")

#

nuc.4 <- update(nuc.3,.~.-ct)
drop1(nuc.4,test="F")

#

summary(nuc.4)


## 22.2.4 ##

car.null <- lm(mpg~1,data=mtcars)

#

car.step <- step(car.null,scope=.~.+wt*hp*factor(cyl)*disp+am+factor(gear)+drat+vs+qsec+carb)

#

summary(car.step)


## 22.2.5 ##

# no R code



################
# Section 22.3 #
################

## 22.3.1 ##

car.step <- step(car.null,scope=.~.+wt*hp*factor(cyl)*disp+am+factor(gear)+drat+vs+qsec+carb)
plot(car.step,which=1)

#

plot(car.step,which=3,add.smooth=FALSE,id.n=2)

#

gal <- data.frame(d=c(573,534,495,451,395,337,253),h=c(1,0.8,0.6,0.45,0.3,0.2,0.1))
gal.mod1 <- lm(d~h,data=gal)
gal.mod2 <- lm(d~h+I(h^2),data=gal)

#

plot(gal$d~gal$h,xlab="Height",ylab="Distance")
abline(gal.mod1)
plot(gal.mod1,which=1,id.n=0)
plot(gal.mod2,which=1,id.n=0)


## 22.3.2 ##

plot(car.step,which=2)

#

shapiro.test(rstandard(car.step))


## 22.3.3 ##

x <- c(1.1,1.3,2.3,1.6,1.2,0.1,1.8,1.9,0.2,0.75)
y <- c(6.7,7.9,9.8,9.3,8.2,2.9,6.6,11.1,4.7,3)

#

p1x <- 1.2
p1y <- 14
p2x <- 5
p2y <- 19
p3x <- 5
p3y <- 5

#

mod.0 <- lm(y~x)
mod.1 <- lm(c(y,p1y)~c(x,p1x))
mod.2 <- lm(c(y,p2y)~c(x,p2x))
mod.3 <- lm(c(y,p3y)~c(x,p3x))

#

plot(x,y,xlim=c(0,5),ylim=c(0,20))

#

points(p1x,p1y,pch=15,cex=1.5)
abline(mod.0)
abline(mod.1,lty=2)
text(2.5,1,labels="Outlier, low leverage, low influence",cex=1.4)

#

plot(x,y,xlim=c(0,5),ylim=c(0,20))
points(p2x,p2y,pch=15,cex=1.5)
abline(mod.0)
abline(mod.2,lty=2)
text(2.5,1,labels="Not outlier, high leverage, low influence",cex=1.4)

#

plot(x,y,xlim=c(0,5),ylim=c(0,20))
points(p3x,p3y,pch=15,cex=1.5)
abline(mod.0)
abline(mod.3,lty=2)
text(2.5,1,labels="Outlier, high leverage, high influence",cex=1.4)


## 22.3.4 ##

X <- cbind(rep(1,10),x)
hii <- diag(X%*%solve(t(X)%*%X)%*%t(X))
hii
plot(hii~x,ylab="Leverage",main="",pch=4)

#

hatvalues(mod.0)


## 22.3.5 ##

x1 <- c(x,p1x)
y1 <- c(y,p1y)
n <- length(x1)
param <- length(coef(mod.1))
yhat.full <- fitted(mod.1)
sigma <- summary(mod.1)$sigma
cooks <- rep(NA,n)
for(i in 1:n){
  temp.y <- y1[-i]
  temp.x <- x1[-i]
  temp.model <- lm(temp.y~temp.x)
  temp.fitted <- predict(temp.model,newdata=data.frame(temp.x=x1))
  cooks[i] <- sum((yhat.full-temp.fitted)^2)/(param*sigma^2)
}

#

cooks 

#

cooks.distance(mod.1)

#

plot(mod.1,which=4)
plot(mod.2,which=4)
plot(mod.3,which=4)
abline(h=c(1,4/n),lty=2)

#

plot(car.step,which=4)
abline(h=4/nrow(mtcars),lty=2)


## 22.3.6 ##

plot(mod.1,which=5,add.smooth=FALSE,cook.levels=c(4/11,0.5,1))
plot(mod.2,which=5,add.smooth=FALSE,cook.levels=c(4/11,0.5,1))
plot(mod.3,which=5,add.smooth=FALSE,cook.levels=c(4/11,0.5,1))

#

plot(mod.1,which=6,add.smooth=FALSE)
plot(mod.2,which=6,add.smooth=FALSE)
plot(mod.3,which=6,add.smooth=FALSE)

#

plot(car.step,which=5,cook.levels=c(4/nrow(mtcars),0.5,1))
plot(car.step,which=6,cook.levels=c(4/nrow(mtcars),0.5,1))



################
# Section 22.4 #
################

## 22.4.1 ##

# no R code


## 22.4.2 ##

cor(survey$Wr.Hnd,survey$NW.Hnd,use="complete.obs")

#

summary(lm(Height~Wr.Hnd,data=survey))

#

summary(lm(Height~NW.Hnd,data=survey))

#

summary(lm(Height~Wr.Hnd+NW.Hnd,data=survey))


######################################
#### SUGGESTED EXERCISE SOLUTIONS ####
######################################

##########
## 17.1 ##
##########
xbar <- 41.1
#(a)
se <- 11.3/sqrt(6)
se
#(b)
pnorm(55,mean=41.1,sd=se)-pnorm(45,mean=41.1,sd=se)
#(c)
pnorm(32.5,41.1,se)
#(d)
140*0.35
140*(1-0.35) # Both are > 5 so rule-of-thumb says using the normal distribution to represent the sampling distribution of the sample proportion is okay.
#(e)
1-pnorm(0.4,mean=0.35,sd=sqrt(0.35*0.65/140))
#(f)
qnorm(0.9,0.35,sqrt(0.35*0.65/140)) # Upper limit - tail above this value has probability 0.1
qnorm(0.1,0.35,sqrt(0.35*0.65/140)) # Lower limit - tail below this value has probability 0.1. Together, these two limits therefore mark of a central area under the curve of exactly 0.8.
#(g)
# Even though raw data are not normal, sample size is large (n>30 by rule-of-thumb). Standard deviation estimated from sample, so sampling distribution for the mean should be a t-distribution with 63-1=62 df.
se <- 34.51/sqrt(63)
se
#(h)
##(i)
1-pt((40-37.8)/se,df=62)
##(ii)
pt((30-37.8)/se,df=62)
##(iii)
pt((40-37.8)/se,df=62)-0.5

##########
## 17.2 ##
##########
x.bar <- 14.22
sigma <- 2.9
#(a)
x.bar+c(-1,1)*qnorm(0.95)*sigma/sqrt(34) # When sd is known, sampling distribution of sample mean is normal. 90% confident that the true mean run time lies somewhere between 13.40 and 15.04 seconds (rounded 2 d.p.)
#(b)
x.bar+c(-1,1)*qt(0.95,df=33)*sigma/sqrt(34) # When sd is estimated from sample, sampling distribution is t with n-1 df. This means more extreme critical values and wider 90% CI when compared to normal version.
#(c)
p.hat <- 37/400
p.hat
p.hat+c(-1,1)*qnorm(0.995)*sqrt(p.hat*(1-p.hat)/400) # 99% confident that the true proportion of left-handedness only is somwhere between 0.055 and 0.130 (rounded 3 d.p.)
#(d)
p.hat <- (37+11)/400
p.hat
p.hat+c(-1,1)*qnorm(0.995)*sqrt(p.hat*(1-p.hat)/400) # 99% confident that the true proportion of left-handed or ambidextrous citizens is somwhere between 0.078 and 0.162 (rounded 3 d.p.)
#(e)
ci.mat <- matrix(NA,5000,3)
n <- 300
lambda.e <- 0.1
mu <- 1/lambda.e
for(i in 1:5000){
  samp <- rexp(n,rate=lambda.e)
  samp.ci <- mean(samp)+c(-1,1)*qt(0.975,n-1)*sd(samp)/sqrt(n)
  ci.mat[i,1:2] <- samp.ci
  ci.mat[i,3] <- mu>=samp.ci[1] && mu<=samp.ci[2]
}
mean(ci.mat[,3])
#(f)
plot(ci.mat[1,1:2],c(1,1),xlim=range(ci.mat[,1:2]),ylim=c(1,100),type="l",xlab="",ylab="")
for(i in 2:100){
  lines(ci.mat[i,1:2],c(i,i))
}
abline(v=mu,col=2)

##########
## 18.1 ##
##########
#(a)
### H0: mu = 3.5; HA: mu != 3.5 (two-sided test)
tstat <- (3.97-3.5)/(2.21/sqrt(73))
tstat
pt(-tstat,df=72)+(1-pt(tstat,df=72))
### p-value is around 0.073, this is > than alpha=0.05, therefore insufficient evidence to reject the null hypothesis. There is no evidence that the true mean cat weight is different to 3.5.
#(b)
# H0: mu = 4.3; HA: mu > 4.3 (one-sided test)
t.test(quakes$mag,mu=4.3,alternative="greater",conf.level=0.99)
# p-value very small; strong evidence to reject the null. There is evidence to suggest that the true mean magnitude is greater than 4.3.
#(c)
mean(quakes$mag)+c(-1,1)*qt(0.995,df=999)*sd(quakes$mag)/sqrt(1000)

##########
## 18.2 ##
##########
#(a)
library("MASS")
?anorexia
t.test(anorexia[,3],anorexia[,2],alternative="greater",paired=TRUE)
# p-value ~0.0023. Less than 0.05; some evidence to reject H0. There is evidence to suggest that the mean post-weight is greater than the mean pre-weight.
#(b)
t.test(anorexia[anorexia$Treat=="Cont",3],anorexia[anorexia$Treat=="Cont",2],alternative="greater",paired=TRUE)
t.test(anorexia[anorexia$Treat=="CBT",3],anorexia[anorexia$Treat=="CBT",2],alternative="greater",paired=TRUE)
t.test(anorexia[anorexia$Treat=="FT",3],anorexia[anorexia$Treat=="FT",2],alternative="greater",paired=TRUE)
# There is no statistical evidence to reject the claim that there is no difference between the pre- and post-weight means in the control group, mild evidence to reject at the 5% level for the CBL treatment, and strong evidence to reject in favor of HA for the FT treatment. It would seem that the FT treatment is the most effective based on these data.
#(c)
?PlantGrowth
control <- PlantGrowth$weight[PlantGrowth$group=="ctrl"]
treated <- PlantGrowth$weight[PlantGrowth$group!="ctrl"]
# H0: mu_control - mu_treated = 0; HA: mu_control - mu_treated < 0
max(c(sd(control),sd(treated)))/min(c(sd(control),sd(treated)))
# Ratio of (large sd) / (small sd) is less than 2 so use pooled variance according to rule-of-thumb.
#(d)
t.test(x=control,y=treated,alternative="less",var.equal=TRUE)
# Large p-value ~0.41. There is no evidence to reject H0. There is insufficient evidence to conclude that the mean treated weight is more than the mean untreated weight.
#(e)
myt.test <- function(x,y,paired=FALSE,var.equal=FALSE,...){	
  if(!paired){
    sdx <- sd(x)
    sdy <- sd(y)
    sd.big <- max(c(sdx,sdy))
    sd.small <- min(c(sdx,sdy))
    var.equal <- (sd.big/sd.small)<2
  }
  return(t.test(x=x,y=y,paired=paired,var.equal=var.equal,...))
}
#(f)
# Example 1
snacks <- c(87.7,80.01,77.28,78.76,81.52,74.2,80.71,79.5,77.87,81.94,80.7,82.32,
            75.78,80.19,83.91,79.4,77.52,77.62,81.4,74.89,82.95,73.59,77.92,77.18,
            79.83,81.23,79.28,78.44,79.01,80.47,76.23,78.89,77.14,69.94,78.54,79.7,
            82.45,77.29,75.52,77.21,75.99,81.94,80.41,77.7)
snacks2 <- c(80.22,79.73,81.1,78.76,82.03,81.66,80.97,81.32,80.12,78.98,79.21,
             81.48,79.86,81.06,77.96,80.73,80.34,80.01,81.82,79.3,79.08,79.47,
             78.98,80.87,82.24,77.22,80.03,79.2,80.95,79.17,81)
myt.test(x=snacks2,y=snacks,alternative="greater",conf.level=0.9)
# Example 2
men <- c(102,87,101,96,107,101,91,85,108,67,85,82)
women <- c(73,81,111,109,143,95,92,120,93,89,119,79,90,126,62,92,77,106,105,111)
myt.test(x=men,y=women,alternative="two.sided",conf.level=0.95)
# Example 3
rate.before <- c(52,66,89,87,89,72,66,65,49,62,70,52,75,63,65,61) 
rate.after <- c(51,66,71,73,70,68,60,51,40,57,65,53,64,56,60,59) 
myt.test(x=rate.after,y=rate.before,alternative="less",paired=TRUE,conf.level=0.95)

##########
## 18.3 ##
##########
#(a)
# H0: p=0.9; HA: p<0.9
n <- 89
n*0.9
n*0.1
# Both np and n(p-1) > 5 so OK to continue with normal distribution according to the rule-of-thumb.
#(b)
p.hat <- 71/n
p.hat
Z <- (p.hat-0.9)/(sqrt(0.9*0.1/n))
Z
pnorm(Z)
prop.test(x=71,n=n,p=0.9,alternative="less",conf.level=0.9,correct=FALSE)
# p-value very small; less than 0.1. There is evidence to reject H0 and conclude the true proportion of women who would recommend in samples of size 89 is less than 0.9.
#(c)
p.hat+c(-1,1)*qnorm(0.95)*sqrt(p.hat*(1-p.hat)/n)
#(d)
x1 <- 97
n1 <- 445
p.hat1 <- x1/n1
p.hat1
x2 <- 90
n2 <- 419
p.hat2 <- x2/n2
p.hat2
p.star <- (x1+x2)/(n1+n2)
p.star
Z <- (p.hat2-p.hat1)/sqrt(p.star*(1-p.star)*(1/n1+1/n2))
Z
2*pnorm(Z)
# p-value is very large; much greater than 0.05. No evidence to reject H0. Retain H0 and conclude there is no evidence to suggest the proportion of citizens who support decriminalization varies between the two states.
#(e)
(p.hat2-p.hat1)+c(-1,1)*qnorm(0.975)*sqrt(p.star*(1-p.star)*(1/n1+1/n2))
# We are 95% confident that the true difference in the proportion of support between State 2 and State 1 lies somewhere between -0.058 and 0.051. CI includes zero; reflects the same result as the hypothesis (no evidence of a difference).
#(f)
Z.test <- function(p1,n1,p2=NULL,n2=NULL,p0,alternative="two.sided",conf.level=0.95){
  if(is.null(p2)||is.null(n2)){
    cat("One-sample Z-test.\n")
    if(p1*n1<=5||n1*(1-p1)<=5){
      warning("Normal distribution may not be valid; np or n(1-p) <= 5 detected.")
    }
    Z <- (p1-p0)/sqrt(p0*(1-p0)/n1)
    CI <- (p1)+c(-1,1)*qnorm(conf.level+(1-conf.level)/2)*sqrt(p0*(1-p0)/n1)
  } else {
    cat("Two-sample Z-test.\n")
    if(p1*n1<=5||n1*(1-p1)<=5||p2*n2<=5||n2*(1-p2)<=5){
      warning("Normal distribution may not be valid; np or n(1-p) <= 5 detected.")
    }
    p.star <- (p1*n1+p2*n2)/(n1+n2)
    Z <- (p1-p2-p0)/sqrt(p.star*(1-p.star)*(1/n1+1/n2))
    CI <- sort((p1-p2)+c(-1,1)*qnorm(conf.level+(1-conf.level)/2)*sqrt(p.star*(1-p.star)*(1/n1+1/n2)))
  }
  
  P <- pnorm(Z)
  if(alternative=="greater"){
    P <- 1-P
  } else if(alternative=="two.sided"){
    if(Z<0){
      P <- 2*P
    } else {
      P <- 2*(1-P)
    }
  }
  return(list(Z=Z,P=P,CI=CI))
}
#(g)
# Example 1
sick <- c(0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1)
Z.test(p1=mean(sick),n1=length(sick),p0=0.2,alternative="two.sided",conf.level=0.95)
# Example 2
x1 <- 180
n1 <- 233
p.hat1 <- x1/n1
x2 <- 175
n2 <- 197
p.hat2 <- x2/n2
Z.test(p.hat2,n2,p.hat1,n1,p0=0,alternative="greater",conf.level=0.95) # ...or you could flip the order of differencing and use alternative="less"
#(h)
Z.test(p1=0.11,n1=10,p0=0.1)

##########
## 18.4 ##
##########
#(a) H0: No relationship between hair and eye color; HA: There is a relationship.
?HairEyeColor
chisq.test(x=HairEyeColor[,,1]+HairEyeColor[,,2])
# Very small P-value. Very strong evidence against the null. Reject H0 and conclude there does appear to be a relationship between hair and eye color of statistics students.
#(b) H0: p1=p2=p3=1/3; HA: H0 is incorrect
library("car")
data(Duncan)
?Duncan
jobtype <- Duncan$type
jobtype.tab <- table(jobtype)
jobtype.tab
chisq.test(jobtype.tab)
##(i) With a significance level of 0.05 and a p-value of 0.015, there is weak evidence to reject H0 and we therefore conclude that the three job types do not appear to be unifomly represented in the data set.
##(ii) With a significance level of 0.01 and a P-value of 0.015, there is no evidence to reject H0 and we therefore conclude that the three job types appear to be unifomly represented in the data set. 

##########
## 18.5 ##
##########
#(a)
typeI.mean <- function(mu0,sigma,n,alpha,test="two.sided",ITERATIONS=10000){
  tstats <- rep(NA,ITERATIONS)
  for(i in 1:ITERATIONS){
    temporary.sample <- rnorm(n=n,mean=mu0,sd=sigma)
    temporary.mean <- mean(temporary.sample)
    temporary.sd <- sd(temporary.sample)
    tstats[i] <- (temporary.mean-mu0)/(temporary.sd/sqrt(n))
  }
  pvals <- pt(tstats,df=n-1)
  if(test=="less"){
    return(mean(pvals<alpha))
  } else if(test=="greater"){
    return(mean((1-pvals)<alpha))
  } else if(test=="two.sided"){
    result <- pvals
    result[tstats>0] <- 1-pvals[tstats>0]
    return(mean(result<alpha/2))
  } else {
    stop("'test' argument not recognised")
  }
}
##(i)
typeI.mean(mu0=0,sigma=1,n=40,alpha=0.05,test="less")
typeI.mean(mu0=0,sigma=1,n=40,alpha=0.05,test="greater")
typeI.mean(mu0=0,sigma=1,n=40,alpha=0.05,test="two.sided")
##(ii)
typeI.mean(mu0=-4,sigma=0.3,n=60,alpha=0.01,test="less")
typeI.mean(mu0=-4,sigma=0.3,n=60,alpha=0.01,test="greater")
typeI.mean(mu0=-4,sigma=0.3,n=60,alpha=0.01,test="two.sided")
#(b)
typeII.mean <- function(mu0,muA,sigma,n,alpha,test="two.sided",ITERATIONS=10000){
  tstats <- rep(NA,ITERATIONS)
  for(i in 1:ITERATIONS){
    temporary.sample <- rnorm(n=n,mean=muA,sd=sigma)
    temporary.mean <- mean(temporary.sample)
    temporary.sd <- sd(temporary.sample)
    tstats[i] <- (temporary.mean-mu0)/(temporary.sd/sqrt(n))
  }
  pvals <- pt(tstats,df=n-1)
  if(test=="less"){
    return(mean(pvals>=alpha))
  } else if(test=="greater"){
    return(mean((1-pvals)>=alpha))
  } else if(test=="two.sided"){
    result <- pvals
    result[tstats>0] <- 1-pvals[tstats>0]
    return(mean(result>=alpha/2))
  } else {
    stop("'test' argument not recognised")
  }
}
##(i)
typeII.mean(mu0=-3.2,muA=-3.3,sigma=0.1,n=25,alpha=0.05,test="two.sided")
##(ii)
typeII.mean(mu0=8994,muA=5600,sigma=3888,n=9,alpha=0.01,test="less")
##(iii)
typeII.mean(mu0=0.44,muA=0.4,sigma=2.4,n=68,alpha=0.05,test="greater")

##########
## 18.6 ##
##########
#(a)
power.mean <- function(nvec,...){
  nlen <- length(nvec)
  result <- rep(NA,nlen)
  pbar <- txtProgressBar(min=0,max=nlen,style=3)
  for(i in 1:nlen){
    result[i] <- 1-typeII.mean(n=nvec[i],...)
    setTxtProgressBar(pbar,i)
  }
  close(pbar)
  return(result)
}
##(i)
power.mean(nvec=50,mu0=10,muA=10.5,sigma=0.9,alpha=0.01,test="two.sided")
##(ii)
power.mean(nvec=44,mu0=80,muA=78.5,sigma=3.1,alpha=0.05,test="less") #Yes, seems statistically powerful.
power.mean(nvec=44,mu0=80,muA=78.5,sigma=3.1,alpha=0.01,test="less") #No, power appears less than 80%, but only just...
#(b)
sample.sizes <- 5:100
pow <- power.mean(nvec=sample.sizes,mu0=80,muA=78.5,sigma=3.1,alpha=0.05,test="less")
minimum.n <- sample.sizes[min(which(pow>=0.8))]
minimum.n
pow2 <- power.mean(nvec=sample.sizes,mu0=80,muA=78.5,sigma=3.1,alpha=0.01,test="less")
minimum.n2 <- sample.sizes[min(which(pow2>=0.8))]
minimum.n2
plot(sample.sizes,pow,xlab="sample size n",ylab="simulated power")
points(sample.sizes,pow2,col="grey")
abline(h=0.8,lty=2)
abline(v=c(minimum.n,minimum.n2),lty=3,col=c("black","grey"))
legend("bottomright",legend=c("alpha=0.05","alpha=0.01"),col=c("black","grey"),pch=1)

##########
## 19.1 ##
##########
mim <- c(93,120,65,105,115,82,99,87,100,90,78,95,93,88,110,85,45,80,28,75,70,65,55,50,40,100,75,65,40,73,65,50,30,45,50,45,55,96,58,95,90,65,80,85,95,82)
site <- c(rep("I",15),rep("II",10),rep("III",12),rep("IV",9))
#(a)
boxplot(mim~site)
mim.means <- tapply(mim,INDEX=site,FUN=mean)
points(1:4,mim.means,pch=4)
#(b)
mim.meancen <- mim-rep(mim.means,table(site))
qqnorm(mim.meancen)
qqline(mim.meancen) # Seems consistent with normality based on a visual inspection.
sds <- tapply(mim,INDEX=site,FUN=sd)
max(sds)/min(sds) # Less than 2 so variances can be assumed equal according to the rule-of-thumb.
#(c)
summary(aov(mim~site)) # Small p-value. Very strong evidence against H0. There is evidence to reject the null and conclude there is a difference in the mean depths of the finds across the four sites.
#(d)
m <- tapply(iris$Sepal.Length,iris$Species,mean)
mc <- iris$Sepal.Length-m[as.numeric(iris$Species)]
qqnorm(mc)
qqline(mc) # Looks approximately normal.
tapply(iris$Sepal.Length,iris$Species,sd) # max(sd)/min(sd) < 2 so variances may be assumed equal.
m <- tapply(iris$Sepal.Width,iris$Species,mean)
mc <- iris$Sepal.Width-m[as.numeric(iris$Species)]
qqnorm(mc)
qqline(mc) # Looks approximately normal.
tapply(iris$Sepal.Width,iris$Species,FUN=sd) # max(sd)/min(sd) < 2 so variances may be assumed equal.
m <- tapply(iris$Petal.Length,iris$Species,mean)
mc <- iris$Petal.Length-m[as.numeric(iris$Species)]
qqnorm(mc)
qqline(mc) # Looks approximately normal; some deviation though.
tapply(iris$Petal.Length,iris$Species,sd) # max(sd)/min(sd) > 2 so variances may not be assumed equal.
m <- tapply(iris$Petal.Width,iris$Species,mean)
mc <- iris$Petal.Width-m[as.numeric(iris$Species)]
qqnorm(mc)
qqline(mc) # Looks approximately normal.
tapply(iris$Petal.Width,iris$Species,FUN=sd) # max(sd)/min(sd) > 2 so variances may not be assumed equal.
#(e)
summary(aov(Sepal.Length~Species,data=iris))
summary(aov(Sepal.Width~Species,data=iris))
# Both p-values are very small. Strong evidence to reject H0 and conclude that the mean sepal lengths and widths do vary according to species.

##########
## 19.2 ##
##########
#(a)
depth.fac <- cut(quakes$depth,breaks=c(0,200,400,680))
table(depth.fac)
#(b)
m <- tapply(quakes$stations,depth.fac,mean)
mc <- quakes$stations-m[as.numeric(depth.fac)]
hist(mc)
qqnorm(mc)
qqline(mc) # Data appear non-normal... Kruskal-Wallis preferred over parametric one-way ANOVA
#(c)
kruskal.test(quakes$stations~depth.fac) # P-value > 0.01 (just barely); retain null. Minimal evidence, at best, of a difference in median number of detecting stations according to depth.fac categories.
#(d)
library("MASS")
cars.means <- aggregate(Cars93$Length,by=list(Cars93$AirBags,Cars93$Man.trans.avail),FUN=mean)
cars.means
#(e)
interaction.plot(x.factor=cars.means[,1],trace.factor=cars.means[,2],respons=cars.means$x,trace.label="Manual avail.",xlab="Airbags",ylab="Mean length") # There is some visual indication of interactive behavior owing to the non-parallel nature of the two lines; but rememeber there are no measures of variability of the various group means displayed on this plot...
#(f)
summary(aov(Length~AirBags+Man.trans.avail+AirBags:Man.trans.avail,data=Cars93)) # No formal statistical evidence of an interactive effect. Strong evidence of main effects, however -- both 'Airbags' and 'Man.trans.avail' appear to be related to car length.

##########
## 20.1 ##
##########
library("MASS")
survfit <- lm(Height~Wr.Hnd,data=survey)
#(a)
predict(survfit,newdata=data.frame(Wr.Hnd=c(12,15.2,17,19.9)),interval="confidence",level=0.99)
#(b)
incomplete.obs <- which(is.na(survey$Height)|is.na(survey$Wr.Hnd))
rho.xy <- cor(survey$Wr.Hnd,survey$Height,use="complete.obs")
b1 <- sd(survey$Height[-incomplete.obs])/sd(survey$Wr.Hnd[-incomplete.obs])*rho.xy
b1
b0 <- mean(survey$Height[-incomplete.obs])-b1*mean(survey$Wr.Hnd[-incomplete.obs])
b0
#(c)
##(i)
plot(survey$Height~survey$Pulse,xlab="Pulse rate (bpm)",ylab="Height (cm)")
survfit <- lm(Height~Pulse,data=survey)
survfit # Model equation is  y = 177.86 - 0.072x
abline(survfit,lwd=2)
##(ii)
summary(survfit) # For each additional bpm, the mean student height is estimated to decrease by 0.072cm; p-value for slope is 0.275. No evidence to reject H0. Insufficient evidence to conclude that pulse rate affects mean student height.
confint(survfit,level=0.9)
##(iii)
xseq <- data.frame(Pulse=seq(30,110,length=100))
survfit.ci <- predict(survfit,newdata=xseq,interval="confidence",level=0.9)
survfit.pi <- predict(survfit,newdata=xseq,interval="prediction",level=0.9)
lines(xseq[,1],survfit.ci[,2],lty=2)
lines(xseq[,1],survfit.ci[,3],lty=2)
lines(xseq[,1],survfit.pi[,2],lty=2,col="grey")
lines(xseq[,1],survfit.pi[,3],lty=2,col="grey")
legend("topleft",legend=c("90% CI","90% PI"),lty=2,col=c("black","grey"))
##(iv)
incomplete.obs <- which(is.na(survey$Height)|is.na(survey$Pulse))
abline(h=mean(survey$Height[-incomplete.obs]),col=2,lty=3,lwd=3) # The line sits in the middle of the CI bands without breaching them. This supports the conclusion that pulse rate is not significantly related to mean student height.
#(d)
?mtcars
plot(mtcars$mpg~mtcars$wt,xlab="Weight (lbs/1000)",ylab="MPG")
#(e)
carfit <- lm(mpg~wt,data=mtcars)
carfit
abline(carfit,lwd=2)
#(f)
summary(carfit) # mean MPG = 37.28 - 5.34*weight # For each extra 1000lbs of weight, the mean MPG decreases by 5.34; p-value for slope is very small---result is statistically significant---strong evidence to suggest that mean MPG changes according to weight (for cars of this era).
#(g)
predict(carfit,newdata=data.frame(wt=6),interval="prediction",level=0.95) # Predicting at 6000lbs seems untrustworthy. Extrapolation is far enough outside the range of the observed data such that the associated PI has a lower limit that is negative, which makes no sense in terms of the response variable of MPG.

##########
## 20.2 ##
##########
library("MASS")
#(a)
table(survey$Exer)
boxplot(survey$Height~survey$Exer)
#(b)
survfit <- lm(Height~Exer,data=survey)
summary(survfit) # The reference level of the predictor defaults to the first level of the factor, which in this case (as is the default in R) is alphabetically arranged to be 'Freq'.
#(c)
# We observe that both of the levels for which we've obtained coefficient estimates yield p-values that suggest evidence the coefficients are different to zero.
# The coefficient corresponding to 'some' has the smallest p-value of the two additive dummy levels.
# The negative point estimates of both estimates tell us that the model predicts the effect on height of being in either the 'none' or the 'some' categories, when compared to the 'frequent' category, is one of a decrease. In other words, it appears those who exercise less than 'frequent' are shorter on average.
# The shortest mean height is reserved for those in the 'none' exercise category; the estimated coefficient (-5.58) is more extreme than the coefficient for 'some' (-4.21).
# Overall statistical significance of the predictor (in terms of the effect of exercise on height) is supported by the global (omnibus F-test) P-value of 0.0035.
#(d)
one.of.each <- factor(levels(survey$Exer))
one.of.each
predict(survfit,newdata=data.frame(Exer=one.of.each),interval="prediction")
#(e)
summary(aov(Height~Exer,data=survey)) # Same 'global' P-value as the lm model summary. There is evidence to suggest mean student height differs according to exercise frequency.
#(f)
ExerReordered <- relevel(survey$Exer,ref="None")
levels(ExerReordered)
summary(aov(Height~ExerReordered,data=survey)) # There is no change to the omnibus F-test if we reorder the the reference level of exercise to be 'none', and nor should we expect there to be. The global test of a difference between the means doesn't depend on what we set the baseline value of the response to be.
#(g)
carfit <- lm(qsec~gear,data=mtcars)
summary(carfit) # The effect of 'gear' when treated as a continuous variable is interpreted as a decrease in quarter-mile time of around 0.5 seconds, on average, for each additional forward gear. However, there is no evidence that this effect is different to zero, with a high P-value of 0.243 (and a similar global p-value)
#(h)
carfit2 <- lm(qsec~factor(gear),data=mtcars)
summary(carfit2) # The effect of 'gear' when treated as a categorical variable now appears to be statistically significant. Having 4 gears (as opposed to the reference level of 3) seems to increase the mean quarter-mile time by 1.27 seconds, having 5 gears appears to decrease the mean quarter-mile time by 2.05 seconds. The global (omnibus F-test) p-value is also quite small, yielding evidence in support of an effect of 'gear' on 'qsec'.
#(i)
plot(mtcars$qsec~mtcars$gear,xlab="No. of Gears",ylab="Quarter-mile Time")
abline(carfit,lwd=2) # The plot indicates clearly that the difference between the two models is due to the fact that the relationship cannot be well explained by a continuous straight line. The first model therefore is incapable of realistically capturing the effect of changing categories in 'gear' on 'qsec'.

##########
## 21.1 ##
##########
library("MASS")
?cats
#(a)
plot(cats$Bwt,cats$Hwt,col=cats$Sex,xlab="Body weight (kg)",ylab="Heart weight (g)")
legend("topleft",legend=c("female","male"),col=c(1,2),pch=c(1,1)) # Females are black, since the levels of the factor vector cats$Sex are in the alphabetical order of "F" "M", which R interprets as 1 and 2 when this factor vector is passed to the 'col' argument of 'plot'.
#(b)
cats.fit <- lm(Hwt~Bwt+Sex,data=cats)
summary(cats.fit)
##(i)
# "Mean heart weight" = -0.415 + 4.076*"Body weight" - 0.082*"is male"
# For cats of the same sex, the effect of an additional kg of body weight is, on average, an extra 4.076 grams of heart weight. For cats of the same body weight, the heart weight of a male is, on average, 0.082 grams lighter than that of a female.
# The model states the effect of body weight is highly statistically significant -- there is evidence to believe body weight does indeed affect heart weight. However, the effect of sex is not significant. There is no statistical evidence to suggest the coefficient for "is male" is any different to zero (when also adjusting for body weight).
# The above notes imply that the inclusion of "sex" as a predictor is statistically unnecessary when it comes to modeling the response for these data.
##(ii)
names(summary(cats.fit))
summary(cats.fit)$r.squared
# The coefficient of determination, 'R-squared', shows that for your fitted model, 64.5% of the variation in heart weight is able to be captured by the regression.
summary(cats.fit)$fstatistic
1-pf(129.1056,df1=2,df2=141)
# Reading from the summary output, or running the line above, the result of the omnibus F-test is a tiny p-value; effectively zero. This suggests very strong evidence against the null hypothesis (H0 is that modeling heart weight isn't improved by taking body weight and sex into account).
#(c)
predict(cats.fit,newdata=data.frame(Bwt=3.4,Sex="F"),interval="prediction",level=0.95)
#(d)
plot(cats$Bwt,cats$Hwt,col=cats$Sex,xlab="Body weight (kg)",ylab="Heart weight (g)")
legend("topleft",legend=c("female","male"),col=c(1,2),pch=c(1,1))
Bwt.seq <- seq(min(cats$Bwt)-0.5,max(cats$Bwt)+0.5,length=30)
n <- length(Bwt.seq)
cats.pred <- predict(cats.fit,newdata=data.frame(Bwt=rep(Bwt.seq,2),Sex=rep(c("M","F"),each=n)))
lines(Bwt.seq,cats.pred[1:n],col=2) 
lines(Bwt.seq,cats.pred[(n+1):(2*n)])
# The two superimposed lines are positively sloped according to the coefficient for "Bwt", but extremely close together, mirroring the minimal impact (and lack of statistical significance) of "Sex".
#(e)
library("boot")
?nuclear
pairs(nuclear)
#(f)
nuc.fit1 <- lm(cost~t1+t2,data=nuclear)
summary(nuc.fit1)
#(g)
nuc.fit2 <- lm(cost~t1+t2+date,data=nuclear)
summary(nuc.fit2)
# By including "date" in the linear model, this completely removes the statistical significance of "t1" as seen in the previous model. In fact, the coefficient for "t1" changes signs! What this implies is the previous positive, significant relationship between "t1" and "cost" is actually explained more by "date", owing to the positive correlation of "t1" with "date", and it is "date" that should probably be used in a fitted model. The "t2" predictor remains non-significant, albeit with a reduced p-value in this latest model.
#(h)
nuc.fit3 <- lm(cost~date+cap+ne,data=nuclear)
summary(nuc.fit3) # Fitted model is "cost" = -6458 + 95.4*"date of permit issue" + 0.42*"capacity" + 126.1*"constructed in north-east"
confint(nuc.fit3) # All intervals exclude null value of zero, reflecting their significance in the model summary.
#(i)
detroit <- data.frame(Murder=c(8.6,8.9,8.52,8.89,13.07,14.57,21.36,28.03,31.49,37.39,46.26,47.24,52.33),Police=c(260.35,269.8,272.04,272.96,272.51,261.34,268.89,295.99,319.87,341.43,356.59,376.69,390.19),Unemploy=c(11,7,5.2,4.3,3.5,3.2,4.1,3.9,3.6,7.1,8.4,7.7,6.3),Guns=c(178.15,156.41,198.02,222.1,301.92,391.22,665.56,1131.21,837.6,794.9,817.74,583.17,709.59))
detroit
pairs(detroit)
# The number of police seems to be the single most telling variable for prediction of murder numbers.
#(j)
murd.fit <- lm(Murder~Police+Unemploy+Guns,data=detroit)
summary(murd.fit)
summary(murd.fit)$r.squared
# "Mean murders" = -68.85 + 0.281*"no. of police" + 0.147*"unemployment" + 0.014*"no. of gun licenses".
# After adjusting for "no. of gun licenses", and "unemployment", each additional police per 100000 population is related to a mean increase of 0.28 murders per 100000 population.
# After adjusting for "no. of police", and "unemployment", each additional gun license per 100000 population results in a mean increase of 0.014 murders per 100000 population.
# After adjusting for "no. of gun licenses", and "no. of police", each additional percentage of unemployment results in a mean increase of 0.147 murders per 100000 population.
# No, it doesn't make sense to claim that *any* of the relationships are causal, particularly based only one data set and analysis. Causality is extremely difficult to prove in general. In this case, the idea that having a larger police force 'causes' more murders, for example, is rather difficult to justify.
#(k)
summary(murd.fit)$r.squared
# Approx. 97.67% of the variability in mean murder numbers is explained by the three-predictor model (this can also be seen from the model summary).
murd.fit2 <- lm(Murder~Police+Guns,data=detroit)
summary(murd.fit2)
summary(murd.fit2)$r.squared
## The coefficient of determination has barely changed from before; approx. 97.63% of the variability in mean murder numbers is now explained by the two-predictor model. Removing the unemployment predictor variable has had little discernable impact on the model's ability to explain the variation in the murder numbers. This mirrors the non-significant nature of "unemployment rate" when it is included in the model---non-significance implies there is no evidence to suggest varying the unemployment rate will change the mean response.
#(l)
predict(murd.fit2,newdata=data.frame(Police=c(300,300),Guns=c(500,0)),interval="confidence",level=0.99)

##########
## 21.2 ##
##########
#(a)
gal <- data.frame(d=c(573,534,495,451,395,337,253),h=c(1,0.8,0.6,0.45,0.3,0.2,0.1))
plot(gal$d~gal$h,pch=19,xlab="Initial height",ylab="Distance traveled")
#(b)
##(i)
gal.fit.order2 <- lm(d~h+I(h^2),data=gal)
summary(gal.fit.order2)
##(ii)
gal.fit.order3 <- lm(d~h+I(h^2)+I(h^3),data=gal)
summary(gal.fit.order3)
gal.fit.order4 <- lm(d~h+I(h^2)+I(h^3)+I(h^4),data=gal)
summary(gal.fit.order4)
# These models reveal that the order 3 model is significant in it's highest-order term, and the fit is improved in terms of the coefficient of determination. The same cannot be said for the order 4 model. 
#(c)
# Based on the above, out of the three fitted models, the cubic function in height seems preferable -- in other words, the relationship between "distance traveled" and "initial height" therefore appears cubic -- the quadratic model seems too simple, and the quartic model seems uneccessarily complex.
hseq <- seq(0.05,1.05,length=30)
gal.pred <- predict(gal.fit.order3,newdata=data.frame(h=hseq),interval="confidence",level=0.9)
lines(hseq,gal.pred[,1])
lines(hseq,gal.pred[,2],lty=2)
lines(hseq,gal.pred[,3],lty=2)
#(d)
library("faraway")
?trees
plot(trees$Volume~trees$Girth,pch=19,xlab="Girth",ylab="Volume")
#(e)
tree.fit1 <- lm(Volume~Girth+I(Girth^2),trees)
summary(tree.fit1) ## "Mean volume" = 10.79 - 2.09*"girth" + 0.254*"girth^2"
tree.fit2 <- lm(log(Volume)~log(Girth),trees)
summary(tree.fit2) ## "Mean log(volume)" = -2.35 + 2.20*"log(girth)"
# Coefficients of determination are similar; the quadratic model is slightly higher. Both indicate a statistically significant positive effect of girth on volume. Based on the F-test, both models are clearly better than fitting an intercept alone (i.e. girth does appear to be able to explain (mean) volume very well).
#(f)
gseq <- seq(7,21,length=30)
tree.pred1 <- predict(tree.fit1,newdata=data.frame(Girth=gseq),interval="prediction")
tree.pred2 <- predict(tree.fit2,newdata=data.frame(Girth=gseq),interval="prediction")
plot(trees$Volume~trees$Girth,pch=19,xlab="Girth",ylab="Volume")
lines(gseq,tree.pred1[,1],lwd=2)
lines(gseq,tree.pred1[,2])
lines(gseq,tree.pred1[,3])
lines(gseq,exp(tree.pred2[,1]),lwd=2,lty=2)
lines(gseq,exp(tree.pred2[,2]),lty=2)
lines(gseq,exp(tree.pred2[,3]),lty=2)
legend("topleft",legend=c("quadratic","logged"),lty=1:2)
# The fitted values of the models themselves are extremely similar. However, the prediction intervals tell a different story. Notably, the quadratic model has far wider limits for small girth values than for larger ones. On the other hand, the limits for the logged model are substantially wider than those of the quadratic model at larger girth values. Which model is 'better'...? It's very difficult to answer that without further information...
#(g)
library("MASS")
car.fit <- lm(mpg~wt+hp+disp,data=mtcars)
summary(car.fit)
#(h)
car.fit <- lm(I(1/mpg)~wt+hp+disp,data=mtcars)
summary(car.fit)
# Both fits to the mtcars data here provide similar levels of significance for the three predictors; though there is a mild yet noticeable improvement in the coefficient of determination for the latter model based on a response of GPM = 1/MPG.

##########
## 21.3 ##
##########
library("MASS")
#(a)
cat.fit <- lm(Hwt~Bwt*Sex,data=cats)
summary(cat.fit)
# The main-effects-only version of the model had a mild negative effect of "sex male", and it was not significant. In this version, the effect of being male is more extreme, and the p-value is far smaller, now providing weak evidence of significance. The interactive term for the slope of 'Bwt' for a male is reduced somewhat compared to females (estimated parameter is negative).
#(b)
plot(cats$Bwt,cats$Hwt,col=cats$Sex,ylab="Heart weight (g)",xlab="Body weight (kg)")
legend("topleft",legend=c("Female","Male"),col=1:2,pch=1)
cat.coefs <- coef(cat.fit)
abline(coef=cat.coefs[1:2])
abline(coef=c(sum(cat.coefs[c(1,3)]),sum(cat.coefs[c(2,4)])),col=2)
# Lines of the fitted model are no longer parallel; the effect of the weakly significant interaction is apparent.
#(c)
predict(cat.fit,newdata=data.frame(Bwt=3.4,Sex="F"),interval="prediction",level=0.95)
# Sigma's heart weight predicted from the new model is around 1.5 grams lighter than predicted from the main-effects-only model in the earlier exercise. The prediction interval is set accordingly lower as well, but is also wider than the interval from earlier.
#(d)
library("faraway")
tree.fit1 <- lm(Volume~Girth+Height,data=trees)
summary(tree.fit1)
tree.fit2 <- lm(Volume~Girth*Height,data=trees)
summary(tree.fit2)
#(e)
tree.fit3 <- lm(log(Volume)~log(Girth)+log(Height),data=trees)
summary(tree.fit3)
tree.fit4 <- lm(log(Volume)~log(Girth)*log(Height),data=trees)
summary(tree.fit4)
# The interactive effect is highly significant in the untransformed model from (d), but completely non-significant after log-transformation of all present variables. This suggests that 'straight-line' relationships are not the best way to model these data (you can experiment with plots if you wish), and that we must account for curvature in the response surface by either working with transformed data or including the two-way interaction between the two untransformed continuous predictors. Once more, it is difficult to decide on which approach ought to be preferred -- we need to know more about the nature of the data themselves in context, as well as the ultimate purpose of the fitted model.
#(f)
car.fit <- lm(mpg~factor(cyl)*hp+wt,data=mtcars)
summary(car.fit)
#(g)
coef(car.fit)
# The interactive effect is between a continuous (hp) and a categorical (factor(cyl)) predictor. As such, each of the two estimated coefficents can be interpreted as the change in the slope of hp for each of the non-reference levels of factor(cyl).
coef(car.fit)[4] # When the car has 4 cylinders (reference level), the slope for hp is -0.0995 (to 4 decimal places)
coef(car.fit)[4] + coef(car.fit)[6] # When the car has 6 cylinders, the slope for hp is -0.0995 + 0.0781 = -0.0214 (to 4 decimal places)
coef(car.fit)[4] + coef(car.fit)[7] # When the car has 8 cylinders, the slope for hp is -0.0995 + 0.0860 = -0.0135 (to 4 decimal places)
# This model suggests that as hp increases, mean MPG decreases (for a fixed wt). However, in comparison to 4-cylinder cars, mean MPG is estimated to decrease at a slower rate with an increasing hp for both 6- and 8-cylinder cars (since the positive additive terms to the baseline slope of -0.0995 still provide negative slopes in hp, but ones that are closer to zero).
#(h)
##(i)
predict(car.fit,newdata=data.frame(wt=c(2.1,3.9,2.9),hp=c(100,210,200),cyl=c(4,8,6)),interval="confidence",level=0.95)
# The first car is the only car that has a point estimate of mean MPG that is higher than your mother's demand of 25, so this would be the initial choice.
##(ii)
# Although the point estimate for Car 3 is much less than 25, looking at the confidence intervals you can see that the interval for Car 3 includes 25. So, you could argue to your mother that you're 95% confident that the true mean MPG of a car like Car 3 lies somewhere in that interval; in particular, your model suggests no evidence against the hypothesis that the true mean MPG of such a car is equal to 25. Of course, the interval also includes possible true values that are far worse than 25... but you don't need to tell your mother that.

##########
## 22.1 ##
##########
#(a)
library("boot")
nuc.null <- lm(cost~1,data=nuclear)
nuc.step <- step(nuc.null,scope=.~.+date+t1+t2+cap+pr+ne+ct+bw+cum.n+pt,direction="both")
summary(nuc.step)
#(b) The model in (a) is again different to that selected by either forward or backward elimination. It's most similar to the model chosen via backward selection, differing only by an additional, non-significant, effect for 'ct'.
#(c)
gal <- data.frame(d=c(573,534,495,451,395,337,253),h=c(1,0.8,0.6,0.45,0.3,0.2,0.1))
gal.mod0 <- lm(d~1,data=gal)
gal.mod1 <- lm(d~h,data=gal)
gal.mod2 <- lm(d~h+I(h^2),data=gal)
gal.mod3 <- lm(d~h+I(h^2)+I(h^3),data=gal)
gal.mod4 <- lm(d~h+I(h^2)+I(h^3)+I(h^4),data=gal)
#(d)
anova(gal.mod0,gal.mod1,gal.mod2,gal.mod3,gal.mod4)
# Yes, the model chosen via nested partial F-tests matches the cubic model you would have likely identified in Exercise 21.2. The improvement to goodness-of-fit is minimal (and not statistically significant) if moving from the order 3 to the order 4 model.
#(e)
library("faraway")
diab <- na.omit(diabetes[,c("chol","age","gender","height","weight","frame","waist","hip","location")])
#(f)
dia.null <- lm(chol~1,data=diab)
dia.full <- lm(chol~age*gender*weight*frame+waist*height*hip+location,data=diab)
#(g)
dia.step <- step(dia.null,scope=.~.+age*gender*weight*frame+waist*height*hip+location)
summary(dia.step)
#(h)
add1(dia.null,scope=.~.+age*gender*weight*frame+waist*height*hip+location,test="F")
dia.1 <- update(dia.null,.~.+age)
add1(dia.1,scope=.~.+age*gender*weight*frame+waist*height*hip+location,test="F")
dia.2 <- update(dia.1,.~.+frame)
add1(dia.2,scope=.~.+age*gender*weight*frame+waist*height*hip+location,test="F")
summary(dia.2)
# The two models in (g) and (h) are quite different. Step-wise AIC has chosen a model that includes the two-way interaction between 'age' and 'frame', as well as a main effect for 'waist'. Performing forward selection halts after main effects for 'age' and 'frame' are added only, with no further additions offering a statistically significant improvement to goodness-of-fit.
#(i)
dia.step2 <- step(dia.full)
summary(dia.step2)
# The model chosen via step-wise AIC when starting from 'dia.full' is far, far more complex than the AIC-selected model in (g). The model contains a three-way interaction between 'age', 'gender' and 'weight', and all lower-order effects, as well as main effects and a two-way interaction between 'waist' and 'hip'. By starting from the very complex model, it seems we've accessed potential candidate models which have AIC values that are lower than those accessible by starting the selection algorithm at the null (intercept-only) model. For this reason, the final model reached in this way refuses to move to the simpler model selected in (g) since it's already found a model with a lower value for that criterion. Whether or not we'd be willing to settle on this rather complex model for predictive purposes is a more difficult question to answer...
#(j)
library("MASS")
car.null <- lm(I(1/mpg)~1,data=mtcars)
car.step <- step(car.null,scope=.~.+wt*hp*factor(cyl)*disp+am+factor(gear)+drat+vs+qsec+carb)
summary(car.step)
# Modeling the response as GPM instead of MPG, step-wise AIC selection offers up a far simpler model, with main effects for 'wt' and 'hp' only. It would seem that this transformation of the response greatly simplifies the relationships in the data.

##########
## 22.2 ##
##########
#(a)
library("boot")
nuc.fit <- lm(cost~date+cap+pt+ne,data=nuclear)
summary(nuc.fit)
#(b)
plot(nuc.fit,which=1)
plot(nuc.fit,which=2)
# The residuals vs. fitted plot doesn't suggest any major causes of concern: the points appear randomly scattered around zero with no indication of heteroscedasticity. However, there is one clear extreme observation, labeled 19, that should be kept in mind. The QQ plot also shows little cause for concern with fairly modest deviation of the residuals from the normal quantiles, though again, observation 19 departs considerably. All in all, it seems reasonable to assume satisfaction of the error assumptions for now (at the very least, there is no clearly obvious violation thereof).
#(c)
cutoff <- 4/nrow(nuclear)
cutoff
plot(nuc.fit,which=4)
abline(h=cutoff,lty=2)
# The Cook's distances demonstrate what you might already have expected: that observation 19 is highly influential. It breaches the rule-of-thumb cut-off of 0.125 by a considerable amount. Observation 10 is the next most influential, albeit one that registers a Cook's distance less than 0.125.
#(d)
plot(nuc.fit,which=5,cook.levels=c(cutoff,0.5,1))
# Observation 19 sits in a low-to-medium leverage position; it's the large residual for that point that pushes it past the rule-of-thumb contour. The same characterization can be applied to the other two labeled points that have more or less the same leverage (observations 10 and 12), though their smaller residuals mean a correspondingly lower influence.
#(e)
nuc.fit2 <- lm(cost~date+cap+pt+ne,data=nuclear[-19,])
summary(nuc.fit2)
plot(nuc.fit2,which=1)
plot(nuc.fit2,which=2)
# Assuming there's a good reason to remove observation 19, once done so, the plots of the residuals are improved in terms of satisfaction of the assumptions of the error component. Independence and homgeneity are illustrated with randomness around zero, and normality also appears reasonable. Points 10 and 12 are now the top two extreme points, joined by observation 7.
#(f)
library("faraway")
dia.fit <- lm(chol~age*frame+waist,data=diabetes)
summary(dia.fit)
#(g)
plot(dia.fit,which=1)
plot(dia.fit,which=2)
# Very similar to the interpretation of the plots in the earlier exercise (b), the error assumptions seem satisfied here, with randomness around zero and homoscedasticity seemingly valid. There are some extreme positive residuals, namely point 63. The QQ plot does appear to exhibit some departure from normality in the upper tail, with the three extreme points 63, 295 and 148 not helping in that respect. Estimation will remain valid if normality is questionable, though the estmates might not be 'optimal' in a theoretical sense.
#(h)
nrow(diabetes)-15
cutoff <- 4/(nrow(diabetes)-15)
cutoff
plot(dia.fit,which=5,cook.levels=c(1,3,5)*cutoff)
# The size of the data set results in a relatively small rule-of-thumb influence cut-off of around 0.01. A number of points breach this mark, though the overall pattern of the standardized residuals vs. leverage appears consistent with the contours, so we shouldn't be especially concerned here. There are three points that breach the '3 times cut-off' contour. Point 63 is in a very low leverage position, but its extreme residual leads it to be highly influential. Points 4 and 148 have much smaller residuals, but their position in much higher leverage areas when compared with 63 which means they too are considered highly influential here.
#(i)
dia.url <- "http://www.amstat.org/publications/jse/v9n2/4cdata.txt"
diamonds <- read.table(dia.url)
names(diamonds) <- c("Carat","Color","Clarity","Cert","Price")
plot(diamonds$Carat,diamonds$Price,pch=19,col=as.numeric(diamonds$Clarity))
legend("topleft",legend=levels(diamonds$Clarity),col=1:length(levels(diamonds$Clarity)),pch=19)
plot(diamonds$Carat,diamonds$Price,pch=19,col=as.numeric(diamonds$Color))
legend("topleft",legend=levels(diamonds$Color),col=1:length(levels(diamonds$Color)),pch=19)
plot(diamonds$Carat,diamonds$Price,pch=19,col=as.numeric(diamonds$Cert))
legend("topleft",legend=levels(diamonds$Cert),col=1:length(levels(diamonds$Cert)),pch=19)
#(j)
sparkly.fit <- lm(Price~Carat+Color+Clarity+Cert,data=diamonds)
summary(sparkly.fit)
plot(sparkly.fit,which=1)
plot(sparkly.fit,which=2)
plot(sparkly.fit,which=3)
# There seems to be clear, systematic non-linearity in the estimated residuals. This violates the linearity assumption of the trends in our data, and suggests our current model is inadequate in terms of representation of the data at hand. Based on the plots from (i), should we try modeling a log-transformation of the response?
# Normality is also affected---there is rather obvious deviation, from the distribution, particularly in the upper tail.
# Barring three extreme points, though, plots 1 and 3 don't show that the assumptions of homoscedasticity is violated---the variability of the residuals remains more or less constant.
#(k)
sparkly.fit2 <- lm(log(Price)~Carat+Color+Clarity+Cert,data=diamonds)
summary(sparkly.fit2)
plot(sparkly.fit2,which=1)
plot(sparkly.fit2,which=2)
plot(sparkly.fit2,which=3)
# Log-transformation of Price has done nothing to curb the clear, systematic, non-linear appearance of the residuals. A non-linear trend still appears very promenent.
# However, the log-transformation has reigned in the extreme points and the severity of the non-normality to a certain extent.
#(l)
sparkly.fit3 <- lm(log(Price)~Carat+I(Carat^2)+Color+Clarity+Cert,data=diamonds)
summary(sparkly.fit3)
plot(sparkly.fit3,which=1)
plot(sparkly.fit3,which=2)
plot(sparkly.fit3,which=3)
# Including an order 2 polynomial term in Carat has eliminated much of the concern of the non-linear curvature.
# Normality is more apparent now, and there's little if any concern of heteroscedasticity.
# This third model is by far the most appropriate of the three (with respect to the residual diagnostics) when it comes to modeling the cost of diamonds in light of the available data.

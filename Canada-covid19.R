rm(list=ls())
setwd("D:/R codes/CANADA_Covid_19")
library(tidyverse)
Covid19_data <- read.csv("Canada June_3_2020_covid19.csv")
str(Covid19_data)
#non-parametric smoother function, s()
#nb to estimate the theta parameter of the negative binomial (gam only).
#bs="cc" specifies a cyclic cubic regression splines 
#Smooth terms are specified in a gam formula using s, te, ti and t2 terms.
names(Covid19_data)
str(Covid19_data)
fix(Covid19_data)
view(Covid19_data)
plot(Covid19_data$numdeaths)
hist(Covid19_data$numdeaths)
ca <- Covid19_data %>% 
  filter(prname == "Canada" & prname != "NA")
ca$numdeaths
ca1 <- c(ca$numdeaths)
length(ca1)
can <-diff(ca1)
length(can)
can <- can[11:length(ca1)-1]
length(can)
mean(can); var(can); var(can)/mean(can)
QU <- Covid19_data %>% 
  filter(prname == "Quebec" & prname != "NA")
QU$numdeaths
qu1 <- c(QU$numdeaths)
que <- diff(qu1)
length(que); length(can)
mean(que); var(que) ;var(que)/mean(que)

decomp <- decompose(ca2)
par(mfrow=c(1,1))
plot(decomp, col='darkcyan')
plot(decomp$figure+25,
     type = 'b',
     xlab = 'Weeks',
     ylab = 'Seasonality Index',
     col  = 2)
mean(can); var(can);var(can)/mean(can) #Overdespersion 
ps <- FALSE
if (ps) postscript("covid-deaths.eps",width=11,height=4.5)
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(can,type="b",xlab="Days between Jan. 31 2020 and June 3rd. 2020",
     ylab="Reported Canada Deaths", col='blue')
plot(que,type="b",xlab="Days between Jan. 31 2020 and June 3rd. 2020",
     ylab="Reported Quebec Deaths", col='blue')
if (ps) dev.off()

#Not normal distribution
par(mfrow=c(1,2))
hist(can,
     xlab = 'Canada',
     main = 'Histogram of Deaths',
     freq = F, col='green')
lines(density(can), col='blue', lw=2)
hist(que,
     xlab = 'Quebec',
     main = 'Histogram of Deaths',
     freq = F, col='green')
lines(density(que), col='blue', lw=2)
#########################################
set.seed(7)
day <- 1:length(can)
dow <- rep(1:7,100)[1:length(can)]
CA <- data.frame(deaths = can, day=day,dow=dow) ## Canada data
QU <- data.frame(deaths = que, day=day,dow=dow) ## Quebec data
## Fit the basic death profile models...
library(car)
scatterplotMatrix(CA,pch=19,cex=.5,reg.line=F, lwd.smooth=1.25,
                  spread=F,ellipse=T, col=c('blue','#2957FF','#FF8000'),
                  col.axis='gray50')

ggplot(aes(x=day,y=deaths), data=CA) +
  geom_point(color='#FF8000',alpha=.75) +
  geom_smooth(se=F, method='gam', formula=y~s(x), color='#2957FF') 
ggplot(aes(x=dow,y=deaths), data=CA) +
  geom_point(color='#FF8000',alpha=.75) +
  geom_smooth(se=F, method='gam', formula=y~s(dow,bs="cc",k=7), color='#2957FF') 
# rnbinom(n, size, prob, mu)
#n	number of observations.
# size	= target for number of successful trials= tetha
#var = (mu +mu^2)/tetha
#dnbinom(x, size, prob, mu, log = FALSE) 
#size=target for number of successful trials
#prob=probability of success in each trial. 0 < prob <= 1.
#prob = size/(size+mu)
#variance = mu + mu^2/size 
library(mgcv)
#te=tensor product
#s(spline) smooth 
#bs =   based smooths
#"tp" for thin plate regression spline, "cr" for cubic regression spline)
#bs="ds"  Duchon.spline
# bs="cr". cubic regression spline 
# bs="cs" specifies a shrinkage version of "cr".
#bs="cc" specifies a cyclic cubic regression splines
#(see cyclic.cubic.spline). i.e. a penalized cubic regression 
#splines whose ends match, up to second derivative.
## cyclic spline example...
set.seed(7)
day <- 1:length(can)
dow <- rep(1:7,100)[1:length(can)]
CA <- data.frame(deaths = can, day=day,dow=dow) ## Canada data
QU <- data.frame(deaths = que, day=day,dow=dow) ## Quebec data
## Fit the basic death profile models...
bc <- gam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=nb(link="log"),data=CA,
          method="REML",knots=list(dow=c(0,150))) 
## much lower AIC than identity link
summary(bc)
par(mfrow=c(2,2))
gam.check(bc)
# 150   0.0819 . 
bq <- gam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=nb(link="log"),data=QU,
          method="REML",knots=list(dow=c(0,7)))
summary(bq)
par(mfrow=c(2,2))
gam.check(bq)
# shows the 19 is for day and 5 for dow
#1. Separate cubic polynomials are fit at each section,
#   and then joined at the knots to create a continuous curve
#2. Theorem = cubic polynomial are the optimized curves for each section (piece)
#3. REML for getting lambda
#4. the rank r of the covariance matrix
#   for the coefficients for a particular smooth, so here conceptually it
#   is the p-value associated with the F(r,n- ed f )
#5. vcov(bc)= covariance
#6. The GCV, or generalized cross validation score can be taken as an
#   estimate of the mean square prediction error based on a leave-one-out
#   cross validation estimation process
#7.The intervals are Bayesian credible
#8. identifiability constraints, the smooths must sum
#   to zero, and thus are presented in a mean-centered fashion.
par(mfrow=c(1,1),mar=c(4,4,1,1))
vis.gam(bc, type='response', plot.type='persp',
        phi=25, theta=30,n.grid=500, border=NA)
vcov(bc)
coef(bc)
model.matrix(bc) # 91 data points and  25 vaiables
dim(model.matrix(bc))
vcov.gam(bc) #shows 19 and 5 variables
plot(bc, pages=1, residuals=T, pch=19, cex=0.25,
     scheme=1, col='blue', shade=T,shade.col='gray90')
vis.gam(bc, type='response', plot.type='contour')
#1. with day on the x axis, dow on the
#   y, with values on the response= deaths
#   scale given by the contours, with lighter
#   color indicating higher values.
#2.Technically we could specify the smoothing parameters
#  (Wood,2006, p. 128)
#   explicitly, and the appendix has some 'by-hand' code taken directly
#   from Wood (2006) with only slight modifications
#   Smoothing parameters are selected which minimize the GCV score
model_matrix <- predict(bc, type = "lpmatrix")
dim(model_matrix) 
model_matrix[1:4,2:8]
model_matrix[1,2:20]
model_matrix[1,21:25]
plot(model_matrix[,10],model_matrix[,20])
plot(model_matrix[,10],model_matrix[,25])
plot(model_matrix[,21],model_matrix[,25])
plot(model_matrix[,18],model_matrix[,19])
# The scale estimate is the scaled deviance, which
# here is equivalent to the residual sums
# Dev.exp= R.sq in this case
# bij=model matrix
## Now plot the estimated effects and basic checking plots...
X <- model.matrix(bc)
dim(X)
X[,21:25] <- 0            ## Canada model matrix, weekly removed
Xq <- model.matrix(bq)
Xq[,21:25] <- 0           ## Quebec model matrix, weekly removed

ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
c1 <- 1; c2=1.3
plot(bc,scale=0,select=1,xlab="day, t",ylab="f(t)", main='Smooth functions',cex.lab=c2, col='red');text(45,-3.07,"a-Canada",cex=c1, col='blue')
plot(bc,scale=0,select=2,xlab="day of week, d",ylab=expression(f[h](d)), main='Smooth functions',cex.lab=c2, col='red');text(4,0.05,"b-Canada",cex=c1, col='blue')
acf(residuals(bc),cex.lab=c2);text(10,.8,"Canada",cex=c1, col='blue')
plot(day[1:length(residuals(bc))],residuals(bc),xlab="day",ylab="residuals",cex.lab=c2, col='blue');text(46, -2.8,"Canada",cex=c1, col='blue')
if (ps) dev.off()

ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(bq,scale=0,select=1,xlab="day, t",ylab="f(t)",cex.lab=c2, col='red');text(45,-3.7,"a-QC",cex=c1, col='blue')
plot(bq,scale=0,select=2,xlab="day of week, d",ylab=expression(f[h](d)),cex.lab=c2, col='red');text(6,.15,"a-QC",cex=c1, col='blue')
acf(residuals(bq),cex.lab=c2);text(10,.8,"Quebec",cex=c1)
plot(day[1:length(residuals(bq))],residuals(bq),xlab="day",ylab="residuals",cex.lab=c2, col='blue');text(49, -4,"Quebec",cex=c1)
if (ps) dev.off()
## Underlying death rate plots on response scale...
##  vcov(b)
if (ps) postscript("death-rates.eps",width=11,height=4.5)
par(mfrow=c(2, 1),mar=c(5,5,1,1))
yla <- c("Canada death rate","Quebec death rate")
for (j in 1:2) { ## 1, Canada; 2, Quebec}
  if (j==1) {
    beta <- coef(bc); Vbc <- vcov(bc); Xj <- X 
  } else {
    beta <- coef(bq); Vbc <- vcov(bq); Xj <- Xq
  }
  fv <- Xj%*%beta            ## death rates - link scale
  se <- rowSums(Xj*(Xj%*%Vbc))**(.5) ## corresponding s.e.
  ilink <- bc$family$linkinv ## inverse link function
  ## get and plot death rate profile and CI...
  mu <- ilink(fv);ll <- ilink(fv-2*se);ul <- ilink(fv+2*se)
  ylim <- if (j==1) range(can) else range(que)
  plot(mu,type="l",xlab="day",ylab='yla[j]',cex.lab=1.1,ylim=ylim, col='blue')
  lines(ll,lty=2, col='blue');lines(ul,lty=2, col='blue')
  ## simulate for distn. of peak location.
  n.rep <- 1000
  bcp <- rmvn(n.rep,beta,Vbc) ## simulate 1000 parameter vectors from posterior
  ## for multidimensional normal deviance
  peak <- apply(Xj%*%t(bcp),2,function(x) which(x==max(x))) ## find peak location for each
  pt <- tabulate(peak) / if (j==1) 1 else 7 ## tabulate peak locations and scale for plotting
  for (i in 1:length(pt)) if (pt[i]>0) lines(c(i,i),c(0,pt[i]),lwd=3,col="red") ## plot
}
if (ps) dev.off()
#########
## Now the infection profile modelling...
library(rjags)
load.module("glm") 
## improved samplers for GLMs often worth loading
#can1 <- c(rep(0,20),can)    # 20 days without deaths in paper
## append run in time before anyone dies
cand1 <- c(rep(0,30),can) ## append run in time before anyone dies
nc <- length(cand1)
day <- 1:nc
dow <- rep(1:7,100)[1:nc]
jdat <- data.frame(deaths = cand1, day=day,dow=dow)
d <- dgamma(1:nc,shape=17.8/4,scale=4) ##symptom onset to death
#dist from Verity et al. with mean 17.8 and variance 71.2 (s.d. 8.44).
#With a shape parameter k and a scale parameter ?? and mean parameter ?? = k??
d[1:9]
par(mfrow=c(1,1))
plot(d, ylab='Canada deaths', main= 'Gamma distribution', col='red')
if (ps) dev.off()
B <- matrix(0,nc,nc)
for (i in 1:nc) { ## map case rate day i-1 to death rate day i ...
  B[,i] <- c(rep(0,i-1),d[1:(nc-i+1)])
}
B[,1]
nc
det(B)
dim(B)
diag(B)
B[1:5, 1:5]
det(B[1:84, 1:84])
det(B[1:85, 1:85])
##########################################################################
#que1 <- c(rep(0,20),que)    # 20 days without deaths in paper
## append run in time before anyone dies
nc <- length(can)
day <- 1:nc
dow <- rep(1:7,100)[1:nc]
jdat <- data.frame(deaths = can, day=day,dow=dow)
d <- dgamma(1:nc,shape=17.8/4,scale=4) ## s-to-d dist from Verity et al.
d[1:9]
par(mfrow=c(1,1))
plot(d, ylab="Quebec deaths", main= 'Gamma distribution', col='red')
if (ps) dev.off()
B <- matrix(0,nc,nc)
for (i in 1:nc) { ## map case rate day i-1 to death rate day i ...
  B[,i] <- c(rep(0,i-1),d[1:(nc-i+1)])
}
det(B)
##############################################################
## following sets up a model template and associated data, suitable for editing
## into desired model...
jad <- jagam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=poisson,data=jdat,
             knots=list(dow=c(0,7)),file="deaths.jags",diagonalize=TRUE)
names(jad)
jad$pregam
names(jad$pregam)
jad$pregam$intercept
jad$pregam$rank
jad$pregam$var.summary
jad$pregam$n
jad$jags.data
jad$jags.ini  # $b (Intercept) = 25 coefficients
names(jad$jags.ini)
jad$jags.ini$b
jad$jags.ini$lambda
names(jad$jags.data)
jad$jags.data$y
jad$jags.data$n
jad$jags.data$X
dim(jad$jags.data$X)
dim(X)
#the diagonal elements of F are where the effective degrees of freedom
#for each covariate come from.
######################################################################
names(jad$jags.data)
jad$jags.data$y
jad$jags.data$n
jad$jags.data$X
jad$jags.data$B <- B ## adding in the forward mapping matrix
jad$jags.data$Z <- -6 ## a bit of extra regularization
names(jad$jags.data)
jad$jags.data$B
jad$jags.data$Z
######################$$$$$$$$$$$$$$
model="model{
  theta~dunif(1,1000)
  eta <- X[,1:20] %*% b[1:20] ## linear predictor for case
  fw <- X[,21:25] %*% b[21:25] ## lp for week cycle
  for (i in 1:n) { muc[i] <- exp(eta[i]) } ## expected cases
  Z~dnorm(eta[1],.1) ## allow slight regularization at start
  mud <- B %*% muc ## expected deaths
  for (i in 1:n) { mu[i] <- exp(log(mud[i])+fw[i]) } ## model expected deaths
  for (i in 1:n) { y[i]~dnegbin(theta/(theta+mu[i]),theta) } ## response
  ## Parametric effect priors tau=1/38^2
  for (i in 1:1) { b[i]~dnorm(0,0.00068) }
  ## prior for s(day)...
  for (i in c(2:19)) { b[i]~dnorm(0, lambda[1]) }
  for (i in c(20)) { b[i]~dnorm(0, lambda[2]) }
  ## prior for s(dow)...
  for (i in c(21:25)) { b[i]~dnorm(0, lambda[3]) }
  ## smoothing parameter priors...
  for (i in 1:3) {
    lambda[i]~dgamma(.05,.005)
    rho[i] <- log(lambda[i]) ## for monitoring
  }
}"
model1.spec<-textConnection(model)
jm <-jags.model(model1.spec,data=jad$jags.data,inits=jad$jags.ini,n.chains=1)
list.samplers(jm)
names(list.samplers(jm))
sam <- jags.samples(jm,c("b","rho","theta"),n.iter=3000000,thin=300)
###################################################################
names(sam)
sam$theta
sam$rho
sam$b; length(sam$b); dim(sam$b); sam$b[,,1]; 
jam <- sim2jam(sam,jad$pregam)
plot(jam,pages=1)
summary(jam)
names(jam)
jam$coefficients 
jam$family
jam$formula
jam$smooth
jam$var.summary
save(sam,file="sam2.rda")
load("sam2.rda") # loads in Global Environment
## check chains...
effectiveSize(as.mcmc.list(sam$theta))  #mcmc= Markof chain Monte carlo
effectiveSize(as.mcmc.list(sam$rho))
effectiveSize(as.mcmc.list(sam$b))
sam$b[1,2,]
sam$b[1,,]
typeof(sam$b[1,,])
class(sam$b[1,,])
sam$b[1,,][1]
typeof(sam$b[1,,][1])
class(sam$b[1,,][1])
dim(sam$b)
##############################################
## check the slow mixers...
par(mfrow=c(5,5),mar=c(4,4,1,1))
for (i in 1:25) plot(sam$b[i,,],type="l")
if (ps) dev.off()
par(mfrow=c(1,1),mar=c(4,4,1,1))
plot(sam$b[1,,])
if (ps) postscript("infections.eps",width=11,height=4.5)

dim(jad$jags.data$X)
typeof(jad$jags.data$X)
class(jad$jags.data$X)
typeof(sam$b[1,,])
class(sam$b[1,,])
dim(sam$b[1,,])
sam$b[,,1]
typeof(sam$b[,,1])
dim(sam$b[,,1])
class(sam$b[,,1])
##############################################################
X <- jad$jags.data$X
X[,21:25] <- 0
bcs <- sam$b[,,1]
bcs[1:5,1:10]
dim(bcs)
dim(X)
dim(X %*% bcs)
f <- exp(X %*% bcs) ## the simulated fatal infection profiles
dim(f)
f[1:5, 1:5]

A <- matrix(c(3,5,-6,7,41,8,5,8,-5,7,3,6), 3)
pk <- apply(A,2,function(x) which(x[1:3]==max(x[1:3]))) 
pt1 <- tabulate(pk)
A1 <- apply(A,1,median)
fm <- apply(f,1,median) ## median profile
length(fm)
fq <- apply(f,1,quantile,probs=c(.025,.1,.9,.975)) ## get profile CIs
day <- 1:length(cand1)-36 ## account for run in + 1 day shift deaths to onsets + 5 day incubation
c1 <- 1.0
par(mfrow=c(1,1),mar=c(2,2,1,1))
plot(day,fm,type="l",ylim=c(0,max(fq)),ylab="Canada fatal infections",xlab="day",cex.lab=c1,lwd=2)
abline(v=30,col="grey",lwd=3) ## lock down day
lines(day,fm, col='red')
lines(day,fq[1,],lty=3,lwd=2);lines(day,fq[2,],lty=2,lwd=2);
lines(day,fq[3,],lty=2,lwd=2);lines(day,fq[4,],lty=3,lwd=2);
## Some smooth diagnostics - looking for high grad or curvature around lock down that could
## indicate model trying to match abrupt change...
df <- abs(diff(fm))
lines(1:length(df)+.5-36,df,col="grey",lty=1)
df2 <- abs(diff(log(fm),difference=2))^2*100000
lines(1:length(df2)+1-36,df2,col="green",lty=2)
peak <- apply(f,2,function(x) which(x[1:91]==max(x[1:91]))) ## locate the peak (ignore any wild boundary highs)
pt <- tabulate(peak) ## tabulate and plot as bar chart...
for (i in 1:length(pt)) if (pt[i]>0) lines(c(i-6,i-6),c(0,pt[i]/4),lwd=3)
if (ps) dev.off()
length(cand1)
## Sanity check the median infection profile.
fi <- round(fm) ## median profile
n <- length(fi)
n.rep <- 100
death <- matrix(0,n.rep,n)
for (j in 1:n.rep) for (i in 1:n){
  t <- round(rgamma(fi[i],shape=17.8/4,scale=4)) ## simulate death times for fi[i] fatal cases
  t <- t[i+t-1<=n] ## discard deaths beyond end of data
  dd <- tabulate(t)
  ii <- 1:length(dd)+i-1
  death[j,ii] <- death[j,ii] + dd
}
par(mfrow=c(1,1),mar=c(5,5,4,4))
plot(day+6,death[1,],type="l",ylim=range(can),col="gray",xlab="day",ylab="Canada deaths",cex.lab=c1)
for (j in 2:n.rep) lines(day+6,death[j,],col="gray")
X <- model.matrix(bc);X[,21:25] <- 0 ## Canada model matrix, weekly removed
beta <- coef(bc); Vb <- vcov(bc); Xj <- X
fv <- Xj%*%beta ## death rates - link scale
se <- rowSums(Xj*(Xj%*%Vb))^.5 ## corresponding s.e.
lines(1:length(fv),exp(fv))
lines(1:length(fv),exp(fv+2*se),lty=2)
lines(1:length(fv),exp(fv-2*se),lty=2)
if (ps) dev.off()
#############################################################
on <- Covid19_data %>% 
  filter(prname == "Ontario" & prname != "NA")
on$numdeaths
plot(on$numdeaths)
on1 <- c(on$numdeaths)
length(on1)
ont <-diff(on1)
length(ont)
ont <- ont[16:length(on1)-1]
mean(ont); var(ont);var(ont)/mean(ont)
ps <- FALSE
if (ps) postscript("covid-deaths.eps",width=11,height=4.5)
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(ont,type="b",xlab="Days between Jan. 31 2020 and June 3rd. 2020",
     ylab="Reported Ontario Deaths", col='blue')
if (ps) dev.off()
#Not normal distribution
hist(ont,
     xlab = 'Ontario',
     main = 'Histogram of Deaths',
     freq = F, col='green')
lines(density(ont),col='blue', lw=2)
if (ps) dev.off()
AL <- Covid19_data %>% 
  filter(prname == "Alberta" & prname != "NA")
AL$numdeaths
al1 <- c(AL$numdeaths)
alb <- diff(al1)
##############################################
par(mfrow=c(1,1))
boxplot(can, que, ont, alb,
        main = "All four cases",
        at = c(1,2,3,4),
        names = c("Canada", "Quebec", "Ontario", "Alberta"),
        las = 2,
        col = c("red","orange","blue","green"),
        border = "brown",
        horizontal = F,
        notch = F)
boxplot(alb, col='green',main = "Alberta")
length(ont);length(alb)
sum(alb)
mean(alb); var(alb) ;var(alb)/mean(alb)
ps <- FALSE
if (ps) postscript("covid-deaths.eps",width=11,height=4.5)
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(alb,type="b",xlab="Days between Jan. 31 2020 and June 3rd. 2020",
     ylab="Reported Alberta Deaths", col='blue')
if (ps) dev.off()
ps <- FALSE
if (ps) postscript("covid-deaths.eps",width=11,height=4.5)
par(mfrow=c(1,1),mar=c(5,5,1,1))
hist(alb,
     xlab = 'Alberta',
     main = 'Histogram of Deaths',
     freq = FALSE, col='green')
lines(density(alb), col='blue', lw=2)
if (ps) dev.off()
length(ont);length(alb)
day <- 1:length(ont)
dow <- rep(1:7,100)[1:length(ont)]
ON <- data.frame(deaths = ont, day=day,dow=dow) ## Ontario data
AL <- data.frame(deaths = alb, day=day,dow=dow) ## Alberta data
## Fit the basic death profile models...
model="model{
  theta~dunif(1,1000)
  eta <- X[,1:20] %*% b[1:20] ## linear predictor for case
  fw <- X[,21:25] %*% b[21:25] ## lp for week cycle
  for (i in 1:n) { muc[i] <- exp(eta[i]) } ## expected cases
  Z~dnorm(eta[1],.1) ## allow slight regularization at start
  mud <- B %*% muc ## expected deaths
  for (i in 1:n) { mu[i] <- exp(log(mud[i])+fw[i]) } ## model expected deaths
  for (i in 1:n) { y[i]~dnegbin(theta/(theta+mu[i]),theta) } ## response
  ## Parametric effect priors tau=1/56^2
  for (i in 1:1) { b[i]~dnorm(0,0.00032) }
  ## prior for s(day)...
  for (i in c(2:19)) { b[i]~dnorm(0, lambda[1]) }
  for (i in c(20)) { b[i]~dnorm(0, lambda[2]) }
  ## prior for s(dow)...
  for (i in c(21:25)) { b[i]~dnorm(0, lambda[3]) }
  ## smoothing parameter priors...
  for (i in 1:3) {
    lambda[i]~dgamma(.05,.005)
    rho[i] <- log(lambda[i]) ## for monitoring
  }
}"
library(mgcv)
b <- gam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=nb(link="log"),data=ON,
         method="REML",knots=list(dow=c(0,7))) ## much lower AIC than identity link
bs <- gam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=nb(link="log"),data=AL,
          method="REML",knots=list(dow=c(0,7)))

summary(b)
summary(bs)
ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
gam.check(b)
if (ps) dev.off()
ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
gam.check(bs)

ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
c1 <- 1.7; c2=1.3
plot(b,scale=0,select=1,xlab="day, t",ylab="f(t)", main='Smoothe functions',cex.lab=c2, col='red')
text(53,-3.07,"f-Ontario",cex=c1, col='blue')
plot(b,scale=0,select=2,xlab="day of week, d",ylab=expression(f[h](d)), main='Smoother functions',cex.lab=c2, col='red')
text(4,.18,"f(w)-Ontario",cex=c1, col='blue')
acf(residuals(b),cex.lab=c2)
plot(day[1:length(residuals(b))],residuals(b),xlab="day",ylab="residuals",cex.lab=c2, col='blue')
if (ps) dev.off()

ps <- FALSE
if (ps) postscript("death-fits.eps",height=5)
par(mfrow=c(2,2),mar=c(4,4,1,1))
c1 <- 1.7; c2=1.3
plot(bs,scale=0,select=1,xlab="day, t",ylab="f(t)",cex.lab=c2, col='red')
text(45,-3.7,"f-Alberta",cex=c1, col='blue')
plot(bs,scale=0,select=2,xlab="day of week, d",ylab=expression(f[h](d)),cex.lab=c2, col='red')
text(4,.5,"f(w)-Alberta",cex=c1, col='blue')
acf(residuals(bs),cex.lab=c2)
plot(day[1:length(residuals(bs))],residuals(bs),xlab="day",ylab="residuals",cex.lab=c2, col='blue')
text(40, 2,"Alberta",cex=c1)
if (ps) dev.off()

if (ps) postscript("death-rates.eps",width=11,height=4.5)
par(mfrow=c(2, 1),mar=c(5,5,1,1))

X <- model.matrix(b)
X[,21:25] <- 0 ## Ontario model matrix, weekly removed
Xs <- model.matrix(bs)
Xs[,21:25] <- 0 ## Alberta model matrix, weekly removed

yla <- c("Ontario death rate","Alberta death rate")
for (j in 1:2) { ## 1, Ontario; 2, Alberta
  if (j==1) {
    beta <- coef(b); Vb <- vcov(b); Xj <- X 
  } else {
    beta <- coef(bs); Vb <- vcov(bs); Xj <- Xs
  }
  fv <- Xj%*%beta ## death rates - link scale
  se <- rowSums(Xj*(Xj%*%Vb))**(.5) ## corresponding s.e.
  ilink <- b$family$linkinv ## inverse link function
  ## get and plot death rate profile and CI...
  mu <- ilink(fv);ll <- ilink(fv-2*se);ul <- ilink(fv+2*se)
  ylim <- if (j==1) range(ont) else range(alb)
  plot(mu,type="l",xlab="day",ylab='yla[j]',cex.lab=1.1,ylim=ylim, col='blue')
  lines(ll,lty=2, col='blue');lines(ul,lty=2, col='blue')
  ## simulate for distn. of peak location.
  n.rep <- 1000
  bp <- rmvn(n.rep,beta,Vb) ## simulate 1000 parameter vectors from posterior
  peak <- apply(Xj%*%t(bp),2,function(x) which(x==max(x))) ## find peak location for each
  pt <- tabulate(peak) / if (j==1) 1 else 7 ## tabulate peak locations and scale for plotting
  for (i in 1:length(pt)) if (pt[i]>0) lines(c(i,i),c(0,pt[i]),lwd=3,col="red") ## plot
}   # bar charts
if (ps) dev.off()
## Now the infection profile modelling...
##########################################
library(rjags)
library(mgcv)
load.module("glm") 
## improved samplers for GLMs often worth loading
#can1 <- c(rep(0,30),can)    # 30 days without deaths in paper
## append run in time before anyone dies
nc <- length(ont)
day <- 1:nc
dow <- rep(1:7,100)[1:nc]
jdat <- data.frame(deaths = ont, day=day,dow=dow)
d <- dgamma(1:nc,shape=17.8/4,scale=4) ##symptom onset to death
#dist from Verity et al. with mean 17.8 and variance 71.2 (s.d. 8.44).
#With a shape parameter k and a scale parameter ?? and mean parameter ?? = k??
d[1:9]
par(mfrow=c(1,1))
plot(d, ylab='Ontario deaths', main= 'Gamma distribution', col='red')
if (ps) dev.off()

B <- matrix(0,nc,nc)
for (i in 1:nc) { ## map case rate day i-1 to death rate day i ...
  B[,i] <- c(rep(0,i-1),d[1:(nc-i+1)])
}
det(B)
## following sets up a model template and associated data, suitable for editing
## into desired model...
jad <- jagam(deaths~s(day,k=20)+s(dow,bs="cc",k=7),family=poisson,data=jdat,
             knots=list(dow=c(0,7)),file="death-Canada.jags",diagonalize=TRUE)
jad$jags.data$B <- B ## adding in the forward mapping matrix
jad$jags.data$Z <- -6 ## a bit of extra regularization
model1.spec<-textConnection(model)
jm <-jags.model(model1.spec,data=jad$jags.data,inits=jad$jags.ini,n.chains=1)
#jm <-jags.model(death-Canada.jags,data=jad$jags.data,inits=jad$jags.ini,n.chains=1)
sam <- jags.samples(jm,c("b","rho","theta"),n.iter=3000000,thin=300)
######################
jam <- sim2jam(sam,jd$pregam)
plot(jam,pages=1)
summary(jam)
jam$coefficients 
jam$family
jam$formula
######################$$$$$$$$$$
save(sam,file="sam2.rda")
load("sam2.rda") # loads in Global Environment
## check chains...
effectiveSize(as.mcmc.list(sam$theta))  #mcmc= Markof chain Monte carlo
effectiveSize(as.mcmc.list(sam$rho))
effectiveSize(as.mcmc.list(sam$b))
## check the slow mixers...
par(mfrow=c(5,5),mar=c(4,4,1,1))
for (i in 1:25) plot(sam$b[i,,],type="l")
if (ps) dev.off()
par(mfrow=c(1,1),mar=c(4,4,1,1))
if (ps) postscript("infections.eps",width=11,height=4.5)
c1 <- 1.0
X <- jad$jags.data$X
X[,21:25] <- 0
bcs <- sam$b[,,1]
f <- exp(X %*% bcs) ## the simulated fatal infection profiles
fm <- apply(f,1,median) ## median profile
fq <- apply(f,1,quantile,probs=c(.025,.1,.9,.975)) ## get profile CIs
day <- 1:length(ont)-16 ## account for run in + 1 day shift deaths to onsets + 5 day incubation
par(mfrow=c(1,1),mar=c(2,2,1,1))
plot(day,fm,type="l",ylim=c(0,max(fq)),ylab="Ontario fatal infections",xlab="day",cex.lab=c1,lwd=2)
abline(v=22,col="grey",lwd=3) ## lock down day
lines(day,fm, col='red')
lines(day,fq[1,],lty=3,lwd=2);lines(day,fq[2,],lty=2,lwd=2, col='blue');
lines(day,fq[3,],lty=2,lwd=2);lines(day,fq[4,],lty=3,lwd=2, col='blue');
## Some smooth diagnostics - looking for high grad or curvature around lock down that could
## indicate model trying to match abrupt change...
df <- abs(diff(fm))
lines(1:length(df)+.5-20,df,col="grey",lty=1)
df2 <- abs(diff(log(fm),difference=2))^2*100000
lines(1:length(df2)+1-20,df2,col="green",lty=2)
require(rjags)
peak <- apply(f,2,function(x) which(x[1:86]==max(x[1:86]))) ## locate the peak (ignore any wild boundary highs)
pt <- tabulate(peak) ## tabulate and plot as bar chart...
for (i in 1:length(pt)) if (pt[i]>0) lines(c(i,i),c(0,pt[i]/4),lwd=3)
if (ps) dev.off()
length(ont);length(alb)
## Sanity check the median infection profile.
fi <- round(fm) ## median profile
n <- length(fi)
n.rep <- 100
death <- matrix(0,n.rep,n)
for (j in 1:n.rep) for (i in 1:n){
  t <- round(rgamma(fi[i],shape=17.8/4,scale=4)) ## simulate death times for fi[i] fatal cases
  t <- t[i+t-1<=n] ## discard deaths beyond end of data
  dd <- tabulate(t)
  ii <- 1:length(dd)+i-1
  death[j,ii] <- death[j,ii] + dd
}
par(mfrow=c(1,1),mar=c(5,5,4,4))
plot(day+15,death[1,],type="l",ylim=range(can),col="gray",xlab="day",ylab="Ontario deaths",cex.lab=c1)
for (j in 2:n.rep) lines(day+15,death[j,],col="gray")
X <- model.matrix(bc);X[,21:25] <- 0 ## Ontario model matrix, weekly removed
beta <- coef(bc); Vb <- vcov(bc); Xj <- X
fv <- Xj%*%beta ## death rates - link scale
se <- rowSums(Xj*(Xj%*%Vb))^.5 ## corresponding s.e.
lines(1:length(fv),exp(fv))
lines(1:length(fv),exp(fv+2*se),lty=2)
lines(1:length(fv),exp(fv-2*se),lty=2)
if (ps) dev.off()
###########&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


## ----eval=FALSE---------------------------------------------------------------
#  fr=function(x){
#    z=rnorm(100)
#    d=sort(abs(z-x))
#    y=numeric(100)
#    for(k in 1:100){
#      y[k]=k/2/100/d[k]
#    }
#    dis=abs(y-dnorm(x))
#    re=which.min(dis)
#    return(re)
#  }

## ----eval=FALSE---------------------------------------------------------------
#  library(Rcpp) # Attach R package "Rcpp"
#  # Define function "fc"
#  cppFunction('int fc(int x, int y, int z) {
#  int sum = x + y + z;
#  return sum;
#  }')

## -----------------------------------------------------------------------------
a=rnorm(100,mean=0,sd=1)
summary(a)

## -----------------------------------------------------------------------------
hist(a, freq = F, breaks = 30)
lines(density(a, bw=.5), col="red", lwd=2)

## -----------------------------------------------------------------------------
b=rnorm(1000,mean=0,sd=1)
summary(b)

## -----------------------------------------------------------------------------
hist(b, freq = F, breaks = 30)
lines(density(b, bw=.5), col="red", lwd=2)

## -----------------------------------------------------------------------------
x1=rnorm(1000,1,1)
x2=rnorm(1000,3,1)
x3=x1+x2
mean(x3)
hist(x3, freq = F, breaks = 30)
lines(density(x3, bw=.5), col="red", lwd=2)
x4=cbind(x1,x2,x3)
library(knitr)
kable(head(x4),format="markdown")

## -----------------------------------------------------------------------------
library(knitr)
y1=runif(100)
y2=rnorm(100)
y3=rexp(100)
y4=cbind(y1,y2,y3,y1^2,y2^2,y3^2)
kable(head(y4),format="markdown")
par(mar=c(1,1,1,1))
plot(y1,type="l",col="pink")
plot(y2,type="l",col="green")
plot(y3,type="l",col="yellow")

## -----------------------------------------------------------------------------
n=10000                        #size of samples
f_sigma=function(sigma){            
  u=runif(n)
  x=sqrt(-2*sigma^2*log(1-u))
  x
}

## -----------------------------------------------------------------------------
par=(mar=c(1,1,1,1,1,1,1,1,1))
sigma=c(1,2,3,4,0.5,0.7,10.3,6.1,20)
for(i in 1:9){
x1=f_sigma(sigma[i])
hist(x1,prob=TRUE,main="Rayleigh density")
y <- seq(0, 100, .001)
lines(y, y/sigma[i]^2*exp(-y^2/(2*sigma[i]^2)),col="red")
}

## -----------------------------------------------------------------------------
n=1000
mixture_p1=function(p1){
r=sample(c(0,1),size=n,replace=TRUE,prob=c(1-p1,p1))
#p1 for N(0,1),p2=1-p1 for N(3,1)
x1=rnorm(n,0,1)
x2=rnorm(n,3,1)
z=r*x1+(1-r)*x2
hist(z,,breaks=50,prob=TRUE)
lines(density(z),col="red")
}

## -----------------------------------------------------------------------------
mixture_p1(0.75)

## -----------------------------------------------------------------------------
mixture_p1(0.05)
mixture_p1(0.15)
mixture_p1(0.25)
mixture_p1(0.35)
mixture_p1(0.45)
mixture_p1(0.55)
mixture_p1(0.65)
mixture_p1(0.85)
mixture_p1(0.95)

## -----------------------------------------------------------------------------
mixture_p1(0.37)
mixture_p1(0.4)
mixture_p1(0.43)
mixture_p1(0.47)
mixture_p1(0.5)
mixture_p1(0.53)
mixture_p1(0.57)
mixture_p1(0.6)
mixture_p1(0.63)

## -----------------------------------------------------------------------------
wishart=function(sigma,n,d){
  if(n<d+2||d<0) print("please input integers n>d+1 and d>=0")
  else{
    A=matrix(numeric(d^2),d,d)
    for(i in 1:d){
      A[i,i]=sqrt(rchisq(1,n+1-i))
      for(j in 1:d){
        if(j<i) A[i,j]=rnorm(1)
        else if(j>i) A[i,j]=0
        else A[i,j]=A[i,i]
      }
    }#definate matrix A
    L=t(chol(sigma))  #Compute the Choleski factorization of sigma.
    X=L%*%A%*%t(A)%*%t(L)   #based on Bartlett's decompsition
  }
  X
}

## -----------------------------------------------------------------------------
sigma=matrix(c(5,1,1,3),2,2)
d=2
n=4
wishart(sigma,n,d)

## -----------------------------------------------------------------------------
set.seed(1)
m=1e4
x=runif(m,min=0,max=pi/3)
y=mean(pi/3*sin(x))
print(c(y,-cos(pi/3)+cos(0)))

## -----------------------------------------------------------------------------
set.seed(12)
m=1e4
x=runif(m)
x1=x[1:m/2]
u1=exp(-x1)/(1+x1^2)
u2=exp(-(1-x1))/(1+(1-x1)^2)
y1=(u1+u2)/2
theta.hat=mean(y1)#with variance reduction
y2=exp(-x)/(1+x^2)
theta.hat2=mean(y2) #without variance reduction
c(theta.hat,theta.hat2,var(y1),var(y2),1-var(y1)/var(y2))

## -----------------------------------------------------------------------------
sub=function(a,b){#write a function for every subinterval to obtain importance sampling
m <- 10000
theta.hat <- se <- numeric(5)
g <- function(x) {
5*exp(-x - log(1+x^2))
}
u <- runif(m/5) #f3, inverse transform method
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))*(x>a)*(x<b)
theta.hat <- mean(fg)
return(theta.hat)
}

## -----------------------------------------------------------------------------
set.seed(54321)
alpha=c(0,0,0,0,0,1)#(0,1)分成五个子区间，第i个子区间的端点是(j-1)/5分位数和j/5分位数
for(j in 1:5){
  f=function(x)(1-exp(-x))/(1-exp(-1))-(j-1)/5
alpha[j]=uniroot(f,c(0,1))$root
}#求分位数
est=c(rep(0,50))
T=numeric(0)
for(i in 1:50){
for(j in 1:5)
  T[j]=sub(alpha[j],alpha[j+1])
  est[i]=mean(T)
}
c(mean(est),sd(est))

## -----------------------------------------------------------------------------
set.seed(11)
alpha=0.05
n=20
x=rchisq(n,2)
s=sd(x)
l1=mean(x)-s/sqrt(n)*qt(1-alpha/2,n-1)
l2=mean(x)+s/sqrt(n)*qt(1-alpha/2,n-1)
c(l1,l2)

## -----------------------------------------------------------------------------
set.seed(12)
m=1000
l1=l2=k=numeric(m)
for(i in 1:m){
x=rchisq(n,2)
s=sd(x)
l1[i]=mean(x)-s/sqrt(n)*qt(1-alpha/2,n-1)#卡方分布期望等于自由度
l2[i]=mean(x)+s/sqrt(n)*qt(1-alpha/2,n-1)
k[i]=(l1[i]<=2)*(2<=l2[i])
}
mean(k)

## -----------------------------------------------------------------------------
g=function(x){
  a=mean((x-mean(x))^3)
  b=(mean((x-mean(x))^2))^(3/2)
  return(a/b)
}

## -----------------------------------------------------------------------------
set.seed(3)
m=10000
n=200
b=numeric(m)
for(i in 1:m){
  x=rnorm(n)
  b[i]=g(x)
}

## -----------------------------------------------------------------------------
qb=numeric(4)
qb=quantile(b,probs=c(0.025,0.05,0.95,0.975))
qb

## -----------------------------------------------------------------------------
varb=function(q){
  xq=qnorm(q,0,sqrt(6/n))
  return(q*(1-q)/n/dnorm(xq,0,sqrt(6*(n-2)/(n+1)/(n+3)))^2)
}
q=c(0.025,0.05,0.95,0.975)
varb1=numeric(4)
for(i in 1:4)
  varb1[i]=varb(q[i])
sqrt(varb1)

## -----------------------------------------------------------------------------
qnorm(q,0,sqrt(6/n))
qb

## -----------------------------------------------------------------------------
library("bootstrap")
 plot(scor$mec,scor$vec)
 plot(scor$mec,scor$alg)
 plot(scor$mec,scor$ana)
 plot(scor$mec,scor$sta)
 plot(scor$vec,scor$alg)
 plot(scor$vec,scor$ana)
 plot(scor$vec,scor$sta)
 plot(scor$alg,scor$ana)
 plot(scor$alg,scor$sta)
 plot(scor$ana,scor$sta)

## -----------------------------------------------------------------------------
cor(scor)

## -----------------------------------------------------------------------------
B=1e3
set.seed(123)
rho12=cor(scor$mec,scor$vec)
rho12star=rho34star=rho35star=rho45star=numeric(B)
for(b in 1:B){
  x1star=sample(1:nrow(scor),replace=TRUE)
  rho12star[b]=cor(scor$mec[x1star],scor$vec[x1star])
  rho34star[b]=cor(scor$alg[x1star],scor$ana[x1star])
  rho35star[b]=cor(scor$alg[x1star],scor$sta[x1star])
  rho45star[b]=cor(scor$ana[x1star],scor$sta[x1star])
}
se12.boot=sd(rho12star)
se34.boot=sd(rho34star)
se35.boot=sd(rho35star)
se45.boot=sd(rho45star)
print(c(se12.boot,se34.boot,se35.boot,se45.boot))

## -----------------------------------------------------------------------------
g=function(x){
  a=mean((x-mean(x))^3)
  b=(mean((x-mean(x))^2))^(3/2)
  return(a/b)
}
m=100
alpha=0.05
CI=function(sample1,sample2,dist){ #The parameter dist is for selecting distribution and CI
b1=b2=numeric(m)
for(j in 1:m){
  x1=sample(sample1,replace=TRUE)
  x2=sample(sample2,replace=TRUE)
  b1[j]=g(x1)
  b2[j]=g(x2)
}
se1=sd(b1)
se2=sd(b2)
thetahat1=g(sample1)
thetahat2=g(sample2)
if(dist==11){
  CI1_stand=c(thetahat1-qnorm(1-alpha/2)*se1,thetahat1-qnorm(alpha/2)*se1)
  return(CI1_stand)
}
else if(dist==12){
  CI2_stand=c(thetahat2-qnorm(1-alpha/2)*se2,thetahat2-qnorm(alpha/2)*se2)
  return(CI2_stand)
}
else if(dist==21){
  CI1_basic=c(2*thetahat1-quantile(b1,probs=1-alpha/2),2*thetahat1-quantile(b1,probs=alpha/2))
  return(CI1_basic)
}
else if(dist==22){
  CI2_basic=c(2*thetahat2-quantile(b2,probs=1-alpha/2),2*thetahat2-quantile(b2,probs=alpha/2))
  return(CI2_basic)
  }
else if(dist==31){
  CI1_perc=c(quantile(b1,probs=alpha/2),quantile(b1,probs=1-alpha/2))
  return(CI1_perc)
}
else if(dist==32){
  CI2_perc=c(quantile(b2,probs=alpha/2),quantile(b2,probs=1-alpha/2))
  return(CI2_perc)
}
else cat("Error,dist should be included in{11,12,21,22,31,32}!","\n")
}

## -----------------------------------------------------------------------------
theta1=0
f=function(x) ((x-5)/sqrt(2*5))^3*dchisq(x,5)
theta2=integrate(f,0,Inf)$value
theta2#the true skewness of chisq

## -----------------------------------------------------------------------------
set.seed(2333)
ita=1000 #times of MC
test1=test2=matrix(numeric(9),ncol=3)#test1 for normal,test2 for chisq
colnames(test1)=colnames(test2)=c("cover","left_miss","right_miss")
rownames(test1)=rownames(test2)=c("normal_CI","basic_CI","perc_CI")
test11=test12=test21=test22=test31=test32=testl11=testl12=testl21=testl22=testl31=testl32=testr11=testr12=testr21=testr22=testr31=testr32=numeric(ita)
for(i in 1:ita){
  s1=rnorm(20)
  s2=rchisq(20,5)
  test11[i]=(CI(s1,s2,11)[1]<=theta1)*(CI(s1,s2,11)[2]>=theta1)
  test12[i]=(CI(s1,s2,12)[1]<=theta2)*(CI(s1,s2,12)[2]>=theta2)
  test21[i]=(CI(s1,s2,21)[1]<=theta1)*(CI(s1,s2,21)[2]>=theta1)
  test22[i]=(CI(s1,s2,22)[1]<=theta2)*(CI(s1,s2,22)[2]>=theta2)
  test31[i]=(CI(s1,s2,31)[1]<=theta1)*(CI(s1,s2,31)[2]>=theta1)
  test32[i]=(CI(s1,s2,32)[1]<=theta2)*(CI(s1,s2,32)[2]>=theta2)
  testl11[i]=(CI(s1,s2,11)[1]>=theta1);testl12[i]=(CI(s1,s2,12)[1]>=theta2);
  testl21[i]=(CI(s1,s2,21)[1]>=theta1);testl22[i]=(CI(s1,s2,22)[1]>=theta2);
  testl31[i]=(CI(s1,s2,31)[1]>=theta1);testl32[i]=(CI(s1,s2,32)[1]>=theta2);
  testr11[i]=(CI(s1,s2,11)[2]<=theta1);testr12[i]=(CI(s1,s2,12)[2]<=theta2);
  testr21[i]=(CI(s1,s2,21)[2]<=theta1);testr22[i]=(CI(s1,s2,22)[2]<=theta2);
  testr31[i]=(CI(s1,s2,31)[2]<=theta1);testr32[i]=(CI(s1,s2,32)[2]<=theta2)
}
test1[1,1]=mean(test11);test1[1,2]=mean(testl11);test1[1,3]=mean(testr11);
test1[2,1]=mean(test21);test1[2,2]=mean(testl21);test1[2,3]=mean(testr21);
test1[3,1]=mean(test31);test1[3,2]=mean(testl31);test1[3,3]=mean(testr31);
test2[1,1]=mean(test12);test2[1,2]=mean(testl12);test2[1,3]=mean(testr12);
test2[2,1]=mean(test22);test2[2,2]=mean(testl22);test2[2,3]=mean(testr22);
test2[3,1]=mean(test32);test2[3,2]=mean(testl32);test2[3,3]=mean(testr32)
#For normal population:
test1
#For $\chi^2\left(5\right)$ population:
test2

## -----------------------------------------------------------------------------
g=function(x){
  a=mean((x-mean(x))^3)
  b=(mean((x-mean(x))^2))^(3/2)
  return(a/b)
}

## -----------------------------------------------------------------------------
f=function(alpha,a){
n=30
m=1000
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
  x=rbeta(n,a,a)
  sktests[i] <- as.integer(abs(g(x)) >= cv)
}
pwr <- mean(sktests)
return(pwr)
}

## -----------------------------------------------------------------------------
set.seed(1)
a=seq(1,20,1)
Pwr=matrix(c(rep(0,60)),nrow=3)
alpha=c(0.05,0.07,0.1)
for(i in 1:3)
for(j in 1:20)
  Pwr[i,j]=f(alpha[i],a[j])
Pwr

## -----------------------------------------------------------------------------
f2=function(alpha,v){
n=30
m=1000
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
  x=rt(n,v)
  sktests[i] <- as.integer(abs(g(x)) >= cv)
}
pwr <- mean(sktests)
return(pwr)
}

## -----------------------------------------------------------------------------
set.seed(2)
v=seq(1,20,1)
Pwr2=matrix(c(rep(0,60)),nrow=3)
alpha=c(0.05,0.07,0.1)
for(i in 1:3)
for(j in 1:20)
  Pwr2[i,j]=f2(alpha[i],v[j])
Pwr2

## -----------------------------------------------------------------------------
set.seed(21)
alpha=0.025
m <- 10000
n=20
test1 <- numeric(m) #test decisions
for (j in 1:m) {
x1 <- rchisq(n,1)
s1=sd(x1)
#test decision is 1 (reject) or 0
l1=1-s1/sqrt(n)*qt(1-alpha/2,n-1)
l2=1+s1/sqrt(n)*qt(1-alpha/2,n-1)
test1[j] <- 1-(mean(x1) >=l1 )*(mean(x1) <=l2)
}
p.reject1 <- mean(test1) #proportion rejected
p.reject1

## -----------------------------------------------------------------------------
set.seed(21)
alpha=0.05
m <- 10000
n=20
test1 <- numeric(m) #test decisions
for (j in 1:m) {
x1 <- rchisq(n,1)
s1=sd(x1)
#test decision is 1 (reject) or 0
l1=1-s1/sqrt(n)*qt(1-alpha/2,n-1)
l2=1+s1/sqrt(n)*qt(1-alpha/2,n-1)
test1[j] <- 1-(mean(x1) >=l1 )*(mean(x1) <=l2)
}
p.reject1 <- mean(test1) #proportion rejected
p.reject1

## -----------------------------------------------------------------------------
set.seed(21)
alpha=0.1
m <- 10000
n=20
test1 <- numeric(m) #test decisions
for (j in 1:m) {
x1 <- rchisq(n,1)
s1=sd(x1)
#test decision is 1 (reject) or 0
l1=1-s1/sqrt(n)*qt(1-alpha/2,n-1)
l2=1+s1/sqrt(n)*qt(1-alpha/2,n-1)
test1[j] <- 1-(mean(x1) >=l1 )*(mean(x1) <=l2)
}
p.reject1 <- mean(test1) #proportion rejected
p.reject1

## -----------------------------------------------------------------------------
set.seed(22)
alpha=0.025
m <- 10000 
n=20
test2 <- numeric(m) #test decisions
for (j in 1:m) {
x2 <- runif(n,min=0,max=2)
s2=sd(x2)
l21=1-s2/sqrt(n)*qt(1-alpha/2,n-1)
l22=1+s2/sqrt(n)*qt(1-alpha/2,n-1)
test2[j] <- 1-(mean(x2) >=l21 )*(mean(x2) <=l22)
}
p.reject2 <- mean(test2) #proportion rejected
p.reject2

## -----------------------------------------------------------------------------
set.seed(22)
alpha=0.05
m <- 10000 
n=20
test2 <- numeric(m) #test decisions
for (j in 1:m) {
x2 <- runif(n,min=0,max=2)
s2=sd(x2)
l21=1-s2/sqrt(n)*qt(1-alpha/2,n-1)
l22=1+s2/sqrt(n)*qt(1-alpha/2,n-1)
test2[j] <- 1-(mean(x2) >=l21 )*(mean(x2) <=l22)
}
p.reject2 <- mean(test2) #proportion rejected
p.reject2

## -----------------------------------------------------------------------------
set.seed(22)
alpha=0.1
m <- 10000 
n=20
test2 <- numeric(m) #test decisions
for (j in 1:m) {
x2 <- runif(n,min=0,max=2)
s2=sd(x2)
l21=1-s2/sqrt(n)*qt(1-alpha/2,n-1)
l22=1+s2/sqrt(n)*qt(1-alpha/2,n-1)
test2[j] <- 1-(mean(x2) >=l21 )*(mean(x2) <=l22)
}
p.reject2 <- mean(test2) #proportion rejected
p.reject2

## -----------------------------------------------------------------------------
set.seed(23)
alpha=0.025
m <- 10000 
n=20
test3 <- numeric(m) #test decisions
for (j in 1:m) {
x3 <- rexp(n,1)
s3=sd(x3)
l31=1-s3/sqrt(n)*qt(1-alpha/2,n-1)
l32=1+s3/sqrt(n)*qt(1-alpha/2,n-1)
test3[j] <- 1-(mean(x3) >=l31 )*(mean(x3) <=l32)
}
p.reject3 <- mean(test3) #proportion rejected
p.reject3

## -----------------------------------------------------------------------------
set.seed(23)
alpha=0.05
m <- 10000 
n=20
test3 <- numeric(m) #test decisions
for (j in 1:m) {
x3 <- rexp(n,1)
s3=sd(x3)
l31=1-s3/sqrt(n)*qt(1-alpha/2,n-1)
l32=1+s3/sqrt(n)*qt(1-alpha/2,n-1)
test3[j] <- 1-(mean(x3) >=l31 )*(mean(x3) <=l32)
}
p.reject3 <- mean(test3) #proportion rejected
p.reject3

## -----------------------------------------------------------------------------
set.seed(23)
alpha=0.1
m <- 10000 
n=20
test3 <- numeric(m) #test decisions
for (j in 1:m) {
x3 <- rexp(n,1)
s3=sd(x3)
l31=1-s3/sqrt(n)*qt(1-alpha/2,n-1)
l32=1+s3/sqrt(n)*qt(1-alpha/2,n-1)
test3[j] <- 1-(mean(x3) >=l31 )*(mean(x3) <=l32)
}
p.reject3 <- mean(test3) #proportion rejected
p.reject3

## -----------------------------------------------------------------------------
library("bootstrap")
n=nrow(scor)
sigma.hat=(n-1)/n*cor(scor)
lambda.hat=eigen(sigma.hat)$values
theta.hat=lambda.hat[1]/sum(lambda.hat)
theta.jack=numeric(n)
for(i in 1:n){
  sigma=(n-2)/(n-1)*cor(scor[-i,])
  lambda=eigen(sigma)$values
  theta.jack[i]=lambda[1]/sum(lambda)
}
#An unbiased estimate of bias
esti.bias=(n-1)*(mean(theta.jack)-theta.hat)
#An estimate of standard error
esti.se=sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
cat("The jackknife estimate of bias is",esti.bias,"\n",
    "The jackknife estimate of standard error is",esti.se,"\n")

## -----------------------------------------------------------------------------
library("DAAG")
attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits
L1 <- lm(magnetic ~ chemical)
plot(chemical, magnetic, main="Linear", pch=16)
yhat1 <- L1$coef[1] + L1$coef[2] * a
lines(a, yhat1, lwd=2)
L2 <- lm(magnetic ~ chemical + I(chemical^2))
plot(chemical, magnetic, main="Quadratic", pch=16)
yhat2 <- L2$coef[1] + L2$coef[2] * a + L2$coef[3] * a^2
lines(a, yhat2, lwd=2)
L3 <- lm(log(magnetic) ~ chemical)
plot(chemical, magnetic, main="Exponential", pch=16)
logyhat3 <- L3$coef[1] + L3$coef[2] * a
yhat3 <- exp(logyhat3)
lines(a, yhat3, lwd=2)
L4 <- lm(magnetic ~ chemical + I(chemical^2)+I(chemical^3))
plot(chemical, magnetic, main="cubic polynomial", pch=16)
yhat4 <- L4$coef[1] + L4$coef[2] * a+ L4$coef[3] * a^2++ L4$coef[4] * a^3
lines(a, yhat4, lwd=2)

## -----------------------------------------------------------------------------
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +J2$coef[3] * chemical[k]^2
e2[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3[k] <- magnetic[k] - yhat3
J4 <- lm(y ~ x + I(x^2)+I(x^3))
yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] +J4$coef[3] * chemical[k]^2+J4$coef[4] * chemical[k]^3
e4[k] <- magnetic[k] - yhat4
}

## -----------------------------------------------------------------------------
c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))

## -----------------------------------------------------------------------------
L2

## -----------------------------------------------------------------------------
adj.r.squa=numeric(4)
adj.r.squa[1]=summary(L1)$adj.r.squared
adj.r.squa[2]=summary(L2)$adj.r.squared
adj.r.squa[3]=summary(L3)$adj.r.squared
adj.r.squa[4]=summary(L4)$adj.r.squared
c(adj.r.squa[1],adj.r.squa[2],adj.r.squa[3],adj.r.squa[4])

## -----------------------------------------------------------------------------
L4

## -----------------------------------------------------------------------------
count5test_new = function(z) {
n = length(z)
x = z[1:(n/2)]
y = z[-(1:(n/2))]
X = x - mean(x)
Y = y - mean(y)
outx = sum(X > max(Y)) + sum(X < min(Y)) 
outy = sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > 5))
}
perm = function(z,R) {
  n = length(z)
  out = numeric(R)
  for (r in 1: R){
      p = sample(1:n ,n ,replace = FALSE)
      out[r] = count5test_new(z[p])
  }
  sum(out)/R
}              
n1 = 20
n2 = 50
mu1 = mu2 = 0
sigma1 = sigma2 = 1
m = 1e3
set.seed(1122)
alphahat = mean(replicate(m, expr={
x = rnorm(n1, mu1, sigma1)
y = rnorm(n2, mu2, sigma2)
x = x - mean(x) #centered by sample mean 
y = y - mean(y)
z = c(x,y)
perm(z,1000) 
})<0.05)
round(alphahat,4)

## -----------------------------------------------------------------------------
library(MASS)
library(Ball)
library(boot)
dCov <- function(x, y) {
x <- as.matrix(x); y <- as.matrix(y)
n <- nrow(x); m <- nrow(y)
if (n != m || n < 2) stop("Sample sizes must agree")
if (! (all(is.finite(c(x, y)))))
stop("Data contains missing or infinite values")
Akl <- function(x) {
d <- as.matrix(dist(x))
m <- rowMeans(d); M <- mean(d)
a <- sweep(d, 1, m); b <- sweep(a, 2, m)
b + M
}
A <- Akl(x); B <- Akl(y)
sqrt(mean(A * B))
}
ndCov2 <- function(z, ix, dims) {
#dims contains dimensions of x and y
p <- dims[1]
q <- dims[2]
d <- p + q
x <- z[ , 1:p] #leave x as is
y <- z[ix, -(1:p)] #permute rows of y
return(nrow(z) * dCov(x, y)^2)
}
n = seq(20,200,20)
m = 100

p.cor = matrix(0,length(n),2)
p.ball = matrix(0,length(n),2)

for (i in 1:length(n)){
 
  p.cor_sim = matrix(0, m, 2)
  p.ball_sim = matrix(0, m, 2)

  for (j in 1:m){
    set.seed(12345+j)
    x = mvrnorm(n[i],mu=rep(0,2),Sigma=diag(rep(1,2)))
    e = mvrnorm(n[i],mu=rep(0,2),Sigma=diag(rep(1,2)))
    y1 = x/4 + e
    y2 = x/4*e
    #model1
    boot.obj1 = boot(data = cbind(x,y1), statistic = ndCov2, R = 99, sim = "permutation", dims = c(2, 2))
    #permutatin: resampling without replacement
    tb1 = c(boot.obj1$t0, boot.obj1$t)
    p.cor_sim[j,1] = mean(tb1>=tb1[1])
    p.ball_sim[j,1] = bcov.test(x,y1,num.permutations=99,seed=j)$p.value
    #model2
    boot.obj2 = boot(data = cbind(x,y2), statistic = ndCov2, R = 99, sim = "permutation", dims = c(2, 2))
    #permutatin: resampling without replacement
    tb2 = c(boot.obj2$t0, boot.obj2$t)
    p.cor_sim[j,2] = mean(tb2>=tb2[1])
    p.ball_sim[j,2] = bcov.test(x,y2,num.permutations =99,seed=j)$p.value
  }
  p.cor[i,] = colMeans(p.cor_sim < 0.05)
  p.ball[i,] = colMeans(p.ball_sim < 0.05 )
}

# plot
plot(x=n, y=p.cor[,1], type="b", pch=20, main="Y = X/4 + e", ylab='power', col='blue', ylim = c(0,1))
legend(x=160, y=0.4, legend = c("distance","ball"), lty=1:2, col = c('blue','red'))
lines(x = n, y=p.ball[,1], type="b", lty=2, col='red')
plot(x=n, y=p.cor[,2], type="b", pch=20, main="Y = X/4 * e", ylab='power', col='blue', ylim = c(0,1))
legend(x=160,y=0.4,legend = c("distance","ball"),lty=1:2, col = c('blue','red'))
lines(x = n, y=p.ball[,2], type="b", lty=2, col='red')

## -----------------------------------------------------------------------------
f=function(x) 0.5*exp(-abs(x))

## -----------------------------------------------------------------------------
rw.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (f(y) / f(x[i-1]))) x[i] <- y 
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 25
rw1 <- rw.Metropolis( sigma[1], x0, N)
rw2 <- rw.Metropolis( sigma[2], x0, N)
rw3 <- rw.Metropolis( sigma[3], x0, N)
rw4 <- rw.Metropolis( sigma[4], x0, N)
#number of candidate points rejected
print(c(rw1$k, rw2$k, rw3$k, rw4$k))

## -----------------------------------------------------------------------------
library(GeneralizedHyperbolic)
refline <- qskewlap(c(.025, .975),param=c(0,1,1))
rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
for (j in 1:4) {
plot(rw[,j], type="l",
xlab=bquote(sigma == .(round(sigma[j],3))),
ylab="X", ylim=range(rw[,j]))
abline(h=refline)
}

## -----------------------------------------------------------------------------
f1 = function(x) log(exp(x))
f2 = function(x) exp(log(x))
x=1:10
f1(x) == f2(x)
f1(x) == x
f2(x) == x

## -----------------------------------------------------------------------------
c(isTRUE(all.equal(f1(x),f2(x))),isTRUE(all.equal(f1(x),x)),isTRUE(all.equal(f2(x),x)))

## -----------------------------------------------------------------------------
f = function(a,k){
  c_k = sqrt(a^2*k/(k+1-a^2))
  c1 = 2*gamma((k+1)/2)
  c2 = sqrt(pi*k)*gamma(k/2)
  f_in = function(x) (1+x^2/k)^(-(k+1)/2)
  c3 = integrate(f=f_in,lower=0,upper=c_k)$value
  return(c1/c2*c3)
}

## -----------------------------------------------------------------------------
k = c(4:25,100);n=length(k);
solution = numeric(n)
for(i in 1:n){
solution[i] = uniroot(function(a){f(a,k=k[i]-1)-f(a,k=k[i])},interval=c(-2,1))$root
}
library(knitr)
kable(cbind(k,solution),format="markdown")

## -----------------------------------------------------------------------------
s = function(a,k){
  q = sqrt(a^2*k/(k+1-a^2))
  s = pt(q,k)
  return(s)
}

## -----------------------------------------------------------------------------
A_k = numeric(n)
for(i in 1:n){
  g2 = function(a) s(a,k[i]-1)-s(a,k[i])
  A_k[i] = uniroot(g2,interval=c(-2,1))$root
}
kable(cbind(k,A_k),format="markdown")

## -----------------------------------------------------------------------------
kable(cbind(k,solution,A_k),format="markdown")

## -----------------------------------------------------------------------------
nA=28;nB=24;nOO=41;nAB=70;
E_log=function(x0,x){p=x[1];q=x[2];p0=x0[1];q0=x0[2];
s=(2*p0/(p0+2*(1-p0-q0))*nA*log(p)+2*q0/(q0+2*(1-p0-q0))*nB*log(q)+2*nOO*log(1-p-q)+2*(1-p0-q0)/(p0+2*(1-p0-q0))*nA*log(2*p*(1-p-q))+2*(1-p0-q0)/(q0+2*(1-p0-q0))*nB*log(2*q*(1-p-q))+nAB*log(2*p*q))
return(-s)}
optim(c(0.3,0.4),E_log,x0=c(0.2,0.3))

## -----------------------------------------------------------------------------
M=numeric(51)
for(i in 1:51){
  M[i]=-optim(c(0.3,0.4),E_log,x0=c(0.2,0.3),control=list(maxit=i))$value
}
plot(M,type="l")

## -----------------------------------------------------------------------------
x=mtcars;mpg=x$mpg;disp=x$disp;wt=x$wt
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

## -----------------------------------------------------------------------------
out1=lapply(formulas,lm)
out1

## -----------------------------------------------------------------------------
out2=vector("list",length(formulas))
for(i in seq_along(formulas)){
  out2[[i]]=lm(formulas[[i]])
}
out2

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
out22=vector("list",length(bootstraps))
for(i in seq_along(bootstraps)){
mpg=bootstraps[[i]]$mpg;disp=bootstraps[[i]]$disp
formula2 <- list(mpg ~ disp)
out22[i]=lapply(formula2,lm)
}
out22

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
lapply(out1,rsq)
lapply(out22,rsq)

## -----------------------------------------------------------------------------
set.seed(12144)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
p <- function(x) x$p.value
out41=sapply(trials,p)

## -----------------------------------------------------------------------------
out42=numeric(length(trials))
for(i in seq_along(trials)){
  out42[[i]]=trials[[i]]$p.value
}

## -----------------------------------------------------------------------------
library(knitr)
kable(cbind(out41[1:25],out42[1:25],out41[26:50],out42[26:50],out41[51:75],out42[51:75],out41[76:100],out42[76:100]),col.names=c("out41_1to25","out42_1to25","out41_25to50","out42_25to50","out41_51to75","out42_51to75","out41_76to100","out42_76to100"))

## -----------------------------------------------------------------------------
set.seed(121451)
fun <- function(x){
  return (x+1);
  }
system.time({
res <- sapply(1:5000000, fun);
})

## -----------------------------------------------------------------------------
library(parallel)
set.seed(12145)
mcsapply=function(x,f){
cl <- makeCluster(getOption("cl.cores", 4))
  time=system.time({
  res <- parLapply(cl, x, f)
  })
 stopCluster(cl)
 return(list(res=res,time=time))
}
mcsapply(1:5000000,fun)$time

## -----------------------------------------------------------------------------
f_r=function(x) 0.5*exp(-abs(x))
rw.Metropolis <- function(sigma, x0, N) {
x <- numeric(N)
x[1] <- x0
u <- runif(N)
k <- 0
for (i in 2:N) {
y <- rnorm(1, x[i-1], sigma)
if (u[i] <= (f_r(y) / f_r(x[i-1]))) x[i] <- y 
else {
x[i] <- x[i-1]
k <- k + 1
}
}
return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('List rwM(double sigma, double x0, int N) {
NumericVector x(N);
x[0]=x0;
int k=0;
for(int i=1;i<N;++i){
double u=runif(1)[0];
double y=rnorm(1,x[i-1],sigma)[0];
if (u<=exp(abs(x[i-1])-abs(y))) {x[i]=y;k=k+1;}
else x[i]=x[i-1];
}
List s=List::create(x,k);
return s;
}')

## -----------------------------------------------------------------------------
set.seed(12252)
N <- 2000
sigma <- c(0.5, 2, 7,16)
x0 <- 25
rw11 <- rw.Metropolis( sigma[1], x0, N)
rw21 <- rw.Metropolis( sigma[2], x0, N)
rw31 <- rw.Metropolis( sigma[3], x0, N)
rw41 <- rw.Metropolis( sigma[4], x0, N)
rw12 <- rwM( sigma[1], x0, N)
rw22 <- rwM( sigma[2], x0, N)
rw32 <- rwM( sigma[3], x0, N)
rw42 <- rwM( sigma[4], x0, N)

## -----------------------------------------------------------------------------
rw1 <- cbind(rw11$x, rw21$x, rw31$x, rw41$x)
rw2 <- cbind(rw12[[1]], rw22[[1]], rw32[[1]], rw42[[1]])
for(i in 1:4){
qqplot(rw1[,i],rw2[,i],xlab="R",ylab="cpp")
abline(0,1,col="red")
}

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(rw11 <- rw.Metropolis( sigma[1], x0, N),
rw21 <- rw.Metropolis( sigma[2], x0, N),
rw31 <- rw.Metropolis( sigma[3], x0, N),
rw41 <- rw.Metropolis( sigma[4], x0, N),
rw12 <- rwM( sigma[1], x0, N),
rw22 <- rwM( sigma[2], x0, N),
rw32 <- rwM( sigma[3], x0, N),
rw42 <- rwM( sigma[4], x0, N))
summary(ts)[,c(1,3,5,6)]


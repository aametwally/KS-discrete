## Simulate data from Zero Inflated Poisson, NB, Beta Binomial, Beta Negative Binomial, Half Normal##
sample.zero.inflated <- function(n, phi, theta, dist="Poisson") {
  fsample <- function(m) { rpois(m, lambda=theta[1]); };
  if(dist=="Negative Binomial") fsample <- function(m) {rnbinom(m, size=theta[1], prob=1-theta[2])};#size=r, p=p(failure), so our p is (1-p)
  if(dist=="Beta Binomial") fsample <- function(m) {rbbinom(m, theta[1], theta[2], theta[3])} #{extraDistr}
  if(dist=="Beta Negative Binomial") fsample <- function(m) {rbnbinom(m, theta[1], theta[2], theta[3])} #{extraDistr}
  if(dist=="Half Normal") fsample <- function(m) {rhnorm(m, theta[1])}
  if(dist=="Log Normal") fsample <- function(m) {rlnorm(m, theta[1], theta[2])} #theta[2]=sigma, not sig squar
  ans=rbinom(n, size=1, prob=1-phi);
  m=sum(ans==1);
  if(m>0) ans[ans==1]=fsample(m);
  ans;
}
# For examle:
x= sample.zero.inflated(n=300, phi=0.03, theta=0.1, dist="Poission")
x= sample.zero.inflated(n=10000, phi=0.2, theta=c(15,0.3), dist="Negative Binomial")
x= sample.zero.inflated(10000, phi=0.50, theta=c(20, 2,5), dist="Beta Binomial")
x= sample.zero.inflated(50000, phi=0.50, theta=c(5, 1,10), dist="Beta Negative Binomial")
x= sample.zero.inflated(500, phi=0.50, theta=4, dist="Half Normal")
x= sample.zero.inflated(400, phi=0.50, theta=c(1,4), dist="Log Normal")

######### Next: Simulate from Hurdle ditributions####

##### Simulate from Hurdle Poisson###########
rhpois <- function(n, phi, theta) {
  ans=rbinom(n, size=1, prob=1-phi);
  m=sum(ans==1);
  p0=exp(-theta) 
  M=ceiling(m+m*p0+3*sqrt(m*p0*(1-p0)))
  z=rpois(M, theta)
  u=z[z>0]
  t=length(u)
  if(t < m ) {
    u1=rep(0, m-t);
    itemp=0;
    while(itemp<(m-t)) {
      temp=rpois(1,theta);
      if(temp>0) {itemp=itemp+1; u1[itemp]=temp; };
    };
    u=c(u,u1);
  };
  ans[ans==1]=u[1:m];
  ans;
}

##### Simulate from Hurdle NB###########
rhnb <- function(n, phi, theta=c(theta[1], theta[2])) { # function(m) {rnbinom(m, size=theta[1], prob=1-theta[2])}
  ans=rbinom(n, size=1, prob=1-phi);
  m=sum(ans==1);
  p0=(1-theta[2])^(theta[1])          #(1-p) ^r
  M=ceiling(m+m*p0+3*sqrt(m*p0*(1-p0)))
  z=rnbinom(M, size=theta[1], prob=1-theta[2])
  u=z[z>0]
  t=length(u)
  if(t < m ) {
    u1=rep(0, m-t);
    itemp=0;
    while(itemp<(m-t)) {
      temp=rnbinom(1,size=theta[1], prob=1-theta[2]);
      if(temp>0) {itemp=itemp+1; u1[itemp]=temp; };
    };
    u=c(u,u1);
  };
  ans[ans==1]=u[1:m];
  ans;
} 

##### Simulate from Hurdle Beta Binomial###########
rhbb <- function(n, phi, theta=c(theta[1], theta[2], theta[3])) { # function(m) {rnbinom(m, size=theta[1], prob=1-theta[2])}
  ans=rbinom(n, size=1, prob=1-phi);
  m=sum(ans==1);
  p0=gamma(theta[1]+theta[3])*gamma(theta[2]+theta[3])/(gamma(theta[1]+theta[2]+theta[3])*gamma(theta[3]))
  M=ceiling(m+m*p0+3*sqrt(m*p0*(1-p0)))
  z=rbbinom(M, theta[1], theta[2], theta[3])  
  u=z[z>0]
  t=length(u)
  if(t < m ) {
    u1=rep(0, m-t);
    itemp=0;
    while(itemp<(m-t)) {
      temp=rbbinom(1,theta[1], theta[2], theta[3]);
      if(temp>0) {itemp=itemp+1; u1[itemp]=temp; };
    };
    u=c(u,u1);
  };
  ans[ans==1]=u[1:m];
  ans;
} 


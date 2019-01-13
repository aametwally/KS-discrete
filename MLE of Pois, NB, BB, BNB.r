##MLE of Poission#####
pois.mle=function(x){
  n=length(x)
  S.x= sum(x)
  lambda=S.x/n
  lambda
}
#x=rpois(1000,5) # To check
#pois.mle(x)
#########MLE of Negative Binomial, Initial parameters is required###
nb.mle=function(x,z){
  N=length(x)
  r<- z[1]
  p<- z[2]
  xsum=sum(x)
  neg.log.lik=function(y) {
    r=y[1];
    p=y[2];
    if((p<=0)||(p>1)||(r<0)) ans=-Inf else {
    ans= N*r*log(1-p)- N*lgamma(r) +xsum*log(p) + sum(lgamma(x+r))-sum(lgamma(x+1));
    }
    -ans;
  }
  gp<- function(g) {
    r=g[1];
    p=g[2];
    dr<- N*log(1-p)-N*digamma(r)+sum(digamma(x+r));
    dp<- (xsum/p)-(N*r/(1-p));
    -c(dr,dp)
  }
  
  estimate=optim(par=c(r,p),fn=neg.log.lik, gr=gp, method = "L-BFGS-B", lower = c(1, 1e-8), upper = c(1e9,1-1e-8))
  ans=estimate$par
  list(r=ans[1], p=1-ans[2])
}
#x=rnbinom(1000, 7, 0.3)
#nb.mle(x,c(5, 0.5))

######## MLE of Beta Binomial####
bb.mle=function(x, b0){
  n=b0[1]; a=b0[2]; b=b0[3]
  N=length(x)
  neg.log.lik=function(y) {
    n=y[1];
    a=y[2];
    b=y[3];
    if((n<max(x))||(a<=0)||(b<=0)) ans=-Inf else {
      ans= N*lgamma(n+1)+N*lgamma(a+b)-N*lgamma(a)-N*lgamma(b)-N*lgamma(a+n+b)+sum(lgamma(x+a))+sum(lgamma(n-x+b))-sum(lgamma(x+1))-sum(lgamma(n-x+1));
    };
    -ans;
  }
  gp<- function(z) {
    n=z[1];
    a=z[2];
    b=z[3];
    dn<- N*digamma(n+1)-N*digamma(a+n+b)+sum(digamma(n-x+b))-sum(digamma(n-x+1));
    da<- N*digamma(a+b)-N*digamma(a)-N*digamma(a+n+b)+sum(digamma(x+a));
    db= N*digamma(a+b)-N*digamma(b)-N*digamma(a+n+b)+sum(digamma(n-x+b));
    -c(dn,da,db)
  }
  estimate=optim(par=c(n,a,b),fn=neg.log.lik, gr=gp, method = "BFGS")
  ans=estimate$par
  list(n=ans[1], a=ans[2], b=ans[3])
}
#x=rbbinom(100,10,1,2)
#bb.mle(x,c(max(x)+1, 2,5))

######## MLE of Beta Negative Binomial####
bnb.mle=function(x, b0){  #initial value is required
  r=b0[1]; a=b0[2]; b=b0[3]
  N=length(x)
  neg.log.lik=function(y) {
    r=y[1];
    a=y[2];
    b=y[3];
    if((r<=0)||(a<=0)||(b<=0)) ans=-Inf else {
      ans= N*lgamma(a+r)+N*lgamma(a+b)-N*lgamma(r)-N*lgamma(a)-N*lgamma(b)+sum(lgamma(x+r))+sum(lgamma(x+b))-sum(lgamma(x+1))-sum(lgamma(a+r+b+x));
    };
    -ans;
  }
  gp<- function(z) {
    r=z[1];
    a=z[2];
    b=z[3];
    dr<- N*digamma(a+r)-N*digamma(r)+sum(digamma(r+x))-sum(digamma(a+r+b+x));
    da<- N*digamma(a+r)+N*digamma(a+b)-N*digamma(a)-sum(digamma(a+r+b+x));
    db= N*digamma(a+b)-N*digamma(b)+sum(digamma(b+x))-sum(digamma(a+r+b+x));
    -c(dr,da,db)
  }
  estimate=optim(par=c(r,a,b),fn=neg.log.lik, gr=gp, method = "BFGS")
  ans=estimate$par
  list(r=ans[1], a=ans[2], b=ans[3])
}
## Remark: This code require large sample size. i.e, N.100000
#x=rbnbinom(1000, 2,1,5)
#bnb.mle(x,c(4,2,5))

# x=rbnbinom(100000, 2,1,5)
# bnb.mle(x,c(4,2,5))

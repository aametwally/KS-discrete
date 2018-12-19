
#################################### MLE of ZIBB################

zibb.mle=function(x, b0){ #z=1/a, w=1/b
  n=b0[1]; a=b0[2]; b=b0[3]; #z=b0[4]; w=b0[5];
  N=length(x)
  t= x[x > 0]
  m=length(t)
  tmax=max(t)
  neg.log.lik=function(y) {
    n=y[1];
    a=y[2];
    b=y[3];
    if((n<tmax)||(a<=0)||(b<=0)) ans=-Inf else {
      logA=lgamma(a+n+b)+lgamma(b);
      logB=lgamma(a+b)+lgamma(n+b);
      ans=-m*log(1-exp(logB-logA))-m*logA +m*lgamma(n+1)+sum(lgamma(t+a))+sum(lgamma(n-t+b))+m*lgamma(a+b)-sum(lgamma(t+1))-sum(lgamma(n-t+1))-m*lgamma(a);
    };
    -ans;
  }
  gp<- function(x) {
    n=x[1];
    a=x[2];
    b=x[3];
    logA=lgamma(a+n+b)+lgamma(b);
    logB=lgamma(a+b)+lgamma(n+b);
    dn=m*exp(logB-logA)*(digamma(n+b)-digamma(a+n+b))/(1-exp(logB-logA))+m*digamma(n+1)+sum(digamma(n+b-t))-sum(digamma(n-t+1))-m*digamma(a+n+b);
    da=m*exp(logB-logA)*(digamma(a+b)-digamma(a+n+b))/(1-exp(logB-logA))+sum(digamma(t+a))+m*digamma(a+b)-m*digamma(a+n+b)-m*digamma(a);       
    db=m*exp(logB-logA)*(digamma(a+b)+digamma(n+b)-digamma(a+n+b)-digamma(b))/(1-exp(logB-logA))-m*digamma(b)+sum(digamma(n-t+b))+m*digamma(a+b)-m*digamma(a+n+b);
    -c(dn, da, db);
  }
  estimate1=optim(par=c(n,a,b),fn=neg.log.lik, gr=gp, method = "BFGS")
  ans=estimate1$par
  fvalue=estimate1$value
  if (m/N <= 1-(beta(ans[2],ans[1]+ans[3])/beta(ans[2],ans[3]))) {
    epsilon=m/N
    phi=1-epsilon/(1-(beta(ans[2],ans[1]+ans[3])/beta(ans[2],ans[3])));
    LRT.fvalue=-1*fvalue+(N-m)*log(1-epsilon)+m*log(epsilon)  # fvalue is neg log likelihood
  } else{
    neg.log.lik1=function(y) {
        n=y[1];
        a=y[2];
        b=y[3];
        if((n<tmax)||(a<=0)||(b<=0)) ans=-Inf else {
          if (m/N <= 1-(beta(a,n+b)/beta(a,b))) {
            ans=-(N-m)*log((N-m)/N)-m*log(m/N)+neg.log.lik(y);
          } else {
            ans=-(N-m)*log(beta(a,n+b)/beta(a,b))-m*log(1-beta(a,n+b)/beta(a,b)) + neg.log.lik(y);
          };
        };
        ans;
    }
    gp1<- function(x) {
        n=x[1];
        a=x[2];
        b=x[3];
        if (m/N <= 1-(beta(a,n+b)/beta(a,b))) {
          ans=gp(x);
        } else {
          dn=(N-m)*(digamma(n+b)-digamma(a+n+b))-m*digamma(a+n+b)+ m*digamma(n+1)+sum(digamma(n-t+b))-sum(digamma(n-t+1));
          da=(N-m)*(digamma(a+b)-digamma(a+n+b))-m*digamma(a+n+b)+sum(digamma(t+a))+m*digamma(a+b)-m*digamma(a); 
          db=(N-m)*(digamma(n+b)+digamma(a+b)-digamma(a+n+b)-digamma(b))-m*digamma(a+n+b)-m*digamma(b)+sum(digamma(n-t+b))+m*digamma(a+b);
          ans=-c(dn, da, db);
        };
        ans;
      }
    estimate1=optim(par=c(n,a,b),fn=neg.log.lik1, gr=gp1, method = "BFGS")
    ans=estimate1$par
    fvalue=estimate1$value  #The same for LRT
    LRT.fvalue=-1*fvalue
    epsilon=min(m/N, 1-(beta(ans[2],ans[1]+ans[3])/beta(ans[2],ans[3])));
    phi=1-epsilon/(1-(beta(ans[2],ans[1]+ans[3])/beta(ans[2],ans[3])));
  }    
  list(phi=phi, n=ans[1], alpha=ans[2], beta=ans[3], fvalue=fvalue, LRT.fvalue=LRT.fvalue)
}

##########KS discrete test by Bootstrapping Simulation #########

ks.zibb.boot <- function(x, nsim=500, boot=F) {
  n.hat.sim=rep(NA,nsim)
  a.hat.sim=rep(NA,nsim)
  b.hat.sim=rep(NA,nsim)
  phi.sim=rep(NA,nsim)
  D2.sim=rep(NA,nsim)
  fvalue=rep(NA,nsim)
  N=length(x)
  
  ans.x=zibb.mle(x, c(max(x)+1, 1, 10))
  y= stepfun(0:max(x), c(0, ans.x$phi+(1-ans.x$phi)*pbbinom(0:max(x), size=round(ans.x$n), alpha=ans.x$alpha, beta=ans.x$beta)))
  Dn0=disc_ks_test(x=x, y=y, exact=T, tol=1e-08)$statistic
  
  for(k in 1:nsim){
    if(boot) {
      xboot=sample(x, size=N, replace=T);
      ans=zibb.mle(xboot, c(max(xboot)+1, 1, 10));
      n.hat.sim[k]=ans$n
      a.hat.sim[k]=ans$alpha
      b.hat.sim[k]=ans$beta
      phi.sim[k]=ans$phi
      Xtemp=sample.zero.inflated(N, phi=phi.sim[k], theta=c(ceiling(n.hat.sim[k]),a.hat.sim[k], b.hat.sim[k]), dist="Beta Binomial") 
    } else Xtemp=sample.zero.inflated(N, phi=ans.x$phi, theta=c(ceiling(ans.x$n),ans.x$a, ans.x$b), dist="Beta Binomial"); 
    ans=zibb.mle(Xtemp, c(max(Xtemp)+1, 1, 10))
    y1= stepfun(0:max(Xtemp), c(0, ans$phi+(1-ans$phi)*pbbinom(0:max(Xtemp), size=round(ans$n), alpha=ans$alpha, beta=ans$beta)))
    D2.sim[k]=disc_ks_test(x=Xtemp, y=y1, exact=T, tol=1e-08)$statistic;
  }
  pvalue=(sum(D2.sim > Dn0) + 1)/(nsim + 1);
  pvalue
}
#ks.zibb.boot(x)
#ks.zibb.boot(x, boot=T)

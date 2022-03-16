

symarfima.sim <-  function(d,phi=NULL,theta=NULL, alpha=0, beta=NULL, xreg=NULL, n, index1, varphi=1, burn.in = 5000)  
{ 
 
  nef=n+burn.in
  model = list(dfrac=d)
  if(!is.null(theta)) model$theta=-theta
  if(!is.null(phi)) model$phi=phi
  resp <- rt(n=nef,df=index1)
  scalevar <- index1/(index1-2)
  serie <- arfima::arfima.sim(n=nef,model=model, rand.gen = function(nef, ...) sqrt(scalevar*varphi)*resp)
  
  if(is.null(phi)){
    serie = serie + alpha
  }
  else{
    if(abs(1-sum(phi))>0.05) serie = serie + alpha/(1-sum(phi))  
  }
  if(length(beta)>0){
    serie=serie+xreg*beta
  }
  return(serie[burn.in+(1:n)])
}

##############################################################################
loglik_fixed <- function(par, p, q, b, y, m, family="Gstudent", index1, index2=NULL, xreg=NULL , fixed.ma, fixed.ar, fit.d, fit.alpha )
{
  
  npar <- length(par)
  n = length(y)
  # the variable names are somewhat misleading below
  # these are the coefficients we want to estimate.
  #  we should have used len.fit.ar etc. but it does what it is supposed to.
  len.fixed.ar = sum(is.na(fixed.ar))  
  len.fixed.ma = sum(is.na(fixed.ma))
  
  
  if (p == 0) {phi = 0
  } else { #  AR present in the model
    if(len.fixed.ar==0) { phi <- par[1:p] 
    } else{ 
      phi <- fixed.ar
      where.fix.ar = which(is.na(fixed.ar))
      phi[where.fix.ar] = par[1:len.fixed.ar]
    }
  }
  
  if (q == 0) {
    theta = 0
  } else { # MA present in the model
    if(len.fixed.ma==0){ # if MA
      if(len.fixed.ar==0) {theta <- par[p+(1:q)] # MA and AR 
      }else{theta <- par[len.fixed.ar+(1:q)]  } # MA but AR is fixed
    }else{
      theta <- fixed.ma
      where.fix.ma = which(is.na(fixed.ma))
      theta[where.fix.ma] = par[len.fixed.ar+(1:len.fixed.ma)]
    }
  }
  
  d <- par[npar-fit.alpha-b-1]
  
  if (b == 0) {
    beta = 0
  } else {    
    beta <- par[(npar-fit.alpha-b):(npar-fit.alpha-1)] 
  }
  alpha <- ifelse(fit.alpha,par[npar-1],0)
  varphi <- par[npar] 
  
  vc <- vc.f(d,theta,m)
  if(is.null(xreg)){ 
    vx = c(rep(0, n))
  }else{
    xreg <- matrix(xreg,ncol = b) 
    vx = xreg%*%beta
  }
  mut=rep(0,n)
  
  
  for( t in 1:n ){
    mut[t]=alpha + vx[t] + phi%*%const(t,p,y-vx)+ vc%*%const(t,m,y-mut)
  }
  
  z <- (y-mut[1:n])/sqrt(varphi)
  g0 <- log(index1^(index1/2)*gamma(0.5+index1/2)/(gamma(0.5)*gamma(index1/2))*(index1+z^2)^(-0.5*(index1+1)))  
  loglik <- -sum(-0.5*log(varphi) + g0)
  return(loglik)
}


##############################################################################

symarfima.fit <- function(y, p = NULL, q = NULL,  family = "Gstudent", start.Sym = NULL, index1 = NULL, index2 = NULL, xreg = NULL, 
                          xreg_hat=NULL, h = NULL, m = 100, fixed.ar = NULL, fixed.ma = NULL, fit.d = T, fit.alpha = T, cons = F, max.fun = 1000, 
                          h.method=1) 
{
  
  p=ifelse(is.null(p),0,p)
  q=ifelse(is.null(q),0,q)
  
  if(!is.null(xreg)){if("numeric"%in%class(xreg)) xreg=as.matrix(xreg)}
  b=ifelse(is.null(xreg),0,ncol(xreg))
  
  if(is.null(start.Sym)){
    if(is.null(fixed.ar) & is.null(fixed.ma) & fit.d & fit.alpha){
      
      if (is.null(start.Sym))  start.Sym = c(rep(0,p+q), 0.1, rep(0,b), mean(y), 1) 
    }else{
      if(is.null(fixed.ar)) {start.Sym = rep(0,p)
      } else{ 
        start.Sym = rep(0 , sum(is.na(fixed.ar)))
      } 
      if(is.null(fixed.ma)) {start.Sym = c(start.Sym,rep(0,q))
      } else{
        start.Sym = c(start.Sym, rep(0 , sum(is.na(fixed.ma))))
      }
      if(fit.d) start.Sym = c(start.Sym,0.1) 
      if(!is.null(xreg)) start.Sym = c(start.Sym, rep(0 , b))
      if(fit.alpha) start.Sym = c(start.Sym,mean(y)) 
      start.Sym=c(start.Sym,1)
    }
  }
  
  
  len.fixed.ar = len.fixed.ma = 0
  ss.len = length(start.Sym)
  
  
  if(is.null(fixed.ar)){
    if(p!=0){names(start.Sym)[1:p] <-  paste(rep("phi", p), 1:p, sep = "")
    len.fixed.ar = p}
  }else{
    where.fix.ar = which(is.na(fixed.ar))
    len.fixed.ar = sum(is.na(fixed.ar))
    names(start.Sym)[1:len.fixed.ar] <-  paste(rep("phi", len.fixed.ar), where.fix.ar , sep = "")
  } 
  
  if(is.null(fixed.ma)) {
    if(q!=0){ names(start.Sym)[len.fixed.ar+(1:q)] <- paste(rep("theta", q), 1:q, sep = "")
    len.fixed.ma = q}
  }else{
    where.fix.ma = which(is.na(fixed.ma))
    len.fixed.ma = sum(is.na(fixed.ma))
    names(start.Sym)[len.fixed.ar +(1:len.fixed.ma)] <-  paste(rep("theta", len.fixed.ma), where.fix.ma , sep = "")
  }
  
  if(fit.d) names(start.Sym)[len.fixed.ar+len.fixed.ma+1] <- "d"
  if (!is.null(xreg)) names(start.Sym)[(len.fixed.ar+len.fixed.ma+1+(1:b))] <- paste(rep("X", b), 1:b, sep = "")
  if(fit.alpha) names(start.Sym)[ss.len-1] <- "alpha"
  names(start.Sym)[ss.len] <- "varphi"
  
  loglik.temp = function(par) loglik_fixed(par, p=p,q=q, b=b, y=y, m=m,index1=index1 , 
                                           index2=NULL, xreg=xreg , family=family,fixed.ma=fixed.ma, 
                                           fixed.ar=fixed.ar, fit.d=fit.d, fit.alpha=fit.alpha)
  
  # checking if the log-likelihood at the starting point is finite
  s = try(loglik.temp(start.Sym), silent = T)
  if(class(s)!="numeric"){
    start.Sym = rep(0, ss.len)
    start.Sym[ss.len] = 1
  }
  
  flagwrn = F
  
  if(!cons){
    fit <- optim(par = start.Sym, fn=loglik.temp,method = "Nelder-Mead",hessian = FALSE) 
    if(fit.d == T & abs(fit$par["d"])>0.5) {
      cons = T 
      flagwrn=T
    }
    if(fit$par["varphi"]<0){
      cons = T 
      flagwrn=T
    }
  }
  if(flagwrn) print("Unconstrained optimization failed, restarting optimization using a constrained approach", sep="\n \n ")
  if(cons){
    low_construct = c(rep(-0.99,len.fixed.ar+len.fixed.ma),rep(-0.49, fit.d),rep(-Inf,b),rep(min(y), fit.alpha) ,0.0001)
    upper_construct = c(rep(0.99,len.fixed.ar+len.fixed.ma),rep(0.49, fit.d),rep(Inf,b),rep(max(y), fit.alpha),Inf)
    
    fit = try(lme4::Nelder_Mead(fn=loglik.temp, par = start.Sym, lower = low_construct, upper = upper_construct, 
                                control = list(warnOnly = TRUE, maxfun = max.fun)), silent = TRUE)
    names(fit$par)=names(start.Sym)
  }
  
  if(h.method == 1) {fit$hessian = numDeriv::hessian(func=loglik.temp, x=fit$par, method="Richardson")}
  if(h.method == 2){fit$hessian = numDeriv::hessian(func=loglik.temp, x=fit$par, method="complex")}
  
  
  
  print(fit$par, sep="\n")
  print("------", sep="\n")
  
  par = fit$par
  n = length(y)
  
  if (p == 0) {phi = 0
  } else { #  AR present in the model
    if(is.null(fixed.ar)) { phi <- par[1:p] # AR not is fix
    } else{ 
      phi <- fixed.ar
      phi[where.fix.ar] = par[1:len.fixed.ar]
    }
  }
  
  if (q == 0) {
    theta = 0
  } else { # MA present in the model
    if(is.null(fixed.ma)){ # if MA 
      if(is.null(fixed.ar)) {theta <- par[p+(1:q)] # MA  and AR 
      }else{theta <- par[len.fixed.ar+(1:q)]  } # MA but AR is fixed
    }else{
      theta <- fixed.ma
      where.fix.ma = which(is.na(fixed.ma))
      theta[where.fix.ma] = par[len.fixed.ar+(1:len.fixed.ma)]
    }
  }
  
  d <- par[ss.len-fit.alpha-b-1]
  
  if (b == 0) {
    beta = 0
  } else {    
    beta <- par[(ss.len-fit.alpha-b):(ss.len-fit.alpha-1)] 
  }
  alpha <- ifelse(fit.alpha,par[ss.len-1],0)
  varphi <- par[ss.len] 
  
  vc <- vc.f(d,theta,m)
  
  if(is.null(xreg)){ 
    vx = c(rep(0, n))
  }else{
    xreg <- matrix(xreg,ncol = b) 
    vx = xreg%*%beta
  }
  mut=rep(0,n)
  
  for( t in 1:n ){
    mut[t]=alpha + vx[t] + phi%*%const(t,p,y-vx)+ vc%*%const(t,m,y-mut)
  }
  
  
  if(!is.null(h)){
    y.ar = c(y[n+((1-p):0)],rep(0,h)) 
    r.ck = c(y[n+((1-m):0)]-mut[n+((1-m):0)],rep(0,h))
    y.prev = rep(0,h)
    if(is.null(xreg)){ 
      hvx = c(rep(0, p+h))
    }else{
      hxreg <- matrix(rbind(xreg[n+((1-p):0)], xreg_hat),ncol = b) 
      hvx = hxreg%*%beta
    }
    
    for(t in 1:h){
      y.prev[t] = alpha + hvx[p+t] + phi%*%const(t,p,y.ar-hvx) + vc%*%const(t,m,r.ck)
      y.ar[p+t] = y.prev[t]
    }
  }# close h
  
  #--------------------------------------------
  vcov1 <- try(solve(fit$hessian),silent = T)
  if("try-error"%in%class(vcov1)){
    sd1=666
    zstat1 <- 666
    pvalues1 = 2
    coef1 = fit$par[1:ss.len]# coefficients
  }
  else{
    sd1 <- sqrt(abs(diag(vcov1)))
    coef1 = fit$par[1:ss.len]# coefficients
    zstat1 <- coef1/sd1 # Wald's
    pvalues1 <- 2*(1-pnorm(abs(zstat1))) #
  }
  
  #AIC/BIC
  Negloglik <- round(loglik.temp(fit$par),4) 
  AIC <- round(2*Negloglik + 2*(ss.len) ,4)
  BIC <- round(2*Negloglik + (ss.len)*log(n),4)
  
  model_presentation1 <-cbind(round(coef1,5),round(sd1,5),round(zstat1,5),pvalues1) #return
  colnames(model_presentation1)<-c("Estimate","Std. Error","z value","Pr(>|z|)") #return names
  
  
  coef = printCoefmat(model_presentation1)
  
  print(matrix(nrow=4, ncol=0, byrow=TRUE,dimnames =list(c(
    paste(c("AIC: ",AIC),collapse=""),
    paste(c("BIC: ",BIC),collapse=""),
    paste(c("Log-likelihood: ",-Negloglik),collapse=""),""))))
  
  
  if(!is.null(h)){return(
    list(coef = fit$par, residual=y-mut, mut_hat=mut, prev = y.prev, p=p, q=q, b=b ,xreg=xreg, k=ss.len, 
         AIC=AIC, BIC=BIC, sd.coef=sd1, Loglik=-Negloglik, family=family, Y=y, index1=index1, h=h))
  }else{return(list(coef = fit$par,residual=y-mut, mut_hat=mut, prev=NULL, p=p, q=q, b=b, xreg=xreg, k=ss.len, 
                    AIC=AIC, BIC=BIC, sd.coef=sd1, Loglik=-Negloglik, family=family, Y=y, index1=index1, h=0))  
  }
}




##############################################################################

vc.f <- function(d, theta = NA, infinite = 100)
{
  if(!any(is.na(theta)) )	
    vtheta <- c(1,theta,rep(0,infinite-length(theta)))
  else 
    vtheta = c(1, rep(0, infinite))		
  
  pik<-c(1,rep(0,infinite))  #  pi_k's are the coefficients in the Laurent's expansion in 3.1
  if(d != 0){ 
    pik[2]= d
    for(j in 3:(infinite+1)){
      pik[j]<- ((j-2+d)/(j-1) )*pik[j-1] 
    }
  }
  
  ck= c(1, rep(0,infinite))
  for(i in 2:(infinite+1))
    ck[i]=sum(rev(pik[1:i])*vtheta[1:i])  #ck[i] is c_(i-1)
  
  return(ck[2:(infinite+1)])
}

##############################################################################

const = function(t,mk,plk)
{
  if(mk==0) return(0)
  ymx = rep(0,mk)
  indice = which(t-(1:mk)>0) 
  ymx[indice] = plk[(t-(1:mk))[indice]]
  return(ymx) #retorn  vector y_{t-1}, \cdots, t_{t-k} 
}
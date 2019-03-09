# A wrap-up R program for the block-coordinate-descent algorithm. 
rbcd <- function(y,x,lambda,z=numeric(0),weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000) {
  #bcd_C(y,x,z,n,p,q,lambda,weight,XTol,maxIter,b,gamma)
  nm <- dim(x)
  n <- nm[1]
  p <- nm[2]
  q <- dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  out <- .C("bcd_C",
            y=as.double(y),
            x=as.double(x),
            z=as.double(z),
            n=as.integer(n),
            p=as.integer(p),
            q=as.integer(q),
            lambda=as.double(lambda),
            weight=as.double(weight),
            XTol=as.double(XTol),
            maxIter=as.integer(maxIter),
            theta=as.double(vector("double",n*p)),
            gamma=as.double(vector("double",q))
            )
  tht=out[[11]]
  dim(tht)<-c(n,p)
  gam=out[[12]]
  dim(gam)<-c(q,1)
  return(list(theta=as.matrix(tht),gamma=gam))
}

findbreaks <- function(theta,h=1,minseg=dim(theta)[2]){
  n=dim(theta)[1]
  p=dim(theta)[2]
  thn = apply(theta,1,norm,"2")
  ts = rep(0,n-1)
  for (i in 1:(n-1) ){
    for (j in 1:(n-1) ){
      ts[i]=ts[i]+dnorm((j-i)/h)*thn[j+1]
    }
  }
#  pks<-findpeaks(ts,minpeakdistance = minseg)
  pks<-findpeaks(ts)
  index=seq(2,n-1)
  return(c(1,index[pks[,2]],n+1))
}

fpostest <-function(y,x,regime=c(1,length(y)+1),z=numeric(0)){
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  m=length(regime)-2
  Q=rep(0,n*((m+1)*p+q))
  dim(Q)<-c(n,(m+1)*p+q)
  for(j in 1:(m+1)){
    Q[regime[j]:(regime[j+1]-1),(p*(j-1)+1):(p*j)]<-x[regime[j]:(regime[j+1]-1),]
  }
  if(q>0){
    Q[,((m+1)*p+1):((m+1)*p+q)]<-z
  }
  Sqq=t(Q) %*% Q
  if(cond(Sqq)>1e10){
    Sqq<-Sqq+1e-10*eye((m+1)*p+q)
  }
  rq=t(Q) %*% y
  alpha=mldivide(Sqq,rq)
  resid=y-Q %*% alpha
  ssr=t(resid) %*% resid
  ydm=y-mean(y)
  R2=1-ssr/(t(ydm)%*% ydm)
  s2 = ssr/(n-(m+1)*p-q);
  Sqqinv=mldivide(Sqq,eye((m+1)*p+q))
  Sigma=as.double(s2)*Sqqinv
  return(list(alpha=alpha,ssr=ssr,R2=R2,resid=resid,Sigma=Sigma))
}

# The the tuning parameter that yields no breaks
get_max_lambda <- function(y,x,z=numeric(0)){
  n=length(y)
  p=dim(x)[2]
  q=dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  if(q==0){
    M=eye(n)
  }else{
    M = eye(n)-z %*% mldivide(t(z)%*%z,eye(q))%*%t(z)
  }
  b0=mldivide(t(x) %*% M %*% x,t(x) %*% M %*% y)
  e=M %*% (y - x %*% b0)
  ne=rep(0,n-1)
  for(i in (1:(n-1))){
    xi2 = x[(i+1):n,]
    dim(xi2)<-c(n-i,p)
    xi=rbind(zeros(i,p),xi2)
    ne[i]=norm(t(xi) %*% e,"2")
  }
  return(max(ne))
}

# Estimate the regression with structural changes, using information criterion (IC) to determine the tuning parameter on the group-fused-Lasso penalty.  
regsc_ic <-function(y,x,z=numeric(0),S=20,h=1,weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000){
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  
  lam_max = get_max_lambda(y,x,z)
  lam_min = 0.1*lam_max
  R = (1.1*lam_max-lam_min)/(S-1)
  L = lam_min+R*(0:(S-1))
  IC = rep(0,S)
  K = rep(0,S)
#  THT = matrix(0,n,S*p)
  listTheta = vector("list",S)
  listRegime = vector("list",S)
  for(i in (1:S)){
    tht=rbcd(y,x,L[i],z,weight,XTol,maxIter)$theta
    listTheta[[i]] = tht
#    THT[,(p*(i-1)+1):(p*i)]=tht
    regime=findbreaks(tht,h)
    listRegime[[i]] = regime
    res=fpostest(y,x,regime,z)
    k=length(regime)-2
    IC[i]=log(res$ssr / n) + ((k+1)*p+q)/sqrt(n)
    K[i]=k
  }
  i=which.min(IC)
  theta=as.matrix(listTheta[[i]])
#  theta=as.matrix(THT[,(p*(i-1)+1):(p*i)])
  regime=listRegime[[i]]
  res=fpostest(y,x,regime,z)
  lambda = L[i]
  uqk = uniq(K)
  K = uqk$b
  L = L[uqk$m]
  IC = IC[uqk$m]
  listRegime = listRegime[uqk$m]
  return(list(regime=regime,alpha=res$alpha,Sigma=res$Sigma,R2=res$R2,ssr=res$ssr,resid=res$resid,lambda=lambda,IC=IC,K=K,L=L,listRegime=listRegime))
}

# Estimate the regression with structural changes, using a rule of thumb to determine the tuning parameter on the group-fused-Lasso penalty.  
regsc_rt <-function(y,x,z=numeric(0),h=1,weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000){
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  lambda = 0.618*get_max_lambda(y,x,z)
  tht=rbcd(y,x,lambda,z,weight,XTol,maxIter)$theta
  regime=findbreaks(tht,h)
  res=fpostest(y,x,regime,z)
  return(list(regime=regime,alpha=res$alpha,Sigma=res$Sigma,R2=res$R2,ssr=res$ssr,resid=res$resid,lambda=lambda))
}

# Estimate the regression with structural changes, using a user-supplied tuning parameter on the group-fused-Lasso penalty.  
regsc_us <-function(y,x,z=numeric(0),lambda,h=1,weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000){
  n=dim(x)[1]
  p=dim(x)[2]
  q=dim(z)[2]
  if(is.null(q)){
    q=0
    z=numeric(0)
  }
  tht=rbcd(y,x,lambda,z,weight,XTol,maxIter)$theta
  regime=findbreaks(tht,h)
  res=fpostest(y,x,regime,z)
  return(list(regime=regime,alpha=res$alpha,Sigma=res$Sigma,R2=res$R2,ssr=res$ssr,resid=res$resid,lambda=lambda))
}


# Estimate the regression with structural changes.
regsc <-function(formula,data,lambda=NULL,method="ic",date=seq(1,dim(data)[1])){
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  f <- Formula(formula)
  mf[[1]] <- as.name("model.frame")
  mf$formula <- f
  mf <- eval(mf, parent.frame())
  y <- model.response(mf)
  mt <- terms(f, data = data)
  mtx <- terms(f, data = data, rhs = 1)
  x <- model.matrix(mtx, mf)
  if(length(f)[2] < 2L) {
    mtz <- NULL
    z <- NULL
  } else {
    mtz <- delete.response(terms(f, data = data, rhs = 2))
    xic = attr(mtx, "intercept")
    zic = attr(mtz, "intercept")
    z <- as.matrix(model.matrix(mtz, mf))
    if(xic*zic==1){#There are two intercepts in the formula
      z = as.matrix(z[,-1])
      attr(mtz,"intercept")<-0L
    }
  }
  if(is.numeric(lambda)){
    res <- regsc_us(y,x,z,lambda)
    method = "user supplied"
  }
  if(method=="ic" | method==0 | method=="information criterion"){
    res <- regsc_ic(y,x,z)
    method = "information criterion"
  }
  if(method=="rot" | method==1 | method=="rule of thumb"){
    res <- regsc_rt(y,x,z)
    method = "rule of thumb"
  }
  res$call <- mf
  res$formula <- formula(f)
  res$termsx <- mtx
  res$termsz <- mtz
  res$date <- as.vector(date)
  res$y <- y
  res$x <- as.matrix(x)
  if(!is.null(z)){
    res$z <- as.matrix(z)
  }else{
    res$z=z
  }
  res$method = method
  class(res) <- "regsc"
  return(res)
}

print.regsc <- function(object){
  regime=object$regime
  m=length(regime)-2
  cat("\n")
  print(object$formula)
  cat("\n")
  cat("The method of choosing tuning parameter (lambda): ", object$method,"\n")
  cat("The value of the tuning parameter (lambda): ", object$lambda,"\n")
  if(m==0){
    cat("No structural change has been found. The usual OLS is performed.\n")
  }else{
    cat("We find ", m, " structural change(s) at:\n")
    cat("\t",object$date[regime[2:(m+1)]],"\n")
  }
  cat("\n")
  
}

plot.regsc <- function(object){
  y = object$y
  resid = object$resid
  fitted = y-resid
  par(mfrow=c(1,2))
  plot(fitted,y,main="fitted v.s. actual values")
  plot(fitted,resid,main="fitted v.s. errors")
  
}

summary.regsc <- function(object){
  regime=object$regime
  m=length(regime)-2
  p=dim(object$x)[2]
  if(!is.null(object$z)){
    q=dim(object$z)[2]
    gamma=object$alpha[(p*(m+1)+1):(p*(m+1)+q)]
  }else{
    q=0
    gamma=NULL
  }
  res <- object$resid
  y <- object$y
  n <- NROW(res)
  df.residual = n-(p*(m+1)+q)
    
  rss <- sum(res^2)
  sigma = sqrt(rss/df.residual)
  xic = attr(object$termsx, "intercept")
  zic = attr(object$termsz, "intercept")
  xnames = attr(object$termsx,"term.labels")
  znames = attr(object$termsz,"term.labels")
  if(is.null(zic)){
    zic = 0
  }
  if(xic==1 | zic == 1) {
    tss <- sum((y - mean(y))^2)
    dfi <- 1    
  } else {
    tss <- sum(y^2)
    dfi <- 0
  }
  r.squared <- 1 - rss/tss
  adj.r.squared <- 1 - (1 - r.squared) * ((n - dfi)/df.residual)

  alpha=object$alpha
  Sigma=object$Sigma
  se=sqrt(diag(Sigma))
  if(q==0){
    listCoef = vector("list",m+1)
  }else{
    listCoef = vector("list",m+2)
  }
  for(i in (1:(m+1))){
    a=alpha[((i-1)*p+1):(i*p)]
    ase = se[((i-1)*p+1):(i*p)]
    az = a/ase
    apv = 2*(1-pnorm(abs(az)))
    C=cbind(a,ase,az,apv)
    colnames(C) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    if(xic==1){
      varnames=c("Intercept",xnames)
    }else{
      varnames=xnames
    }
    rownames(C) <- varnames
    listCoef[[i]] = C
  }
  if(q>0){
    g = alpha[((m+1)*p+1):((m+1)*p+q)]
    gse = se[((m+1)*p+1):((m+1)*p+q)]
    gz = g/gse
    gpv = 2*(1-pnorm(abs(gz)))
    C = cbind(g,gse,gz,gpv)
    colnames(C) <- c("Estimate", "Std.Err", "Z value", "Pr(>z)")
    if(zic==1){
      varnames=c("Intercept",znames)
    }else{
      varnames=znames
    }
    rownames(C) <- varnames
    listCoef[[m+2]] = C
  }
    #  listRegime = vector
  
  rval <- list(
    call = object$call,
    formula = object$formula,
    termsx = object$termsx,
    termsz = object$termsz,
    date = object$date,
    x = object$x,
    y = object$y,
    z = object$z,
    method=object$method,
    lambda=object$lambda,
    p = p,
    q = q,
    sigma = sigma,
    residuals = res,
    df.residual = df.residual,
    regime = regime,
    listCoef = listCoef,
#    coefficients = cf,
#    sigma = object$sigma,
    r.squared = r.squared,
    adj.r.squared = adj.r.squared
#    waldtest = waldtest,
#    vcov = vc,
#    diagnostics = diag
  )
  
  class(rval) <- "summary.regsc"
#  invisible(rval)
  return(rval)
}


print.summary.regsc <- function(object){
  regime=object$regime
  m=length(regime)-2
  cat("\n")
  cat("The method of choosing tuning parameter: ", object$method,"\n")
  cat("The value of the tuning parameter (lambda): ", object$lambda,"\n")
  if(m==0){
    cat("No structural change has been found. The usual OLS is performed.\n")
  }else{
    cat("We find ", m, " structural change(s) at:\n")
    cat("\t",object$date[regime[2:(m+1)]],"\n")
  }
  cat("\n")
  for(i in (1:(m+1))){
    cat("Coefficients with Time-Varying Effects, Regime ", i, ":\n")
    cat("From ",object$date[regime[i]]," to ", object$date[regime[i+1]-1],"\n")
    C=object$listCoef[[i]]
    printCoefmat(C,signif.legend=FALSE)
    cat("\n")
  }
  if(object$q>0){
    cat("Coefficients with Time-Invariant Effects :\n")
    C=object$listCoef[[m+2]]
    printCoefmat(C,signif.legend=FALSE)
    cat("\n")
  }
  cat("Residual standard error: ",object$sigma," on ",object$df.residual," degrees of freedom\n")
  cat("R-squared: ", object$r.squared,"\n")
  invisible(object)
}

require(pracma)
require(Formula)
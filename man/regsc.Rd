\name{regsc}
\alias{regsc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Regressions with Structural Changes}
\description{
regsc is used to estimate time-series regressions with an unknown number of structural changes. To estimate structural changes, a penalized least squares problem with the group-fused-Lasso (GFL) penalty is solved. The tuning parameter on the GFL penalty can be selected either by information criterion (default) or by a rule of thumb.
}
\usage{
regsc(formula,data,lambda=NULL,method="ic",date=seq(1,dim(data)[1]))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
  An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
  }
  \item{data}{
  An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).
  }
  \item{lambda}{
  An optional positive numeric, the tuning parameter on the group-fused-Lasso penalty. If lambda is not supplied, it takes \code{NULL} and \code{regsc} would choose a value for \code{lambda} by "information criterion" or "rule of thumb", depending on the value of \code{method}. If a non-positive numeric is supplied, then \code{regsc} would print an error message and go on to estimate the model using a lambda that is selected by information criterion.  
  }
  \item{method}{
  An optional string or integer indicating the method of selecting tuning parameter. "ic" or "information criterion" or 0 indicates the method of information criterion. "rot" or "rule of thumb" or 1 indiciates the method of rule of thumb. The rule of thumb is simply 0.618*maxlambda, where maxlambda is obtained by get_max_lambda(y,x,z). If no method is specified, then "ic" would be used.
  }
  \item{date}{
  An optional vector of dates. If it is not provided, the procedure would use 1:n, where n is the sample size.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of the class "regsc". The function summary can be used to obtain and print a summary of results. An object of class "regsc" is a list containing the following components:
  \item{regime }{a \code{(m+2)}-by-1 numeric vector, containing \code{{1,T1,...,Tm,n+1}}, where \code{T1,...,Tm} are m breaks}
  \item{alpha }{a \code{(p*(m+1)+q)}-by-1 vector of estimated coefficients}
  \item{Sigma}{the estimated covariance matrix for \code{alpha}}
  \item{ssr}{the sum of squared residuals}
  \item{R2}{the overall goodness-of-fit}
  \item{resid}{an \code{n}-by-1 vector of the fitted residual}
  \item{lambda}{the tuning parameter selected either by the information criterion (Qian & Su, 2016) or the rule of thumb}

If the tuning parameter \code{lambda} is selected by the information criterion, "regsc" further returns:
  \item{L}{a numeric vector in ascending order, each element of which is a candidate for the tuning parameter \code{lambda}}
  \item{IC}{a numeric vector of calculated information criteria corresponding to different values of the tuning parameter in \code{L}}
  \item{K}{a vector of integers, each of which is the number of breaks corresponding to different values of the tuning parameter in \code{L}}
  \item{listRegime}{a list of vectors, each of which is the estimated \code{regime} corresponding different values of the tuning parameter in \code{L}}
  
  Note that the repetitive elements of \code{L}, \code{IC}, \code{K}, \code{listRegime} are deleted.
}
  
\references{
Qian, J., L. Su, 2016, "Shrinkage estimation of regression models with multiple structural changes", \emph{Econometric Theory}, 32 (6), 1376-1433.
}
\author{
Junhui Qian and Liangjun Su
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{print}}, \code{\link{summary}}, \code{\link{plot}}
}
\examples{
n=120;p=2;q=3;sigma=0.5

x=rnorm(n*p)
dim(x)<-c(n,p)
z=rnorm(n*q)
dim(z)<-c(n,q)

beta0=c(rep(1,n/2),rep(0.5,n/2))
beta0=rep(beta0,p)
dim(beta0)<-c(n,p)
gamma0=rep(1,q)
dim(gamma0)<-c(q,1)

y = rowSums(x*beta0) + z \%*\% gamma0 + sigma*rnorm(n)

data = as.data.frame(cbind(y,x,z))
colnames(data) <- c("y","x1","x2","z1","z2","z3")
  
res=regsc(y~x1+x2|z1+z2+z3,data)
summary(res)

res=regsc(y~x1+x2|z1+z2+z3,data,method="rule of thumb")
summary(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{ts}

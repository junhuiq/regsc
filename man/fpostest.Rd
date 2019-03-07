\name{fpostest}
\alias{fpostest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Post-Lasso estimation of the regression with structural changes
}
\description{
After estimating the set of break dates, \code{fpostest} estimates the regression using the least squares. The procedure also works for the case where there is no structural change in the model.
}
\usage{
fpostest(y,x,regime,z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An \code{n}-by-1 numeric matrix (vector), the dependent variable}
  \item{x}{An \code{n}-by-\code{p} numeric matrix, the regressors with possibly time-varying effect on \code{y}}
  \item{regime}{a \code{(m+2)}-by-1 numeric vector, containing \code{{1,T1,...,Tm,n+1}}, where \code{T1,...,Tm} are m breaks}
  \item{z}{An optional \code{n}-by-\code{q} numeric matrix, the regressors with time-invariant effect on \code{y}}
}

\details{

}

\value{
A list containing the following components:
  \item{alpha }{a \code{(p*(m+1)+q)}-by-1 vector of estimated coefficients.The first \code{p*(m+1)} elements of \code{alpha} correspond to 
  the time-varying effects of \code{x} on \code{y}. For example, \code{alpha[1:p]} is the estimated coefficient on \code{x} in the first regime.} 
  \item{Sigma}{the estimated covariance matrix for \code{alpha}}
  \item{ssr}{the sum of squared residuals}
  \item{R2}{the overall goodness-of-fit}
  \item{resid}{an \code{n}-by-1 vector of the fitted residual}
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
\code{\link{regsc}}
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

regime=c(1,61,121)
res=fpostest(y,x,regime,z)
}
\keyword{ ~regression }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~ts }% __ONLY ONE__ keyword per line

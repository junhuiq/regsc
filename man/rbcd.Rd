\name{rbcd}
\alias{rbcd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The block-coordinate-descent algorithm
}
\description{
A wrap-up R program for the block-coordinate-descent algorithm}
\usage{
rbcd(y,x,lambda,z=numeric(0),weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An \code{n}-by-1 numeric matrix (vector), the dependent variable}
  \item{x}{An \code{n}-by-\code{p} numeric matrix, the regressors with possibly time-varying effect on \code{y}}
  \item{lambda}{A positive number, the tuning parameter on the group-fused-Lasso penalty}
  \item{z}{An optional n-by-q numeric matrix, the regressors with time-invariant effect on \code{y}}
  \item{weight}{An optional \code{(n-1)}-by-1 numeric vector}
  \item{XTol}{An optional small number, the level of error tolerance}
  \item{maxIter}{An optional integer, the maximum number of iterations allowed}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of two elements: {\code{theta}, \code{gamma}}.     
  \item{theta}{an \code{n}-by-\code{p} numeric matrix. \code{theta[1,]=beta[1,]}, \code{theta[t,]=beta[t,]-beta[t-1,]} for \code{t>1}.}
  \item{gamma}{a \code{q}-by-1 vector (if \code{z} is non-empty). }
%% ...
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
\code{\link{regsc}}}

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

lambda = 10
res=rbcd(y,x,lambda,z)

beta = cumsum(res$theta)
dim(beta)<-c(n,p)
time = 1:n
plot(time,beta[,1],type="l")
}
\keyword{block-coordinate descent}

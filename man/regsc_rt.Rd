\name{regsc_rt}
\alias{regsc_rt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate the regression with structural changes by the group fused Lasso}
\description{
Estimate the regression with structural changes, using a rule of thumb to determine the tuning parameter on the group-fused-Lasso penalty.}
\usage{
regsc_rt(y,x,z=numeric(0),h=1,weight=rep(1,length(y)-1),XTol=1e-6,maxIter=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An \code{n}-by-1 numeric matrix (vector), the dependent variable}
  \item{x}{An \code{n}-by-\code{p} numeric matrix, the regressors with possibly time-varying effect on \code{y}}
  \item{z}{An optional n-by-q numeric matrix, the regressors with time-invariant effect on \code{y}}
  \item{h}{An optional positive number, used in the determination of break dates. }
  \item{weight}{An optional \code{(n-1)}-by-1 numeric vector}
  \item{XTol}{An optional small number, the level of error tolerance}
  \item{maxIter}{An optional integer, the maximum number of iterations allowed}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list containing the following components:
  \item{regime }{a \code{(m+2)}-by-1 numeric vector, containing \code{{1,T1,...,Tm,n+1}}, where \code{T1,...,Tm} are m breaks}
  \item{alpha }{a \code{(p*(m+1)+q)}-by-1 vector of estimated coefficients}
  \item{Sigma}{the estimated covariance matrix for \code{alpha}}
  \item{ssr}{the sum of squared residuals}
  \item{R2}{the overall goodness-of-fit}
  \item{resid}{an \code{n}-by-1 vector of the fitted residual}
  \item{lambda}{the tuning parameter selected either by the information criterion (Qian & Su, 2016) or the rule of thumb}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{ts}

\name{get_max_lambda}
\alias{get_max_lambda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The the tuning parameter that yields no breaks}

\description{
Calculate the the tuning parameter (\code{lmax}) that yields no breaks. To ensure no breaks, one is advised to use 1.1*\code{lmax}.
}
\usage{
get_max_lambda(y,x,z=numeric(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An \code{n}-by-1 numeric matrix (vector), the dependent variable}
  \item{x}{An \code{n}-by-\code{p} numeric matrix, the regressors with possibly time-varying effect on \code{y}}
  \item{z}{An optional n-by-q numeric matrix, the regressors with time-invariant effect on \code{y}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A positive number, the the tuning parameter (\code{lmax}) that yields no breaks. To ensure no breaks, one is advised to use 1.1*\code{lmax}.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ lasso }% use one of  RShowDoc("KEYWORDS")

\name{policyrule}
\alias{policyrule}
\docType{data}
\title{
The Data for Estimating the Taylor's Rule 
}
\description{
Data for estimating the Taylor's Rule for the US monetary policy 
}
\usage{data("policyrule")}
\format{
  A data frame with 200 observations on the following 4 variables.
  \describe{
    \item{\code{quarter}}{a factor of dates, in the form of YYYYQQ.}
    \item{\code{fedrate}}{a numeric vector, the federal funds rate.}
    \item{\code{inflation}}{a numeric vector, the inflation rate.}
    \item{\code{gap}}{a numeric vector, the GDP gap.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
http://dx.doi.org/10.1016/j.jedc.2014.09.034
}

\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(policyrule)
res = regsc(fedrate~inflation+gap,policyrule,date=policyrule$quarter)
summary(res)
}
\keyword{datasets}

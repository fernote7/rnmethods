% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heston_ea.R
\name{hestonea_di}
\alias{hestonea_di}
\title{Drift Interpolated Heston European Call}
\usage{
hestonea_di(S, X, r, v, theta, rho, k, sigma, t = 0, dt = NULL, tau = 1,
  N)
}
\arguments{
\item{S}{Spot price}

\item{X}{Strike price}

\item{r}{Asset's rate of return}

\item{v}{Instantaneous variance}

\item{theta}{Long variance}

\item{rho}{Processes' correlation}

\item{k}{Rate at which v returns to theta}

\item{sigma}{Vol of vol}

\item{t}{Starting time, Default: 0}

\item{dt}{Stepsize, Default: NULL}

\item{tau}{Ending time, Default: 1}

\item{N}{Number of simulations}
}
\value{
List with call price, values used to compute the call and all simulated paths.
}
\description{
Simulate a drift interpolated exact algorithm Heston European call
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[stats]{rchisq}},\code{\link[stats]{rnorm}}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/philrob.R
\name{philrob_jacobian}
\alias{philrob_jacobian}
\title{Jacobian of the Phillips and Robinson's model equations}
\usage{
philrob_jacobian(y, parms = philrob_default_parms())
}
\arguments{
\item{y}{The state}

\item{parms}{Model's parameters (optional, see \code{\link{philrob_default_parms}})}
}
\value{
The jacobian of the model's equations (for constant forcing)
}
\description{
Jacobian of the Phillips and Robinson's model equations
}
\examples{
y <- c(Vv = 1, Vm = 1, H = 1)
J <- philrob_jacobian(y)
}
\references{
Phillips AJK, Robinson PA.
A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
}
\seealso{
\code{\link{philrob}, \link{dPhilrob}}
}
\author{
Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
}

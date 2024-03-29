% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DELTD.R
\name{Gamma}
\alias{Gamma}
\title{Estimate Density Values by Gamma kernel}
\usage{
Gamma(x = NULL, y, k = NULL, h = NULL)
}
\arguments{
\item{x}{scheme for generating grid points}

\item{y}{a numeric vector of positive values}

\item{k}{number of gird points}

\item{h}{the bandwidth}
}
\value{
\item{x}{grid points}
        \item{y}{estimated values of density}
}
\description{
This function provide the estimated Kernel density values by using Gamma Kernel.The Gamma kernel is developed by Chen (2000). He was first to introduce asymetrical kernels to control boundary Bias.
Gamma Kernel is
\deqn{K_{Gam1( \frac{x}{h}+1, h)}(y) = \frac{y^ \frac{x}{h} exp(-\frac{y}{h})}{ \Gamma (\frac{x}{h}+1)h^{ \frac{x}{h}+1}}}
}
\details{
see the details in the \code{\link{BS}}.
}
\examples{
##Number of grid points "k" should be at least equal to the data size.
###If user defines the generating scheme of grid points then length
####of grid points should be equal or greater than "k". Otherwise NA will be produced.
y <- rexp(100, 1)
xx <- seq(min(y) + 0.05, max(y), length = 500)
h <- 2
den <- Gamma(x = xx, y = y, k = 200, h = h)

##If scheme for generating grid points is unknown
y <- rexp(200, 1)
h <- 3
Gamma(y = y, k = 90, h = h)

\dontrun{
y <- data(TUNA)
xx <- seq(min(y) + 0.05, max(y), length = 500)
h <- 2
den <- Gamma(x = xx, y = y, k = 200, h = h)
}

\dontrun{
##If user do not mention the number of grid points
y <- rexp(1000, 1)
xx <- seq(0.001, 1000, length = 1000)

#any bandwidth can be used
require(KernSmooth)
h <- dpik(y)
Gamma(x = xx, y = y, h = h)
}

\dontrun{
#if generating scheme and number of grid points are missing then function generate NA
y <- rexp(1000, 1)
band = 3
Gamma(y = y, h = band)
}

#if bandwidth is missing
y <- rexp(100,1)
xx <- seq(0.001, max(y), length = 100)
Gamma(x = xx, y = y, k = 90)

}
\references{
Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
}
\seealso{
For further kernels see \code{\link{Erlang}}, \code{\link{BS}}, \code{\link{Beta}}and \code{\link{LogN}}. To plot its density see \code{\link{plot.Gamma}} and to calculate MSE \code{\link{mse}}.
}
\author{
Javaria Ahmad Khan, Atif Akbar.
}

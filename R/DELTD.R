#' DELTD
#'
#' Kernel Density Estimation using Lifetime Distributions
#'
#' @description  A collection of asymmetrical kernels belong to lifetime distributions for kernel density estimation is presented. i.e. \code{\link{plot.BS}}, \code{\link{plot.Beta}}, \code{\link{plot.Erlang}},
#' \code{\link{plot.Gamma}} and \code{\link{plot.LogN}}. Estimated values can also observed by using
#' \code{\link{Beta}}, \code{\link{BS}}, \code{\link{Erlang}}, \code{\link{Gamma}} and \code{\link{LogN}}. For calculating mean squared error by using different kernel functions are \code{\link{mse}} can be used.
#'@author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#'\item Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#'\item Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation using Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#'\item Chen, S. X. 2000. Probability density function estimation using Gamma kernels. \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#'\item Chen, S. X. 2000. Beta kernel smothers for regression curves. \emph{Statistica Sinica} \strong{10}, 73-91.
#'}
"_PACKAGE"

#' Estimate Density Values by Gamma kernel
#'
#' This function provide the estimated Kernel density values by using Gamma Kernel.The Gamma kernel is developed by Chen (2000). He was first to introduce asymetrical kernels to control boundary Bias.
#' Gamma Kernel is
#' \deqn{K_{Gam1( \frac{x}{h}+1, h)}(y) = \frac{y^ \frac{x}{h} exp(-\frac{y}{h})}{ \Gamma (\frac{x}{h}+1)h^{ \frac{x}{h}+1}}}

#' @details see the details in the \code{\link{BS}}.
#' @param x scheme for generating grid points
#' @param y a vector of positive values
#' @param k number of grid points
#' @param h the bandwidth
#' @import stats
#' @examples
#' ##Number of grid points "k" should be at least equal to the data size.
#' ###If user defines the generating scheme of grid points then length
#' ####of grid points should be equal or greater than "k". Otherwise NA will be produced.
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 500)
#' h <- 2
#' den <- Gamma(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating grid points is unknown
#' y <- rexp(200, 1)
#' h <- 3
#' Gamma(y = y, k = 90, h = h)
#'
#' \dontrun{
#' ##If user do not mention the number of grid points
#' y <- rexp(1000, 1)
#' xx <- seq(0.001, 1000, length = 1000)
#'
#' #any bandwidth can be used
#' require(KernSmooth)
#' h <- dpik(y)
#' Gamma(x = xx, y = y, h = h)
#' }
#'
#' \dontrun{
#' #if generating scheme and number of grid points are missing then function generate NA
#' y <- rexp(1000, 1)
#' band = 3
#' Gamma(y = y, h = band)
#' }
#'
#' #if bandwidth is missing
#' y <- rexp(100,1)
#' xx <- seq(0.001, max(y), length = 100)
#' Gamma(x = xx, y = y, k = 90)
#'
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references \itemize{
#'\item Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' \item Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
#' }
#' @seealso For further kernels see  \code{\link{Beta}}, \code{\link{BS}}, \code{\link{Erlang}} and \code{\link{LogN}}. To plot its density see \code{\link{plot.Gamma}} and to calculate MSE \code{\link{mse}}.
#' @export
Gamma <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length = k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1/5)
  Kgamma <- matrix(rep(0, k * n), ncol = k)
  ###########gamma###########
  for(j in 1:k) {
    for(i in 1:n) {
      fn <- gamma((x[j] / h) + 1)
      Kgamma[i, j] <- ((y[i] ^ (x[j] / h)) * (exp( - (y[i] / h))))/(h ^ ((x[j] / h) + 1) * fn)
    }}
  fhat <- colMeans(Kgamma)
  results <- list(x = x,
                  y = fhat)
  class ( results ) <- c('list',
                         'Gamma')
  results
}
#' Estimate Density Values by Beta kernel
#'
#' This function provide the estimated Kernel density values by using Beta Kernel. The Beta kernel is developed by Chen (2000) by using Beta distribution of first kind. He was first to introduce asymetrical kernels to control boundary Bias.
#' Beta Kernel is
#' \deqn{K_{Beta( \frac{x}{h}+1, \frac{1-x}{h}+1)}(y)=\frac{y^ \frac{x}{h} (1-y)^{\frac{1-x}{b}}} { B \{ \frac{x}{h}+1, \frac{(1-x)}{h}+1 \}}}
#' @details In this function, choice of bandwidth, number of grid points and scheme that how these grid points are generated are user based. If any parameter(s) is missing then function used default parameters.
#' But at least \code{x} or \code{k} should be specified otherwise \code{NA} will be produced. If \code{x} is missing then function will generate \code{k} grid points by using uniform distribution. Similarly, if
#' \code{k} is missing then function consider it same to length of main vector \code{y}. In case if \code{h} is missing then function used normal scale rule bandwidth for non-normal data and described in Silverman (1986). This function can be onlt used if
#' data is between (0, 1). Similarly, \code{x} should be also lies between (0, 1).
#' @param x scheme for generating grid points
#' @param y a vector of positive values
#' @param k number of grid points
#' @param h the bandwidth
#' @import stats
#' @examples
#' ## Data: Simulated or real data can be used
#' ## Number of grid points "k" should be at least equal to the data size.
#' ## If user defines the generating scheme of grid points then length
#' ## of grid points should be equal or greater than "k", Otherwise NA will be produced.
#' y <- runif(50)
#' xx <- sample(0.00001:900, 500, replace = FALSE)/1000
#' h <- 0.9
#' Beta(x = xx, y = y, k = 500, h = h)
#'
#' ## If scheme for generating grid points is unknown
#' y <- runif(500)
#' h <- 0.9
#' Beta(x = xx, y = y, k = 500, h = h)
#'
#' \dontrun{
#' ## If user do not mention the number of grid points
#' y <- runif(1000)
#' xx <- seq(0.001, 1000, length = 2000)
#'
#' ## any bandwidth can be used
#' require(kedd)
#' h <- h.bcv(y) ## Biased cross validation
#' Beta(x = xx, y = y, h = h)
#' }
#'
#' \dontrun{
#' ##if both generating scheme and number of grid points are missing then function generate NA
#' y <- runif(1000)
#' band = 0.8
#' Beta(y = y, h = band)
#' }
#'
#' ## if bandwidth is missing
#' y <- runif(100)
#' xx <- seq(0.001, 100, length = 300)
#' Beta(x = xx, y = y, k = 200)
#'
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Beta kernel smothers for regression curves. \emph{Statistica Sinica} \strong{10}, 73-91.
#' Silverman, B. W. 1986. \emph{Density Estimation}. Chapman & Hall/ CRC, London.
#' @seealso For further kernels see \code{\link{BS}}, \code{\link{Gamma}}, \code{\link{Erlang}} and \code{\link{LogN}}. To plot its density see \code{\link{plot.Beta}} and to calculate MSE \code{\link{mse}}.
#' @export
Beta <- function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x <- runif(k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h <- 0.79 * IQR(y) * length(y) ^ (-1/5)

  Kbeta <- matrix(rep(0, k * n), ncol = k)
  ###########Beta###########
  for(j in 1:k) {
    for(i in 1:n) {
      a <- (x[j] / h)
      b <- ((1 - x[j]) / h)
      Kbeta[i, j] <- ((y[i] ^ a) * (1 - y[i]) ^ b) / (beta(a + 1, b + 1))
    }
  }
  fhat <- colMeans(Kbeta)
  results <- list(x = x,
                  y = fhat)

  class ( results ) <- c('list',
                         'Beta')
  results
}

#' Density Plot by Beta kernel
#'
#' Plot density by using Beta Kernel.
#' @param x an object of class "Beta"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- runif(100)
#' h <- 0.5
#' xx <- sample(0.00001:900, 50, replace = FALSE)/1000
#' den <- Beta(x = xx, y = y, k = 50, h = h)
#' plot(den, type = "p")

#' ##other details can also be added
#' y <- runif(100)
#' h <- 0.7
#' xx <- sample(0.00001:900, 50, replace = FALSE)/1000
#' den <- Beta(x = xx, y = y, k = 50, h = h)
#' plot(den, type = "l", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Beta kernel smothers for regression curves. \emph{Statistica Sinica} \strong{10}, 73-91.
#' @seealso For further kernels see \code{\link{plot.BS}}, \code{\link{plot.Erlang}}, \code{\link{plot.Gamma}} and \code{\link{plot.LogN}}. To calculate its estimated values see \code{\link{Beta}} and for
#' MSE see \code{\link{mse}}.
#' @export
plot.Beta <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Density Plot by Gamma kernel
#'
#' Plot density by using Gamma Kernel.
#' @param x an object of class "Gamma"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100, 1)
#' h <- 1.5
#' xx <- seq(min(y) + 0.05, max(y), length =200)
#' den <- Gamma(x=xx, y=y, k=200, h=h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- rexp(100, 2)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' gr <- Gamma(x=xx, y=y, k=200, h=h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw=h)
#' lines(d1, type="p", col="red")
#' legend("topright", c("Real Density", "Density by Gamma Kernel"),
#' col=c("red", "black"), lty=c(1,2))
#'
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' @seealso For further kernels see \code{\link{plot.Beta}},  \code{\link{plot.BS}}, \code{\link{plot.Erlang}} and \code{\link{plot.LogN}}. To calculate its estimated values see \code{\link{Gamma}} and for
#' MSE \code{\link{mse}}.
#' @export
plot.Gamma <- function(x,...) {
  plot(x$x, x$y,...)
}
#' Calculate Mean Squared Error(MSE) by using different Kernels
#'
#' This function calculates the mean squared error (MSE) by using user specified kernel. But distribution of vector should be Exponential, Gamma or Weibull. Any other choice of distribution will result \code{NaN}.
#' @param kernel type of kernel which is to be used
#' @param type mention distribution of vector. If exponential distribution then use \code{"Exp"}.
#'     If use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references  \itemize{
#'\item Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103-124.
#'\item Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation using Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141-152.
#'\item Chen, S. X. 2000. Probability density function estimation using Gamma kernels. \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#'\item Chen, S. X. 2000. Beta kernel smothers for regression curves. \emph{Statistica Sinica} \strong{10}, 73-91.
#'}
#' @examples
#' y <- rexp(100, 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 500)
#' h <- 2
#' gr <- Gamma(x = xx, y = y, k = 200, h = h)
#' mse(kernel = gr, type = "Exp")
#' ## if distribution is other than mentioned \code{type} is used then NaN will be produced.
#' \dontrun{
#' mse(kernel = gr, type ="Beta")
#' }
#' @return Mean Squared Error (MSE)
#' @export
mse<-function(kernel,type){
  ftrue<-switch(type,
                Exp = dexp(kernel$x, (1 / mean(kernel$x))),
                Gamma = dgamma(kernel$x, (mean(kernel$x) / (var(kernel$x) / mean(kernel$x))), (var(kernel$x) / mean(kernel$x))),
                Weibull = dweibull(kernel$x, ((sd(kernel$x) / mean(kernel$x)) ^ - 1.806), scale = 1, log = FALSE)
  )
  return(mean((ftrue-kernel$y)^2))#mse of fhat w.r.t. the true density
}

#' Estimate Density Values by Lognormal kernel
#'
#' The \code{LogN} estimate Values of density by using Lognormal Kernel.The Lognomal kernel is developed by Jin and Kawczak (2003). For this too, they claimed that performance of their developed kernel is better near the
#' boundary points in terms of boundary reduction.
#' Lognormal Kernel is
#' \deqn{K_{LN(\ln(x),4\ln(1+h))}=\frac{1}{\sqrt{( 8\pi \ln(1+h))} y)} exp\left[-\frac{(\ln(y)-\ln(x))^2}{(8\ln(1+h))}\right]}
#' @details see the details in the \code{\link{BS}}.
#' @param x scheme for generating grid points
#' @param y a vector of positive values.
#' @param k grid points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' ## Data: Simulated or real data can be used
#' ## Number of grid points "k" should be at least equal to the data size.
#' ## If user defines the generating scheme of grid points then length
#' ## of grid points should be equal or greater than "k", Otherwise NA will be produced.
#' y <- rweibull(350, 1)
#' xx <- seq(0.001, max(y), length = 500)
#' h <- 2
#' den <- LogN(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating grid points is unknown
#' n <- 1000
#' y <- abs(rlogis(n, location = 0, scale = 1))
#' h <- 3
#' LogN(y = y, k = 90, h = h)
#'
#' \dontrun{
#' ##If user do not mention the number of grid points
#' y <- rweibull(350, 1)
#' xx <- seq(0.00001, max(y), 500)
#'
#' #any bandwidth can be used
#' require(ks)
#' h <- hscv(y)   #Smooth cross validation bandwidth
#' LogN(x = xx, y = y, h = h)
#' }
#'
#' \dontrun{
#' #if both scheme and number of grid points are missing then function generate NA
#' n <- 1000
#' y <- abs(rlogis(n, location = 0, scale = 1))
#' band = 3
#' LogN(y = y, h = band)
#' }
#'
#' #if bandwidth is missing
#' y <- rweibull(350, 1)
#' xx <- seq(0.001, 100, length = 500)
#' LogN(x = xx, y = y, k = 90)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103-124.
#' @seealso For further kernels see \code{\link{Beta}}, \code{\link{BS}}, \code{\link{Erlang}} and \code{\link{Gamma}}. To plot its density see \code{\link{plot.LogN}} and to calculate MSE use \code{\link{mse}}.
#' @export
LogN<-function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length =k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1/5)
  KLNormal <- matrix(rep(0, k * n), ncol = k)
  ###########Lognormal###########
  for(j in 1:k) {
    for(i in 1:n) {

      KLNormal [i, j] <-(1 / (y[i] * sqrt(8 * pi * log(1 + h)))) * exp( - (((log(y[i]) - log(x[j])) ^ 2) / (8 * log(1 + h))))
    }}
  fhat<- colMeans(KLNormal)
  results <- list(x = x,
                  y = fhat)
  class ( results ) <- c('list',
                         'LogN')
  results
}

#' Density Plot by Lognormal kernel
#'
#' Plot Kernel density by using Lognormal Kernel.
#' @param x An object of class "LogN"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' n <- 1000
#' y <- abs(rlogis(n, location = 0, scale = 1))
#' xx <- seq(min(y) + 0.05, max(y), length =90)
#' h <- 0.00003
#' den <- LogN(x = xx, y = y, k = 90, h = h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- abs(rlogis(n, location = 0, scale = 1))
#' h <- 3
#' gr <- LogN(x = xx, y = y, k = 90, h = h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "green")
#' legend("topleft", c("Real Density", "Density by Lognormal Kernel"),
#' col = c("green", "black"), lty = c(1,2))
#' @return Nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103-124.
#' @seealso For further kernels see \code{\link{plot.Beta}}, \code{\link{plot.BS}}, \code{\link{plot.Erlang}} and \code{\link{plot.Gamma}}. To calculate MSE use \code{\link{mse}} and for estimated values for density
#' estimation see \code{\link{LogN}}.
#' @export
plot.LogN <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Estimate Density Values by Erlang kernel
#'
#' This function provide the estimated values for density by using Erlang Kernel. Erlang kernel is developed by Salha et al. (2014). They developed this asymmetrical kernal with its hazard function and also
#' proved its asymtotic normality.
#' \deqn{K_{E(x,\frac{1}{h})}  (y)=\frac{1}{\Gamma (1+\frac{1}{h})} \left[\frac{1}{x} (1+\frac{1}{h}) \right]^\frac{h+1}{h} y^\frac{1}{h} exp\left(-\frac{y}{x} (1+\frac{1}{h}) \right)}
#' @details see the details in the \code{\link{BS}}.
#' @param x scheme for generating grid points
#' @param y a vector of positive values.
#' @param k grid points.
#' @param h the bandwidth
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation using Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141-152.
#' @seealso For further MSE by using other kernels see \code{\link{Beta}}, \code{\link{BS}}, \code{\link{Gamma}} and \code{\link{LogN}}. For plotting these estimated values \code{\link{plot.Erlang}} and for calculating MSE use \code{\link{mse}}.
#' @examples
#' ## Data: Simulated or real data can be used
#' ## Number of grid points "k" should be at least equal to the data size.
#' ## If user defines the generating scheme of grid points then length
#' ## of grid points should be equal or greater than "k", Otherwise NA will be produced.
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' xx <- seq(min(y) + 0.05, max(y), length = 500)
#' h <-2
#' den <- Erlang(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating grid points is unknown
#' y <- rlnorm(1000, meanlog = 0, sdlog = 1)
#' h <- 3
#' Erlang(y = y, k = 90, h = h)
#'
#' \dontrun{
#' ##If user do not mention the number of grid points
#'  y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' xx <- seq(0.001, 1000, length = 1000)
#'
#' #any bandwidth can be used
#' require(kedd)
#' h <- h.ucv(y)     #Unbaised cross validation bandwidth
#' Erlang(x = xx, y = y, h = h)
#' }
#'
#' \dontrun{
#' #if generating scheme and number of grid points are missing then function generate NA
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' band = 3
#' Erlang(y = y, h = band)
#' }
#'
#' #if bandwidth is missing
#'  y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' xx <- seq(0.001, 100, length = 100)
#' Erlang(x = xx, y = y, k = 90)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @export
Erlang<-function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length =k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1/5)
  KErlang <- matrix(rep(0, k * n), ncol = k)
  ###########erlang###########
  for(j in 1:k) {
    for(i in 1:n) {
      KErlang[i, j] <- (1 / (gamma(1 + (1 / h)))) * ((1 / x[j]) * (1 + (1 / h))) ^ ((h + 1) / h) * y[i] ^ (1 / h) * exp( - y[i] / x[j] * (1 + (1 / h)))
    }
  }
  fhat<- colMeans(KErlang)
  results <- list(x = x,
                  y = fhat)
  class ( results ) <-c('list',
                        'Erlang')
  results
}

#' Density Plot by Erlang kernel
#'
#' Plot Kernel density by using Erlang Kernel.
#' @param x An object of class "Erlang"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation using Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141-152.
#' @seealso For further MSE by using other kernels see \code{\link{plot.Beta}}, \code{\link{plot.BS}}, \code{\link{plot.Gamma}} and \code{\link{plot.LogN}}. For estimated values \code{\link{Erlang}} and for calculating MSE see \code{\link{mse}}.
#' @examples
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' h <- 1.5
#' xx <- seq(min(y) + 0.05, max(y), length = 200)
#' den <- Erlang(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- rlnorm(100, meanlog = 0, sdlog = 1)
#' grid <- seq(min(y) + 0.05, max(y), length = 200)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' gr <- Erlang(x = grid, y = y, k = 200, h = h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "red")
#' legend("topright", c("Real Density", "Density by Erlang Kernel"),
#' col=c("red", "black"), lty=c(1,2))
#' @return Nothing
#' @export
plot.Erlang <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Estimate Density Values by Birnbaum-Saunders kernel
#'
#' This function calculates the estimated Values by using Birnbaum-Saunders Kernel. The Birnbaum-Saunders kernel is developed by Jin and Kawczak (2003). They claimed that performance of their developed kernel is better near the
#' boundary points in terms of boundary reduction.
#' \deqn{K_{BS(h^\frac{1}{2},x)} (y)=\frac{1}{2\sqrt 2 \pi h} \left(\sqrt \frac{1}{xy} +\sqrt\frac{x}{y^3}\right)exp\left(-\frac{1}{2h}\left(\frac{y}{x}-2+\frac{x}{y}\right)\right)}
#'
#' @details In this function, choice of bandwidth, number of grid points and scheme that how these grid points are generated are user based. If any parameter(s) is missing then function used default parameters.
#' But at least \code{x} or \code{k} should be specified otherwise \code{NA} will be produced. If \code{x} is missing then function will generate \code{k} grid points between minimum and maximum values of vector. Similarly, if
#' \code{k} is missing then function consider it same to length of main vector \code{y}. In case if \code{h} is missing then function used normal scale rule bandwidth for non-normal data and described in Silverman (1986).

#' @param x scheme for generating grid points
#' @param y a vector of positive values.
#' @param k grid points
#' @param h the bandwidth
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103-124.
#' @seealso For further kernels see \code{\link{Beta}}, \code{\link{Erlang}}, \code{\link{Gamma}} and \code{\link{LogN}}. To plot the density by using BS kernel \code{\link{plot.BS}} and to calculate MSE by \code{\link{mse}}.
#'
#' @examples
#' ## Data: Simulated or real data can be used
#' ## Number of grid points "k" should be at least equal to the data size.
#' ## If user defines the generating scheme of grid points then length
#' ## of grid points should be equal or greater than "k", Otherwise NA will be produced.
#' alpha = 10
#' theta = 15 / 60
#'
#'y <- rgamma(n = 1000, shape = alpha, scale = theta)
#'xx <- seq(min(y) + 0.05, max(y), length =200)
#'h <- 1.1
#'den <- BS(x = xx, y = y, k = 200, h = h)
#'
#' ##If scheme for generating grid points is unknown
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' h <- 3
#' BS(y = y, k = 90, h = h)
#'
#' \dontrun{
#' ##If user do not mention the number of grid points
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' xx <- seq(0.001, 1000, length = 1000)
#'
#' #any bandwidth can be used
#' require(KernSmooth)
#' h <- dpik(y)     #Direct Plug-In Bandwidth
#' BS(x = xx, y = y, h = h)
#' }
#'
#' \dontrun{
#' #if both generating scheme and number of grid points are missing then function generate NA
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' band = 3
#' BS(y = y, h = band)
#' }
#'
#' #if bandwidth is missing
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' xx <- seq(0.001, 100, length = 1000)
#' BS(x = xx, y = y, k = 900)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @export

BS<-function(x = NULL, y, k = NULL, h = NULL){
  n <- length(y)
  if(is.null(x))
    x = seq(min(y) + 0.05, max(y), length = k)
  if(is.null(k))
    k = length(y)
  if(is.null(h))
    h = 0.79 * IQR(y) * length(y) ^ (-1/5)

  Kbs <- matrix(rep(0, k * n), ncol = k)
  ###########BS###########
  for(j in 1:k) {
    for(i in 1:n) {
      Kbs[i, j] <- (1 / (2 * sqrt(2 * h * pi))) * ((sqrt(1 / (x[j] * y[i]))) + (sqrt(x[j] / (y[i] ^ 3)))) * exp( - (y[i] / (2 * h * x[j])) + (1 / h) - (x[j] / (2 * h * y[i])))
    }}
  fhat<- colMeans(Kbs)
  results <- list(x = x,
                  y = fhat)

  class ( results ) <- c('list',
                         'BS')
  results
}

#' Density Plot by Birnbaum-Saunders kernel
#'
#' Plot Kernel density by using Birnbaum-Saunders Kernel.
#' @param x An object of class "BS"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103-124.
#' @seealso For further kernels see  \code{\link{plot.Beta}}, \code{\link{plot.Erlang}}, \code{\link{plot.Gamma}} and \code{\link{plot.LogN}}. For estimated values \code{\link{BS}} and for MSE \code{\link{mse}}.
#' @examples
#' alpha = 10
#' theta = 15 / 60
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' h <- 1.5
#' xx <- seq(min(y) + 0.05, max(y), length = 200)
#' den <- BS(x = xx, y = y, k = 200, h = h)
#' plot(den, type = "l")
#'
#' ##other details can also be added
#' y <- rgamma(n = 1000, shape = alpha, scale = theta)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)  #Normal Scale Rule Bandwidth
#' gr <- BS(x = xx, y = y, k = 200, h = h)
#' plot(gr, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#'
#' ## To add true density along with estimated
#' d1 <- density(y, bw = h)
#' lines(d1, type = "p", col = "red")
#' legend("topright", c("Real Density", "Density by Birnbaum-Saunders Kernel"),
#' col=c("red", "black"), lty = c(1,2))
#' @return Nothing
#' @export
plot.BS <- function(x,...) {
  plot(x$x, x$y,...)
}

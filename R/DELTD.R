#' DELTD
#'
#' Kernel Density Estimation using Lifetime Distributions
#'
#' @description  A collection of asymmetrical kernels belong to lifetime distributions for kernel density estimation is presented. i.e. \code{\link{plot.BS}}, \code{\link{plot.Erlang}},
#' \code{\link{plot.Gamma}}, \code{\link{plot.LN}} and for plotting all densities at same time \code{\link{dencomb}}. Estimated values can also observed by using
#' \code{\link{BS}}, \code{\link{Gamma}}, \code{\link{Erlang}} and \code{\link{LN}}. For calculating mean squared error by using different kernels functions are \code{\link{mseBS}}, \code{\link{mseEr}},
#' \code{\link{mseGamma}} and \code{\link{mseLN}}.
#'@author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#'\item Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#'\item Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation ksing Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#'\item Chen, S. X. 2000. Probability density function estimation using Gamma kernels. \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#'}
"_PACKAGE"


#' Estimated Density Values by Gamma kernel
#'
#' Estimated Kernel density values by using Gamma Kernel.
#' @details The Gamma kernel is developed by Chen (2000). He was first to introduce asymetrical kernels to control boundary Bias.
#' Gamma Kernel is
#' \deqn{K_{Gam1( \frac{x}{h+1}, h)}(y) = \frac{y^ \frac{x}{h} exp(-\frac{y}{h})}{ \Gamma \frac{x}{(h+1)}h^{ \frac{x}{h+1}}}}

#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' Gamma(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' @seealso For further kernels see \code{\link{Erlang}}, \code{\link{BS}} and \code{\link{LN}}. To plot its density see \code{\link{plot.Gamma}} and to calculate MSE by using Gamma Kernel \code{\link{mseGamma}}.
#' @export
Gamma<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  Kgamma <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########gamma###########
  for(j in 1:k) {
    for(i in 1:n) {
      fn<-gamma((x[j]/h)+1)
      Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
    }
    fhat[j] <- 1/n * (sum(Kgamma[, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'Gamma')
  results
}

#' Density Plot by Gamma kernel
#'
#' Plot density by using Gamma Kernel.
#' @param x an object of class "BS"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- Gamma(y,200,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' @return nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' @seealso For further kernels see \code{\link{plot.Erlang}}, \code{\link{plot.BS}} and \code{\link{plot.LN}}. To calculate its estimated values see \code{\link{Gamma}} and for
#' MSE by using Gamma Kernel \code{\link{mseGamma}}.
#' @export
plot.Gamma <- function(x,...) {
  plot(x$x, x$y,...)
}
#' Calculate Mean Squared Error( MSE) when Gamma kernel is used
#'
#' Calculate MSE by using Gamma Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     If use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#' @seealso For further MSE by using other kernels see \code{\link{mseBS}}, \code{\link{mseEr}} and \code{\link{mseLN}}. For density estimation by using Gamma Kernel \code{\link{plot.Gamma}} and for estimated values
#' of density \code{\link{Gamma}}.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' mseGamma(y,200,h,"Exp")
#' @return MSE
#' @export
mseGamma<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,

                Exp = dexp(x,(1/mean(x))),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

  )
  Kgamma <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########gamma###########
  for(j in 1:k) {
    for(i in 1:n) {
      fn<-gamma((x[j]/h)+1)
      Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
    }
    fhat[j] <- 1/n * (sum(Kgamma[, j]))
  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end

#' Estimated Density Values by Lognormal kernel
#'
#' Estimated Values of Density estimation by using Lognormal Kernel.
#' @details The Lognomal kernel is also developed by Jin and Kawczak (2003). For this too, they claimed that performance of their developed kernel is better near the
#' boundary points in terms of boundary reduction.
#' Lognormal Kernel is
#' \deqn{K_{LN(\ln(x),4\ln(1+h))}=\frac{1}{\sqrt{( 8\pi \ln(1+h))} y)} exp\left[-\frac{(\ln(y)-\ln(x))^2}{(8\ln(1+h))}\right]}
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' LN(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further kernels see \code{\link{Erlang}}, \code{\link{Gamma}} and \code{\link{BS}}. To plot its density see \code{\link{plot.LN}} and to calculate MSE by using Lognormal Kernel \code{\link{mseLN}}.
#' @export
LN<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  KLNormal <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########Lognormal###########
  for(j in 1:k) {
    for(i in 1:n) {

      KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
    }
    fhat[j] <- 1/n * (sum(KLNormal [, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'LN')
  results
}

#' Density Plot by Lognormal kernel
#'
#' Plot Kernel density by using Lognormal Kernel.
#' @param x An object of class "LN"
#' @param \dots Not presently used in this implementation
#' @import graphics
#' @import stats
#' @examples
#' y <- rexp(23,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- LN(y,90,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' ## To add true density along with estimated
#' d1 <- density(y,bw=h)
#' lines(d1,type="p",col="red")
#' @return Nothing
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further kernels see \code{\link{plot.Erlang}}, \code{\link{plot.Gamma}} and \code{\link{plot.BS}}. To calculate MSE by using Lognormal Kernel \code{\link{mseLN}} and for estimated values for density
#' estimation see \code{\link{LN}}.
#' @export
plot.LN <- function(x,...) {
  plot(x$x, x$y,...)
}
#' Calculate Mean Squared Error( MSE) when Lognormal Kernel is used
#'
#' Calculate MSE by using Lognormal Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     if use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further MSE by using other kernels see \code{\link{mseBS}}, \code{\link{mseEr}} and \code{\link{mseGamma}}. For estimated values and for density estimation by using Lognormal Kernel see \code{\link{LN}} and \code{\link{plot.LN}}, respectively.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' mseLN(y,200,h,"Exp")
#' @return MSE
#' @export
mseLN<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,

                Exp = dexp(x,(1/mean(x))),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

  )
  KLNormal <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########Lognormal###########
  for(j in 1:k) {
    for(i in 1:n) {

      KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
    }
    fhat[j] <- 1/n * (sum(KLNormal [, j]))
  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end

#' Estimated Density Values by Erlang kernel
#'
#' Estimated values for density by using Erlang Kernel.
#' @details Erlang kernel is developed by Salha et al. (2014). They developed this asymmetrical kernal with its hazard function and also
#' proved its asymtotic normality.
#' \deqn{K_{E(x,\frac{1}{h})}  (y)=\frac{1}{\Gamma (1+\frac{1}{h})} \left[\frac{1}{x} (1+\frac{1}{h}) \right]^\frac{h+1}{h} y^\frac{1}{h} exp\left(-\frac{y}{x} (1+\frac{1}{h}) \right)}
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation ksing Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#' @seealso For further MSE by using other kernels see \code{\link{BS}}, \code{\link{Gamma}} and \code{\link{LN}}. For plotting these estimated values \code{\link{plot.Erlang}} and for calculating MSE by using Erlang Kernel \code{\link{mseEr}}.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' Erlang(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @export
Erlang<-function(y,k,h){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)

  KErlang <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########erlang###########
  for(j in 1:k) {
    for(i in 1:n) {
      KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
    }
    fhat[j] <- 1/n * (sum( KErlang[, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'Erlang')
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
#' @references Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation ksing Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#' @seealso For further MSE by using other kernels see \code{\link{plot.BS}}, \code{\link{plot.Gamma}} and \code{\link{plot.LN}}. For estimated values \code{\link{Erlang}} and for calculating MSE by using Erlang Kernel \code{\link{mseEr}}.
#' @examples
#' y <- rexp(23,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' ans <- Erlang(y,90,h)
#' plot(ans, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' ## To add true density along with estimated
#' d1<-density(y,bw=h)
#' lines(d1,type="p",col="red")
#' legend("topright", c("Real Density", "Density by Erlang Kernel"), col=c("red", "black"), lty=c(1,2))
#' @return Nothing
#' @export
plot.Erlang <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Calculate Mean Squared Error( MSE) when Erlang kernel is used
#'
#' Calculate MSE by using Erlang Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     if use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation ksing Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#' @seealso For further MSE by using other kernels see \code{\link{mseBS}}, \code{\link{mseLN}} and \code{\link{mseGamma}}. For estimated values for density estimation \code{\link{Erlang}} and for density estimation by using Erlang Kernel \code{\link{plot.Erlang}}.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' mseEr(y,200,h,"Exp")
#' @return MSE
#' @export
mseEr<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,

                Exp = dexp(x,(1/mean(x))),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

  )
  KErlang <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########erlang###########
  for(j in 1:k) {
    for(i in 1:n) {
      KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
    }
    fhat[j] <- 1/n * (sum( KErlang[, j]))
  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end

#' Estimated Density Values by Birnbaum-Saunders kernel
#'
#' Estimated Values by using Birnbaum-Saunders Kernel.
#' @details  The Birnbaum-Saunders kernel is developed by Jin and Kawczak (2003). They claimed that performance of their developed kernel is better near the
#' boundary points in terms of boundary reduction.
#' \deqn{K_{BS(h^\frac{1}{2},x)} (y)=\frac{1}{2\sqrt 2 \pi h} \left(\sqrt \frac{1}{xy} +\sqrt\frac{x}{y^3}\right)exp\left(-\frac{1}{2h}\left(\frac{y}{x}-2+\frac{x}{y}\right)\right)}
#' @param y a numeric vector of positive values.
#' @param k gird points
#' @param h the bandwidth
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further kernels see \code{\link{Erlang}}, \code{\link{Gamma}} and \code{\link{LN}}. To plot the density by using BS kernel \code{\link{plot.BS}} and to calculate MSE by using Birnbaum-Saunders Kernel \code{\link{mseBS}}.
#' @examples
#'  y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' BS(y,200,h)
#' @return \item{x}{grid points}
#'         \item{y}{estimated values of density}
#' @export

BS<-function(y,k,h){
  n<-length(y)
  x <- seq(min(y) + 0.05, max(y), length =k)
  Kbs <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########BS###########
  for(j in 1:k) {
    for(i in 1:n) {
      Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
    }
    fhat[j] <- 1/n * (sum(Kbs[, j]))
  }
  results <- list(x=x, y=fhat)
  class ( results ) <-c('list', 'BS')
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
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further kernels see \code{\link{plot.Erlang}}, \code{\link{plot.Gamma}} and \code{\link{plot.LN}}. For estimated values \code{\link{BS}} and for MSE \code{\link{mseBS}}.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' den <- BS(y,200,h)
#' plot(den, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
#' ## To add true density along with estimated
#' d1<-density(y,bw=h)
#' lines(d1,type="p",col="green")
#' @return Nothing
#' @export
plot.BS <- function(x,...) {
  plot(x$x, x$y,...)
}

#' Calculate Mean Squared Error( MSE) when Birnbaum-Saunders kernel is used
#'
#' Calculate MSE by using Birnbaum-Saunders Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use \code{"Exp"}.
#'     if use gamma distribution then use \code{"Gamma"}.If Weibull distribution then use \code{"Weibull"}.
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#' @references Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#' @seealso For further MSE by using other kernels see \code{\link{mseLN}}, \code{\link{mseEr}} and \code{\link{mseGamma}}. For density estimation by using Birnbaum-Saunders Kernel \code{\link{plot.BS}} and to examine estimated values with grid points see \code{\link{BS}}.
#' @examples
#' y <- rexp(100,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' mseBS(y,200,h,"Exp")
#' @return MSE
#' @export
mseBS<-function(y,k,h,type){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  ftrue<-switch(type,

                Exp = dexp(x,(1/mean(x))),
                Weibull = dweibull(x, ((sd(x)/mean(x))^-1.806), scale = 1, log = FALSE),
                Gamma = dgamma(x,(mean(x)/(var(x)/mean(x))),(var(x)/mean(x)))

  )
  Kbs <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########BS###########
  for(j in 1:k) {
    for(i in 1:n) {
      Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
      #Kbs[is.nan(Kbs)] <- 0
    }
    fhat[j] <- 1/n * (sum(Kbs[, j]))
  }#2nd loop end
  return(mean((ftrue-fhat)^2))#mse of fhat w.r.t. the true density
}#function end

#' Plot the Densities for Comparison
#'
#' Plot the Estimated Densities with Real Density for Comparison .
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param comb mention the combination which kernel estimated densities are to be compared. If Lognormal and Birnbaum-Saunders kernel densities
#'           are to be compared along with real density then use \code{"TLB"}. If Lognormal and Erlang then use \code{"TLE"}. If Lognormal and Gamma to be compared
#'           then use \code{"TLG"}. For Birnbaum-Saunders and Erlang use \code{"TBE"}. For Birnbaum-Saunders and Gamma then use \code{"TBG"}. For Erlang and Gamma use \code{"TEG"}.
#'           For Lognormal, Birnbaum-Saunders and Erlang use \code{"TLBE"}. For Lognormal, Birnbaum-Saunders and Gamma use \code{"TLBG"}. For Birnbaum-Saunders, Erlang and Gamma
#'           use \code{"TBEG"}. To compare all densities in one graph use \code{"TLBEG"}.
#' @import graphics
#' @import stats
#' @author Javaria Ahmad Khan, Atif Akbar.
#'@references \itemize{
#'\item Jin, X.; Kawczak, J. 2003. Birnbaum-Saunders & Lognormal kernel estimators for modeling durations in high frequency financial data. \emph{Annals of Economics and Finance} \strong{4}, 103–124.
#'\item Salha, R. B.; Ahmed, E. S.; Alhoubi, I. M. 2014. Hazard rate function estimation ksing Erlang Kernel. \emph{Pure Mathematical Sciences} \strong{3} (4), 141–152.
#'\item Chen, S. X. 2000. Probability density function estimation using Gamma kernels.  \emph{Annals of the Institute of Statistical Mathematics} \strong{52} (3), 471-480.
#'}
#' @examples
#' \dontrun{
#' y <- rexp(10,1)
#' h <- 0.79 * IQR(y) * length(y) ^ (-1/5)
#' dencomb(y,20,h,"TLB")}
#' @return Plot of Estimated Densities with Real Data Density.
#' @seealso For indivisual densities of each kernels see \code{\link{plot.Erlang}}, \code{\link{plot.BS}}, \code{\link{plot.LN}}, and \code{\link{plot.Gamma}}

#' @export
dencomb<-function(y,k,h,comb){
  n <- length(y)
  x <- seq(min(y) + 0.05, max(y), length = k)
  KLNormal <- matrix(rep(0, k * n), ncol = k)
  Kbs <- matrix(rep(0, k * n), ncol = k)
  KErlang <- matrix(rep(0, k * n), ncol = k)
  Kgamma <- matrix(rep(0, k * n), ncol = k)
  fhat1 <- rep(0, k)
  fhat2<- rep(0, k)
  fhat3 <- rep(0, k)
  fhat4 <- rep(0, k)
  new<-switch(comb,

              TLB = 	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat2[j] <- 1/n * (sum(Kbs[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat2,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Lognormal Kernel","Density by Birnbaum-Saunders Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TLE = 	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat3,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Lognormal Kernel","Density by Erlang Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TLG = 	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Lognormal Kernel","Density by Gamma Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TBE = 	for(j in 1:k) {
                for(i in 1:n) {
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                }
                fhat2[j] <- 1/n * (sum(Kbs [, j]))
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat2, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat3,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Birnbaum-Saunders Kernel","Density by Erlang Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TBG = 	for(j in 1:k) {
                for(i in 1:n) {
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat2[j] <- 1/n * (sum(Kbs [, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat2, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Birnbaum-Saunders Kernel","Density by Gamma Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TEG = 	for(j in 1:k) {
                for(i in 1:n) {
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat3, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Erlang Kernel","Density by Gamma Kernel"),
                       col=c("blue", "black", "red"), lty=c(1,2,6))
              },
              TLBE =	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat2[j] <- 1/n * (sum(Kbs[, j]))
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                par(new=T)
                plot(x,fhat2, type = "s",pch="+",col="green",xaxt="n", yaxt="n", xlab="", ylab="")
                points(x,fhat3,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Log Normal Kernel","Density by Birnbaum-Saunders Kernel", "Density by Erlang Kernel"),
                       col=c("blue", "black", "green","red"), lty=c(1,2,3,6))
              },
              TLBG = 	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat2[j] <- 1/n * (sum(Kbs[, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                par(new=T)
                plot(x,fhat2, type = "s",pch="+",col="green",xaxt="n", yaxt="n", xlab="", ylab="")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Log Normal Kernel","Density by Birnbaum-Saunders Kernel", "Density by Gamma Kernel"),
                       col=c("blue", "black", "green","red"), lty=c(1,2,3,6))
              },
              TBEG = 	for(j in 1:k) {
                for(i in 1:n) {
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat2[j] <- 1/n * (sum(Kbs[, j]))
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat2, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                par(new=T)
                plot(x,fhat3, type = "s",pch="+",col="green",xaxt="n", yaxt="n", xlab="", ylab="")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Birnbaum-Saunders Kernel", "Density by Erlang Kernel","Density by Gamma Kernel"),
                       col=c("blue", "black", "green","red"), lty=c(1,2,3,6))
              },
              TLBEG = 	for(j in 1:k) {
                for(i in 1:n) {
                  KLNormal [i, j] <-(1/(y[i]*sqrt(8*pi*log(1+h))))*exp(-(((log(y[i])-log(x[j]))^2)/(8*log(1+h))))
                  Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
                  KErlang[i, j] <-(1/(gamma(1+(1/h))))*((1/x[j])*(1+(1/h)))^((h+1)/h)*y[i]^(1/h)*exp(-y[i]/x[j]*(1+(1/h)))
                  fn<-gamma((x[j]/h)+1)
                  Kgamma[i, j] <-((y[i]^(x[j]/h))*(exp(-(y[i]/h))))/(h^((x[j]/h)+1)*fn)
                }
                fhat1[j] <- 1/n * (sum(KLNormal [, j]))
                fhat2[j] <- 1/n * (sum(Kbs[, j]))
                fhat3[j] <- 1/n * (sum(KErlang[, j]))
                fhat4[j] <- 1/n * (sum(Kgamma[, j]))
                d1<-density(y,bw=h)
                plot(x,fhat1, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
                par(new=T)
                plot(x,fhat2, type = "s",pch="+",col="darkorange",xaxt="n", yaxt="n", xlab="", ylab="")
                par(new=T)
                plot(x,fhat3, type = "s",pch="+",col="green",xaxt="n", yaxt="n", xlab="", ylab="")
                points(x,fhat4,type="p",col="red")
                lines(d1,type="S",col="blue",lty=3,lwd=2)
                legend("topright", c("Real Density", "Density by Log Normal Kernel","Density by Birnbaum-Saunders Kernel", "Density by Erlang Kernel","Density by Gamma Kernel"),
                       col=c("blue", "black", "darkorange","green","red"), lty=c(1,2,3,4,6))
              }
  )

}





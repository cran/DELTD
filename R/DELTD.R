#' Density Plot by Gamma kernel.
#'
#' Plot Kernel density by using Gamma Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' denGamma(y,80,h)
#' @return Plot the density by using Gamma kernel
#' @export
denGamma<-function(y,k,h){
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

    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by Gamma Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}

#' Calculate Mean Square Error( MSE) when Gamma kernel is used.
#'
#' Calculate MSE by using Gamma Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use "Exp".
#'     If use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(100,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
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

#' Plot Density by Lognormal kernel.
#'
#' Plot Kernel density by using Lognormal Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' denLN(y,80,h)
#' @return Plot the density by using Lognormal kernel
#' @export
denLN<-function(y,k,h){
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
    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1,ylim=c(0,1), xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by Lognormal Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}

#' Calculate Mean Square Error( MSE) when LogNormal Kernel is used.
#'
#' Calculate MSE by using Lognormal Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use "Exp".
#'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(100,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
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

#' Plot Density by Erlang kernel.
#'
#' Plot Kernel density by using Erlang Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' denEr(y,80,h)
#' @return Plot the density by using Erlang kernel
#' @export
denEr<-function(y,k,h){
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
    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by Erlang Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}
#' Calculate Mean Square Error( MSE) when Erlang kernel is used.
#'
#' Calculate MSE by using Erlang Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use "Exp".
#'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(100,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
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

#' Plot Density by Birnbaum-Saunders kernel.
#'
#' Plot Kernel density by using Birnbaum-Saunders Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(23,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' denBS(y,80,h)
#' @return Plot the density by using Birnbaum-Saunders kernel
#' @export

denBS<-function(y,k,h){
  n<-length(y)
  x <- seq(min(y) + 0.05, max(y), length =k)
  Kbs <- matrix(rep(0, k * n), ncol = k)
  fhat <- rep(0, k)
  ###########BS###########
  for(j in 1:k) {
    for(i in 1:n) {
      Kbs[i, j] <-(1/(2*sqrt(2*h*pi)))*((sqrt(1/(x[j]*y[i])))+(sqrt(x[j]/(y[i]^3))))*exp(-(y[i]/(2*h*x[j]))+(1/h)-(x[j]/(2*h*y[i])))
      #Kbs[is.nan(Kbs)] <- 0
    }
    fhat[j] <- 1/n * (sum(Kbs[, j]))

    d1<-density(y,bw=h)
    plot(x,fhat, type = "s", ylab = "Density Function", lty = 1, xlab = "Time")
    lines(d1,type="p",col="red")
    legend("topright", c("Real Density", "Density by Birnbaum-Saunders Kernel"),
           col=c("red", "black"), lty=c(1,2))
  }}

#' Calculate Mean Square Error( MSE) when Birnbaum-Saunders kernel is used.
#'
#' Calculate MSE by using Birnbaum-Saunders Kernel.
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param type mention distribution of vector.If exponential distribution then use "Exp".
#'     if use gamma distribution then use "Gamma".If Weibull distribution then use "Weibull".
#' @import graphics
#' @import stats
#' @examples
#' y<-rexp(100,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
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

#' Plot the densities for comparison.
#'
#' Plot the Estimated Densities with Real Density for Comparison .
#' @param y a numeric vector of positive values.
#' @param k gird points.
#' @param h the bandwidth
#' @param comb mention the combination which kernel estimated densities are to be compared. If Lognormal and Birnbaum-Saunders kernel densities
#'           are to be compared along with real density then use "TLB". If Lognormal and Erlang then use "TLE". If Lognormal and Gamma to be compared
#'           then use "TLG". For Birnbaum-Saunders and Erlang use "TBE". For Birnbaum-Saunders and Gamma then use "TBG". For Erlang and Gamma use "TEG".
#'           For Lognormal, Birnbaum-Saunders and Erlang use "TLBE". For Lognormal, Birnbaum-Saunders and Gamma use "TLBG". For Birnbaum-Saunders, Erlang and Gamma
#'           use "TBEG". To compare all densities in one graph use "TLBEG".
#' @import graphics
#' @import stats
#' @examples
#' \dontshow{
#' y<-rexp(10,1)
#' h<-0.79 * IQR(y) * length(y) ^ (-1/5)
#' dencomb(y,20,h,"TLB")}
#' @return Plot of Estimated Densities with Real Data Density.
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





#' European Call Function
#'
#'#' This function allows you to simulate a call price using the Euler method.
#' @param X0 Initial value.
#' @param mu Drift coeficient.
#' @param sigma Diffusion coefficient.
#' @param Dt Step size.
#' @param t Initial time period.
#' @param T Final time period.
#' @param N Number of simulations.
#' @param K Strike price.
#' @param plt Plot? Defaults to FALSE.
#' @keywords Euler method
#' @importFrom graphics abline axis box grid matplot
#' @export
#' @author Fernando Teixeira
#' @examples
#' euler_call2(307.65, 0.75, 0.3, 0.001, 0, 1, 10000)
#'


euler_call2 <- function(X0, mu, sigma, Dt = 0.1, t=0, T=1, N = 1000, K = 0, plt=FALSE){
    set.seed(1)
    linhas = floor((T-t)/Dt + 1)            # número de linhas da matriz de resultados
    out <- matrix(stats::rnorm(linhas*N), nrow=linhas, ncol = N)
    sqrdt = sqrt(Dt)

    w=sim(out=out, X0=X0, mu=mu, Dt=Dt, sqrdt = sqrdt, sigma=sigma)


    Result = w[linhas,] - K
    Result[Result <= 0] = 0
    call = exp(-mu*T)*mean(Result)

    if (plt == TRUE){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 6){
            for (i in seq_along(w[1,])){
                points(w[,i], col='red', pch=20)
            }
        }
        if (K != 0){
            abline(h=K, col='red', lwd = 3)
        }
        box()
        axis(side = 1, at= seq(1,linhas,linhas-1), labels=c('0', as.character(linhas-1)))
        axis(side = 2, las = 1)
    }
    lista = list('realizations' = as.data.frame(w), 'mean' = mean(w[linhas,]),
                 'solution' = w[linhas,], 'eurocall' = call)
    return(invisible(lista))
}

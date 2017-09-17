#' European Call Function
#'
#' This function allows you to simulate a call price using the Milstein method.
#' @param X0 Initial value.
#' @param mu Drift coeficient.
#' @param sigma Diffusion coefficient.
#' @param Dt Step size.
#' @param t Initial time period.
#' @param T Final time period.
#' @param N Number of simulations.
#' @param K Strike price.
#' @param plt Plot? Defaults to FALSE.
#' @keywords milstein method
#' @export
#' @importFrom graphics abline axis box grid matplot
#' @author Fernando Teixeira
#' @examples
#' milstein_call(307.65, 0.75, 0.3, 0.001, 0, 1, 10000, 300)
#'
#'


milstein_call <- function(X0, mu, sigma, Dt, t, T, N, K = 0, plt=F){
    set.seed(1)
    X0 = rep(X0,N)
    linhas = floor((T-t)/Dt + 1)
    w = matrix(nrow=linhas, ncol = N)
    n = seq(from=Dt, to=T, by=Dt)
    w[1,] = X0
    l2 = linhas-1                           # linhas da matriz de números aleatórios
    bw = sqrt(Dt)*stats::rnorm(N*l2,0,1)           # números aleatórios criados
    bw = matrix(bw, nrow = l2, ncol = N)    # coloca em forma matricial
    b <- function(x) (sigma*t)


    for (i in seq_along(n)){
        w[i+1,] = w[i,] + mu *w[i,] *Dt + sigma * w[i,] * bw[i,] +
            (1/2) * sigma^2 * w[i,]  *(bw[i,]^2 - Dt)
    }

    Result = w[linhas,] - K
    Result[Result <= 0] = 0
    call = exp(-mu*T)*mean(Result)


    if (plt == T){
        matplot(w, ylim = c(min(w),max(w)), axes=F, lty=1, pch = 20 , type="l")
        grid()
        if (N <= 10 && linhas <= 5){
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

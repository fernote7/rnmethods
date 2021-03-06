#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param S Spot price
#' @param X Strike price
#' @param r Asset's rate of return
#' @param v Instantaneous variance
#' @param theta Long variance
#' @param rho Processes' correlation
#' @param k Rate at which v returns to theta
#' @param sigma Vol of vol
#' @param t Starting time, Default: 0
#' @param dt Stepsize, Default: NULL
#' @param tau Ending time, Default: 1
#' @param N Number of simulations
#' @return List with call price, values used to compute the call and all simulated paths.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{rnorm}}
#' @rdname hestoncallkj
#' @export
#' @importFrom stats rnorm


hestoncallkj <- function(S, X, r, v, theta, rho, k, sigma, t = 0, dt = NULL, tau = 1, N){

    cont = 1
    set.seed(103)
    if(is.null(dt)){ dt <- (T-t)/1000}
    v <- rep(v,N)
    theta<- rep(theta,N)
    sequencia <- seq(t,tau,dt)
    ST <- matrix(NA, length(sequencia), N) #transformar em matrix
    S <- log(S)
    j <-1

    for(i in seq(t,tau,dt)){
        Zv <- stats::rnorm(N)
        Zt <- stats::rnorm(N)
        Zs <- rho * Zv + sqrt(1 - rho^2) * Zt


        vt <- (v + k * theta * dt + sigma * sqrt(v) * Zv * sqrt(dt) +
                  (1/4) * sigma^2 * dt * ((Zv)^2 - 1))/(1 + k * dt)
        vt[vt <= 0] <- v[vt <= 0] + k * dt * (theta[vt <= 0] - v[vt <= 0]) +
                       sigma * sqrt(v[vt <= 0]) * Zv[vt <= 0] * sqrt(dt)
        v <- vt
        v[v<=0] <- 0
        vt[vt<=0] <- 0
        S <- S + (r - (v+vt)/4) * dt + rho * sqrt(v) * Zv * sqrt(dt) +
             (1/2) * (sqrt(v) + sqrt(vt)) * (Zs + rho * Zv) * sqrt(dt) +
             ((rho * sigma * dt)/2) * ((Zv)^2 - 1)

        S[S <= 0] = 0
        ST[j,] <- S
        j <- j + 1
        cont = cont + 1
        if(cont %% 50 == 0){print(cont)}
    }

    ST <- as.matrix(ST, ncol=N)
    Result <- exp(ST[nrow(ST),]) - X
    Result[Result <= 0] = 0
    call = mean(exp(-r*tau)*Result)

    lista = list('call' = call, 'Result' = Result, 'Spot' = ST)
    return(lista)

}

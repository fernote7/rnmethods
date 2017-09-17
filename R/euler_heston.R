hestoneuler <- function(S, X, r, v, theta, rho, k, sigma, t = 0, dt = NULL, tau = 1, N){

    cont <- 0
    if(is.null(dt)){ dt <- (tau-t)/1000}
    sequencia <- seq(t,tau,dt)
    ST <- matrix(NA, length(sequencia), N) #transformar em matrix
    aux <- NULL
    sqrt_dt <- sqrt(dt)
    j <-1

    for(i in sequencia){
        Zv <- stats::rnorm(N)
        Zt <- stats::rnorm(N)
        Zs <- rho * Zv + (sqrt(1 - (rho^2)) * Zt)
        aux <- v
        aux[v < 0] <- 0
        sqrt_aux <- sqrt(aux)
        S <- S * (1 + r * dt + sqrt_aux * Zs * sqrt_dt)
        S[S <= 0] = 0
        v <- v + k * dt * (theta - aux) + sigma * sqrt_aux * Zv * sqrt_dt
        ST[j,] <- S
        j <- j + 1
        cont <- cont + 1
        if(cont %% 50 == 0) {print(cont)}
    }

    rm(aux, v, Zv, Zt, Zs, S, j)
    ST <- as.matrix(ST, ncol=N)
    Result <- ST[nrow(ST),] - X
    Result[Result <= 0] = 0
    call = mean(exp(-r*(tau-t))*Result)
    lista = list('call' = call, 'Result' = Result, 'Spot' = ST)
    return(lista)
}

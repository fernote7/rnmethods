hestonea <- function(S, X, r, v, theta, rho, k, sigma, t = 0, dt = NULL, tau = 1, N){

    cont = 1
    if(is.null(dt)){ dt <- (tau-t)/1000}
    sequencia <- seq(t,tau,dt)
    ST <- matrix(NA, length(sequencia), N)
    j <-1
    d <- (4 * k * theta)/(sigma)^2
    c0 <- (sigma^2 * (1 - exp(-k*dt)))/(4*k)

    for(i in seq(t,tau,dt)){

        # sampling V
        lambda <- (4*k*exp(-k*dt)*v)/(sigma^2 * (1-exp(-k*dt)))
        vt <- c0 * stats::rchisq(n = N, df = d, ncp = lambda)

        # Sampling int{V}
        int_v <- dt * ((1/2) * v + (1/2) * vt)

        # Sampling int{v}dw
        int_vdw <- (1/sigma) * (vt - v - k * theta * dt + k  * int_v)

        # Sampling S
        m <- log(S) + (r * dt - (1/2) * int_v + rho * int_vdw)
        std <- sqrt((1 - rho^2)) * sqrt(int_v)
        S <- exp(m + std * stats::rnorm(N))
        v <- vt
        ST[j,] <- S
        j <- j + 1
        cont = cont + 1
        if(cont %% 50 == 0){print(cont)}
    }

    rm(v, S, j)
    ST <- as.matrix(ST, ncol=N)
    Result <- ST[nrow(ST),] - X
    Result[Result <= 0] = 0
    call = mean(exp(-r*tau)*Result)
    lista = list('call' = call, 'Result' = Result, 'Spot' = ST)
    return(lista)
}
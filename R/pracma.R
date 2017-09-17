phi_heston <- function(a, v0, v_t, d){

    gamma_a <- sqrt(k^2 - 2 * sigma^2 * 1i*a)
    gammadt <- gamma_a * (tau-t)
    sqrtv0vt <- sqrt(v0*v_t)
    delta <- -k * (tau - t)

    part1 <- (gamma_a * exp(-(gamma_a - k)/2 * (tau - t)) * (1 - exp(delta)))/
        (k * (1- exp(- gammadt)))

    part2 <- exp((v0+v_t)/(sigma^2) *
                     ((k * (1 + exp(delta)))/(1-exp(delta)) -
                        (gamma_a * (1 + exp(- gammadt)))/(1-exp(- gammadt))))


    part3 <- Bessel::BesselI(z = ((4 * gamma_a * sqrtv0vt)/(sigma^2) *
                                      exp(- gammadt/2)/
                                      (1 - exp(- gammadt))), nu = 0.5*d - 1) /
        Bessel::BesselI(z = ((4 * k * sqrtv0vt)/(sigma^2) * (exp(delta/2))/
                                 (1-exp(delta))), nu = 0.5*d - 1)

    result <- part1 * part2 * part3
    return (result)
}

integrate_gauss_laguerre <- function(f, M=30){

    points_weights <- gaussquad::laguerre.quadrature.rules(M)[[M]]

    int <- 0

    for(i in 1:M){

        int <- int + points_weights[i, 'w'] *
            f(points_weights[i, 'x']) * exp(points_weights[i, 'x'])
    }
    return (int)
}


intv <- function(n, cf, v_t){

    integrand <- function(x, phi = cf){

        f2 <- function(u){
            Im(phi(u) * exp(-1i * u * x)) /u
        }
        return(f2)
    }

    ## integrate to "cdf"

    F_x <- function (x) {

        y <- 0.5 - 1/pi * integrate_gauss_laguerre(integrand(x))

        return (y)

    }

    ## endsign

    endsign <- function(f, sign = 1, n) {
        browser()
        b <- sign
        while (all(sign * f(b) < 0)) b <- 10 * b
        return(b)
    }


    ## inversion
    spdf.lower = -Inf
    spdf.upper = Inf

    invcdf <- function(u) {
        subcdf <- function(t) {
            F_x(t) - u}
        cat("calculando root \n")
        return(pracma::fzero(subcdf, v, maxiter = 1000, tol = 0.01))
    }
    U <- stats::runif(n)
    sapply(U, invcdf)
}


cont <- 0

hestonea_mod <- function(S, X, r, v, theta, rho, k, sigma, t = 0, tau = 1, N){

    ST <- NULL
    d1 <- (4 * k * theta)/(sigma)^2
    c0 <- (sigma^2 * (1 - exp(-k*(tau-t))))/(4*k)

    # sampling V
    lambda <- (4*k*exp(-k*(tau-t))*v)/(sigma^2 * (1-exp(-k*(tau-t))))
    vt <- c0 * stats::rchisq(n = N, df = d1, ncp = lambda)

    # Sampling int{V}

    phi <- function(a, v0=v, v_t=vt, d=d1){phi_heston(a, v0=v, v_t=vt, d=d1)}
    int_v <- intv(N, cf = phi, v_t=vt)
    int_v <- unlist(int_v[1,])

    # Sampling int{v}dw
    int_vdw <- (1/sigma) * (vt - v - k * theta * (tau-t) + k  * int_v)


    # Sampling S
    if( int_v >= 0){
        m <- log(S) + (r * (tau - t) - (1/2) * int_v + rho * int_vdw)
        std <- sqrt((1 - rho^2)) * sqrt(int_v)
        S <- exp(m + std * stats::rnorm(N))
        v <- vt
        ST <- rbind(ST,S)
    } else {
        v <- vt
        ST <- rbind(ST,S)
    }

    ST <- as.matrix(ST, ncol=N)
    Result <- ST[nrow(ST),] - X
    Result[Result <= 0] = 0
    call = exp(-r*tau)*Result
    cont <<- cont + 1
    print(cont)
    lista = list('call' = call, 'Result' = Result, 'Spot' = ST)
    return(lista)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param a PARAM_DESCRIPTION
#' @param v0 PARAM_DESCRIPTION
#' @param v_t PARAM_DESCRIPTION
#' @param d PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[Bessel]{BesselI}}
#' @rdname phi_heston
#' @importFrom Bessel BesselI


phi_heston <- function(a, v0, v_t, d){

    gamma_a <- sqrt(k^2 - 2 * sigma^2 * 1i*a)
    gammadt <- gamma_a * (tau-t)
    sqrtv0vt <- sqrt(v0*v_t)
    delta <- -k * (tau-t)

    part1 <- (gamma_a * exp(-(gamma_a - k)/2 * (tau-t)) * (1 - exp(delta)))/
        (k * (1- exp(- gammadt)))

    part2 <- exp((v0+v_t)/(sigma^2) * ( (k * (1 + exp(delta)))/(1-exp(delta)) -
                        (gamma_a * (1 + exp(- gammadt)))/(1-exp(- gammadt))))

    part3 <- Bessel::BesselI(z = ((4 * gamma_a * sqrtv0vt)/(sigma^2) *
                                      exp(- gammadt/2)/
                                      (1 - exp(- gammadt))), nu = 0.5*d - 1) /
        Bessel::BesselI(z = ((4 * k * sqrtv0vt)/(sigma^2) * (exp(delta/2))/
                                 (1-exp(delta))), nu = 0.5*d - 1)

    return (part1 * part2 * part3)
}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @param cf PARAM_DESCRIPTION
#' @param v_t PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{integrate}},\code{\link[stats]{uniroot}},\code{\link[stats]{runif}}
#' @rdname intv
#' @importFrom stats integrate uniroot runif

intv <- function(n, cf, v_t){
    integrand <- function(x, phi = cf){

        f2 <- function(u){
            Im(phi(u) * exp(-1i * u * x)) /u
        }
        return(f2)
    }

    ## integrate to "cdf"
    F_x <- function (x) {
        y <- 0.5 - 1/pi * stats::integrate(integrand(x),  lower= 0, upper= 1000,
                                    rel.tol = 0.001, stop.on.error = FALSE)$value
        return(y)
    }


    ## endsign

    endsign <- function(f, sign = 1) {
        b <- sign
        while (sign * f(b) < 0) b <- 10 * b
        return(b)
    }

    ## inversion

    spdf.lower = -Inf
    spdf.upper = Inf
    invcdf <- function(u) {
        subcdf <- function(t) F_x(t) - u
        if (spdf.lower == -Inf)
            spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf)
            spdf.upper <- endsign(subcdf)
        return(stats::uniroot(subcdf, lower=spdf.lower, upper=spdf.upper,
                              tol = 0.001220703)$root)
    }
    U <- stats::runif(n)
    sapply(U, invcdf)
}


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
#' @param tau Ending time, Default: 1
#' @return List with call price, values used to compute the call and all simulated paths.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[stats]{rchisq}},\code{\link[stats]{rnorm}}
#' @rdname hestonea
#' @export
#' @importFrom stats rchisq rnorm


hestonea <- function(S, X, r, v, theta, rho, k, sigma, t = 0, tau = 1){

    d1 <- (4 * k * theta)/(sigma)^2
    c0 <- (sigma^2 * (1 - exp(-k*tau)))/(4*k)
    dt <- (tau-t)
    ST <- NULL

    # sampling V

    lambda <- (4*k*exp(-k*dt)*v)/(sigma^2 * (1-exp(-k*dt)))
    vt <- c0 * stats::rchisq(n = 1, df = d1, ncp = lambda)

    # Sampling int{V}

    phi <- function(a, v0=v, v_t=vt, d=d1){phi_heston(a, v0=v, v_t=vt, d=d1)}
    int_v <- intv(1, cf = phi, v_t=vt)

    # Sampling int{v}dw
    int_vdw <- (1/sigma) * (vt - v - k * theta * dt + k  * int_v)

    # Sampling S
    if( int_v >= 0){
        m <- log(S) + (r * (tau - t) - (1/2) * int_v + rho * int_vdw)
        std <- sqrt((1 - rho^2)) * sqrt(int_v)
        S <- exp(m + std * stats::rnorm(1))
        v <- vt
        ST <- S
    } else {
        v <- vt
        ST <- rbind(ST,NA)}

    Result <- ST - X
    Result[Result <= 0] = 0
    call = exp(-r*tau)*Result
    lista = list('call' = call, 'Result' = Result, 'Spot' = ST)
    return(lista)
}

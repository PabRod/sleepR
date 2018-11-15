#' Phillips and Robinson's model flow
#'
#' Generates the flow (right hand side of the differential equation) for the Phillips and Robinson's model
#'
#' @param time The time (in h)
#' @param y  The state
#' @param parms Model's parameters (optional, see \code{\link{philrob_default_parms}})
#' @param C External forcing as a function of t (in h) (optional, see \code{\link{philrob_default_forcing}})
#'
#' @return The flow (right hand side of the differential equation)
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}, \link{philrob_default_parms}, \link{saturating_function}}
#'
#' @examples
#' \dontrun{
#' t <- 0
#' y <- c(Vv = 1, Vm = 1, H = 1)
#' dy <- dPhilrob(t, y)
#' }
dPhilrob <- function(time, y, parms = philrob_default_parms(), C = philrob_default_forcing(parms = philrob_default_parms())) {
  with(as.list(c(y, parms)), {

    S <- function(V) saturating_function(V, parms)

    # Dynamics
    dVv <- (            -vvm*S(Vm) + vvh*H - vvc*C(time) - Vv)/tauv # Ventro-lateral preoptic area activity
    dVm <- (-vmv*S(Vv)                     + vmaSa             - Vm)/taum # Mono-aminergic group activity
    dH <-  (              mu*S(Vm)                             -  H)/Xi # Homeostatic pressure

    return(list(c(dVv, dVm, dH)))
  })
}

#' Saturating function
#'
#' @param V voltage (in mV)
#' @param parms Model's parameters (optional, see \code{\link{philrob_default_parms}})
#'
#' @return The saturated effect of the voltage
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}, \link{philrob_default_parms}}
#'
#' @examples
#' \dontrun{
#' saturating_function(1)
#' }
saturating_function <- function(V, parms = philrob_default_parms()) {
  with(as.list(parms), {
    Qmax/(1 + exp((theta-V)/sigma))
  })
}

#' Derivative of the saturating function
#'
#' @param V voltage (in mV)
#' @param parms Model's parameters (optional, see \code{\link{philrob_default_parms}})
#'
#' @return The derivative of the saturation function
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}, \link{saturating_function}}
#'
#' @examples
#' \dontrun{
#' dSaturating_function(1)
#' }
dSaturating_function <- function(V, parms = philrob_default_parms()) {
  S <- function(V) saturating_function(V, parms)

  with(as.list(parms), {
    exp((theta-V)/sigma)/(Qmax * sigma) * S(V)^2 # Derived analitically
  })
}

#' Jacobian of the Phillips and Robinson's model equations
#'
#' @param y  The state
#' @param parms Model's parameters (optional, see \code{\link{philrob_default_parms}})
#'
#' @return The jacobian of the model's equations (for constant forcing)
#'
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}, \link{dPhilrob}}
#'
#' @examples
#' y <- c(Vv = 1, Vm = 1, H = 1)
#' J <- philrob_jacobian(y)
philrob_jacobian <- function(y, parms = philrob_default_parms()) {

  # Auxiliary functions
  dS <- function (V) dSaturating_function(V, parms)

  with(as.list(c(y, parms)), {
    # Analytically derived jacobian
    vec <- c(-1/tauv            , -vvm/tauv * dS(Vm), vvh/tauv,
             -vmv/taum * dS(Vv) , -1/taum           , 0       ,
             0                  , mu/Xi * dS(Vm)    , -1/Xi)

    # Reshape as matrix
    matrix(vec, nrow = 3, ncol = 3, byrow = TRUE)
  })

}

#' Solve Phillip and Robinson's model
#'
#' Solves the Phillip and Robinson's model for the given times, initial condition and parameters
#'
#' @param ts Vector of times (in h)
#' @param y0 Initial condition
#' @param parms Model parameters (optional, see \code{\link{philrob_default_parms}})
#' @param tStabil Stabilization time (in h)
#' @param ... Additional arguments passed to the \code{\link{ode}} integrator
#'
#' @return Results of the simulation, including times, states and asleep/awake status
#'
#' @export
#' @importFrom deSolve ode
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{dPhilrob}, \link{philrob_default_parms}, \link{ode}}
#'
#' @examples
#' y0 <- c(Vv = -13, Vm = 1, H = 10)
#' nDays <- 30
#' ts <- seq(0, nDays*24, length.out = nDays*24*20)
#' sol <- philrob(ts, y0)
philrob <- function(ts, y0, parms = philrob_default_parms(), tStabil = 0, ...) {

  # Stabilize the solution if required...
  if(tStabil != 0) {
    sol_aux <- ode(y = y0,
                func = dPhilrob,
                times = seq(ts[1], ts[1] + tStabil, length.out = 100),
                parms = parms,
                ...)

    last <- tail(sol_aux, 1) # ... by using the last simulated point...
    y0 <- c(Vv = last[,'Vv'], Vm = last[,'Vm'], H = last[,'H']) # ... as new initial condition
  }

  # Solve
  sol <- ode(y = y0,
             func = dPhilrob,
             times = ts,
             parms = parms,
             ...)

  # Transform into data frame
  sol <- as.data.frame(sol)

  # Implement sleep-awake criterion
  for (i in 1:length(ts)-1) {
    dH <- sol$H[i+1] - sol$H[i]
    sol$asleep[i] <- (dH < 0)
  }

  return(sol)

}

#' Default parameters of Phillips and Robinson's model
#'
#' Loads the parameters used in Phillips and Robinson's model
#'
#' @return The default parameters for Phillips and Robinson's model
#'
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob_default_forcing}, \link{philrob}, \link{dPhilrob}}
#'
#' @examples
#' parms <- philrob_default_parms()
philrob_default_parms <- function() {

  T <- 24 # h

  parms = c(Qmax = 100*3600, # h^-1
            theta = 10, # mV
            sigma = 3, # mV
            vmaSa = 1, # mV
            vvm = 1.9/3600, # mV h (negative effect)
            vmv = 1.9/3600, # mV h (negative effect)
            vvc = 6.3, # mV (negative effect)
            vvh = 0.19, # mV nM^-1
            Xi = 10.8, # h
            mu = 1e-3, # nM h
            taum = 10/3600, # h
            tauv = 10/3600, # h
            w = 2*pi/T, # h^-1
            alpha = 0.0) # rad

}

#' Default forcing function
#'
#' @param parms Parameters (optional)
#'
#' @return The standard offset + cosine forcing function
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}, \link{dPhilrob}}
#'
#' @examples
#' \dontrun{
#' C <- function(t) { philrob_default_forcing(t) }
#' t <- seq(0, 3, length.out = 10)
#' }
philrob_default_forcing <- function(parms = philrob_default_parms()) {
    C <- function(t) { 0.5*(1 + cos(parms['w']*(t - parms['alpha']))) }
}

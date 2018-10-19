#' Phillips and Robinson's model flow
#'
#' Generates the flow (right hand side of the differential equation) for the Phillips and Robinson's model
#'
#' @param time The time (in h)
#' @param y  The state
#' @param parms Model's parameters (optional, see \code{\link{philrob_default_parms}})
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
#' @seealso \code{\link{philrob}, \link{philrob_default_parms}}
#'
#' @examples
#' \dontrun{
#' t <- 0
#' y <- c(Vv = 1, Vm = 1, H = 1)
#' dy <- dPhilrob(t, y)
#' }
dPhilrob <- function(time, y, parms = philrob_default_parms()) {
  with(as.list(c(y, parms)), {

    # Auxiliary functions

    ## Sigmoid growth
    S <- function(V) {
      Qmax/(1 + exp((theta-V)/sigma))
    }

    ## Forcing (approximation)
    C <- function(t) {
      0.5*(1 + cos(w*(t - alpha)))
    }

    # Dynamics
    dVv <- (            -vvm*S(Vm) + vvh*H - vvc*C(time) - Vv)/tauv # Ventro-lateral preoptic area activity
    dVm <- (-vmv*S(Vv)                     + vmaSa       - Vm)/taum # Mono-aminergic group activity
    dH <-  (              mu*S(Vm)                       -  H)/Xi # Homeostatic pressure

    return(list(c(dVv, dVm, dH)))
  })
}

#' Solve Phillip and Robinson's model
#'
#' Solves the Phillip and Robinson's model for the given times, initial condition and parameters
#'
#'
#' @param ts Vector of times (in h)
#' @param y0 Initial condition
#' @param parms Model parameters (optional, see \code{\link{philrob_default_parms}})
#' @param method Integration method (optional, by default 'lsoda')
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
#' @seealso \code{\link{dPhilrob}, \link{philrob_default_parms}}
#'
#' @examples
#' y0 <- c(Vv = -13, Vm = 1, H = 10)
#' nDays <- 30
#' ts <- seq(0, nDays*24, length.out = nDays*24*20)
#' sol <- philrob(ts, y0)
philrob <- function(ts, y0, parms = philrob_default_parms(), method = 'lsoda') {

  # Solve
  sol <- ode(y = y0,
             func = dPhilrob,
             times = ts,
             parms = parms,
             method = method)

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
#' Loads the parameters used in Phillips and Robinson's model model
#'
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @return The default parameters for Strogatz's model
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
            D = 0.77, # mV
            Da = 0.42, # mV
            alpha = 0.0) # rad

}

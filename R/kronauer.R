#' Kronauer's 1990 model flow
#'
#' Generates the flow (right hand side of the differential equation) for the Kronauer's 1990 model
#'
#' @param time The time (in h)
#' @param y  The state
#' @param I Illumination function (in lux)
#' @param parms Model's parameters (optional, see \code{\link{kronauer_default_parms}})
#'
#' @return The flow (right hand side of the differential equation)
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Jewett ME, Kronauer RE.
#' Refinement of Limit Cycle Oscillator Model of the Effects of Light on the Human Circadian Pacemaker.
#' J Theor Biol. 1998 Jun 21;192(4):455–65. \url{https://www.sciencedirect.com/science/article/pii/S0022519398906671}
#'
#' @seealso \code{\link{kronauer}, \link{kronauer_default_parms}}
#'
#' @examples
#' \dontrun{
#' t <- 0
#' y <- c(x = 1, xc = 1)
#' dy <- dKronauer(t, y)
#' }
dKronauer <- function(time, y, I, parms = kronauer_default_parms()) {
  with(as.list(c(y, parms)), {

    # Auxiliary functions

    ## Drive of light
    B <- function(x, I) {
      (1 - m*x)*c*I^n
    }

    # Dynamics
    dx <-  (pi/12)*(xc + mu*(x - 4/3*x^3) + B(x, I(time))) # Temperature
    dxc <- (pi/12)*(q*xc*B(x, I(time)) -(24/taux)^2*x) # Auxiliary

    return(list(c(dx, dxc)))
  })
}

#' Solve Kronauer's 1990 model
#'
#' Solves the Kronauer's 1990 model for the given times, initial condition and parameters
#'
#' @param ts Vector of times (in h)
#' @param y0 Initial condition
#' @param parms Model parameters (optional, see \code{\link{kronauer_default_parms}})
#' @param ... Additional arguments passed to the \code{\link{ode}} integrator
#'
#' @return Results of the simulation, including times and states
#'
#' @export
#' @importFrom deSolve ode
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Jewett ME, Kronauer RE.
#' Refinement of Limit Cycle Oscillator Model of the Effects of Light on the Human Circadian Pacemaker.
#' J Theor Biol. 1998 Jun 21;192(4):455–65. \url{https://www.sciencedirect.com/science/article/pii/S0022519398906671}
#'
#' @seealso \code{\link{dKronauer}, \link{kronauer_default_parms}, \link{ode}}
#'
#' @examples
#' y0 <- c(x = 1, xc = 0)
#' nDays <- 5
#' ts <- seq(0, nDays*24, length.out = nDays*24*20)
#' sol <- kronauer(ts, y0)
kronauer <- function(ts, y0, parms = kronauer_default_parms(), ...) {

  # Solve
  sol <- ode(y = y0,
             func = dKronauer,
             times = ts,
             parms = parms,
             ...)

  # Transform into data frame
  sol <- as.data.frame(sol)

  # Implement sleep-awake criterion

  return(sol)

}

#' Default parameters of Kronauer's 1990 model
#'
#' Loads the parameters used in Kronauer's 1990 model
#'
#' @return The default parameters for Kronauer's 1990 model
#'
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Jewett ME, Kronauer RE.
#' Refinement of Limit Cycle Oscillator Model of the Effects of Light on the Human Circadian Pacemaker.
#' J Theor Biol. 1998 Jun 21;192(4):455–65. \url{https://www.sciencedirect.com/science/article/pii/S0022519398906671}
#'
#' @seealso \code{\link{kronauer}, \link{dKronauer}}
#'
#' @examples
#' parms <- kronauer_default_parms()
kronauer_default_parms <- function() {

  parms = c(mu = 0.13,
            taux = 24.2, # h
            q = 1/3,
            m = 1/3,
            c = 0.018, # lux^-1/n
            n = 1/3)

}

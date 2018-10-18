#' Strogatz's model flow
#'
#' Generates the flow (right hand side of the differential equation) for the Strogatz's model
#'
#' @param time The time (in h)
#' @param y  The state
#' @param parms Model's parameters (optional, see \code{\link{strogatz_default_parms}})
#'
#' @return The flow (right hand side of the differential equation)
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Strogatz, S. H. (1987).
#' Human sleep and circadian rhythms: a simple model based on two coupled oscillators.
#' Journal of Mathematical Biology, 25(3), 327–347. \url{http://doi.org/10.1007/BF00276440}
#'
#' @seealso \code{\link{strogatz}, \link{strogatz_default_parms}}
#'
#' @examples
#' t <- 0
#' y <- c(th1 = 0.1, th2 = 0.2)
#' dy <- dStrogatz(t, y)
dStrogatz <- function(time, y, parms = strogatz_default_parms()) {

  with(as.list(c(y, parms)), {

    dth1 <- w1 - C1*cos(2*pi*(th1-th2)) # Phase of circadian oscillator
    dth2 <- w2 + C2*cos(2*pi*(th1-th2)) # Phase of sleep-wake oscillator

    return(list(c(dth1, dth2)))
  })

}

#' Solve Strogatz's model
#'
#' Solves the Strogatz's model for the given times, initial condition and parameters
#'
#'
#' @param ts Vector of times (in h)
#' @param y0 Initial condition
#' @param parms Model parameters (optional, see \code{\link{strogatz_default_parms}})
#'
#' @return Results of the simulation, including times, states and asleep/awake status
#' @export
#' @importFrom deSolve ode
#' @importFrom dplyr mutate
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Strogatz, S. H. (1987).
#' Human sleep and circadian rhythms: a simple model based on two coupled oscillators.
#' Journal of Mathematical Biology, 25(3), 327–347. \url{http://doi.org/10.1007/BF00276440}
#'
#' @seealso \code{\link{dStrogatz}, \link{strogatz_default_parms}}
#'
#' @examples
#' y0 <- c(th1 = 0.1, th2 = 0.05)
#' nDays <- 30
#' ts <- seq(0, nDays*24, length.out=nDays*24*10)
#' ys <- strogatz(ts, y0)
strogatz <- function(ts, y0, parms = strogatz_default_parms()) {

  # Solve
  sol <- ode(y = y0,
            func = dStrogatz,
            times = ts,
            parms = parms)

  # Transform into data frame
  sol <- as.data.frame(sol)

  # Implement sleep-awake criterion
  sol <- mutate(sol, th1 = th1%%1, th2 = th2%%1)
  sol <- mutate(sol, asleep = (th2 >= 0.0) & (th2 <= 1/3))

}

#' Default parameters of Strogatz's model
#'
#' Loads the parameters used in Strogatz's model
#'
#' Strogatz, S. H. (1987).
#' Human sleep and circadian rhythms: a simple model based on two coupled oscillators.
#' Journal of Mathematical Biology, 25(3), 327–347. \url{http://doi.org/10.1007/BF00276440}
#'
#' @return The default parameters for Strogatz's model
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Strogatz, S. H. (1987).
#' Human sleep and circadian rhythms: a simple model based on two coupled oscillators.
#' Journal of Mathematical Biology, 25(3), 327–347. \url{http://doi.org/10.1007/BF00276440}
#'
#' @seealso \code{\link{strogatz}, \link{dStrogatz}}
#'
#' @examples
#' parms <- strogatz_default_parms()
strogatz_default_parms <- function() {

  parms = c(w1 = 1/24, # h^-1
            w2 = 0.84/24, # h^-1
            C1 = 0/24, # h^-1
            C2 = 0.16/24) # h^-1

}

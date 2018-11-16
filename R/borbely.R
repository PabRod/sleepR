#' Solve Borbely's model
#'
#' Solves the Borbely's sleep-wake model for the given times, initial condition and parameters
#'
#' @param ts Vector of times (in h)
#' @param y0 Initial condition
#' @param awake0 Initial state (TRUE for awake, FALSE for asleep)
#' @param parms Model parameters (optional, see \code{\link{borbely_default_parms}})
#'
#' @return Results of the simulation, including times, states and asleep/awake status
#'
#' @export
#' @importFrom deSolve ode
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Skeldon AC, Dijk D-J, Derks G. Mathematical Models for Sleep-Wake Dynamics:
#' Comparison of the Two-Process Model and a Mutual Inhibition Neuronal Model.
#' PLoS One. 2014 Aug 1;9(8):e103877. Available from: http://dx.plos.org/10.1371/journal.pone.0103877
#'
#' @seealso \code{\link{borbely_default_parms}}
#'
#' @examples
#' ts <- seq(0, 10, length.out = 2000)
#' sol <- borbely(seq(0, 10, length.out = 2000), 1)
borbely <- function(ts, y0, awake0 = FALSE, parms = borbely_default_parms()) {

  # Auxiliary functions
  C <- function(t) as.numeric(sin(parms['w']*t - parms['alpha'])) # Circadian process

  # Upper and lower turning bounds
  upper <- as.numeric(parms['Hu0'] + parms['a']*C(ts))
  lower <- as.numeric(parms['Hl0'] + parms['a']*C(ts))

  # Pieces to construct the piecewise function
  Hs <- function(t, H0) as.numeric(H0*exp(-t/parms['Xis'])) # Sleeping
  Ha <- function(t, H0) as.numeric(parms['mu'] + (H0 - parms['mu'])*exp(-t/parms['Xiw'])) # Awake
  Hc <- function(t, H0, awake) awake*Ha(t, H0) + (!awake)*Hs(t, H0); # Auxiliary function combining both previous

  # Initialize containers
  H <- as.vector(matrix(NaN, length(ts), 1))
  awake <- as.logical(matrix(NaN, length(ts), 1))
  H[1] <- y0
  awake[1] <- awake0

  # Simulate
  for(i in 2:length(ts)) {
    # Update homeostatic pressure
    H[i] <- Hc(ts[i] - ts[i-1], H[i-1], awake[i-1])

    # Set triggers
    wake_trigger <- ( (!awake[i-1]) & (H[i] <= lower[i]) ) # Sleeping and recovering
    sleep_trigger <- ( (awake[i-1]) & (H[i] >= upper[i]) ) # Awake and tired

    # Update awake/sleeping status
    if (wake_trigger) {
      awake[i] <- TRUE # Good morning!
    } else if (sleep_trigger) {
      awake[i] <- FALSE # Good night!
    } else {
      awake[i] <- awake[i-1] # Keep status
    }
  }

  # Transform into data frame
  sol <- data.frame(time = ts, H, awake, upper, lower)

  return(sol)

}

#' Plots the time series generated with borbely
#'
#' @param sol Solution of a borbely simulation
#'
#' @return The simulation's plot
#'
#' @export
#' @importFrom graphics plot lines axis
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Skeldon AC, Dijk D-J, Derks G. Mathematical Models for Sleep-Wake Dynamics:
#' Comparison of the Two-Process Model and a Mutual Inhibition Neuronal Model.
#' PLoS One. 2014 Aug 1;9(8):e103877. Available from: http://dx.plos.org/10.1371/journal.pone.0103877
#'
#' @seealso \code{\link{borbely}}
#'
#' @examples
#' # Generate a solution
#' y0 <- 0.4
#' nDays <- 3
#' ts <- seq(0, nDays*24, length.out = nDays*24*20)
#' sol <- borbely(ts, y0)
#' # Plot it
#' borbelyPlot(sol)
borbelyPlot <- function(sol) {
  ts_days <- sol$time / 24 # Use days

  plot(ts_days, sol$H, type = 'l', col = 'blue',
       ylim = c(min(sol$lower), max(sol$upper)),
       xlab = 'time (d)', ylab = 'sleep pressure')
  lines(ts_days, sol$upper, lty = 'dashed')
  lines(ts_days, sol$lower, lty = 'dashed')
}

#' Default parameters of Borbely's model
#'
#' Loads the parameters used in Borbely's model (see figure 1.b in the reference)
#'
#' @return The default parameters for Borbely's model
#'
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Skeldon AC, Dijk D-J, Derks G. Mathematical Models for Sleep-Wake Dynamics:
#' Comparison of the Two-Process Model and a Mutual Inhibition Neuronal Model.
#' PLoS One. 2014 Aug 1;9(8):e103877. Available from: http://dx.plos.org/10.1371/journal.pone.0103877
#'
#' @seealso \code{\link{borbely}}
#'
#' @examples
#' parms <- borbely_default_parms()
borbely_default_parms <- function() {

  parms = c(Xis = 4.2, # h
            Xiw = 18.2, # h
            mu = 1, # 1
            a = 0.1, # 1
            w = 2*pi/24, # h^-1
            alpha = 0, # h
            Hu0 = 0.6, # 1
            Hl0 = 0.17) # 1

}

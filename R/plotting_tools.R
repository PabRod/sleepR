#' Raster plot / somnogram
#'
#' @param sol Simulation (output of \code{\link{strogatz}})
#' @param dailySamples Number of samples per day (optional)
#' @param asleepCol Color for asleep state (optional)
#' @param awakeCol Color for awake state (optional)
#'
#' @return Plots the raster plot / somnogram
#'
#' @export
#' @importFrom graphics image
#' @importFrom utils tail
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @examples
#' # Simulate a solution
#' y0 <- c(th1 = 0.1, th2 = 0.05)
#' nDays <- 60.2
#' ts <- seq(0, nDays*24, length.out=nDays*24*20)
#' sol <- strogatz(ts, y0)
#' # Create raster plot
#' rasterPlot(sol)
rasterPlot <- function(sol, dailySamples = 480, asleepCol = 'grey26', awakeCol = 'cyan' ) {
  # Use day as time unit
  sol$time_days <- sol$time / 24

  t_last <- tail(sol$time_days, 1)
  nDays <- floor(t_last)

  # Get data to plot
  asleepRaster <- reshape_as_raster(sol, dailySamples = dailySamples)
  xs <- seq(0, 48, length.out =  2 * dailySamples)
  ys <- seq(0, nDays-1)

  # And plot it
  image(x = xs, y = ys, z = t(asleepRaster), xlab = 'h', ylab = 'd', col = c(awakeCol, asleepCol), axes = FALSE)
  axis(1, at = c(0, 6, 12, 18, 24, 30, 36, 42, 48), labels = c(0, 6, 12, 18, 24, 6, 12, 18, 24))
  axis(2, at = seq(0, nDays))
}

#' Plots the time series generated with philrob
#'
#' @param sol Solution of a philrob simulation
#'
#' @return The simulation's plot
#'
#' @export
#' @importFrom graphics plot lines axis
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @references
#' Phillips AJK, Robinson PA.
#' A Quantitative Model of Sleep-Wake Dynamics Based on the Physiology of the Brainstem Ascending Arousal System.
#' J Biol Rhythms. 2007 Apr 29;22(2):167–79. Available from: http://journals.sagepub.com/doi/10.1177/0748730406297512
#'
#' @seealso \code{\link{philrob}}
#'
#' @examples
#' # Generate a solution
#' y0 <- c(Vv = -13, Vm = 1, H = 10)
#' nDays <- 30
#' ts <- seq(0, nDays*24, length.out = nDays*24*20)
#' sol <- philrob(ts, y0)
#' # Plot it
#' philrobPlot(sol)
philrobPlot <- function(sol) {
  ts_days <- sol$time / 24 # Use days
  states <- subset(sol, select = c('Vv', 'Vm', 'H')) # Used to define limits

  plot(ts_days, sol$Vm, type = 'l', col = 'blue',
       ylim = c(min(states), max(states)),
       xlab = 'time (d)', ylab = 'states')
  lines(ts_days, sol$Vv, col = 'red')
  lines(ts_days, sol$H, col = 'green4')
}

#' Reshape as raster vector
#'
#' @param sol Simulation (output of \code{\link{strogatz}})
#' @param dailySamples Number of samples per day (optional)
#'
#' @return The asleep vector reshaped as a raster vector
#'
#' @importFrom dplyr filter
#' @importFrom stats approxfun
#' @importFrom utils tail
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @examples
#' \dontrun{
#' # Simulate a solution
#' y0 <- c(th1 = 0.1, th2 = 0.05)
#' nDays <- 60.2
#' ts <- seq(0, nDays*24, length.out=nDays*24*20)
#' sol <- strogatz(ts, y0)
#' # Get raster sleep vector
#' asleep_raster <- reshape_as_raster(sol)
#' }
reshape_as_raster <- function(sol, dailySamples = 480) {

  # Use day as time unit
  sol$time_days <- sol$time / 24

  # Remove last non-integer day
  t_last <- tail(sol$time_days, 1)
  nDays <- floor(t_last)
  sol <- filter(sol, time_days <= nDays)

  # Resample at a higher resolution
  ts_interp <- seq(0, nDays, length.out = 2 * dailySamples * nDays)
  interpolator <- approxfun(sol$time_days, sol$asleep, method = 'constant')

  ## Reshape the vector as a raster plot
  #
  # Row 1: days 1 and 2
  # Row 2: days 2 and 3
  # Row 3: days 3 and 4
  # ...
  # Row N: days N and N+1
  asleep <- matrix(NaN, nrow=nDays-1, ncol = 2 * dailySamples)
  for (i in 1:nDays-1) {
    ts_eval <- seq(i-1, i+1, length.out = 2 * dailySamples) # Evaluate two days per row
    asleep[i, ] <- interpolator(ts_eval)
  }

  return(asleep)
}

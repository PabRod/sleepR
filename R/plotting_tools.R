#' Raster plot / somnogram
#'
#' @param sol Simulation (output of \code{\link{strogatz}})
#' @param dailySamples Number of samples per day (optional)
#' @param asleepCol Color for asleep state (optional)
#' @param awakeCol Color for awake state (optional)
#'
#' @return Plots the raster plot / somnogram
#' @export
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
  image(x = xs, y = ys, z = t(asleepRaster), xlab = 'h', ylab = 'd', col = c(awakeCol, asleepCol))
}

#' Reshape as raster vector
#'
#' @param sol Simulation (output of \code{\link{strogatz}})
#' @param dailySamples Number of samples per day (optional)
#'
#' @return The asleep vector reshaped as a raster vector
#' @importFrom dplyr filter
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

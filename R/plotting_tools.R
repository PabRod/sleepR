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

#' Lissajous figure
#'
#' @param times Times vector
#' @param ys State vector
#' @param wx (Optional) Reference frequency. Default: 2*pi/24
#' @param Ax (Optional) Reference amplitude. Default: 1
#' @param alphax (Optional) Reference phase delay. Default: 0
#' @param ... (Optional) Graphical parameters
#'
#' @return Plots the Lissajous figure
#' @export
#'
#' @author Pablo Rodríguez-Sánchez (\url{https://pabrod.github.io})
#'
#' @examples
#' # Simulate a solution
#' nDays <- 8
#' ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
#' y0 <- c(Vv = -12.6404, Vm = 0.8997, H = 12.5731) # This point is already close to the attractor
#' sol <- philrob(ts, y0, method = 'lsode')
#' # Create the figure
#' lissajous_figure(ts, sol$H)
lissajous_figure <- function(times, ys, wx = 2*pi/24, Ax = 1, alphax = 0, ...) {

  lissajous <- function(times, ys, wx = 2*pi/24, Ax = 1, alphax = 0) {
    xs <- Ax*sin(wx*times + alphax)

    lis <- data.frame(times, xs, ys)
  }

  lis <- lissajous(times, ys, wx = wx, Ax = Ax, alphax = alphax)

  plot(lis$xs, lis$ys, ...)
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

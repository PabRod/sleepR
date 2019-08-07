## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  
  library(sleepR)
)

## ----Problem posing with default pars, include = TRUE--------------------
## Problem setting
y0 <- c(x = 1, xc = 0) # Initial conditions

nDays <- 12
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

I <- function(t) { (t%%24 > 8)*10000 + 
                  !(t%%24 > 8)*150} # Illumination function (in lux)


## ----Simulate, include = TRUE--------------------------------------------
# Simulate
sol <- kronauer(ts, y0, I)

## ----Plot illumination, warning = FALSE, echo = FALSE, fig.show = 'hold'----
ts_days <- sol$time / 24 # Day is a better unit for plots
plot(ts_days, I(sol$time), type = 'l', col = 'orange', xlab = 'Time (days)', ylab = 'lux')
title('Light intensity')


## ----Plot time series, warning = FALSE, echo = FALSE, fig.show = 'hold'----

plot(ts_days, sol$x, type = 'l', col = 'red', xlab = 'Time (days)', ylab = 'states')
lines(ts_days, sol$xc)
title('Time series')

amplitude <- sqrt(sol$x^2 + sol$xc^2)
plot(ts_days, amplitude, type = 'l', xlab = 'Time (days)', ylab = 'amplitude')
title('Time series')

plot(sol$x, sol$xc, type = 'l', xlab = 'x', ylab = 'xc')
title('Phase plane')



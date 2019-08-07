## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(sleepR)

## ---- Plot-saturation, include = TRUE, echo = FALSE, fig.height = 3------
parms <- philrob_default_parms()
S <- function(V) {
  parms['Qmax']/(1 + exp((parms['theta']-V)/parms['sigma']))
}

xs <- seq(-10, 35, length.out = 100)
plot(xs, S(xs), type = 'l', xlab = 'V', ylab = 'S(V)', col = 'red', axes = FALSE)

# Asymptote
lines(c(-10, 35), c(parms['Qmax'], parms['Qmax']), lty = 'dashed')

# Half life
lines(c(-10, parms['theta']), c(parms['Qmax']/2, parms['Qmax']/2), lty = 'dashed')
lines(c(parms['theta'], parms['theta']), c(0, parms['Qmax']/2), lty = 'dashed')

axis(2, at = c(0, parms['Qmax']/2, parms['Qmax']), labels = c('0', 'Qmax/2', 'Qmax'))
axis(1, at = c(-10, 0, parms['theta'], 35), labels = c('', '0', 'theta', ''))

## ----Problem posing with default pars, include = TRUE--------------------
## Problem setting
y0 <- c(Vv = -13, Vm = 1, H = 10) # Initial conditions

nDays <- 5
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

# Simulate
sol <- philrob(ts, y0)


## ----Problem posing with overriden pars, include = TRUE------------------
## Problem setting
y0 <- c(Vv = -13, Vm = 1, H = 10) # Initial conditions

nDays <- 3
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

parms <- philrob_default_parms() # Load default parameters...
parms['vvc'] <- 6 # .. and modify one

# Simulate
sol <- philrob(ts, y0, parms)

## ----Problem posing with custon forcing, include = TRUE------------------
## Problem setting
y0 <- c(Vv = -13, Vm = 1, H = 10) # Initial conditions

nDays <- 3
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

C <- function(t) { 0 }

# Simulate
sol <- philrob(ts, y0, parms = philrob_default_parms(), forcing = C)

## ----Problem posing with stabilization run, include = TRUE---------------
## Problem setting
y0 <- c(Vv = -13, Vm = 1, H = 10) # Initial conditions

nDays <- 5
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

# Simulate
sol <- philrob(ts, y0, tStabil = 3*24)


## ----Display output, echo=FALSE, results='asis'--------------------------
knitr::kable(head(sol, 5))

## ----Plotting, include = TRUE, fig.show='hold'---------------------------
philrobPlot(sol)
rasterPlot(sol)


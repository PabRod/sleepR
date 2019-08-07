## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(sleepR)

## ----Problem posing with default pars, include = TRUE--------------------
## Problem setting
y0 <- c(th1 = 0.1, th2 = 0.05) # Initial conditions

nDays <- 5
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

# Simulate
sol <- strogatz(ts, y0)


## ----Problem posing with overriden pars, include = TRUE------------------
## Problem setting
y0 <- c(th1 = 0.1, th2 = 0.05) # Initial conditions

nDays <- 5
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate
parms <- c(w1=1/24, w2=0.85/24, C1=0/24, C2=0.16/24) # Parameters

# Simulate
sol <- strogatz(ts, y0, parms)

## ----Display output, echo=FALSE, results='asis'--------------------------
knitr::kable(head(sol, 5))

## ----Plotting, include = TRUE, warning = FALSE, fig.show='hold'----------
plot(sol$time/24, sol$asleep, type = 'line', xlab = 'Time (d)', ylab = 'Asleep')
rasterPlot(sol)

## ----Entrained vs not, echo = FALSE, include = TRUE, fig.show='hold'-----
## Problem setting
y0 <- c(th1 = 0.1, th2 = 0.05) # Initial conditions

nDays <- 60
ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

# Parameters
# With this settings, the bifurcation happens at w2 = 0.84/24 h^-1
parms_entrained <- c(w1=1/24, w2=0.845/24, C1=0/24, C2=0.16/24)
parms_non_entrained <- c(w1=1/24, w2=0.835/24, C1=0/24, C2=0.16/24)

# Simulate
sol_entrained <- strogatz(ts, y0, parms_entrained)
sol_non_entrained <- strogatz(ts, y0, parms_non_entrained)

# Plot
rasterPlot(sol_entrained)
title('Entrained')
rasterPlot(sol_non_entrained)
title('Not entrained')


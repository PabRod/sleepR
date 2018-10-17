dStrogatz <- function(time, y, parms = strogatz_default_parms()) {

  with(as.list(c(y, parms)), {
    dth1 <- w1 - C1*cos(2*pi*(th1-th2))
    dth2 <- w2 + C2*cos(2*pi*(th1-th2))

    return(list(c(dth1, dth2)))
  })

}

strogatz <- function(ts, y0, parms = strogatz_default_parms()) {

  # Load required libraries
  library(deSolve)
  library(dplyr)

  # Solve
  ys <- ode(y = y0,
            func = dStrogatz,
            times = ts,
            parms = parms)

  # Transform into data frame
  ys <- as.data.frame(ys)

  # Implement sleep-awake criterion
  ys <- mutate(ys, th1 = th1%%1, th2 = th2%%1)
  ys <- mutate(ys, asleep = (th2 >= 0.0) & (th2 <= 1/3))

}

strogatz_default_parms <- function() {

  parms = c(w1 = 1/24, # h^-1
            w2 = 0.84/24, # h^-1
            C1 = 0/24, # h^-1
            C2 = 0.16/24) # h^-1

}

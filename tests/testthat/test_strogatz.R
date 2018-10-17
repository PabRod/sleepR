context('Model by Strogatz')

test_that('Defaults',
          {
            # Create a simple input
            t <- 0
            y <- c(th1 = 0, th2 = 0)

            # Calculate flow with default parameters
            dy <- dStrogatz(t, y)

            # Calculate expected result
            parms <- strogatz_default_parms()
            dy_expected <- c(as.numeric(parms['w1'] - parms['C1']), as.numeric(parms['w2'] + parms['C2']))

            # Compare
            expect_equal(dy[[1]][1], dy_expected[1])
            expect_equal(dy[[1]][2], dy_expected[2])
          }
)

test_that('Entrained',
          {
            # Problem posing
            nDays <- 60
            y0 <- c(th1 = 0, th2 = 0.1)
            ts <- seq(0, nDays*24, length.out = nDays*24*20)

            parms <- strogatz_default_parms()
            parms['w2'] <- 0.95/24 # Non default frequency (easier to entrain)

            # Problem solving
            sol <- strogatz(ts, y0, parms = parms)

            # Measured difference
            psi <- sol$th1 - sol$th2
            psi_sync <- tail(psi, n = 1)

            # Auxiliary variables
            omega <- as.numeric(parms['w1'] - parms['w2'])
            C <- as.numeric(parms['C1'] + parms['C2'])

            # Expected difference
            psi_sync_expected <- -1/(2*pi)*acos(omega/C)

            # Compare (don't forget our system has periodicity 1)
            expect_equal(psi_sync%%1, psi_sync_expected%%1)

          }
)

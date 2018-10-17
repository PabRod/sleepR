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

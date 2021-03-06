context('Model by Kronauer et al.')

test_that('Parameters',
          {
            # Load the default parameters
            parms <- kronauer_default_parms()

            # Check
            class_expected <- 'numeric'
            names_expected <- c('mu', 'taux', 'q', 'm', 'c', 'n')

            expect_true(class(parms) == class_expected)
            expect_true(length(parms) == length(names_expected))
            expect_true(all(names(parms) == names_expected))
          }
)

test_that('Flow without light',
          {
            # Create a simple input
            t <- 0
            x_test <- 1
            xc_test <- 2

            y <- c(x = x_test, xc = xc_test)

            I_test <- 0
            I <- function(t) { I_test }

            # Calculate flow with default parameters
            dy <- dKronauer(t, y, I)

            # Calculate expected result
            parms <- kronauer_default_parms()
            B <- (1 - parms['m']*x_test) * parms['c'] * I_test^(1/3)
            dy_expected <- c(as.numeric(pi/12*(xc_test + parms['mu']*(x_test - 4/3*x_test^3) + B)),
                             as.numeric(pi/12*(parms['q']*B*xc_test-(24/parms['taux'])^2*x_test)))

            # Compare
            expect_equal(dy[[1]][1], dy_expected[1])
            expect_equal(dy[[1]][2], dy_expected[2])

          }
)

test_that('Flow with light',
          {
            # Create a simple input
            t <- 0
            x_test <- 1
            xc_test <- 2

            y <- c(x = x_test, xc = xc_test)

            I_test <- 9500
            I <- function(t) { I_test }

            # Calculate flow with default parameters
            dy <- dKronauer(t, y, I)

            # Calculate expected result
            parms <- kronauer_default_parms()
            B <- (1 - parms['m']*x_test) * parms['c'] * I_test^(1/3)
            dy_expected <- c(as.numeric(pi/12*(xc_test + parms['mu']*(x_test - 4/3*x_test^3) + B)),
                             as.numeric(pi/12*(parms['q']*B*xc_test-(24/parms['taux'])^2*x_test)))

            # Compare
            expect_equal(dy[[1]][1], dy_expected[1])
            expect_equal(dy[[1]][2], dy_expected[2])

          }
)

test_that('Basic test',
          {
            ## Problem setting
            y0 <- c(x = 1, xc = 0) # Initial conditions

            nDays <- 12
            ts <- seq(0, nDays*24, length.out=nDays*24*20) # Times to simulate

            ## Simple illumination condition
            I <- function(t) { (t%%24 > 8)*10000 + # Daylight after 8:00 AM
                !(t%%24 > 8)*150} # Night at 00:00 AM

            ## Simulate
            sol <- kronauer(ts, y0, I)

            expected_colnames <- c('time', 'x', 'xc')
            expected_class <- 'data.frame'

            expect_true(all(colnames(sol) == expected_colnames))
            expect_true(class(sol) == expected_class)
          }
)

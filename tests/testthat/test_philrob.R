context('Model by Phillips and Robinson')

test_that('Parameters',
          {
            # Load the default parameters
            parms <- philrob_default_parms()

            # Check
            class_expected <- 'numeric'
            names_expected <- c('Qmax', 'theta', 'sigma', 'vmaSa', 'vvm', 'vmv',
                                'vvc', 'vvh', 'Xi', 'mu', 'taum', 'tauv', 'w',
                                'D', 'Da', 'alpha')

            expect_true(class(parms) == class_expected)
            expect_true(length(parms) == length(names_expected))
            expect_true(all(names(parms) == names_expected))
          }
)

test_that('Stable solution',
          {
            # Create a simple input
            nDays <- 8
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -13, Vm = 1, H = 10)

            # Simulate with parameters allowing stable solution
            parms <- philrob_default_parms()
            parms['vmaSa'] <- 0.0
            parms['vvc'] <- 0.0
            sol <- philrob(ts, y0, parms)

            # Extract last point
            ts_last <- tail(sol$time, 1)
            Vv_last <- tail(sol$Vv, 1)
            Vm_last <- tail(sol$Vm, 1)
            H_last <- tail(sol$H, 1)

            y_last <- c(Vv = Vv_last, Vm = Vm_last, H = H_last)

            # Check that it is an equilibrium
            dy_last <- dPhilrob(ts_last, y_last, parms)
            dy_last <- dy_last[[1]]

            dy_tol <- 1e-3*3600 # h^-1
            expect_true(all(abs(dy_last) < dy_tol))

          }
)

test_that('Only V_m and H stable',
          {
            # Create a simple input
            nDays <- 8
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -13, Vm = 1, H = 10)

            # Simulate with parameters allowing the desired behavior
            parms <- philrob_default_parms()
            parms['vmv'] <- 0.0
            parms['vvm'] <- 0.0
            sol <- philrob(ts, y0, parms)

            # Extract last point
            ts_last <- tail(sol$time, 1)
            Vv_last <- tail(sol$Vv, 1)
            Vm_last <- tail(sol$Vm, 1)
            H_last <- tail(sol$H, 1)

            y_last <- c(Vv = Vv_last, Vm = Vm_last, H = H_last)

            # Check that it is an equilibrium
            dy_last <- dPhilrob(ts_last, y_last, parms)
            dy_last <- dy_last[[1]]

            y_tol <- 1e-3
            dy_tol <- 1e-3*3600 # h^-1

            # The solution for Vv should keep oscillating
            Vv_last_day <- tail(sol$Vv, 24 * 20)
            expect_false(abs(max(Vv_last_day) - min(Vv_last_day)) < y_tol )

            # The time series for Vm and H should be stable
            expect_true(all(abs(dy_last[2:3]) < dy_tol))

          }
)


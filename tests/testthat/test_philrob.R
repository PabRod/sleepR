context('Model by Phillips and Robinson')

test_that('Parameters',
          {
            # Load the default parameters
            parms <- philrob_default_parms()

            # Check
            class_expected <- 'numeric'
            names_expected <- c('Qmax', 'theta', 'sigma',
                                'vmaSa', 'vvm', 'vmv',
                                'vvc', 'vvh', 'Xi',
                                'mu', 'taum', 'tauv',
                                'w', 'alpha')

            expect_true(class(parms) == class_expected)
            expect_true(length(parms) == length(names_expected))
            expect_true(all(names(parms) == names_expected))
          }
)

test_that('Saturation function',
          {
            parms <- philrob_default_parms()
            S <- function(V) saturating_function(V, parms)

            numerical_inf <- 1e10

            with(as.list(parms), {
              expect_equal(S(numerical_inf), Qmax) # Horizontal asymptote
              expect_equal(S(theta), Qmax/2) # Half saturation value
              expect_equal(S(0), Qmax/(1 + exp(theta/sigma))) # Intercept
            })
          }
)

test_that('Derivative of saturation function',
          {
            parms <- philrob_default_parms()
            dS <- function(V) dSaturating_function(V, parms)

            numerical_inf <- 1e10

            with(as.list(parms), {
              expect_equal(dS(numerical_inf), 0) # Horizontal asymptote
              expect_equal(dS(theta), Qmax/(4*sigma)) # Half saturation value
            })
          }
)

test_that('Jacobian',
          {
            # Constant case
            parms <- philrob_default_parms()
            parms['vvm'] <- 0
            parms['vvh'] <- 0
            parms['vmv'] <- 0
            parms['mu'] <- 0

            y <- c(Vv = 2, Vm = 3, H = 4)
            J <- philrob_jacobian(y, parms)

            with(as.list(parms),
                 {
                   expected <- diag(-1,3) * c(1/tauv, 1/taum, 1/Xi)
                   expect_true(all(J == expected))
                 }
                 )

            # Comparison with analytical derivation
            # parms <- philrob_default_parms()
            #
            # y <- c(Vv = 2, Vm = 3, H = 4)
            # dS <- function(x) dSaturating_function(y, parms)
            # J <- philrob_jacobian(y, parms)
            #
            # vec <- c(-1/parms['tauv']            , -parms['vvm']/parms['tauv']*dS(y['Vm']), parms['vvh']/parms['tauv'],
            #          -parms['vmv']/parms['taum']*dS(y['Vv']) , -1/parms['taum']           , 0 ,
            #          0                  , parms['mu']/parms['Xi']*dS(y['Vm'])    , -1/parms['Xi'])
            #
            # expected <- matrix(vec, nrow = 3, ncol = 3, byrow = TRUE)
            # expect_equal(J, expected)

            # with(as.list(parms), {
            #   vec <- c(-1/tauv            , -vvm/tauv * dS(Vm), vvh/tauv,
            #            -vmv/taum * dS(Vv) , -1/taum           , 0 ,
            #            0                  , mu/Xi * dS(Vm)    , -1/Xi)
            #
            #   expected <- matrix(vec, nrow = 3, ncol = 3, byrow = TRUE)
            #   expect_equal(J, expected)
            # })

          }
)

test_that('Default forcing function',
          {
            # Load the default forcing
            C <- philrob_default_forcing()

            ts <- seq(0, 24, length.out = 100)
            fs <- C(ts)

            # Check that it is bounded between 0 and 1
            tol <- 0.002
            expect_equal(max(fs), 1, tolerance = tol)
            expect_equal(min(fs), 0, tolerance = tol)

          }
)

test_that('Custom forcing',
          {
            # Create a simple input
            nDays <- 3
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -13, Vm = 1, H = 10)

            # Load parameters
            parms <- philrob_default_parms()

            # Case 1: no forcing through forcing equal to zero
            C <- function(t) { 0.0 }
            sol_custom_forcing <- philrob(ts, y0, parms, C = C)

            # Case 2: no forcing through forcing coupling equal to zero
            parms['vvc'] <- 0
            sol_custom_pars <- philrob(ts, y0, parms)

            # Both examples should be equivalent
            expect_equal(sol_custom_pars, sol_custom_forcing)

          }
)

test_that('Stabilization run',
          {
            # Create a simple input
            nDays <- 3
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -13, Vm = 1, H = 10)

            # Simulate with parameters allowing stable solution
            parms <- philrob_default_parms()
            parms['vmaSa'] <- 0.0
            parms['vvc'] <- 0.0

            # Simulate three days without stabilization
            sol <- philrob(ts, y0, parms)

            # Simulate one day after two days of stabilization
            nDays_stabilized <- 1
            ts_stabilized <- seq(0, nDays_stabilized * 24, length.out = nDays_stabilized * 24 * 20)
            sol_stabilized <- philrob(ts_stabilized, y0, parms, tStabil = 2 * 24)

            # Compare the last day for both simulations
            sol_last_day <- tail(sol, nDays_stabilized * 24 * 20)

            expect_equal(sol_last_day$Vv, sol_stabilized$Vv, tolerance = 1e-4)
            expect_equal(sol_last_day$Vm, sol_stabilized$Vm, tolerance = 1e-4)
            expect_equal(sol_last_day$H, sol_stabilized$H, tolerance = 1e-4)
          }
          )

test_that('Stable solution',
          {
            # Create a simple input
            nDays <- 3
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

test_that('Paper 2007',
          {
            # Simulate the situation in the paper
            nDays <- 8
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -12.6404, Vm = 0.8997, H = 12.5731) # This point is already close to the attractor
            sol <- philrob(ts, y0, method = 'lsode')

            # Extract relevant results
            Vv <- sol$Vv
            Vm <- sol$Vv
            H <- sol$H

            # Approximate limits taken from Figure 5
            expect_true((min(Vv) > -13) & (min(Vv) < -10))
            expect_true((max(Vv) > 1) & (max(Vv) < 2))
            expect_true((min(Vm) > -13) & (min(Vm) < -10))
            expect_true((max(Vm) > 0) & (max(Vm) < 3))
            expect_true((min(H) > 7) & (min(H) < 9))
            expect_true((max(H) > 14) & (max(H) < 15))

            # Check the asleep status
            expect_true(class(sol$asleep) == 'logical')
            expect_false(all(sol$asleep[1:3])) # This simulation begins with an awake subject
            expect_true((mean(sol$asleep) > 0.3) & (mean(sol$asleep) < 0.4)) # Average sleeping time should be around 0.3

          }

)

test_that('Philrob plot',
          {
            # Calculate a solution
            y0 <- c(Vv = 1, Vm = 1, H = 1) # Initial conditions
            nDays <- 6 # Times
            ts <- seq(0, nDays*24, length.out=nDays*24*20)
            sol <- philrob(ts, y0) # Simulate

            philrobPlot(sol)

            # Just check that the code doesn't crash
            expect_true(TRUE)
          }

)


context('Model by Borbely et al.')

test_that('Parameters',
          {
            # Load the default parameters
            parms <- borbely_default_parms()

            # Check
            class_expected <- 'numeric'
            names_expected <- c('Xis', 'Xiw', 'mu', 'a', 'w', 'alpha', 'Hu0', 'Hl0')

            expect_true(class(parms) == class_expected)
            expect_true(length(parms) == length(names_expected))
            expect_true(all(names(parms) == names_expected))
          }
          )

test_that('Figure 1b',
          {
            # Problem posing
            nDays <- 3
            awake0 <- FALSE
            y0 <- 0.4
            ts <- seq(0, nDays*24, length.out = nDays*24*20)

            parms <- borbely_default_parms()

            # Problem solving
            sol <- borbely(ts, y0, awake0, parms) # Stabilization run
            y0_att <- tail(sol$H, 1)
            awake0_att <- tail(sol$awake, 1)
            sol_att <- borbely(ts, y0_att, awake0_att, parms) # Run in the attractor

            # Uncomment to generate figure 1.b
            # borbelyPlot(sol_att)

            # Approximate limits taken from Figure 1.b
            sol_1d <- sol_att[which(sol_att$time <= 24), ] # Keep the first day only
            sol_max <- sol_1d[which(sol_1d$H == max(sol_1d$H)), ]
            sol_min <- sol_1d[which(sol_1d$H == min(sol_1d$H)), ]

            y_max <- sol_max$H
            t_max <- sol_max$time/24
            y_min <- sol_min$H
            t_min <- sol_min$time/24

            # Check that limits are within the boundaries
            expect_true((y_max > 0.55) & (y_max < 0.65))
            expect_true((y_min > 0.05) & (y_min < 0.15))
            expect_true((t_max > 0.45) & (t_max < 0.55))
            expect_true((t_min > 0.8) & (t_min < 0.85))

          }
)

test_that('Figure 1c',
          {
            # Problem posing
            nDays <- 15
            awake0 <- FALSE
            y0 <- 0.4
            ts <- seq(0, nDays*24, length.out = nDays*24*20)

            # Set parameters
            parms <- borbely_default_parms()
            parms['Hu0'] <- 0.85
            parms['alpha'] <- -pi/2

            # Problem solving
            sol <- borbely(ts, y0, awake0, parms) # Stabilization run
            y0_att <- tail(sol$H, 1)
            awake0_att <- tail(sol$awake, 1)
            sol_att <- borbely(ts, y0_att, awake0_att, parms) # Run in the attractor

            # Approximate limits taken from Figure 1.c
            sol_2d <- sol_att[which(sol_att$time <= 48), ] # Keep the first day only
            sol_max <- sol_2d[which(sol_2d$H == max(sol_2d$H)), ]
            sol_min <- sol_2d[which(sol_2d$H == min(sol_2d$H)), ]

            y_max <- sol_max$H
            t_max <- sol_max$time/24
            y_min <- sol_min$H
            t_min <- sol_min$time/24

            # Check that limits are within the boundaries
            expect_true((y_max > 0.85) & (y_max < 0.95))
            expect_true((y_min > 0.05) & (y_min < 0.15))
            expect_true((t_max > 0.15) & (t_max < 0.25))
            expect_true((t_min > 0.55) & (t_min < 0.65))

          }
)

test_that('Borbely plot',
          {
            # Calculate a solution
            y0 <- 0.4 # Initial conditions
            nDays <- 6 # Times
            ts <- seq(0, nDays*24, length.out=nDays*24*20)
            sol <- borbely(ts, y0) # Simulate

            borbelyPlot(sol)

            # Just check that the code doesn't crash
            expect_true(TRUE)
          }

)

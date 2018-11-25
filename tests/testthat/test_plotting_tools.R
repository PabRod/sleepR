context('Plotting tools')

test_that('Raster vector',
          {
            # Calculate a solution
            y0 <- c(th1 = 0.1, th2 = 0.05) # Initial conditions
            nDays <- 60.2 # Times
            ts <- seq(0, nDays*24, length.out=nDays*24*20)
            sol <- strogatz(ts, y0) # Simulate

            # Get raster sleep vector
            dailySamples <- 480
            asleep_raster <- reshape_as_raster(sol, dailySamples = dailySamples)

            # Check
            class_expected <- 'matrix'
            dim_expected <- c(floor(nDays-1), 2 * dailySamples)

            expect_true(class(asleep_raster) == class_expected)
            expect_true(all(dim(asleep_raster) == dim_expected))

            expect_equal(max(asleep_raster, na.rm = TRUE), 1)
            expect_equal(min(asleep_raster, na.rm = TRUE), 0)
          }

)

test_that('Raster plot',
          {
            # Calculate a solution
            y0 <- c(th1 = 0.1, th2 = 0.05) # Initial conditions
            nDays <- 6 # Times
            ts <- seq(0, nDays*24, length.out=nDays*24*20)
            sol <- strogatz(ts, y0) # Simulate

            # Generate the raster plot
            dailySamples <- 480
            rasterPlot(sol, dailySamples = dailySamples)

            # Just check that the code doesn't crash
            expect_true(TRUE)
          }
)

test_that('Lissajous',
          {
            # Calculate a solution
            nDays <- 8
            ts <- seq(0, nDays * 24, length.out = nDays * 24 * 20)
            y0 <- c(Vv = -12.6404, Vm = 0.8997, H = 12.5731) # This point is already close to the attractor
            sol <- philrob(ts, y0, method = 'lsode')

            # Create the figure
            lissajous_figure(ts, sol$H)

            # Just check that the code doesn't crash
            expect_true(TRUE)
          }
)

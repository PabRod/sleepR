context('Raster plots')

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
            dim_expected <- c(floor(nDays-1), dailySamples)

            expect_true(class(asleep_raster) == class_expected)
            expect_true(all(dim(asleep_raster) == dim_expected))

            expect_equal(max(asleep_raster, na.rm = TRUE), 1)
            expect_equal(min(asleep_raster, na.rm = TRUE), 0)
          }

)

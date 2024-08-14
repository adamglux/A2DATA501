#test model:
data1 <- mtcars
model1 <- lm(mpg ~ wt + hp, data = data1)

### generated data
data2 <- randu
model2 <- glm(z ~ x + y, data = data2)
data3 <- data.frame(a = rep(c(NA,runif(1)), 50), b = rnorm(100), c=(runif(100)))
model3 <- lm(b ~ a + c, data = data3)
data4 <- data.frame(y = c(runif(10),Inf), x = 1:11, z = 21:31)
model4 <- lm(z ~ x, data = data4)
data5 <- data.frame(x = c(1:4), y = c(rnorm(4)))
model5 <- lm(y ~ x, data = data5)
data6 <- data.frame(y = c(1:10))
model6 <- lm(y ~ 1, data = data6)
matx7 <- matrix(1:30, nrow = 5)
data7 <- data.frame(matx7)
model7 <- lm(X1 ~ X2, data = data7)



## invalid inputs
test_that("function influenceR gives helpful errors",
          {
            expect_error(influenceR(model1, matrix(data1)),
                         "Input data must be a data.frame.")

            expect_error(influenceR(mtcars, mtcars),
                         "Input model must be an lm or glm object.")

            expect_error(influenceR(model6, data6),
                         "Input data must have at least 2 columns.")

            expect_error(influenceR(model5, data5),
                         "Dataset is too small")

            expect_error(influenceR(model7, data7),
                         "Dataset is the wrong dimension")

            expect_error(influenceR(model4, data4),
                         "Input data contains infinite values.")

            expect_error(influenceR(model3, data3),
                         "Input data contains more than 5% NA values.")

            expect_error(influenceR(model1, data1, plot = "phil collins"),
                         "Input for plot must be either: none, all, cooks, DFFITS, or Hadi")

            output <- influenceR(model1, data1, plot = "none")
            expect_true(is.data.frame(output))
            expect_equal(names(output), c("Cook's D", "DFFITs", "Hadi's Influence"))

            output <- influenceR(model1, data1, plot = "all")
            expect_true(is.list(output))
            expect_equal(names(output$influence_measures), c("Cook's D", "DFFITs", "Hadi's Influence"))


            # Assuming pre-calculated values for a specific row for testing
            expect_equal(round(output[[1]][1, "Cook's D"], 3), 0.016)
            expect_equal(round(output[[1]][1, "DFFITs"], 3), -0.223)
            expect_equal(round(output[[1]][1, "Hadi's Influence"], 3), 1.076)

          })

## valid inputs
test_that("function influenceR works with dataframes, lm/glm objects, and plots", {

  pdf(NULL)
  on.exit(dev.off())

  expect_silent(influenceR(model1, mtcars, plot = "none"))
  expect_silent(influenceR(model1, mtcars, plot = "cooks"))
  expect_silent(influenceR(model1, mtcars, plot = "DFFITS"))
  expect_silent(influenceR(model1, mtcars, plot = "Hadi"))
  expect_silent(influenceR(model1, mtcars, plot = "all"))

})

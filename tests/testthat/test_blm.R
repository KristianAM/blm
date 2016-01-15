context('testing blm function')

test_that('blm correctly estimates weights', {
  x <- runif(500000)
  w1 <- runif(1)
  w0 <- runif(1)
  y <- rnorm(500000, mean = x * w1 + w0)
  df <- data.frame(x,y)
  blmodel <- blm(y ~ x, Data = df, beta = 1, prior = NULL)
  post <- blmodel$posterior$mean
  expect_equal(c(post), c(w0,w1), tolerance = 1e-2)
})

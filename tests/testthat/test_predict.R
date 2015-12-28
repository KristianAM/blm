context('testing predict function')
test_that( 'blm correctly predicts y values', {
  x <- runif(500000)
  w1 <- runif(1)
  w0 <- runif(1)
  y <- rnorm(500000, mean = x * w1 + w0)
  df <- data.frame(x,y)
  model <- blm(y ~ x, df = df, beta = 1)
  newx <- data.frame(runif(20))
  colnames(newx) <- 'x'
  predictedy <- predict.blm(blm = model, x =  newx)
  expectedy <- c(newx * w1 + w0)[[1]]
  expect_equal(c(predictedy), expectedy, tolerance = 1e-2)
})

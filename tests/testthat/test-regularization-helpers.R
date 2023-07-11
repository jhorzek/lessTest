test_that("regularization helper work", {
  
  # The following is adapted from ?lavaan::sem
  library(lessSEM)
  model <- ' 
  # latent variable definitions
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
  dem60 ~ ind60
  dem65 ~ ind60 + dem60

  # residual correlations
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8
'
  
  fit <- sem(model, data = PoliticalDemocracy)
  
  l <- loadings(fit)
  
  testthat::expect_true(
    all(sort(l) == sort(c("ind60=~x2","ind60=~x3", "a", "b", "c")))
  )
  
  r <- regressions(fit)
  
  testthat::expect_true(
    all(sort(r) == sort(c("dem60~ind60","dem65~ind60", "dem65~dem60")))
  )
  
  v <- variances(fit)
  
  testthat::expect_true(
    all(sort(v) == sort(c("x1~~x1","x2~~x2", "x3~~x3",
                          "y1~~y1", "y2~~y2", "y3~~y3",
                          "y4~~y4", "y5~~y5", "y6~~y6",
                          "y7~~y7", "y8~~y8",
                          "ind60~~ind60", "dem60~~dem60",
                          "dem65~~dem65")))
  )
  
  cv <- covariances(fit)
  
  testthat::expect_true(
    all(sort(cv) == sort(c("y1~~y5", 
                           "y2~~y4",
                           "y2~~y6", 
                           "y3~~y7",
                           "y4~~y8",
                           "y6~~y8")))
  )
  
})

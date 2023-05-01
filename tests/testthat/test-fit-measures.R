test_that("fit measures work", {
  library(lessSEM)
  library(testthat)
  
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
  lavaanModel = sem(model, PoliticalDemocracy)
  
  # Regularization:
  lambdas <- 0
  
  lsem <- lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("a"),
    # in case of lasso and adaptive lasso, we can specify the number of lambda
    # values to use. lessSEM will automatically find lambda_max and fit
    # models for nLambda values between 0 and lambda_max. For the other
    # penalty functions, lambdas must be specified explicitly
    lambdas = lambdas)
  
  lfit <- fitIndices(lsem)
  
  testthat::expect_true(abs(lfit$AIC[1] - AIC(lavaanModel)) < 1e-5)
  testthat::expect_true(abs(lfit$BIC[1] - BIC(lavaanModel)) < 1e-5)
  testthat::expect_true(abs(lfit$rmsea[1] - fitMeasures(lavaanModel, "rmsea")) < 1e-5)
  testthat::expect_true(abs(lfit$chisq[1] - fitMeasures(lavaanModel, "chisq")) < 1e-5)
  
})

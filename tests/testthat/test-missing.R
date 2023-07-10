test_that("missing data works", {
  library(lessSEM)
  set.seed(123)
  # Identical to regsem, lessSEM builds on the lavaan
  # package for model specification. The first step
  # therefore is to implement the model in lavaan.
  
  dataset <- simulateExampleData(N = 200, 
                                 percentMissing = 20)
  
  # testing some extreme missingness patterns
 
  dataset[3,-1] <- NA
  dataset[4,-2] <- NA
  dataset[1,] <- NA
  dataset[2,] <- NA
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  lavaanModel <- lavaan::sem(lavaanSyntax,
                             data = dataset,
                             meanstructure = TRUE,
                             std.lv = TRUE,
                             missing = "ml")
  
  # Regularization:
  lsem <- lasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    # in case of lasso and adaptive lasso, we can specify the number of lambda
    # values to use. lessSEM will automatically find lambda_max and fit
    # models for nLambda values between 0 and lambda_max. For the other
    # penalty functions, lambdas must be specified explicitly
    nLambdas = 50)
  
  testthat::expect_equal(abs(-2*logLik(lavaanModel) - lsem@fits$m2LL[50]) < 1e-5,
                         TRUE)
  
  # use the plot-function to plot the regularized parameters:
  plot(lsem)
  
  # the coefficients can be accessed with:
  coef(lsem)
  # if you are only interested in the estimates and not the tuning parameters, use
  coef(lsem)@estimates
  # or
  estimates(lsem)
  
  # elements of lsem can be accessed with the @ operator:
  lsem@parameters[1,]
  
  unregularizedAt <- which(lsem@fits$lambda == 0)
  
  # fit Measures:
  testthat::expect_equal(abs(fitIndices(lsem)$AIC[unregularizedAt] - AIC(lavaanModel)) < 1e-4, 
                         TRUE)
  
  testthat::expect_equal(abs(fitIndices(lsem)$BIC[unregularizedAt] - BIC(lavaanModel)) < 1e-4, 
                         TRUE)
  # The best parameters can also be extracted with:
  coef(lsem, criterion = "AIC")
  # or
  estimates(lsem, criterion = "AIC") 
  
  lsem <- lavaanModel |>
    # create template for regularized model with mixed penalty:
    mixedPenalty() |>
    # add lasso penalty on loadings l6 - l10:
    addLasso(regularized = paste0("l", 6:10), 
             lambdas = seq(0,1,length.out = 4)) |>
    # add scad penalty on loadings l11 - l15:
    addScad(regularized = paste0("l", 11:15), 
            lambdas = seq(0,1,length.out = 3),
            thetas = 3.1) |>
    # fit the model:
    fit()
  
  testthat::expect_equal(abs(-2*logLik(lavaanModel) - lsem@fits$m2LL[1]) < 1e-5,
                         TRUE)
  
  # the coefficients can be accessed with:
  coef(lsem)
  # if you are only interested in the estimates and not the tuning parameters, use
  coef(lsem)@estimates
  # or
  estimates(lsem)
  
  # elements of lsem can be accessed with the @ operator:
  lsem@parameters[1,]
  
  # fit Measures:
  testthat::expect_equal(abs(fitIndices(lsem)$AIC[1] - AIC(lavaanModel)) < 1e-4, 
                         TRUE)
  
  testthat::expect_equal(abs(fitIndices(lsem)$BIC[1] - BIC(lavaanModel)) < 1e-4, 
                         TRUE)
  
  # The best parameters can also be extracted with:
  coef(lsem, criterion = "AIC")
  # or
  estimates(lsem, criterion = "AIC") 
  
  # Regularization:
  lsem <- cvLasso(
    # pass the fitted lavaan model
    lavaanModel = lavaanModel,
    # names of the regularized parameters:
    regularized = paste0("l", 6:15),
    # in case of lasso and adaptive lasso, we can specify the number of lambda
    # values to use. lessSEM will automatically find lambda_max and fit
    # models for nLambda values between 0 and lambda_max. For the other
    # penalty functions, lambdas must be specified explicitly
    lambdas = seq(0,1,.1))
  
  coef(lsem)
  estimates(lsem)
  
  testthat::expect_true(!anyNA(lsem@cvfits))
  
})

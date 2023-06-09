test_that("testing elasticNet-ridge-c", {
  library(lslx)
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  N <- 50
  l1 <- 1; l2 <- .2; l3 <- 0;
  v1 <- .2; v2 <- .8; v3 <- 1
  
  f <- matrix(stats::rnorm(N, 0, 1), ncol = 1)
  L <- matrix(c(rep(l1,5), rep(l2,5), rep(l3,15)), nrow = 1)
  y <- matrix(NA, nrow = N, ncol = ncol(L))
  
  covs <- c(rep(v1,5), rep(v2,5), rep(v3,15))
  
  for(i in 1:N){
    y[i,] <- L*f[i,] +  mvtnorm::rmvnorm(1, sigma = diag(covs))
  }
  
  yNames <- paste0("y", 1:ncol(y))
  colnames(y) <- yNames
  
  modelSyntax <- paste0('f =~ 1*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "))
  modelFit = cfa(modelSyntax, y, meanstructure = TRUE)
  
  # fit model using lslx
  lslxModelSyntax <- paste0('fix(1)*', yNames[1], ' + ', paste0(yNames[2:length(yNames)], collapse = " + "), " <=: f") 
  fitLslx <- lslx$new(model = lslxModelSyntax,
                      sample_cov = stats::cov(y),
                      sample_size = nrow(y)
  )
  
  fitLslx$penalize_coefficient(name = paste0("y", 6:ncol(y)," <- f"))
  
  lambdas <- seq(0,50,10)
  fitLslx$fit(penalty_method = "ridge",lambda_grid = lambdas, loss = "ml")
  
  # extract fits
  lslxParameter <- matrix(NA, 
                          nrow = length(lambdas), 
                          ncol = length(fitLslx$extract_coefficient(lambda = 0)))
  colnames(lslxParameter) <- names(fitLslx$extract_coefficient(lambda = 0))
  
  for(l in 1:length(lambdas)){
    pars <- fitLslx$extract_coefficient(lambda = lambdas[l])
    lslxParameter[l,names(pars)] <- pars
  }
  regularized <- paste0("y", 6:ncol(y),"<-f/g") 
  
  # replicate with regularizedSEM
  
  rsemIsta <- ridge(lavaanModel = modelFit, 
                    regularized = paste0("f=~y",6:ncol(y)),
                    lambdas = lambdas,
                    method = "ista",
                    control = controlIsta()
  )
  
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] - lslxParameter[,regularized]) < .002), TRUE)
  #plot(rsemIsta)
  coef(rsemIsta)
  
  rsemGlmnet <- ridge(lavaanModel = modelFit, 
                      regularized = paste0("f=~y",6:ncol(y)),
                      lambdas = lambdas,
                      method = "glmnet",
                      control = controlGlmnet()
  )
  
  testthat::expect_equal(all(abs(rsemGlmnet@parameters[,rsemGlmnet@regularized] - lslxParameter[,regularized]) < .002), TRUE)
  #plot(rsemGlmnet)
  coef(rsemGlmnet)
  
  rsemBfgs <- ridgeBfgs(lavaanModel = modelFit,
                        regularized = paste0("f=~y",6:ncol(y)),
                        lambdas = lambdas
  )
  testthat::expect_equal(all(abs(rsemIsta@parameters[,rsemIsta@regularized] -
                                   rsemBfgs@parameters[,rsemBfgs@regularized]) < .002), TRUE)
  
})

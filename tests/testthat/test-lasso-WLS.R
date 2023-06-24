test_that("testing elasticNet-lasso-WLS", {
  library(lslx)
  library(lavaan)
  library(lessSEM)
  set.seed(123)
  N <- 1000
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
  
  for(estimator in c("wls", "uls", "dwls")){
    
    modelFit = cfa(modelSyntax, y, 
                   meanstructure = FALSE,
                   estimator = estimator)
    
    # fit model using lslx
    lslxModelSyntax <- paste0(paste0('fix(1)*', yNames[1], ' + ', paste0(yNames[2:5], collapse = " + "), " <=: f"),"\n",
                              paste0(paste0(yNames[6:length(yNames)], collapse = " + "), " <~: f"),"\n",
                              paste0(yNames, collapse = " + "), " <= 1"
    )
    fitLslx <- lslx$new(model = lslxModelSyntax,
                        sample_cov = stats::cov(y),
                        sample_size = nrow(y)
    )
    
    #fitLslx$penalize_coefficient(name = paste0("y", 6:ncol(y)," <- f"))
    
    lambdas <- seq(0,.3,length.out = 15)
    fitLslx$fit(penalty_method = "lasso",lambda_grid = lambdas, loss = estimator)
    
    # extract fits
    lslxParameter <- matrix(NA, 
                            nrow = length(lambdas), 
                            ncol = length(fitLslx$extract_coefficient(lambda = 0)))
    colnames(lslxParameter) <- names(fitLslx$extract_coefficient(lambda = 0))
    
    regWls <- rep(NA, length(lambdas))
    
    for(l in 1:length(lambdas)){
      pars <- fitLslx$extract_coefficient(lambda = lambdas[l])
      lslxParameter[l,names(pars)] <- pars
      regWls[l] <- fitLslx$extract_numerical_condition(lambda = lambdas[l])["objective_value"]
    }
    
    if(abs(regWls[1] - 2*(fitMeasures(modelFit, "fmin")) * (N/(N-1))) > .001){
      warning("Different fit for lslx and lavaan when using ", estimator, ".")
      next
    }
    
    # replicate with regularizedSEM
    
    rsemGlmnet <- lasso(lavaanModel = modelFit, 
                        regularized = paste0("f=~y",6:ncol(y)), 
                        lambdas = lambdas,
                        method = "glmnet",
                        control = controlGlmnet()
    )
    
    testthat::expect_true(all(abs(rsemGlmnet@fits$regObjectiveValue/N -
                                    regWls) < 5e-3))
    
  }
})

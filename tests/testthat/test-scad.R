test_that("testing scad", {
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
  
  #### Regularize ####
  lavaanParameters <- lessSEM::getLavaanParameters(modelFit)
  lambdas <- seq(0,1,length.out = 5)
  thetas <- seq(2.1,3,length.out = 3)
  
  regularizedLavaan <- paste0("f=~y",6:ncol(y))
  
  rsemGlmnet <- scad(lavaanModel = modelFit, 
                   regularized = regularizedLavaan, 
                   lambdas = lambdas, 
                   thetas = thetas, 
                   method = "glmnet",
                   control = controlGlmnet()
  )
  coef(rsemGlmnet)
  coef(rsemGlmnet, criterion = "AIC")
  coef(rsemGlmnet, criterion = "BIC")
  
  rsemIsta <- scad(lavaanModel = modelFit, 
                   regularized = regularizedLavaan, 
                   lambdas = lambdas, 
                   thetas = thetas, 
                   method = "ista",
                   control = controlIsta()
  )
  
  coef(rsemIsta)
  coef(rsemIsta, criterion = "AIC")
  coef(rsemIsta, criterion = "BIC")
  
  testthat::expect_equal(
    all(round(rsemIsta@fits$regM2LL -
                rsemGlmnet@fits$regM2LL
              ,3)==0),TRUE)
  
  testthat::expect_equal(
    all(round(rsemIsta@fits$m2LL[rsemIsta@fits$lambda == 0] -
                -2*logLik(modelFit)
              ,4)==0),TRUE)
  
  # compare to smoothed version
  tuningParameters <- expand.grid(theta = thetas,
                                  lambda = lambdas
  )
  
  penaltyFunctionArguments <- list(
    eps = 1e-10,
    regularizedParameterLabels = paste0("f=~y",6:ncol(y))
  )
  
  #### Now we are ready to optimize! ####
  regsemApprox <- lessSEM:::.regularizeSEMWithCustomPenaltyRsolnp(lavaanModel = modelFit,
                                                                  individualPenaltyFunction = lessSEM:::.smoothScadValue,
                                                                  tuningParameters = tuningParameters,
                                                                  penaltyFunctionArguments = penaltyFunctionArguments)
  
  sel <- regsemApprox@parameters$theta == thetas[1] & 
    regsemApprox@parameters$lambda == lambdas[2]
  pars <- unlist(regsemApprox@parameters[sel,
                                         regsemApprox@parameterLabels]
  )
  weights <- lavaanParameters
  weights[] <- 0
  weights[rsemIsta@inputArguments$weights] <- 1
  lambda <- lambdas[2]
  theta <- thetas[1]
  pen <- 0
  for(p in 1:length(pars)) pen <- pen + lessSEM:::scadPenalty_C(par = pars[p], 
                                                      lambda_p = weights[p]*lambda, 
                                                      theta = theta)
  
  testthat::expect_equal(
    all(round(pen -
                (1/N)*(regsemApprox@fits$regM2LL[sel] - regsemApprox@fits$m2LL[sel])
              ,4)==0),TRUE)
  
  for(i in 1:length(rsemIsta@fits$theta)){
    th <- rsemIsta@fits$theta[i]
    la <- rsemIsta@fits$lambda[i]
    
    testthat::expect_equal(
      round(rsemIsta@fits$regM2LL[rsemIsta@fits$theta == th &
                                    rsemIsta@fits$lambda == la] -
              regsemApprox@fits$regM2LL[regsemApprox@fits$theta == th &
                                          regsemApprox@fits$lambda == la]
      )==0,TRUE)
    testthat::expect_equal(
      all(
        abs(
          rsemIsta@parameters[rsemIsta@fits$theta == th &
                                rsemIsta@fits$lambda == la,rsemIsta@parameterLabels]-
            regsemApprox@parameters[regsemApprox@fits$theta == th &
                                      regsemApprox@fits$lambda == la,rsemIsta@parameterLabels]
        ) < .1
      ),
      TRUE
    )
    
  }
  
})

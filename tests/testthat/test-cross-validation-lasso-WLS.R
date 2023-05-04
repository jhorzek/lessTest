test_that("WLS: testing cross-validation for lasso", {
  library(lessSEM)
  set.seed(123)
  N <- 10000
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
  
  estimators <- c("WLS", "ULS", "GLS", "DWLS")
  
  for(i in 1:length(estimators)){

    modelFit = sem(modelSyntax, 
                   y, 
                   meanstructure = TRUE, 
                   estimator = estimators[i])
    
    regularizedLavaan <- paste0("f=~y",6:ncol(y))
    rsem <- lessSEM::lasso(lavaanModel = modelFit, 
                           regularized = regularizedLavaan,
                           nLambdas = 20)
    
    lambdas <- rsem@fits$lambda
    
    ## Test cross-validation
    
    cv <- cvLasso(lavaanModel = modelFit, 
                  regularized = regularizedLavaan,
                  lambdas = lambdas, 
                  returnSubsetParameters = TRUE,
                  k = 2)
    
    testthat::expect_equal(ncol(cv@subsets), 2)
    
    coef(cv)
    #plot(cv)
    
    subsets <- cv@subsets
    pars <- cv@subsetParameters
    
    parameterLabels <- cv@parameterLabels
    
    # lavaan tends to sort the data differently
    ySorted <- try(lavaan::lavInspect(modelFit, 
                                      "data"))
    
    for(ro in 1:nrow(pars)){
      
      trainSet <- pars$trainSet[ro]
      testSet <- subsets[,trainSet]
      
      modelFitTest = sem(modelSyntax, 
                         ySorted[testSet,], 
                         meanstructure = TRUE, 
                         estimator = estimators[i],
                         do.fit = FALSE)
      
      SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = modelFitTest)
      
      SEM <- lessSEM:::.setParameters(SEM = SEM, 
                                      labels = parameterLabels, 
                                      values = unlist(pars[ro, parameterLabels]),
                                      raw = FALSE)
      
      sel <- cv@cvfits$lambda == pars$lambda[ro]
      if(sum(sel) != 1) stop("Error when selecting cv target")
      cv@cvfits[sel,"cvfit"] <- cv@cvfits[sel,"cvfit"] - SEM$fit()
      
    }
    
    testthat::expect_equal(all(abs(cv@cvfits$cvfit)< 1e-6), TRUE)
    
    # test subset parameters
    
    modelFitTrain = sem(modelSyntax, 
                        ySorted[!testSet,], 
                        meanstructure = TRUE, 
                        estimator = estimators[i])
    subsetPars <- pars[pars$trainSet == pars$trainSet[ro],]
    
    subsetLasso <- lasso(lavaanModel = modelFitTrain, 
                         regularized = regularizedLavaan,
                         lambdas = lambdas)
    
    testthat::expect_equal(all(abs(subsetLasso@parameters - subsetPars[,colnames(subsetLasso@parameters)])< 1e-2), TRUE)
    
    lassoError <- try(cvLasso(lavaanModel = modelFit, 
                              regularized = paste0("f=~y",6:(ncol(y)+1)), 
                              lambdas = lambdas, 
                              k = 3), silent = TRUE)
    testthat::expect_equal(is(lassoError, "try-error"), TRUE)
  }
})
test_that("testing automatic standardization for cross-validation with missing data", {
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
  
  regularizedLavaan <- paste0("l", 6:15)
  
  rsem <- lessSEM::lasso(lavaanModel = lavaanModel, 
                         regularized = regularizedLavaan,
                         nLambdas = 5)
  
  lambdas <- rsem@fits$lambda
  
  # test standardization
  cv <- cvLasso(lavaanModel = lavaanModel, 
                regularized = regularizedLavaan,
                lambdas = lambdas, 
                returnSubsetParameters = TRUE,
                standardize = TRUE,
                k = 2)
  
  subsets <- cv@subsets
  pars <- cv@subsetParameters
  cvfits <- cv@cvfits
  
  parameterLabels <- cv@parameterLabels  
  
  # lavaan tends to sort the data differently
  ySorted <- try(lavaan::lavInspect(lavaanModel, 
                                    "data"))
  # remove rows with all NA
  ySorted <- ySorted[!apply(ySorted,1,function(x)all(is.na(x))),]
  
  SEM <- lessSEM:::.SEMFromLavaan(lavaanModel = lavaanModel)
  
  for(ro in 1:nrow(pars)){
    
    trainSet <- ySorted[!subsets[,pars$trainSet[ro]],,drop = FALSE]
    testSet <- ySorted[subsets[,pars$trainSet[ro]],,drop = FALSE]
    
    means <- apply(trainSet,2,mean,na.rm = TRUE)
    standardDeviations <- apply(trainSet,2,sd, na.rm = TRUE)
    
    testSet <- lessSEM::cvScaler(testSet = testSet, 
                                 means = means, 
                                 standardDeviations = standardDeviations)
    
    SEM <- lessSEM:::.setParameters(SEM = SEM, 
                                    labels = parameterLabels, 
                                    values = unlist(pars[ro, parameterLabels]),
                                    raw = FALSE)
    SEM$fit()
    
    m2LL <- 0
    for(i in 1:nrow(testSet)){
      isMissing <- is.na(testSet[i,])
      if(all(isMissing))
        next
      
      m2LL <- m2LL + (-2)*sum(mvtnorm::dmvnorm(
        x = testSet[i,!isMissing,drop = FALSE], 
        mean = SEM$impliedMeans[!isMissing,,drop = FALSE],
        sigma = SEM$impliedCovariance[!isMissing,!isMissing,drop = FALSE],
        log = TRUE
      ))
    }
    
    sel <- cvfits$lambda == pars$lambda[ro]
    if(sum(sel) != 1) stop("Error when selecting cv target")
    cvfits[sel,"cvfit"] <- cvfits[sel,"cvfit"] - m2LL
    
  }
  
  testthat::expect_equal(all(abs(cvfits$cvfit)< 1e-6), TRUE)
})

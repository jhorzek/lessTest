test_that("wls basics works", {
  # The following is adapted from lavaan::simulateData
  library(lessSEM)
  library(testthat)
  
  population.model <- 'f1 =~ x1 + 0.8*x2 + 1.2*x3
                      f2 =~ x4 + 0.5*x5 + 1.5*x6
                      f3 =~ x7 + 0.1*x8 + 0.9*x9
                      
                      x1 ~~ .3*x2
                      
                      f3 ~ 0.3*f1 + 0*f2
                    '
  
  # generate data
  set.seed(1234)
  N <- 1000L
  myData <- simulateData(population.model, sample.nobs = N)
  
  # fit model
  myModel <- ' f1 =~ x1 + x2 + x3
               f2 =~ x4 + x5 + x6
               f3 =~ x7 + x8 + x9
               x1 ~~ x2
               f3 ~ r1*f1 + r2*f2 '
  
  for(estimator in c("wls", "uls", "dwls", "gls")){
    for(meanstructure in c(TRUE, FALSE)){
      fit <- sem(myModel, 
                 data = myData, 
                 estimator = estimator,
                 meanstructure = meanstructure)
      
      fitLessSEM <- lessSEM::lasso(lavaanModel = fit, 
                                   regularized = c("r1", "r2"), 
                                   nLambdas = 100, 
                                   reverse = FALSE)
      
      lavaanFit <- .5*(fitLessSEM@fits$objectiveValue * (N-1) / N^2)
      
      par <- coef(fitLessSEM)@estimates
      
      for(i in 1:nrow(par)){
        pt <- parameterTable(fit)
        pt$label[pt$label == ""] <- paste0(pt$lhs, pt$op, pt$rhs)[pt$label == ""]
        for(l in pt$label){
          if(l %in% colnames(par)){
            pt$ustart[pt$label == l] <- par[i, l]
            pt$start[pt$label == l] <- par[i, l]
            pt$est[pt$label == l] <- par[i, l]
          }
        }
        fit_i <- lavaan::sem(pt, 
                             data = myData, 
                             estimator = estimator,
                             do.fit = FALSE)
        testthat::expect_true(abs(fit_i@Fit@fx - lavaanFit[i]) < 1e-6)
        
      }
    }
  }
  
})

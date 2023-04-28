test_that("WLS works", {
  library(testthat)
  library(lessSEM)
  
  set.seed(123)
  
  dataset <- simulateExampleData(N = 1000)
  datasetOrdered <- round(dataset)
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  for(estimator in c("WLS", "ULS", "GLS", "DWLS")){
    
    for(meanstructure in c(TRUE, FALSE)){
      
      for(ordered in c(TRUE, FALSE)){
        
        if(ordered){
          if(estimator != "ULS")
            next
          
          dat <- datasetOrdered
        }else{
          dat <- dataset
        }
        
        lavaanModel <- lavaan::sem(lavaanSyntax,
                                   data = dat,
                                   estimator = estimator,
                                   meanstructure = meanstructure,
                                   std.lv = TRUE,
                                   ordered = ordered)
        
        coefLavaan <- coef(lavaanModel)
        coefLavaanStart <- coefLavaan
        coefLavaanStart[] <- runif(length(coefLavaanStart), 
                                   .1, 
                                   max = .3)
        
        fit <- try(bfgs(lavaanModel = lavaanModel, 
                        control = controlBFGS(startingValues = coefLavaanStart)),
                   silent = TRUE)
        
        if(ordered){
          expect_true(is(fit, "try-error"))
        }else{
          
          N <- nrow(dat)
          
          expect_true(
            abs(
              (.5*(N-1)/N^2)*fit@fits$objectiveValue -
                fitMeasures(lavaanModel, "fmin")
            ) < .0001)
          
          # note: we take the absolute value of the parameters, because the model is only
          # locally identified
          expect_true(all(abs(abs(coef(fit)@estimates[,names(coefLavaan)]) - abs(coefLavaan)) < 1e-2))
        }
      }
    }
    
  }
  
})

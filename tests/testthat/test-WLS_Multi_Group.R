test_that("Multigroup WLS works", {
  library(testthat)
  library(lessSEM)
  
  set.seed(123)
  N = 1000
  dataset1 <- simulateExampleData(N = N)
  dataset2 <- simulateExampleData(N = N)
  dataset3 <- simulateExampleData(N = N)
  
  lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"
  
  for(estimator in c("WLS", "ULS", "GLS", "DWLS")){
    
    for(meanstructure in c(TRUE, FALSE)){
      
      lavaanModel <- lavaan::sem(lavaanSyntax,
                                 data = dataset1,
                                 estimator = estimator,
                                 meanstructure = meanstructure,
                                 std.lv = TRUE)
      pt <- parTable(lavaanModel)
      pt$label[pt$label == ""] <- paste0(pt$lhs, pt$op, pt$rhs)[pt$label == ""]
      
      pt1 <- pt2 <- pt3 <- pt
      pt1$label <- paste0(pt1$label, "_g1")
      pt2$label <- paste0(pt2$label, "_g2")
      pt3$label <- paste0(pt3$label, "_g3")
      
      model1 <- lavaan::sem(model = pt1,
                            data = dataset1,
                            estimator = estimator,
                            meanstructure = meanstructure,
                            std.lv = TRUE)
      
      model2 <- lavaan::sem(pt2,
                            data = dataset2,
                            estimator = estimator,
                            meanstructure = meanstructure,
                            std.lv = TRUE)
      
      model3 <- lavaan::sem(pt3,
                            data = dataset3,
                            estimator = estimator,
                            meanstructure = meanstructure,
                            std.lv = TRUE)
      
      
      
      coefLavaan <- c(
        coef(model1),
        coef(model2), 
        coef(model3)
      )
      
      coefLavaanStart <- coefLavaan
      
      mgSEM <- lessSEM:::.multiGroupSEMFromLavaan(lavaanModels = c(model1, model2, model3), addMeans = FALSE)
      
      expect_true(abs(
        (.5*(N-1)/N)*mgSEM$fit()/N - 
          
          (fitMeasures(model1, "fmin") +
             fitMeasures(model2, "fmin") +
             fitMeasures(model3, "fmin"))) < 1e-5)
      
      coefLavaanStart[] <- runif(length(coefLavaanStart), 
                                 .1, 
                                 max = .3)
      
      fit <- try(
        bfgs(
          lavaanModel = c(model1, model2, model3), 
          control = controlBFGS(startingValues = coefLavaanStart, breakOuter = 1e-20)
        ),
        silent = TRUE)
      
      expect_true(abs(
        (.5*(N-1)/N^2)*fit@fits$objectiveValue - 
          
          (fitMeasures(model1, "fmin") +
             fitMeasures(model2, "fmin") +
             fitMeasures(model3, "fmin"))) < 1e-3)
      
      # note: we take the absolute value of the parameters, because the model is only
      # locally identified
      expect_true(all(abs(abs(coef(fit)@estimates[,names(coefLavaan)]) - abs(coefLavaan)) < 1e-2))
      
    }
  }
  
})

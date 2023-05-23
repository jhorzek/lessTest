test_that("AnomAuth works", {
  library(ctsemOMX)
  library(lessSEM)
  set.seed(123)
  
  data("AnomAuth")
  
  # initial time point
  lavaanSyntax <- 
    "eta1_T0 =~ 1 * Y1_T0
eta2_T0 =~ 1 * Y2_T0
Y1_T0 ~~ 0*Y1_T0
Y2_T0 ~~ 0*Y2_T0\n"
  
  # variances
  lavaanSyntax <- c(lavaanSyntax,
                    "eta1_T0 ~~ v0_11 * eta1_T0 + v0_12 * eta2_T0\neta2_T0 ~~ v0_22 * eta2_T0\n"
  )
  
  # means
  lavaanSyntax <- c(lavaanSyntax,
                    "eta1_T0 ~ 1\neta2_T0 ~ 1\nY1_T0~mMean1*1\nY2_T0~mMean2*1\n" 
  )
  
  for(tp in c(0,1,2,3)){
    if(tp < 2) {
      a <- "a1"
      v <- "v1"
    }else{
      a <-"a2"
      v <- "v2"
    }
    
    # autoregressive and cross-lagged
    lavaanSyntax <- c(lavaanSyntax,
                      paste0(
                        "eta1_T", tp+1, " ~ ", a, "_11 * eta1_T", tp, " + ", a, "_12 * eta2_T", tp,"\n",
                        "eta2_T", tp+1, " ~ ", a, "_21 * eta1_T", tp, " + ", a, "_22 * eta2_T", tp, "\n"
                      )
    )
    
    # variances
    lavaanSyntax <- c(lavaanSyntax,
                      paste0(
                        "eta1_T", tp+1, " ~~ ", v, "_11 * eta1_T", tp+1, " + ", v, "_12 * eta2_T", tp+1,"\n",
                        "eta2_T", tp+1, " ~~ ", v, "_22 * eta2_T", tp+1, "\n"
                      )
    )
    
    # loadings
    lavaanSyntax <- c(lavaanSyntax,
                      paste0(
                        "eta1_T", tp+1, " =~ 1 * Y1_T", tp+1,"\n",
                        "eta2_T", tp+1, " =~ 1 * Y2_T", tp+1,"\n"
                      )
    )
    
    # manifest variances
    lavaanSyntax <- c(lavaanSyntax,
                      paste0(
                        "Y1_T", tp+1, " ~~ 0* Y1_T", tp+1, "\n",
                        "Y2_T", tp+1, " ~~ 0* Y2_T", tp+1, "\n"
                      )
    )
    
    # manifest means
    lavaanSyntax <- c(lavaanSyntax,
                      paste0(
                        "Y1_T", tp+1, " ~ mMean1 * 1\n",
                        "Y2_T", tp+1, " ~ mMean2 * 1\n"
                      )
    )
  }
  lavaanSyntax <- paste0(lavaanSyntax, collapse = "")
  cat(lavaanSyntax)
  
  lavaanFit <- sem(model = lavaanSyntax, data = AnomAuth,
                   orthogonal.y = TRUE, 
                   orthogonal.x = TRUE,
                   missing = "ml",
                   do.fit = FALSE)
  
  transformations <- "
// Define all parameters which we want to use:
parameters: a1_11, a1_12, a1_21, a1_22, a2_11, a2_12, a2_21, a2_22, 
ctA_11, ctA_12, ctA_21, ctA_22, 
v1_11, v1_12, v1_22, v2_11, v2_12, v2_22, 
ctV_11, ctV_12, ctV_22

// Define the starting values for the continuous time parameters:
start: ctA_11 = -1, ctA_12 = 0, ctA_21 = 0, ctA_22 = -1, 
ctV_11 = .1, ctV_12 = 0, ctV_22 = .1

// transformations:
arma::mat drift(2,2);
arma::mat ARCL1(2,2);
arma::mat ARCL2(2,2);
arma::mat driftHash(4,4);
drift(0,0) = ctA_11;
drift(1,0) = ctA_21;
drift(0,1) = ctA_12;
drift(1,1) = ctA_22;
ARCL1 = expmat(drift);
ARCL2 = expmat(drift*2.0);

driftHash = kron(drift, arma::eye(2,2)) + kron(arma::eye(2,2), drift);

arma::mat diffusion(2,2);
arma::mat discreteDiff1(2,2);
arma::mat discreteDiff2(2,2);
diffusion(0,0) = ctV_11;
diffusion(1,0) = ctV_12;
diffusion(0,1) = ctV_12;
diffusion(1,1) = ctV_22;
discreteDiff1 = arma::reshape(arma::inv(driftHash) * 
  (expmat(driftHash) - arma::eye(arma::size(expmat(driftHash))))*
  arma::vectorise(diffusion),2,2);
discreteDiff2 = arma::reshape(arma::inv(driftHash) * 
  (expmat(driftHash*2.0) - arma::eye(arma::size(expmat(driftHash*2.0))))*
  arma::vectorise(diffusion),2,2);

// extract parameters

a1_11 = ARCL1(0,0);
a1_12 = ARCL1(0,1);
a1_21 = ARCL1(1,0);
a1_22 = ARCL1(1,1);

a2_11 = ARCL2(0,0);
a2_12 = ARCL2(0,1);
a2_21 = ARCL2(1,0);
a2_22 = ARCL2(1,1);

v1_11 = log(discreteDiff1(0,0)); // we take the log because of the internal 
// transformation in lessSEM
v1_12 = discreteDiff1(0,1);
v1_22 = log(discreteDiff1(1,1)); // we take the log because of the internal 
// transformation in lessSEM

v2_11 = log(discreteDiff2(0,0)); // we take the log because of the internal 
// transformation in lessSEM
v2_12 = discreteDiff2(0,1);
v2_22 = log(discreteDiff2(1,1)); // we take the log because of the internal 
// transformation in lessSEM
"
  
  lessSEMFit <- bfgs(lavaanModel = lavaanFit,
                     # Our model modification must make use of the modifyModel - function:
                     modifyModel = modifyModel(transformations = transformations),
                     control = controlBFGS(breakOuter = 1e-20)
  )
  
  AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2), 
                           Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL) 
  AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
  
  testthat::expect_true(all(
    abs(
      AnomAuthfit$mxobj$fitfunction$result[[1]] - lessSEMFit@fits$m2LL
    ) < 1e-3))
  
  testthat::expect_true(all(
    abs(
      c(
        c(matrix(coef(lessSEMFit)@estimates[,c("ctA_11", "ctA_21", "ctA_12", "ctA_22")],2,2) - AnomAuthfit$mxobj$DRIFT$values),
        c(matrix(coef(lessSEMFit)@estimates[,c("ctV_11", "ctV_12", "ctV_12", "ctV_22")],2,2) - AnomAuthfit$mxobj$DIFFUSION$result)
      )
    ) < 1e-3))
})

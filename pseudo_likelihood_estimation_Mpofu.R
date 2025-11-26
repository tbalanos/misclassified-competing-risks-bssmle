################################################################################
# Author: Adapted by Theofanis Balanos from Philani Brian Mpofu (2020)
# Purpose: R functions for performing pseudo-likelihood estimation of
#          misclassification probabilities under the double-sampling framework.
#
# Reference:
# Mpofu, P. B., Bakoyannis, G., Yiannoutsos, C. T., Mwangi, A., & Mburu, M. (2020).
# "Semiparametric regression analysis of misclassified competing risks data."
# Statistics in Medicine, 39(23–24), 3199–3216.
#
# Description:
# The script defines a series of functions that implement the pseudo-likelihood
# estimation procedure proposed by Mpofu et al. (2020) for studies with
# internal-validation sampling (double sampling).
################################################################################

# Load required packages
library(sandwich)
library(Hmisc)
library(lmtest)

# ------------------------------------------------------------------------------
# Function: data_preper
# Purpose : Prepare a dataset for pseudo-likelihood estimation
# ------------------------------------------------------------------------------
data_preper <- function(dat, true_out, surr_out, out_interest, ds_var) {
  dat$r <- dat[, names(dat) %in% c(ds_var)]
  dat$c <- dat[, names(dat) %in% c(true_out)]
  dat$c_obs <- dat[, names(dat) %in% c(surr_out)]
  dat$d1 <- as.numeric(dat$c == out_interest)
  dat$d1_obs <- as.numeric(dat$c_obs == out_interest)
  return(dat)
}

# ------------------------------------------------------------------------------
# Function: prob.compute
# Purpose : Estimate predictive values for missing true outcomes
# ------------------------------------------------------------------------------
prob.compute <- function(dat, formula_pred_val) {
  data_list <- list()
  dat$w1 <- as.numeric((1 - dat$d1_obs) == 1 & dat$r == 1 & dat$c_obs > 0)
  dat$w2 <- as.numeric(dat$d1_obs == 1 & dat$r == 1 & dat$c_obs > 0)
  
  # Convert formula to character for flexible construction
  formula_pred_val <- as.character(formula_pred_val)
  formula_pred_val <- paste(formula_pred_val[1], formula_pred_val[2])
  pred_val1 <- paste("d1", formula_pred_val)
  pred_val2 <- paste("(1-d1)", formula_pred_val)
  
  # Fit logistic models for predictive values
  model_cause1 <- tryCatch(
    glm(pred_val1, family = binomial(link = "logit"), weights = w1, data = dat),
    warning = function(w) return(99)
  )
  model_cause2 <- tryCatch(
    glm(pred_val2, family = binomial(link = "logit"), weights = w2, data = dat),
    warning = function(w) return(99)
  )
  
  # Predict probabilities if both models converge
  if (typeof(model_cause1) != "double" && typeof(model_cause2) != "double") {
    dat$p_t1_obs2 <- predict(model_cause1, newdata = dat, type = "response")
    dat$p_t2_obs1 <- predict(model_cause2, newdata = dat, type = "response")
  }
  
  # Augment dataset with predictive probabilities
  if (typeof(model_cause1) != "double" &&
      typeof(model_cause2) != "double" &&
      typeof(dat$p_t1_obs2) != "character" &&
      typeof(dat$p_t2_obs1) != "character") {
    data <- dat
    data$y1 <- (data$r * data$d1) +
      (1 - data$r) * ((dat$d1_obs) * (1 - data$p_t2_obs1) +
                        (1 - dat$d1_obs) * data$p_t1_obs2)
    
    score1 <- estfun(model_cause1)
    score2 <- estfun(model_cause2)
    v1 <- vcov(model_cause1)
    v2 <- vcov(model_cause2)
    
    data_list <- list(data = data, score1 = score1, score2 = score2, v1 = v1, v2 = v2)
  }
  return(data_list)
}

# ------------------------------------------------------------------------------
# Function: coef.mat
# Purpose : Fit the pseudo-likelihood misclassification model
# ------------------------------------------------------------------------------
coef.mat <- function(data, formula_mis) {
  formula_mis <- as.character(formula_mis)
  formula_mis <- paste(formula_mis[1], formula_mis[2])
  formula_misclass <- paste("(1-d1_obs)", formula_mis)
  
  mod <- suppressWarnings(glm(formula_misclass,
                              family = binomial(link = "logit"),
                              weights = y1,
                              data = data))
  mat <- mod$coefficients
  s <- estfun(mod)
  V <- vcov(mod)
  p <- predict(mod, newdata = data, type = "response")
  names(mat) <- names(mod$coefficients)
  
  list(pseudo.coef = mat,
       pseudo.score = s,
       pseudo.V = V,
       pseudo.p = p,
       model_miscl = mod)
}

# ------------------------------------------------------------------------------
# Function: R_create
# Purpose : Compute derivative of the score w.r.t. predictive parameters
# ------------------------------------------------------------------------------
R_create <- function(dat, formula_pred_val, formula_mis, ps_list) {
  Z <- model.matrix(formula_pred_val, data = dat)
  X <- model.matrix(formula_mis, data = dat)
  w11 <- with(dat, (1 - r) * (1 - d1_obs) * (1 - p_t1_obs2) * p_t1_obs2)
  w12 <- with(dat, (1 - r) * d1_obs * (-p_t2_obs1) * (1 - p_t2_obs1))
  w2 <- with(dat, d2_obs - ps_list$pseudo.p)
  
  R1 <- (t(Z * w11 * w2) %*% X)
  R2 <- (t(Z * w12 * w2) %*% X)
  list(R1 = R1, R2 = R2)
}

# ------------------------------------------------------------------------------
# Function: sand_compute
# Purpose : Closed-form sandwich variance estimator
# ------------------------------------------------------------------------------
sand_compute <- function(plug_list, ps_list, R) {
  H <- ps_list$pseudo.V
  v1 <- plug_list$v1
  v2 <- plug_list$v2
  var1 <- H %*% (t(ps_list$pseudo.score)
                 + R$R1 %*% v1 %*% t(plug_list$score1)
                 + t(R$R2) %*% v2 %*% t(plug_list$score2))
  var1 <- var1 %*% t(var1)
  var1
}

# ------------------------------------------------------------------------------
# Function: misclass_ps_est
# Purpose : Wrapper to compute pseudo-likelihood estimates and SEs
# ------------------------------------------------------------------------------
misclass_ps_est <- function(dat, true_out, surr_out, out_interest,
                            formula_pred_val, formula_mis, ds_var) {
  dat_set <- data_preper(dat, true_out, surr_out, out_interest, ds_var)
  pred_val_mod <- prob.compute(dat_set, formula_pred_val)
  misMod <- coef.mat(pred_val_mod$data, formula_mis)
  jacob <- R_create(dat = pred_val_mod$data,
                    formula_pred_val = formula_pred_val,
                    formula_mis = formula_mis,
                    ps_list = misMod)
  varEst <- sand_compute(plug_list = pred_val_mod, ps_list = misMod, R = jacob)
  parEst <- misMod$pseudo.coef
  test_summary <- coeftest(misMod$model_miscl, varEst)
  
  # Return model_miscl so we can predict later
  list(parEst = parEst,
       varEst = varEst,
       test_summary = test_summary,
       model_miscl = misMod$model_miscl)
}


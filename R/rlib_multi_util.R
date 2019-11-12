##' @importFrom glmnet cv.glmnet
# approx_susie_no_intercept = function(x, y, w, only_weight = F) {
#   x_susie = mydiag(sqrt(w)) %*% x
#   y_susie = mydiag(sqrt(w)) %*% y
#   if(only_weight == T) {
#     return(list(x = x_susie, y = y_susie))
#   }
#   qr_x = qr(x)
#   x_tilde = qr.Q(qr_x)
#   bhat = solve(x_tilde, y)
#   mod = cv.glmnet(x, y, weights = w, standardize = F, intercept = F)
#   beta = as.vector(mod$glmnet.fit$beta[, which.min(mod$cvm)])
#   sigma = sqrt(mean((y - x %*% beta)^2))
#   list(x = x_susie / sigma, y = y_susie / sigma)
# }

#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef residuals
approx_susie = function(x, y, w = NULL, intercept = T, only_weight = F) {
  if(!is.null(w)) {
    x_susie = mydiag(sqrt(w)) %*% x
    y_susie = mydiag(sqrt(w)) %*% y
    if(intercept == T) {
      intercept_susie = sqrt(w)
      mod = lm(y_susie ~ intercept_susie - 1)
      y_susie = residuals(mod)
    }
    if(only_weight == T) {
      return(list(x = x_susie, y = y_susie))
    }
    mod = glmnet::cv.glmnet(x, y, weights = w, standardize = F, intercept = intercept)
    beta = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
    x_pred = cbind(rep(1, nrow(x)), x)
    # sigma = sqrt(sum((y - x_pred %*% beta)^2) / sum(1 / w))
    sigma = sqrt(mean((y - x_pred %*% beta)^2 * w))

    return(list(x = x_susie / sigma, y = y_susie / sigma, sigma = sigma))
  } else {
    mod = glmnet::cv.glmnet(x, y, standardize = F, intercept = intercept)
    beta = as.vector(coef(mod$glmnet.fit)[, which.min(mod$cvm)])
    x_pred = cbind(rep(1, nrow(x)), x)
    sigma = sqrt(mean((y - x_pred %*% beta)^2))
    if(intercept == T) {
      y = scale(y, T, F)
    }
    return(list(x = x / sigma, y = y / sigma, sigma = sigma))
  }
}

#' @importFrom susieR susie
run_susie_default = function(x, y, ...) {
  mod = tryCatch({
    susie(x, y, ...)
  }, error = function(e) {
    vars = data.frame(variable = 1 : ncol(x), variable_prob = 0, cs = -1)
    list(vars = vars, cs = NULL)
  })
  mod
}

mydiag = function(vals) {
  diag(matrix(vals)[, 1])
}

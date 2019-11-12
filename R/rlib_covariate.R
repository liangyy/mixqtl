#' Regress against covariates for trcQTL
#'
#' Given total read count (trc, N x 1) and library size (lib_size, N x 1) along with covarites (N x C),
#' regress covariates against log(trc / lib_size / 2) by two-step approach
#' 1. Regress all covariates against response and select significant ones (alpha = 0.05)
#' 2. Regress the selected covariates against response and return the predicted response
#'
#' @param trc total read count (dimension = N x 1)
#' @param lib_size library size (dimension = N x 1)
#' @param covariates covariates (dimension = C x N) with the first column indicating the covariate name
#'
#' @return predicted response based on the selected covariates (if no covariates are selected, fill in zeros)
#'
#' @examples
#' regress_against_covariate(
#'   trc = rpois(100, 100),
#'   lib_size = rpois(100, 10000),
#'   covariates = cbind(
#'     data.frame(name = paste0('cov', 1:10)),
#'     as.data.frame(matrix(runif(1000), nrow = 10))
#'   )
#' )
#'
#' @export
regress_against_covariate = function(trc, lib_size, covariates) {
  out = effect_size_trc_with_covariates_get_coef_refactored(trc, lib_size, covariates)
  selected = abs(out[-1, 3]) > 2
  if(sum(selected) > 0) {
    out = effect_size_trc_with_covariates_get_coef_refactored(trc, lib_size, covariates)
    indiv_offset = t(covariates[selected, -1]) %*% out[-1, 1]
  } else {
    indiv_offset = rep(0, ncol(covariates) - 1)  # t(covariates[selected, -1]) %*% out[-1, 1]
  }
  indiv_offset
}

effect_size_trc_with_covariates_get_coef_refactored = function(trc, lib_size, cov_i) {
  lhs = log(trc / lib_size / 2)
  rhs_indiv = 1 : length(trc)
  regress_out_covar(lhs, rhs_indiv, cov_i)
}

#' @import dplyr
#' @importFrom stats as.formula lm
regress_out_covar = function(y, indiv, cov_i) {
  df = data.frame(lhs = y, indiv = indiv)
  df = df %>% filter(rowSums(is.na(df)) == 0)
  is_infinite = rowSums(matrix(unlist(lapply(df, is.infinite)), byrow = F, ncol = ncol(df))) > 0
  df = df %>% filter(!is_infinite)
  row.names(cov_i) = paste0('c_', cov_i[, 1])
  cov_i = cov_i[, -1]
  df = cbind(df, t(cov_i[ match(df$indiv, 1 : ncol(cov_i))]))
  form = paste0('lhs ~  1 + ', paste0(row.names(cov_i), collapse = ' + '))
  model = lm(as.formula(form), data = df)
  out = summary(model)$coefficients
  return(out)
}

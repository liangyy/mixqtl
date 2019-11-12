#' Filling in missing genotype
#'
#' Filling in the missing genotype as the average.
#'
#' @param x genotype (dimension N x P)
#'
#' @return genotype matrix with missing values being imputed
#'
#' @examples
#' impute_geno(
#'   x = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2)
#' )
#'
#' @export
impute_geno = function(x) {
  mx = colMeans(x, na.rm = T)
  if(length(mx) == 1) {
    x[is.na(x)] = mx
    return(x)
  }
  t(apply(x, 1, function(x) {
    x[is.na(x)] = mx[is.na(x)]
    return(x)
  }))
}

harmonic_sum_ = function(x, y) {
  1 / (1 / x + 1 / y)
}

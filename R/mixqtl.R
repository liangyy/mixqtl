#' Combining total and allele-specific reads for QTL mapping
#'
#' Given total and allele-specific read counts along with library size,
#' MixQTL procedure (by meta-analysis) is implemented.
#' Specifically, first, ascQTL and trcQTL estimates are obtained.
#' Then, second, meta-analyze ascQTL and trcQTL estimates to obtain mixQTL estimate.
#' The function outputs the estimates of outcome versus each of the input variants  
#'
#' @param geno1 genotype of haplotype 1 (dimension = N x P)
#' @param geno2 genotype of haplotype 2 (dimension = N x P)
#' @param y1 allele-specific read count for haplotype 1 (dimension = N x 1)
#' @param y2 allele-specific read count for haplotype 2 (dimension = N x 1)
#' @param ytotal total read count (dimension = N x 1)
#' @param lib_size library size (dimension = N x 1)
#' @param cov_offset predicted effect of covariates on total read count response term,
#' log(ytotal / lib_size / 2) (dimension = N x 1)
#' @param trc_cutoff total read count cutoff to exclude observations with ytotal lower than the cutoff
#' @param asc_cutoff allele-specific read count cutoff to exclude observations with y1 or y2 lower than asc_cutoff
#' @param weight_cap the maximum weight difference (in fold) is min(weight_cap, floor(sample_size / 10)). The ones exceeding the cutoff is capped.
#' @param asc_cap exclude observations with y1 or y2 higher than asc_cap
#'
#' @return QTL mapping result
#'
#' @examples
#' mixqtl(
#'   geno1 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   geno2 = matrix(sample(c(0, 0.5, 1), 200, replace = TRUE), ncol = 2),
#'   y1 = rpois(100, 100),
#'   y2 = rpois(100, 100),
#'   ytotal = rpois(100, 1000),
#'   lib_size = rpois(100, 10000),
#'   cov_offset = runif(100),
#'   trc_cutoff = 100,
#'   asc_cutoff = 50,
#'   weight_cap = 100,
#'   asc_cap = 1000
#' )
#'
#' @export
mixqtl = function(geno1, geno2, y1, y2, ytotal, lib_size, cov_offset = NULL, trc_cutoff = 20, asc_cutoff = 5, weight_cap = 100, asc_cap = 5000) {
  if(is.null(cov_offset)) {
    cov_offset = rep(0, length(lib_size))
  }
  # prepare X
  h1 = geno1
  h1[is.na(h1)] = 0.5
  h2 = geno2
  h2[is.na(h2)] = 0.5
  Xasc = h1 - h2
  Xtrc = (h1 + h2) / 2

  trc = matrix_ls_trc(ytotal, lib_size, Xtrc, cov_offset, trc_cutoff = trc_cutoff)
  asc = matrix_ls_asc(y1, y2, Xasc, asc_cutoff = asc_cutoff, weight_cap = weight_cap, asc_cap = asc_cap)
  
  meta_result = meta_analyze(trc, asc)
  meta_list = list()
  for(i in names(meta_result)) {
    d = as.data.frame(lapply(meta_result[[i]], t))
    colnames(d) = paste0(colnames(d), '.', i)
    meta_list[[length(meta_list) + 1]] = d
  }
  merge = do.call(cbind, meta_list)
  merge
}

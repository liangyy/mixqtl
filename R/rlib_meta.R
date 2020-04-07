#' Meta-analyze trcQTL and ascQTL results
#'
#' For each variant/gene pair, meta-analyze summary statistics from trcQTL and ascQTL.
#' If sample size is greater than or equal to n_cutoff,
#' approximate test statistic as z-score,
#' otherwise compute p-value from student-t distribution with d.o.f. = 1.
#' If both trcQTL and ascQTL have no fewer than n_cutoff samples, perform
#' inverse variance-based meta-analysis.
#'
#' @param trc list of effect size estimates and standard deviation (e.g. output of matrix_ls_trc)
#' @param asc list of effect size estimates and standard deviation (e.g. output of matrix_ls_asc)
#' @param n_cutoff threshold to approximate test statistic as z-score.
#'
#' @return a list of summary statistics and p-values obtained from trc, asc, and meta-analysis of the two
#'         the object within list is named by the method (trc, asc, and meta) and each of which includes
#'         pval (p-value), stat (test statistic), stat_type (z-value or t-value),
#'         bhat (effect size estimate), se (standard deviation of effect size),
#'         method (only for 'meta' object showing if the result is from 'trc', 'asc', or 'meta')
#'
#' @examples
#' meta_analyze(
#'   trc = list(
#'     beta_hat = matrix(rnorm(100), ncol = 20),
#'     beta_se = abs(matrix(rnorm(100), ncol = 20)),
#'     sample_size = 100
#'   ),
#'   asc = list(
#'     beta_hat = matrix(rnorm(100), ncol = 20),
#'     beta_se = abs(matrix(rnorm(100), ncol = 20)),
#'     sample_size = 100
#'   )
#' )
#'
#' @export
meta_analyze = function(trc, asc, n_cutoff = 15) {
  df = list(
    beta_trc = trc$beta_hat, beta_se_trc = trc$beta_se, sample_size_trc = trc$sample_size,
    beta_asc = asc$beta_hat, beta_se_asc = asc$beta_se, sample_size_asc = asc$sample_size
  )
  my_meta_fast_(df$beta_trc, df$beta_se_trc, df$sample_size_trc, df$beta_asc, df$beta_se_asc, df$sample_size_asc, n_cutoff = n_cutoff)
}

# get_significance_from_estimate = function(est) {
#   df = list(
#     beta_hat = est$beta_hat, beta_se = est$beta_se, sample_size = est$sample_size
#   )
#   o = get_pval_fast_(df$beta_hat, df$beta_se, df$sample_size)
#   o$bhat = df$beta_hat
#   o$se = df$beta_se
#   o
# }

my_meta_fast_ = function(b1, s1, n1, b2, s2, n2, method = c('trc', 'asc'), n_cutoff = 15) {
  pval1 = get_pval_fast_(b1, s1, n1)
  pval1$bhat = b1
  pval1$se = s1
  pval2 = get_pval_fast_(b2, s2, n2)
  pval2$bhat = b2
  pval2$se = s2

  if(n1 >= n2) {
    pval_meta = pval1
    pval_meta$method = matrix(method[1], ncol = dim(b1)[2], nrow = dim(b1)[1])
    pval_meta$bhat = b1
    pval_meta$se = s1
    na_1_notna_2 = is.na(b1) & !is.na(b2)
    pval_meta$bhat[na_1_notna_2] = b2[na_1_notna_2]
    pval_meta$se[na_1_notna_2] = s2[na_1_notna_2]
    pval_meta$pval[na_1_notna_2] = pval2$pval[na_1_notna_2]
    pval_meta$stat[na_1_notna_2] = pval2$stat[na_1_notna_2]
    pval_meta$method[na_1_notna_2] = method[2]
    pval_meta$stat_type = NULL
  } else {
    pval_meta = pval2
    pval_meta$method = matrix(method[2], ncol = dim(b2)[2], nrow = dim(b2)[1])
    pval_meta$bhat = b2
    pval_meta$se = s2
    na_2_notna_1 = is.na(b2) & !is.na(b1)
    pval_meta$bhat[na_2_notna_1] = b1[na_2_notna_1]
    pval_meta$se[na_2_notna_1] = s1[na_2_notna_1]
    pval_meta$pval[na_2_notna_1] = pval1$pval[na_2_notna_1]
    pval_meta$stat[na_2_notna_1] = pval1$stat[na_2_notna_1]
    pval_meta$method[na_2_notna_1] = method[1]
    pval_meta$stat_type = NULL
  }

  if(n1 >= n_cutoff & n2 >= n_cutoff) {
    o = combine_two_stats(b1, b2, s1, s2)
    pval0 = get_pval_fast_(o$b, o$s, n1 + n2)
    b0 = o$b
    s0 = o$s
    notna_meta = !is.na(pval0$pval)
    pval_meta$bhat[notna_meta] = b0[notna_meta]
    pval_meta$se[notna_meta] = s0[notna_meta]
    pval_meta$pval[notna_meta] = pval0$pval[notna_meta]
    pval_meta$stat[notna_meta] = pval0$stat[notna_meta]
    pval_meta$method[notna_meta] = 'meta'
  }
  out = list()
  out[[method[1]]] = pval1
  out[[method[2]]] = pval2
  out[['meta']] = pval_meta
  return(out)
}

my_meta_ = function(b1, s1, n1, b2, s2, n2, method = c('trc', 'asc')) {
  if(n1 >= 15 & n2 >= 15) {
    o = combine_two_stats(b1, b2, s1, s2)
    pval = get_pval(o$b, o$s, n1 + n2)
    pval$method = 'meta'
    bout = o$b
    sout = o$s
  } else if (n1 > n2) {
    pval = get_pval(b1, s1, n1)
    bout = b1
    sout = s1
  } else if (n1 < n2) {
    pval = get_pval(b2, s2, n2)
    pval$method = method[2]
    bout = b2
    sout = s2
  } else {
    pval1 = get_pval(b1, s1, n1)
    pval2 = get_pval(b2, s2, n2)
    if(pval1$pval < pval2$pval) {
      pval = pval1
      pval$method = method[1]
      bout = b1
      sout = s1
    } else {
      pval = pval2
      pval$method = method[2]
      bout = b2
      sout = s2
    }
  }
  pval$bhat = bout
  pval$se = sout
  return(pval)
}

get_pval_fast_ = function(b, s, n, n_cutoff = 15) {
  if(n > n_cutoff) {
    o = list(pval = z2p_(b / s), stat = b / s)
    o$stat_type = 'z-val'
  } else {
    o = list(pval = t2p_(b / s, n), stat = b / s)
    o$stat_type = 't-val'
  }
  # o$stat_type[is.na(o$stat)] = NA
  o
}

get_pval = function(b, s, n) {
  if(n > 15) {
    z = b / s
    return(data.frame(pval = z2p_(z), stat = z, stat_type = 'z-val'))
  } else {
    t = b / s
    return(data.frame(pval = t2p_(z, n), stat = t, stat_type = 't-val'))
  }
}

#' @importFrom stats pnorm
#' @importFrom stats pt
z2p_ = function(z) {
  exp(pnorm(abs(z), log.p=T, lower.tail=F)) * 2
}

t2p_ = function(t, dof) {
  exp(pt(abs(t), df = dof, log.p=T, lower.tail=F)) * 2
}

combine_two_stats = function(b1, b2, s1, s2) {
  w1 = 1 / s1^2
  w2 = 1 / s2^2
  b = (w1 * b1 + w2 * b2) / (w1 + w2)
  s = sqrt(1 / (w1 + w2))
  return(list(b = b, s = s))
}

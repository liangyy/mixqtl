generate_gene_instance = function(L_gene, library_dist, theta, f, maf, vec_pos) {

  ## This function generates a gene region to work with
  ## For a gene region, the following parameters are required
  ##   1. L_gene: length of gene region
  ##   2. lambda_lib: mean library size, Ti_lib ~ Poisson(lambda_lib)
  ##   3. theta: parameters for thetai, for instance, when theta$type == 'beta', thetai ~ Beta(theta$alpha, theta$beta)
  ##   4. S: position of heterozygous sites, Sj ~ Sample(1 : L_gene) where the length of S is defined by L_max ~ Binomial(L_gene, f)
  ##   5. fj: maf of heterozygous sites, fj ~ Unif(maf[1], maf[2])
  ##   6. vec_pos: probability of having reads by position
  ## Input
  ##   f: frequency of heterozygous sites
  ##   maf: MAF range of the heterozygous sites within the gene region

  L_max = rbinom(1, L_gene, f)  # the total number of heterozygous sites on the transcripts
  S = sample(1 : L_gene, L_max, replace = F)  # positions of heterozygous sites
  fj = runif(length(S), min = maf[1], max = maf[2])  # maf of heterozygous sites
  return(list(L_gene = L_gene, library_dist = library_dist, theta = theta, S = S, fj = fj, vec_pos = vec_pos))

}

guess_alpha_beta_from_mu_sd = function(mub, sdb) {

  ## return alpha, beta such that approximately X ~ Beta(alpha, beta) with E(X) = mub, Var(X) = sdb^2
  ## here assumes that mub is small (< 1/100)
  ## approximate is done by:
  ##   mub ~ alpha / beta
  ##   sdb^2 ~ mub * (1 - mub) / beta

  beta = (mub * (1 - mub)) / (sdb ^ 2)
  alpha = mub * beta
  return(list(type = 'beta', alpha = alpha, beta = beta))

}

guess_k_sigma_from_mu_sd = function(mub, sdb) {

  ## return K, sigma such that approximately log X ~ N(K, sigma) with E(X) = mub, Var(X) = sdb^2
  ## this is done by:
  ##   e^(K + sigma/2) = mub
  ##   (e^sigma - 1) e^(2K) e^(sigma) = sdb

  sigma = log(1 + sdb^2 / mub^2)
  k = log(mub) - 1 / 2 * sigma
  return(list(type = 'lognormal', k = k, sigma = sigma))

}

read2data = function(P1, P2, S, Zij, L_read) {
  Si = S[Zij == 1]
  if(length(P1) > 0) {
    overlap_heter1 = sapply(P1, is_overlap, sites = Si, lread = L_read)
    y1 = sum(overlap_heter1)
    if(y1 > 0) {
      pos1 = P1[overlap_heter1]
      ase_site = sapply(pos1, function(x) {
        min(Si[x <= Si])
      })
      ase_count1 = rep(0, length(Si))
      for(i in 1 : length(Si)) {
        ase_count1[i] = sum(ase_site == Si[i])
      }
    } else {
      ase_count1 = NA
    }
  } else {
    y1 = 0
    ase_count1 = NA
  }
  if(length(P2) > 0) {
    overlap_heter2 = sapply(P2, is_overlap, sites = Si, lread = L_read)
    y2 = sum(overlap_heter2)
    if(y2 > 0) {
      pos2 = P2[overlap_heter2]
      ase_site = sapply(pos2, function(x) {
        min(Si[x <= Si])
      })
      ase_count2 = rep(0, length(Si))
      for(i in 1 : length(Si)) {
        ase_count2[i] = sum(ase_site == Si[i])
      }
    } else {
      ase_count2 = NA
    }
  } else {
    y2 = 0
    ase_count2 = NA
  }
  y = length(P1) + length(P2)
  ystar = y - y1 - y2
  y1star = length(P1) - y1
  y2star = length(P2) - y2
  observed = data.frame(y1 = y1, y2 = y2, ystar = ystar, L = sum(Zij))
  hidden = data.frame(y1star = y1star, y2star = y2star)
  return(list(observed = observed, hidden = hidden, ase_bysite = list(h1 = ase_count1, h2 = ase_count2)))
}

is_overlap = function(pos, sites, lread) {
  pos_start = pos
  pos_end = pos + lread - 1
  temp = (pos_start <= sites) & (sites <= pos_end)
  sum(temp) > 0
}

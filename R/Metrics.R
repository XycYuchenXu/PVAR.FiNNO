#' @keywords internal
objfun = function(GK, XTS, eta, Phi_BL, Phi, W, S, Gamma, rho, M = nrow(W),
                  p = nrow(Phi), TT = sapply(XTS, ncol) - 1){
  dPhi = Phi - Phi_BL
  s0 = 0
  s1 = rho * (.5 * sum(dPhi^2) + sum(dPhi * Gamma))
  
  for (m in 1:M) {
    Gm = GK[[1]][m,,]; Km = GK[[2]][m,,]
    
    Am = W[m,] * Phi + S[[m]]
    s0 = s0 + .5 * sum(crossprod(Am) * Gm) + .5 * sum(XTS[[m]][,2:(TT[m]+1)]^2) / TT[m]-
      sum(t(Am) * Km) + eta * sum(abs(S[[m]]))
  }
  return(c(s0, s1 + s0))
  
}

#' @keywords internal
IC_PVAR = function(XTS, W, S, Phi, C = 1, TT = sapply(XTS, ncol) - 1,
                   M = nrow(W), p = ncol(W)){

  r = which(cumsum(svd(Phi)$d) >= 0.95 * C)[1]#rankMatrix(Phi_L, method = 'qr', tol = 0.05)
  rss = 0; dof = p * (M - 1) + (2 * p - r) * r
  for (m in 1:M) {
    ym = XTS[[m]][,2:(TT[m]+1)]
    xm = XTS[[m]][,1:TT[m]]
    res = ym - tcrossprod(W[m,] * Phi + S[[m]], t(xm))
    # pres_m = chol2inv(cov(t(res)))
    rss = rss + sum(res^2)# * crossprod(pres_m, res))
    dof = dof + sum(S[[m]] != 0)
  }
  ics = c(rss, sum(TT) * p * log(rss / (sum(TT) * p)) + dof * c(2, log(sum(TT * p)), 2 * log(log(sum(TT * p)))))
  ics = c(ics, ics[3] + sum(log((p * M + p^2 - p + M * p^2 - dof + 1:dof) / (1:dof))))
  names(ics) = c('RSS', 'AIC', 'BIC', 'HQC', 'eBIC')
  return(ics)
}

#' @keywords internal
distPhi = function(Phi0, Phi1, Phi_BL, C){
  signs_01 = sign(rowSums(Phi0 * Phi1))
  dist_01 = sqrt(sum((signs_01 * Phi0 - Phi1)^2)) / C
  
  dist_BL = sqrt(sum((Phi1 - Phi_BL)^2)) / C
  return(c(dist_01, dist_BL))
}

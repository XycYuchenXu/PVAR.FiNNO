#' @keywords internal
objfun = function(GK, XTS, eta, Phi_BL, Phi, W, S, Gamma, rho, M = nrow(W),
                  p = nrow(Phi), TT = dim(XTS)[3] - 1){
  dPhi = Phi - Phi_BL
  s = rho * (.5 * sum(dPhi^2) + sum(dPhi * Gamma))
  
  for (m in 1:M) {
    Gm = GK[[1]][m,,]; Km = GK[[2]][m,,]
    
    Am = W[m,] * Phi + S[[m]]
    s = s + .5 * sum(crossprod(Am) * Gm) + .5 * sum(XTS[m,,2:(TT+1)]^2) / TT -
      sum(Am * Km) + eta * sum(abs(S[[m]]))
  }
  return(s)
  
}

#' @keywords internal
IC_PVAR = function(XTS, W, S, Phi, C = 1, TT = dim(XTS)[3] - 1,
                    M = nrow(W), p = ncol(W)){

  r = which(cumsum(svd(Phi)$d) >= 0.95 * C)[1]#rankMatrix(Phi_L, method = 'qr', tol = 0.05)
  rss = 0; dof = p * (M - 1) + (2 * p - r) * r
  for (m in 1:M) {
    ym = XTS[m,,2:(TT+1)]
    xm = XTS[m,,1:TT]
    res = ym - tcrossprod(W[m,] * Phi + S[[m]], t(xm))
    # pres_m = chol2inv(cov(t(res)))
    rss = rss + sum(res^2)# * crossprod(pres_m, res))
    dof = dof + sum(S[[m]] != 0)
  }
  ics = c(rss, TT * M * log(rss / (TT * M)) + dof * c(2, log(TT * M), 2 * log(log(TT * M))))
  names(ics) = c('RSS', 'AIC', 'BIC', 'HQC')
  return(ics)
}

#' @keywords internal
distPhi = function(Phi0, Phi1, C){
  signs = sign(rowSums(Phi0 * Phi1))
  return(sqrt(sum((signs * Phi0 - Phi1)^2)) / C)
}

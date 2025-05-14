#' @keywords internal
compGK = function(XTS, M = length(XTS), p = nrow(XTS[[1]]), TT = sapply(XTS, ncol) - 1){
  G = array(0, dim = c(M, p, p)); K = array(0, dim = c(M, p, p))
  # G = vector('list', M)
  # K = vector('list', M)
  for (m in 1:M) {
    G[m,,] = tcrossprod(XTS[[m]][,1:TT[m]]) / TT[m]
    K[m,,] = tcrossprod(XTS[[m]][,1:TT[m]], XTS[[m]][,2:(TT[m]+1)]) / TT[m]
  }
  return(list(G, K))
}

#' @keywords internal
initPhi = function(GK, C, M = dim(GK$G)[1], p = dim(GK$G)[2]){
  Phi = matrix(0, p, p)
  Phis = array(0, dim = c(M, p, p))

  for (m in 1:M) {
    Phis[m,,] = crossprod(MASS::ginv(GK[[1]][m,,]), GK[[2]][m,,])
  }

  for (i in 1:p) {
    Phi[i,] = svd(Phis[,i,])$v[,1]
  }
  Phi = C * Phi / sum(svd(Phi)$d)
  return(Phi)
}

#' @keywords internal
updateWSm = function(Gm, Km, Phi, Wm, Sm, eta, p = nrow(Phi)){
  wsm = matrix(0, p, p+1)
  gram_mat = matrix(0, p+1, p+1)
  gram_mat[2:(p+1), 2:(p+1)] = Gm
  
  for (i in 1:p) {
    gram_mat[1, 2:(p+1)] = crossprod(Phi[i,], Gm)
    gram_mat[2:(p+1), 1] = gram_mat[1, 2:(p+1)]
    gram_mat[1,1] = crossprod(Phi[i,], gram_mat[1, 2:(p+1)])
    
    Ki = c(crossprod(Phi[i,], Km[,i]), Km[,i])
    beta = c(Wm[i], Sm[i,])
    wsm[i,] = lasso_regression(gram_mat, Ki, beta, lambda = eta, penalty_factors = c(0, rep(1, p)))
  }
  return(wsm)
}

#' @keywords internal
updateWS = function(GK, Phi, eta, M = dim(GK[[1]])[1], p = nrow(Phi), WS = NULL){
  if (is.null(WS)) {
    W = matrix(0, M, p); S = rep(list(matrix(0, p, p)), M)
  } else {
    W = WS$W; S = WS$S
  }
  for (m in 1:M) {
    wsm = updateWSm(GK[[1]][m,,], GK[[2]][m,,], Phi, W[m,], S[[m]], eta, p)
    W[m,] = wsm[,1]
    S[[m]] = wsm[,-1]
  }
  return(list(W=W, S=S))
}

#' @keywords internal
refineWSm = function(Gm, Km, Phi, Sm, p = nrow(Phi)){
  wsm = matrix(0, p, p+1)
  
  for (i in 1:p) {
    Si = Sm[i,]
    supp_ind = which(Si != 0)
    Si_l0 = length(supp_ind)
    
    if (Si_l0 == 0) {
      wsm[i,1] = crossprod(Phi[i,], Km[,i]) / tcrossprod(Phi[i,], crossprod(Phi[i,], Gm))
    } else {
      Gi = matrix(0, 1 + Si_l0, 1 + Si_l0)
      Gi[2:(Si_l0 + 1), 2:(Si_l0 + 1)] = Gm[supp_ind, supp_ind]
      Gi[1, 2:(Si_l0 + 1)] = crossprod(Phi[i,], Gm[,supp_ind])
      Gi[2:(Si_l0 + 1), 1] = Gi[1, 2:(Si_l0 + 1)]
      Gi[1, 1] = tcrossprod(Phi[i,], crossprod(Phi[i,], Gm))
      Ki = c(crossprod(Phi[i,], Km[,i]), Km[supp_ind,i])
      
      wsm[i,c(1, supp_ind+1)] = crossprod(MASS::ginv(Gi), Ki)
    }
    
  }
  return(wsm)
}

#' @keywords internal
refineWS = function(GK, Phi, S, M = dim(GK[[1]])[1], p = nrow(Phi)){
  W = matrix(0, M, p)
  for (m in 1:M) {
    wsm = refineWSm(GK[[1]][m,,], GK[[2]][m,,], Phi, S[[m]], p)
    W[m,] = wsm[,1]
    S[[m]] = wsm[,-1]
  }
  return(list(W=W, S=S))
}

#' @keywords internal
updatePhi = function(W, S, GK, rho, Phi_BL, Gamma, M = nrow(W), p = ncol(W)){
  Phi = matrix(0, p, p)
  
  if (is.null(rho)) {
    rho = max(sqrt(colSums(W^2 * apply(GK[[1]], 1, function(x) sum(diag(x))))))
  }
  
  B = rho * (Phi_BL - Gamma)
  for (m in 1:M) {
    B = B + W[m,] * (t(GK[[2]][m,,]) - tcrossprod(S[[m]], GK[[1]][m,,]))
  }
  
  for (i in 1:p) {
    A = rho * diag(p)
    for (m in 1:M) {
      A = A + W[m,i]^2 * GK[[1]][m,,]
    }
    Phi[i,] = crossprod(B[i,], MASS::ginv(A))
  }
  return(list(Phi = Phi, rho = rho))
}

#' @keywords internal
#' @importFrom irlba irlba
P_L = function(Phi, r, C){
  p = nrow(Phi)
  if (r < p / 4) {svdPhiGam = irlba(Phi, nu = r, nv = r)} else {svdPhiGam = svd(Phi, nu = r, nv = r)}
  U = svdPhiGam$u
  V = svdPhiGam$v
  D = svdPhiGam$d[1:r]
  
  cum_rem_d = (cumsum(D) - C) / (1:r)
  rr = max(which(D - cum_rem_d >= C / r^2))
  
  lambda = D[1:rr] - cum_rem_d[rr]
  if (rr < r) {lambda[(rr+1):r] = 0}
  return(crossprod(t(U), lambda * t(V)))
}

#' @keywords internal
P_B = function(Phi){
  rowL = sqrt(rowSums(Phi^2))
  Phi_B = Phi * mean(rowL) / rowL
  return(Phi_B)
}

#' @keywords internal
updatePhi_BL = function(Phi_n_Gamma, Phi_BL, rho_p, r, C, p = nrow(Phi_n_Gamma),
                        err=1e-4, maxiter=30){
  Phi_0 = (Phi_n_Gamma + rho_p * Phi_BL) / (1 + rho_p)
  Phi_L = Phi_BL; Phi_B = Phi_BL; Gamma_BL = matrix(0, p, p)
  obj_old = sum((Phi_L - Phi_0)^2); diff_L_old = Inf
  for (i in 1:maxiter) {
    temp = P_L((Phi_0 + Phi_B + Gamma_BL)/2, r, C)
    diff_L_new = sum((temp - Phi_L)^2)
    Phi_L = temp
    Phi_B = P_B(Phi_L - Gamma_BL)
    Gamma_BL = Gamma_BL + Phi_B - Phi_L
    obj_new = sum((Phi_L - Phi_0)^2)
    if (abs(diff_L_new - diff_L_old) < err || abs(obj_new - obj_old) < err) {
      break
    }
    obj_old = obj_new; diff_L_old = diff_L_new
  }
  return(Phi_L)
}

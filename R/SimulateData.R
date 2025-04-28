#' Simulate PVAR Parameters
#'
#' Aims to simulate Panel VAR data parameters, including the transition matrices of the form \eqn{A_m = W_m \Phi + S_m} with all these components, and covariance matrices of the innovations being scalar matrices \eqn{\Sigma_m = \sigma_m^2 I_p}.
#'
#' If there exist mixture patterns in the simulated parameters, they are arranged in the following order: 1. Clusters; 2. Isolates; 3. Singular S; Singular W.
#'
#' @param M Size of the panel, i.e., number of entities.
#' @param p Dimension of the PVAR, i.e., number of variables.
#' @param r Rank of the low-rank component.
#' @param s The average fraction of nonzero elements in the sparse components of the coefficient matrices.
#' @param C The specified nuclear norm. Default is \code{C = sqrt(p * r).}
#' @param G.W Number of rescaling groups in the panel, excluding singular ones and isolates. Default is \code{G.W = M - sg_w - sg_s - isolate}.
#' @param G.S Sparse group assignment in the panel, excluding \code{sg_s} singular ones. Default is \code{G.S = 1:(M - sg_s)}.
#' @param isolate Integer or fraction (\code{isolate = round(isolate * M)}), the number of rescaling isolates in the panel model.
#' @param sg_w Integer or fraction (\code{sg_w = round(sg_w * M)}), the number of entities with singular W (i.e., purely sparse entities) in the panel.
#' @param sg_s Integer or fraction (\code{sg_s = round(sg_s * M)}), the number of entities with singular S (i.e., purely low-rank entities) in the panel.
#' @param GW.sd The standard deviation of the rescaling effects \code{W} in the same group when \code{GW.sd >= 0}. Default is \code{GW.sd = 0} meaning exact equality in the group. If \code{GW.sd < 0}, then \code{-GW.sd} is the ratio of standard deviation versus the mean magnitude of \code{W} in the same group.
#' @param GS.sd The standard deviation of the sparse components \code{S} in the same group when \code{GS.sd >= 0}. Default is \code{GS.sd = 0} meaning exact equality in the group. If \code{GS.sd < 0}, then \code{-GS.sd} is the ratio of standard deviation versus the mean magnitude of \code{S} in the same group.
#' @param GS.frac The density ratio of the sparse perturbation versus mean in the same group. Default is \code{GS.frac = 0} meaning no extra non-zero elements.
#'
#' @return A named list:\itemize{
#' \item \code{Coef}: A named list of coefficients:\itemize{
#' \item \code{A}: An \code{M} x \code{p} x \code{p} array of \code{M} coefficient matrices of size \code{p} x \code{p}.
#' \item \code{Sigma}: A length-\code{M} vector of the scalar variances \eqn{\sigma_m}, and the covariance matrices are then constructed as \eqn{\Sigma_m = \sigma_m^2 I_p}.
#' \item \code{W}: An \code{M} x \code{p} matrix of rescaling effects, with rows corresponding to entities.
#' \item \code{Phi}: The shared \code{p} x \code{p} low rank basis.
#' \item \code{S}: A length-\code{M} list of \code{p} x \code{p} sparse matrices.
#' }
#' \item \code{r, M, p, s}: the same as the input, for bookkeeping purpose.
#' }
#'
#' @import foreach
#' @import progressr
#' @import Matrix
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom stats rpois
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples data = simuPar(5, 10, 3, 0.02)
simuPar = function(M, p, r, s, C = sqrt(p * r), G.W = NULL, G.S = NULL, isolate = 0,
                   sg_w = 0, sg_s = 0, GW.sd = 0, GS.sd = 0, GS.frac = 0){
  L = matrix(rnorm(p * p), nrow = p)
  L = expm(L - t(L))[,1:r]
  R = matrix(rnorm(p * p), nrow = p)
  R = expm(R - t(R))[,1:r]

  Lambda = exp(runif(r))
  Phi = as.matrix(crossprod(t(L), t(R) * Lambda))
  Phi = Phi / sqrt(rowSums(Phi^2))
  Phi = C * Phi / sum(svd(Phi)$d)
  maxPhi_p = max(abs(Phi)) * sqrt(p)
  
  if (GW.sd > .4 / sqrt(sum(Phi[1,]^2))) {
    cat('The input "GW.sd" will likely lead to unstationary VAR. Try a smaller value.')
    return()
  }
  if (GS.sd > .6 / sqrt(s * p)) {
    cat('The input "GS.sd" will likely lead to unstationary VAR. Try a smaller value.')
    return()
  }
  
  Am = array(0, dim = c(M, p, p))
  Wm = matrix(0, M, p)
  Sm = vector('list', M)
  Zm = sqrt(1/rgamma(M, shape = 5, scale = 5))

  grs = c()
  lab = c()
  sg_s = max(0, sg_s)
  if (sg_s > 0) {
    if (sg_s < 1 && sg_s > 0) {sg_s = M * sg_s}
    sg_s = round(sg_s)
    if (sg_s > M) {sg_s = 0; cat('Invalid sg_s parameter.')}
    if (sg_s > 0) {grs = c(rep(1, sg_s), grs); lab = c(rep('s', sg_s), lab)}
  }
  sg_w = max(0, sg_w)
  if (sg_w > 0) {
    if (sg_w < 1 && sg_w > 0) {sg_w = M * sg_w}
    sg_w = round(sg_w)
    if (sg_w > M - sg_s) {sg_w = 0; cat('Invalid sg_w parameter.')}
    if (sg_w > 0) {grs = c(sg_w, grs); lab = c('w', lab)}
  }
  isolate = max(0, isolate)
  if (isolate > 0) {
    if (isolate < 1 && isolate > 0) {isolate = M * isolate}
    isolate = round(isolate)
    if (isolate > M - sg_w - sg_s) {isolate = 0; cat('Invalid isolate parameter.')}
    if (isolate > 0) {grs = c(rep(1, isolate), grs); lab = c(rep('o', isolate), lab)}
  }
  M_cur = M - sum(grs)
  if (!is.null(G.W)) {
    grs = c(rep(M_cur %/% G.W, G.W) + c(rep(1, M_cur %% G.W), rep(0, G.W - M_cur %% G.W)), grs)
    lab = c(paste('c', 1:G.W, sep = ''), lab)
  } else {
    grs = c(rep(1, M_cur), grs)
    lab = c(paste('c', 1:M_cur, sep = ''), lab)
  }
  G.W = length(grs)
  
  if (is.null(G.S)) {G.S = 1:(M - sg_s)} else {G.S = as.integer(factor(G.S[1:(M - sg_s)]))}
  Sms = vector('list', length = length(unique(G.S)))
  if (GS.sd < 0) {sd.ws = rep(0, length(Sms))}
  for (ss in 1:length(Sms)) {
    sm = rsparsematrix(p, p, nnz = max(1, min(rpois(1,s*p^2), 1.2*s*p^2))) * maxPhi_p
    dsm = which(abs(diag(sm)) >= .5 * maxPhi_p)
    if (length(dsm) > 0) {
      diag(sm)[dsm] = runif(length(dsm), -0.5, 0.5) * maxPhi_p
    }
    Sms[[ss]] = sm
    if (GS.sd < 0) {sd.ws[ss] = - GS.sd * sqrt(mean(sm@x^2))}
  }
  
  cur = 0
  Wm_g = matrix(exp(runif(G.W * p, -2, 1)) * sample(c(-1, 1), G.W * p, replace = T), p)
  for (g in 1:G.W) {
    gr = grs[g]
    for (m in (cur+1):(cur+gr)) {
      if (lab[g] != 's') {
        sm = Sms[[G.S[m]]]
        if (GS.sd < 0) {
          sm@x = sm@x + rnorm(length(sm@i), sd = sd.ws[G.S[m]])
          sm = sm + rsparsematrix(p, p, density = s * GS.frac) * sd.ws[G.S[m]]
        }
      } else {
        sm = rsparsematrix(p, p, 0)
      }
      Sm[[m]] = sm
      Am[m,,] = as.matrix(sm)

      wm = Wm_g[,g]
      if (lab[g] != 'w') {
        if (GW.sd < 0) {
          wm = wm + rnorm(p, sd = - GW.sd * sqrt(mean(wm^2)))
        }
        Wm[m,] = wm
      }
      Am[m,,] = Am[m,,] + Phi * wm
    }
    cur = cur + gr
  }
  
  omega = apply(Am, 1, function(x){max(abs(eigen(x)$values))})
  rescale = runif(1, .8, .95) / max(omega)
  Wm = Wm * rescale; Am = Am * rescale
  Sm = lapply(Sm, function(x) x * rescale)
  
  if (GW.sd > 0 || GS.sd > 0) {
    repeated_simu = c()
    for (m in 1:M) {
      simu = TRUE
      while (simu) {
        am = Am[m,,]; wm = rep(0,p); sm = rsparsematrix(p,p,0)
        if (GW.sd > 0 && startsWith(lab[which(m <= cumsum(grs))[1]], 'c')) {
          wm = rnorm(p, sd = GW.sd)
          am = am + Phi * wm
        }
        
        if (GS.sd > 0 && m <= M - sg_s && sum(G.S == G.S[m]) > 1) {
          sm_temp = Sm[[m]]
          sm = sm_temp
          sm@x = rnorm(length(sm@i), sd = GS.sd)
          sm = sm + rsparsematrix(p, p, density = s * GS.frac) * GS.sd
          dsm = which(abs(diag(sm)) >= .8)
          if (length(dsm) > 0) {
            diag(sm)[dsm] = runif(length(dsm), -0.5, 0.5) * max(1 - max(abs(diag(am)[dsm])), GS.sd)
          }
          am = am + as.matrix(sm)
        }
        if (max(abs(eigen(am)$values)) < 1) {
          simu = FALSE
          Am[m,,] = am; Wm[m,] = Wm[m,] + wm; Sm[[m]] = Sm[[m]] + sm
        } else {repeated_simu = c(repeated_simu, m)}
      }
    }
    if (length(repeated_simu) > 0) {
      warning(
        paste('The Perturbation for the entity/entities',
              paste(repeated_simu, collapse = ', '),
              'is/are re-simulated to ensure a stationary VAR. The noise distribution may be biased.')
      )
    }
  }
  
  return(list(Coef = list(A = Am, Sigma = Zm, W = Wm, Phi = Phi, S = Sm),
              r = r, M = M, p = p, s = s))
}

#' Simulate PVAR parameters and time series data
#'
#' One can either input all the coefficients \code{(M, p, r, s)} to start from sampling the coefficient matrices, or it also supports inputting sampled coefficient matrices from \code{simuPar} for only time series simulation. In addition, if the input parameters \code{(M, p, r, s)} are of the vector format, then every combination of the setup is simulated.
#'
#' @param N Number of replicates of simulated time series per setting.
#' @param TT The length of time series, default a scalar \code{TT = p * r * 2}. The user can input an integer or a lengths-\code{M} vector of integers. If the length of the input integer vector is less than \code{M}, the first number is recycled.
#' @param Pars The PVAR parameters generated from \code{simuPar}.
#' @param M Size of the panel, i.e., number of entities; can be a vector.
#' @param p Dimension of the PVAR, i.e., number of variables; can be a vector.
#' @param r Rank of the low-rank component; can be a vector.
#' @param s The average fraction of nonzero elements in the sparse components of the coefficient matrices; can be a vector
#' @param C The specified nuclear norm. Default is \code{C = sqrt(p * r)}.
#' @param G.W Number of rescaling groups in the panel, excluding singular ones and isolates. Default is \code{G.W = M - sg_w - sg_s - isolate}.
#' @param G.S Sparse group assignment in the panel, excluding \code{sg_s} singular ones. Default is \code{G.S = 1:(M - sg_s)}.
#' @param isolate Integer or fraction, the number of rescaling isolates in the panel model.
#' @param sg_w Integer or fraction (\code{sg_w = round(sg_w * M)}), the number of entities with singular W (i.e., purely sparse entities) in the panel.
#' @param sg_s Integer or fraction (\code{sg_s = round(sg_s * M)}), the number of entities with singular S (i.e., purely low-rank entities) in the panel.
#' @param GW.sd The standard deviation of the rescaling effects \code{W} in the same group when \code{GW.sd >= 0}. Default is \code{GW.sd = 0} meaning exact equality in the group. If \code{GW.sd < 0}, then \code{-GW.sd} is the ratio of standard deviation versus the mean magnitude of \code{W} in the same group.
#' @param GS.sd The standard deviation of the sparse components \code{S} in the same group when \code{GS.sd >= 0}. Default is \code{GS.sd = 0} meaning exact equality in the group. If \code{GS.sd < 0}, then \code{-GS.sd} is the ratio of standard deviation versus the mean magnitude of \code{S} in the same group.
#' @param GS.frac The density ratio of the sparse perturbation versus mean in the same group. Default is \code{GS.frac = 0} meaning no extra non-zero elements.
#' @param seed Random seed.
#' @param prl If \code{is.numeric(prl)} and \code{prl >= 1}, then its rounded integer is treated as the number of cores for parallel simulation. By default \code{prl = NULL} and simulations are generated sequentially.
#'
#' @return A named list of PVAR parameters and data. In particular denote \code{nM = length(M)}, \code{np = length(p)}, \code{nr = length(r)} and \code{ns = length(s)}, then\itemize{
#' \item \code{Data} A list of panel time series data with dimension \code{nM} x \code{np} x \code{nr} x \code{ns} corresponding to the number of different parameter combinations. Each entry in the list is a length-\code{N} list of named lists\itemize{
#' \item \code{XTS}: A length-\code{M} list of time series, with the \code{m}-th data matrix having size \code{p} x \code{(TT_m + 1)}.
#' \item \code{ind}: The index of the replicate in the simulation, ranging from 1 to \code{N}.
#' \item \code{M, p, r, s, TT}: the same as the input, for bookkeeping purpose.
#' }
#' \item \code{Pars} is a list of simulated PVAR coefficients with the same dimension as \code{Data}, corresponding to all the combinations of the parameters. Each entry is an output object returned by the function \code{simuPar}.
#' }
#'
#' @import foreach
#' @import future
#' @import doFuture
#' @import progressr
#' @importFrom Rdpack reprompt
#' @export
#'
#' @examples simuDP(1, M = 5, p = 10, r = 3, s = 0.02, prl = 2)
simuDP = function(N, TT = NULL, Pars = NULL, M = NULL, p = NULL, r = NULL, s = NULL, C = sqrt(p * r),
                  G.W = NULL, G.S = NULL, isolate = 0, sg_w = 0, sg_s = 0,
                  GW.sd = 0, GS.sd = 0, GS.frac = 0, seed = NULL, prl = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  if (is.null(Pars)) {
    size_par = c(length(M), length(p), length(r), length(s))
    M.i = p.i = r.i = s.i = NULL
    Pars = foreach(M.i = M, p.i = p, r.i = r, s.i = s) %do%
      {
        simuPar(M.i, p.i, r.i, s.i, C = C, G.W = G.W, G.S = G.S,
                isolate = isolate, sg_w = sg_w, sg_s = sg_s,
                GW.sd = GW.sd, GS.sd = GS.sd, GS.frac = GS.frac)
      }
    dim(Pars) = size_par
  }

  if (N == 1) {prl = FALSE}
  Par = n = NULL
  if (is.numeric(prl) && prl >= 1) {
    prl = prl %/% 1
    registerDoFuture()
    plan(multisession, workers = prl)

    pb = progressor(N * length(Pars))
    Data <<- foreach(Par = Pars, .inorder = TRUE, .combine = list, .multicombine = TRUE) %:%
      foreach(n = 1:N, .combine = c, .inorder = FALSE, .options.future = list(scheduling = FALSE)) %dopar% {
        M = Par$M; p = Par$p; r = Par$r; s = Par$s
        result = list(sampleData(Par$Coef, M, p, r, s, TT, n)); pb()
        return(result)
      }

    plan(sequential)

  } else {
    Data = foreach(Par = Pars, .inorder = TRUE, .combine = list, .multicombine = TRUE) %:%
      foreach(n = 1:N, .combine = c, .inorder = FALSE) %do% {
        M = Par$M; p = Par$p; r = Par$r; s = Par$s
        return(list(sampleData(Par$Coef, M, p, r, s, TT, n)))
      }
  }
  if (length(Pars) > 1) {dim(Data) = size_par}
  while (length(Data) == 1) {Data = Data[[1]]}
  while (length(Pars) == 1) {Pars = Pars[[1]]}
  return(list(Data = Data, Pars = Pars))
}

#' Simulate PVAR time series data
#'
#' @param Coef The PVAR coefficients, the first element of the output list generated from \code{simuPar}.
#' @param M Size of the panel, i.e., number of entities.
#' @param p Dimension of the PVAR, i.e., number of variables.
#' @param r Rank of the low-rank component.
#' @param s The average fraction of nonzero elements in the sparse components of the coefficient matrices.
#' @param TT The length of time series, default a scalar \code{TT = p * r * 2}. The user can input an integer or a lengths-\code{M} vector of integers. If the length of the input integer vector is less than \code{M}, the first number is recycled.
#' @param n The index of the replicate of the simulated time series.
#'
#' @return \code{Data} A list of panel time series data with dimension \code{nM} x \code{np} x \code{nr} x \code{ns} corresponding to the number of different parameter combinations. Each entry in the list is a length-\code{N} list of named lists\itemize{
#' \item \code{XTS}: A length-\code{M} list of time series, with the \code{m}-th data matrix having size \code{p} x \code{(TT_m + 1)}.
#' \item \code{ind}: The index of the replicate in the simulation, ranging from 1 to \code{N}.
#' \item \code{M, p, r, s, TT}: the same as the input, for bookkeeping purpose.
#' }
#' @export
#' 
#' @examples sampleData(simuPar(5, 10, 3, 0.02)$Coef, 5, 10, 3, 0.02)
sampleData = function(Coef, M = dim(Coef$A)[1], p = dim(Coef$A)[2], r = NULL, s = NULL, TT = NULL, n = 1){
  XTS = vector('list', length = M)
  if (is.null(TT)) {TT = rep(p*r*2, M)}
  else if (length(TT) < M) {TT = rep(TT[1], M)}
  else {TT = TT[1:M]}
  for (m in 1:M) {
    Zm = Coef$Sigma[m]
    Xm = matrix(0, p, TT[m] + 1)
    x0 = rnorm(p)
    Atemp = Coef$A[m,,]
    for (t in 1:(TT[m]+1)) {
      x0 = tcrossprod(Atemp, t(x0)) + Zm * rnorm(p, sd = 1)
      Xm[,t] = x0
    }
    XTS[[m]] = Xm
  }
  return(list(XTS = XTS, M = M, p = p, r = r, s = s, TT = TT, ind = n))
}

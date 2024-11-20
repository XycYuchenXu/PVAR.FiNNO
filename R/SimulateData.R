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
#' @param C The specified nuclear norm.
#' @param G Number of groups in the panel.
#' @param isolate Integer or fraction (\code{isolate = round(isolate * M)}), the number of isolates in the panel model.
#' @param sg_w Integer or fraction (\code{sg_w = round(sg_w * M)}), the number of entities with singular W (i.e., purely sparse entities) in the panel.
#' @param sg_s Integer or fraction (\code{sg_s = round(sg_s * M)}), the number of entities with singular S (i.e., purely low-rank entities) in the panel.
#' @param G.sd The ratio of standard deviation versus the norm of the rescaling effects \code{W} in the same group. Default is \code{G.sd = 0} meaning exact equality in the group.
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
simuPar = function(M, p, r, s, C = 1, G = NULL, isolate = 0, sg_w = 0, sg_s = 0, G.sd = 0){
  L = matrix(rnorm(p * p), nrow = p)
  L = expm(L - t(L))[,1:r]
  R = matrix(rnorm(p * p), nrow = p)
  R = expm(R - t(R))[,1:r]

  Lambda = exp(runif(r))
  Phi = as.matrix(crossprod(t(L), t(R) * Lambda))
  Phi = Phi / sqrt(rowSums(Phi^2))
  Phi = C * Phi / sum(svd(Phi)$d)

  maxPhi_p = max(Phi) * sqrt(p)

  Am = array(0, dim = c(M, p, p))
  Wm = matrix(0, M, p)
  # Sm = array(0, dim = c(M, p, p))
  # Am = vector('list', M)
  # Wm = vector('list', M)
  Sm = vector('list', M)
  Zm = sqrt(1/rgamma(M, shape = 5, scale = 5))

  grs = c()
  lab = c()
  sg_w = max(0, sg_w)
  if (sg_w > 0) {
    if (sg_w < 1 && sg_w > 0) {sg_w = M * sg_w}
    sg_w = round(sg_w)
    if (sg_w > M) {sg_w = 0; cat('Invalid sg_w parameter.')}
    if (sg_w > 0) {grs = c(sg_w, grs); lab = c('w', lab)}
  }
  sg_s = max(0, sg_s)
  if (sg_s > 0) {
    if (sg_s < 1 && sg_s > 0) {sg_s = M * sg_s}
    sg_s = round(sg_s)
    if (sg_s > M - sg_w) {sg_s = 0; cat('Invalid sg_w parameter.')}
    if (sg_s > 0) {grs = c(sg_s, grs); lab = c('s', lab)}
  }
  isolate = max(0, isolate)
  if (isolate > 0) {
    if (isolate < 1 && isolate > 0) {isolate = M * isolate}
    isolate = round(isolate)
    if (isolate > M - sg_w - sg_s) {isolate = 0; cat('Invalid isolate parameter.')}
    if (isolate > 0) {grs = c(rep(1, isolate), grs); lab = c(rep('o', isolate), lab)}
  }
  M_cur = M - sum(grs)
  if (!is.null(G)) {
    grs = c(rep(M_cur %/% G, G) + c(rep(1, M_cur %% G), rep(0, G - M_cur %% G)), grs)
    lab = c(paste('c', 1:G, sep = ''), lab)
  } else {
    grs = c(rep(1, M_cur), grs)
    lab = c(paste('c', 1:M_cur, sep = ''), lab)
  }
  G = length(grs)

  size_K = max(round(log2(G)), sample(G, 1))
  K = sample(size_K, p, replace = TRUE)

  cur = 0
  for (g in 1:G) {
    gr = grs[g]
    omega = Inf


    wg = runif(p, 1/2, 1)
    sd.wg = G.sd * sqrt(sum(wg^2))
    if (G < M) {wg = wg * sample(c(-1, 1), size_K, replace = TRUE)[K]}

    for (m in (cur+1):(cur+gr)) {
      
      if (lab[g] != 's') {
        sm = rsparsematrix(p, p, nnz = min(rpois(1,s*p*p), 2*s*p*p)) * maxPhi_p
      } else {
        sm = matrix(0, p, p)
      }
      
      if (lab[g] == 'w') {
        Am[m,,] = as.matrix(sm)
        Sm[[m]] = sm
      } else {
        Wm[m,] = wg
        if (sd.wg > 0) {Wm[m,] = Wm[m,] + rnorm(p, sd = sd.wg)}
        Am[m,,] = as.matrix(Phi + sm) * Wm[m,]
        Sm[[m]] = sm * Wm[m,]
      }
      omega = min(omega, 1 / max(abs(eigen(Am[m,,])$values)))
    }

    rescale = runif(1, min = omega/2, omega)# * ifelse(lab[g] == 's', 1/5, 1)
    for (m in (cur+1):(cur+gr)) {
      Wm[m,] = Wm[m,] * ifelse(lab[g] == 'w', 0, rescale)
      Sm[[m]] = Sm[[m]] * rescale
      Am[m,,] = Am[m,,] * rescale
    }
    cur = cur + gr
  }
  return(list(Coef = list(A = Am, Sigma = Zm, W = Wm, Phi = Phi, S = Sm),
              r = r, M = M, p = p, s = s))
}

#' Simulate PVAR parameters and time series data
#'
#' One can either input all the coefficients \code{(M, p, r, s)} to start from sampling the coefficient matrices, or it also supports inputting sampled coefficient matrices from \code{simuPar} for only time series simulation. In addition, if the input parameters \code{(M, p, r, s)} are of the vector format, then every combination of the setup is simulated.
#'
#' @param N Number of replicates of simulated time series per setting.
#' @param TT The length of time series, default a scalar \code{TT = p * r * 2}. The user can input an integer or a lengths-\code{M} vector of integers.
#' @param Pars The PVAR parameters generated from \code{simuPar}.
#' @param M Size of the panel, i.e., number of entities; can be a vector.
#' @param p Dimension of the PVAR, i.e., number of variables; can be a vector.
#' @param r Rank of the low-rank component; can be a vector.
#' @param s The average fraction of nonzero elements in the sparse components of the coefficient matrices; can be a vector
#' @param C The specified nuclear norm.
#' @param G Number of groups in the panel.
#' @param isolate Integer or fraction, the number of isolates in the panel model.
#' @param sg_w Integer or fraction (\code{sg_w = round(sg_w * M)}), the number of entities with singular W (i.e., purely sparse entities) in the panel.
#' @param sg_s Integer or fraction (\code{sg_s = round(sg_s * M)}), the number of entities with singular S (i.e., purely low-rank entities) in the panel.
#' @param G.sd The ratio of standard deviation versus the norm of the rescaling effects \code{W} in the same group. Default is \code{G.sd = 0} meaning exact equality in the group.
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
simuDP = function(N, TT = NULL, Pars = NULL, M = NULL, p = NULL, r = NULL, s = NULL, C = 1,
                  G = NULL, isolate = 0, sg_w = 0, sg_s = 0, G.sd = 0, seed = NULL, prl = NULL){
  if (!is.null(seed)) {set.seed(seed)}
  if (is.null(Pars)) {
    size_par = c(length(M), length(p), length(r), length(s))
    M.i = p.i = r.i = s.i = NULL
    Pars = foreach(M.i = M, p.i = p, r.i = r, s.i = s) %do%
      {
        simuPar(M.i, p.i, r.i, s.i, C = C, G = G,
                isolate = isolate, sg_w = sg_w, sg_s = sg_s, G.sd = G.sd)
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
#' @param TT The length of time series, default a scalar \code{TT = p * r * 2}. The user can input an integer or a lengths-\code{M} vector of integers. If the length of the input integer vector is not \code{M}, the first number is recycled.
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
  else if (length(TT) != M) {TT = rep(TT[1], M)}
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

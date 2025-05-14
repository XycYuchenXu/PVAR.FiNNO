#' ADMM Algorithm for Panel Vector Auto-Regression Model
#'
#' Consider a multivariate Panel VAR model of the form
#' \deqn{X_t^m = A_m X_{t-1}^m + \epsilon_t^m}
#' with the coefficient matrices \eqn{A_m} have the structure
#' \deqn{A_m = W_m \Phi + S_m,}
#' where \eqn{W_m} are diagonal matrices, \eqn{S_m} are sparse matrices, and \eqn{\Phi}
#' is the low-rank basis shared by the whole panel. The optimization problem is constructed under the ADMM framework, with the objective function
#' \deqn{G(W, S, \Phi, \Phi_c, \Gamma; \eta, \rho) = \sum_{m=1}^M (\frac{1}{2T_m} \|Y_m - (W_m \Phi + S_m) X_m\|_F^2 + \eta \|S_m\|_1)\\ + \rho \|\Phi - \Phi_c\|_F^2 + \rho \langle \Gamma, \Phi - \Phi_c \rangle}
#' subject to fixed nuclear-norm constraint \eqn{\|\Phi_c\|_* = \ell} and row-balanced constraint \eqn{\|e_i' \Phi_c\| = \|e_j' \Phi_c\|} for any \eqn{i \ne j}. Note here \eqn{Y_m = (X_1^m, \dots, X_T^m)} and \eqn{X_m = (X_0^m, \dots, X_{T-1}^m)} are stacked time series data of dimension \eqn{p \times T_m}.
#'
#' While constructing the objective function, note that \eqn{\Phi} is set to be flexible, and all the low-rank constraints are imposed through specifying the support domain of the augmented variable \eqn{\Phi_c}. It may be shown that the optimization within the fixed nuclear norm space can automatically induce sparsity in the singular values, as a simplex vertex of the support domain is likely to be hit. For more details, see (CITE).
#'
#' @param XTS A length-\code{M} list of time series, where the \code{m}-th data matrix has size \code{p} x \code{(TT_m + 1)} with potentially unbalanced number of observations \code{TT_m} from \code{TT}.
#' @param r The upper bound for the rank of the estimated low-rank component.
#' @param eta The penalty coefficient for the sparse component.
#' @param TT The effective lengths of the time series, i.e., the notations \eqn{T_m} in description above.
#' @param M The number of subjects in the panel.
#' @param p The number of variables in the multivariate time series data.
#' @param C The specified nuclear norm for the low-rank component.
#' @param rho The step size coefficient in the ADMM setting.
#' @param maxiter Maximum number of iterations.
#' @param miniter Minimum number of iterations.
#' @param err Tolerance error to determine convergence.
#' @param pb The progress bar object defined from the \code{'progressr'} package to track the estimation progress. By default \code{pb = NULL}, no progress bar is inputted and progress is instead tracked by \code{message()} if \code{verbose = TRUE}.
#' @param verbose Logical, whether verbal messages should be printed for progress tracking. Effective only when \code{pb} is \code{NULL}.
#' @param Phi_BL The initializer for the augmented variable \eqn{\Phi_c}.
#' @param Phi The initializer for \eqn{\Phi}.
#' @param Gamma The initializer for the dual variable (Lagrange multiplier) \eqn{\Gamma}.
#' @param WS The list of initializers for the rescaling \code{W} (a matrix of size \code{M} x \code{p}) and the sparse components \code{S} (a length-\code{M} list of \code{p} x \code{p} matrices).
#' @param bulk If \code{verbose = TRUE}, the number of iterations between permanent progress tracking messages.
#' @param perupdate If \code{verbose = TRUE}, the number of iterations between live updates about progress tracking.
#' @param kappa The proximal step size coefficient for the subproblem of \eqn{\Phi_c}, i.e., \code{Phi_BL}.
#' @param normalize Logical, whether the time series data should be normalized. Default is True.
#' @param adap_rho Logical, whether the step size coefficient \code{rho} should be adaptively updated.
#'
#' @return A named list of estimators and some metrics:\itemize{
#' \item \code{Phi}: the \code{p} x \code{p} estimator of \eqn{\Phi} in the objective function above, not necessarily low-rank.
#' \item \code{W}: the \code{M} x \code{p} estimated rescaling effects, with each row corresponding to \eqn{W_m}.
#' \item \code{S}: a length-\code{M} list of sparse components, with the entries being \eqn{S_m} of dimension \code{p} x \code{p}.
#' \item \code{Phi_BL}: the \code{p} x \code{p} estimator of \eqn{\Phi_c} in the objective function above, that is strictly constrained with low-rankness.
#' \item \code{Gamma}: the \code{p} x \code{p} estimator of the dual variable matrix \eqn{\Gamma} in the objective function above.
#' \item \code{status}: 1 if convergence criteria is met within \code{maxiter} iterations, otherwise 0.
#' \item \code{iters}: if \code{status == 1}, it is the number of iterations until reaching convergence; otherwise \code{maxiter} (when \code{status == 0}).
#' \item \code{traj}: a vector of length \code{iters}, the trajectory of the evaluation of the objective function along the sequence of estimators.
#' \item \code{rele}: a vector of length \code{iters}, the trajectory of the relative errors of the estimator \eqn{\Phi} between consecutive iterations.
#' \item \code{ics}: the vector of RSS/AIC/BIC/HQC/eBIC evaluation at the final output estimators.
#' \item \code{time}: computation time for the whole estimation process.
#' \item \code{eta, C, rho}: the same as the input, for bookkeeping purpose.
#' }
#'
#' @import Rcpp
#' @importFrom Rdpack reprompt
#' @importFrom stats sd
#' @export
#'
#' @examples DP = simuDP(1, M = 5, p = 10, r = 3, s = 0.02)
#' est = PVAR_ADMM(DP$Data$XTS, 5, 0.01, maxiter = 10, miniter = 3)
PVAR_ADMM = function(XTS, r, eta, TT = sapply(XTS, ncol) - 1, M = length(XTS), p = nrow(XTS[[1]]),
                     C = sqrt(p * r), rho = M / 10, maxiter = 1e4, miniter = 200, err = 1e-5,
                     pb = NULL, verbose = FALSE, Phi_BL = NULL, Phi = NULL, Gamma = NULL, WS = NULL,
                     bulk = 1, perupdate = 1, kappa = NULL, normalize = T, adap_rho = F){
  tm = proc.time()[3]
  status = 0
  
  if (adap_rho) {
    const_cnt = 0
  }
  
  if (normalize) {
    XTS = lapply(1:M, function(x) XTS[[x]] - rowMeans(XTS[[x]]))
    XTS = lapply(1:M, function(x) XTS[[x]] / sqrt(mean(XTS[[x]]^2)))
  }

  if (is.null(Gamma)) Gamma = matrix(0, p, p)
  if (is.null(kappa)) kappa = M/rho

  traj = Inf; rele = c()

  GK = compGK(XTS, M, p, TT)
  if (is.null(Phi)) {
    Phi = initPhi(GK, C, M, p)
  }

  if (is.null(Phi_BL)) {
    Phi_BL = updatePhi_BL(Phi + Gamma, Phi, kappa, r, C, p)
  }
  
  for (i in 1:maxiter) {
    WS = updateWS(GK, Phi, eta, M, p, WS)

    Phi_BL = updatePhi_BL(Phi + Gamma, Phi_BL, kappa, r, C, p)

    Phi0 = Phi
    Phi = updatePhi(WS$W, WS$S, GK, rho, Phi_BL, Gamma, M, p)

    Gamma = Gamma + Phi - Phi_BL
    traj = c(traj,
             objfun(GK, XTS, eta, Phi_BL, Phi, WS$W, WS$S, Gamma, rho, M, p, TT))
    dist_Phi = distPhi(Phi0, Phi, Phi_BL, C)
    rele = c(rele, max(dist_Phi))
    
    if (adap_rho) {
      if (const_cnt <= 500) {
        if (const_cnt >= 50 && dist_Phi[1] > 10 * dist_Phi[2]) {
          rho = rho / 2
          kappa = kappa * 2
          const_cnt = 0
        } else if (const_cnt >= 50 && dist_Phi[1] < dist_Phi[2] / 10) {
          rho = rho * 2
          kappa = kappa / 2
          const_cnt = 0
        } else {
          const_cnt = const_cnt + 1
        }
      } else {
        adap_rho = F
      }
    }

    if (!is.null(pb)) {
      if (i %% perupdate == 0){
        pb(sprintf('e:%.3f, i:%d, G:%.4f, E:%.3f, R:%.1f', eta, i, traj[i+1], 10000 * rele[i], rho),
           class = if (i %% bulk == 0) "sticky")
      }
    } else if (verbose) {
      if (i %% bulk == 0) {
        message(sprintf('\r%d, %.4f, %.3f', i, traj[i+1], 100 * rele[i]), appendLF = T)
      } else if (i %% perupdate == 0) {
        message(sprintf('\r%d, %.4f, %.3f', i, traj[i+1], 100 * rele[i]), appendLF = F)
      }
    }
    if (i >= miniter) {
      cond1 = (abs(rele[i]) < err) && (abs(rele[i-1]) < err)
      # cond2 = (abs(traj[i+1] - traj[i]) < err) && (abs(traj[i] - traj[i-1]) < err)
      if (cond1) {status = 1; break}
    }
  }
  # W = WS$W
  # S = refineS(GK, WS$S, W, Phi, M, p)
  WS = refineWS(GK, Phi, WS$S, M, p)
  ics = IC_PVAR(XTS, WS$W, WS$S, Phi, C, TT, M, p)
  tm = as.numeric(proc.time()[3] - tm)
  if (!is.null(pb)) {
    pb(sprintf('Done! e:%.3f i:%d, G:%.4f, E:%.3f, R:%.1f', eta, i, traj[i+1], 10000 * rele[i], rho),
       amount = maxiter %/% perupdate - i %/% perupdate, class = 'sticky')
  }
  return(list(Phi = Phi, W = WS$W, S = WS$S, Phi_BL = Phi_BL, Gamma = Gamma,
              eta = eta, traj = traj[-1], rele = rele, C = C, rho = rho,
              ics = ics, iters = i, status = status, time = tm))
}

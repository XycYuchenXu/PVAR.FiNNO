#' Solving Lasso Regression
#' 
#' For a Lasso regression problem of the format
#' \deqn{f(\beta) = \frac{1}{2} \beta' A \beta - \beta' \alpha + \lambda \|\beta\|_1,}
#' this is the function to solve it using the input of \eqn{A} and \eqn{\alpha}.
#'
#' @param xtx The quadratic coefficient \eqn{A}, usually corresponding to the data gram matrix \eqn{X X' / T} in Lasso.
#' @param xty The linear coefficient \eqn{\alpha}, usually corresponding to \eqn{X Y' / T} in Lasso.
#' @param lambda The penalty coefficient.
#' @param penalty_factors Whether the variables is penalized, and its penalty scale with respect to \code{lambda}.
#' @param beta0 The initialization for the estimator.
#' @param max_iter Maximum number of iterations.
#' @param tolerance Convergence criterion.
#'
#' @return The estimated sparse vector \code{beta}.
#' @export
#'
#' @examples X = matrix(rnorm(40), 4, 10); Y = rnorm(10)
#' lasso_gram(tcrossprod(X), tcrossprod(X, t(Y)), 0.05, rep(1, 4))
lasso_gram = function(xtx, xty, lambda, penalty_factors,
                      beta0 = NULL, max_iter = 1000, tolerance = 1e-6){
  p = nrow(xtx)
  if (is.null(beta0)) {
    beta0 = rep(0, p)
  }
  
  return(lasso_regression(xtx, xty, beta0, lambda, penalty_factors, max_iter, tolerance))
}
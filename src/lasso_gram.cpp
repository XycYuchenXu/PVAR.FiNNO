// src/lasso_gram.cpp
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Soft-thresholding function with scaled penalty
double soft_threshold(double z, double gamma, double gram_jj) {
  if (z > gamma) {
    return (z - gamma) / gram_jj;
  } else if (z < -gamma) {
    return (z + gamma) / gram_jj;
  } else {
    return 0.0;
  }
}

// Lasso regression function with a single penalty and scaling factors
// [[Rcpp::export]]
NumericVector lasso_regression(NumericMatrix gram_matrix, NumericVector xy, double lambda, NumericVector penalty_factors, int max_iter = 1000, double tolerance = 1e-6) {
  int p = gram_matrix.ncol();
  NumericVector beta(p);
  NumericVector beta_old(p);
  
  // Coordinate descent
  for (int iter = 0; iter < max_iter; iter++) {
    double max_change = 0.0;
    for (int j = 0; j < p; j++) {
      // Save the old coefficient
      beta_old[j] = beta[j];
      
      // Calculate the partial residual
      double diff = xy[j];
      for (int k = 0; k < p; k++) {
        if (k != j) {
          diff -= gram_matrix(j, k) * beta[k];
        }
      }
      
      // Update beta using the soft-thresholding rule with scaled penalties
      double scaled_penalty = lambda * penalty_factors[j];
      beta[j] = soft_threshold(diff, scaled_penalty, gram_matrix(j, j));
      
      // Update the maximum change
      double change = std::fabs(beta[j] - beta_old[j]);
      if (change > max_change) {
        max_change = change;
      }
    }
    
    // Check for convergence
    if (max_change < tolerance) {
      // Rcpp::Rcout << "Converged after " << iter + 1 << " iterations." << std::endl;
      break;
    }
  }
  return beta;
}

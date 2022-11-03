functions {
/**
  * Compute the normal kernel of a multivariate normal density
  *   on the log scale: `-0.5 * delta' * (LL')^{-1} * delta`,
  *   where L is the lower Cholesky factor of the covariance matrix,
  *   delta = `y - E(y)`, and y is the data vector.
  *   The expression `delta' * (LL')^{-1} * delta` factors into
  *   the inner product of `L^{-1} * delta` with itself. Below,
  *   `mdivide_left_tri_low(cholesky, delta)` computes `L^{-1} * delta`,
  *   and then sum(columns_dot_self()) computes the inner product.
  * @return Real scalar with the density evaluated at delta.
  * @param cholesky Lower-triangular square matrix with one row per
  *   dimension of the multivariate normal. Cholesky factor of the
  *   covariance matrix.
  * @param delta Matrix with one row per dimension of the multivariate normal
  *   and one column per data record.
  */
  real log_normal_kernel(matrix cholesky, matrix delta) {
    real out;
    out = cols(delta) * sum(log(diagonal(cholesky)));
    out += 0.5 * sum(columns_dot_self(mdivide_left_tri_low(cholesky, delta)));
    return -out;
  }
/**
  * Impute missing data as model parameters as described at
  * https://mc-stan.org/docs/stan-users-guide/missing-data.html.
  * missing data points are treated as parameters.
  * @return Matrix of imputed data.
  * @param y Matrix of non-imputed data.
  * @param y_missing Vector of model parameters representing
  *   missing data points to be imputed by the model.
  * @param missing Integer matrix with dimensions equal to those of y.
  *   Element missing[i,j] should be 1 if y[i,j] is missing and 0 otherwise.
  */
  matrix impute(matrix y, vector y_missing, array[,] int missing) {
    int m = rows(y);
    int n = cols(y);
    int index = 0;
    matrix[m,n] out;
    for (j in 1:n) {
      for (i in 1:m) {
        if (missing[i,j] == 1) {
          index += 1;
          out[i,j] = y_missing[index];
        } else {
          out[i,j] = y[i,j];
        }
      }
    }
    return out;
  }
}
data {
  int<lower=0> n_patient;
  int<lower=0> n_visit;
  int<lower=0> n_beta;
  int<lower=0> n_missing;
  array[n_visit, n_patient] int<lower=0, upper=1> missing;
  real<lower=0> s_beta;
  real<lower=0> s_sigma;
  real<lower=0> s_lambda;
 /* The data records (columns) in the model matrix x are ordered
  * by visit within patient, i.e. same as the row order of
  * tidyr::expand_grid(patient = seq_len(n_patient), visit = seq_len(n_visit))
  */
  matrix[n_beta, n_patient * n_visit] x;
  matrix[n_visit, n_patient] y;
}
parameters {
  vector[n_missing] y_missing;
  row_vector[n_beta] beta;
  vector<lower=0,upper=s_sigma>[n_visit] sigma;
  cholesky_factor_corr[n_visit] lambda;
}
model {
  matrix[n_visit, n_patient] y_imputed;
  matrix[n_visit, n_patient] mu;
  matrix[n_visit, n_visit] cholesky = diag_pre_multiply(sigma, lambda);
  y_imputed = impute(y, y_missing, missing);
  mu = to_matrix(beta * x, n_visit, n_patient);
  target += log_normal_kernel(cholesky, y_imputed - mu);
  beta ~ normal(0, s_beta);
  sigma ~ uniform(0, s_sigma);
  lambda ~ lkj_corr_cholesky(s_lambda);
}

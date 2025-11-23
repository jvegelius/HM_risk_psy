library(tidyverse)
library(xtable)
library(mvtnorm)
library(maxLik)
library(lavaan)

############################################################
# Simulation script (N = 100)
#
# Purpose:
#   1) Specify true parameters for the hierarchical model
#   2) Generate simulated data sets
#   3) Fit the hierarchical model (HM) with est_fun()
#   4) Fit a two-step factor-score (IE) approach
#   5) Compare parameter accuracy and partial correlations
#
# Output:
#   - Est_obj_sim_list: list of fitted HM objects for each replication
#   - Multivariate_est_list: list of IE regression models
#   - Summary tables of bias, RMSE, and partial correlations
############################################################


###########################
# 1. Input lists and simulation settings
###########################

Info     <- list()   # Info for generating data
Truepar  <- list()   # True parameter values used to generate data
Par_ini  <- list()   # Initial values of parameters for estimation
Estpat   <- list()   # Parameter restriction patterns
Control  <- list()   # Optimization control parameters

# Number of subjects and number of blocks of trials
Info$n              <- 100
Info$n_block_trials <- 8


###########################
# 2. True parameter values (used for data generation)
###########################

# True parameters are chosen to be similar to those in the empirical example

Truepar <- list()

# Covariance of subject-level decision parameters xi (3×3)
Truepar$phi <- matrix(
  c(0.74, -0.2, 0.47,
    -0.2, 0.11, -0.11,
    0.47, -0.11, 0.49),
  nrow = 3
)

# Covariance of latent factors eta (3×3, diagonal)
Truepar$psi <- diag(c(0.44, 0.65, 0.66))

# Covariance of observed symptom variables y (17×17, diagonal)
Truepar$omega_y <- diag(c(
  0.22, 0.54, 0.83, 0.73, 1.66, 0.84, 0.71, 0.46, 0.57,
  0.79, 0.78, 0.54, 0.2, 0.48, 0.44, 0.09, 0.7
))

# Mean of decision parameters xi
Truepar$kappa <- c(-1.37, 0.65, -0.86)

# Regression matrix gamma (3×3) linking xi to eta
Truepar$gamma <- matrix(
  c(-0.57, -0.47,  0.33,
    -2.15, -1.50,  1.53,
    0,     0,     0),
  nrow = 3
)

# Factor loading matrix lambda_y (17 observed variables × 3 factors)
Truepar$lambda_y <- matrix(0, nrow = 17, ncol = 3)
Truepar$lambda_y[1:11, 1]  <- c(1, 0.7, 0.94, 0.86, 0.8, 0.74, 0.82, 0.92, 0.94, 1.11, 1.03)
Truepar$lambda_y[11:14, 2] <- c(1, 0.67, 0.39, 0.49)
Truepar$lambda_y[16:17, 3] <- c(1, 0.51)

# Means of y
Truepar$mu_y <- c(
  3.03, 2.87, 2.72, 1.85, 2.58, 1.97, 2.43, 2.21, 2.62,
  2.61, 2.53, 2.28, 1.73, 1.6, 1.48, 2.5, 1.72
)

# Example alternative gamma values (first column) – overwriting above
Truepar$gamma[1:3, 1] <- c(-1, -0.9, 0.8)


###########################
# 3. Initial values for estimation (Par_ini)
###########################

Par_ini$phi      <- diag(0.01, 3)
Par_ini$psi      <- diag(0.01, 3)
Par_ini$omega_y  <- diag(rep(0.2, 17))
Par_ini$kappa    <- rep(0, 3)
Par_ini$gamma    <- matrix(rep(0, 9), nrow = 3)

# Initial factor loading pattern (same structure as in empirical example)
Par_ini$lambda_y <- matrix(
  c(
    1,          rep(0, 10),
    rep(0, 6),
    rep(0, 11),
    c(1, 0, 0, 0),
    c(0, 0),
    rep(0, 15),
    c(1, 0)
  ),
  ncol = 3
)

Par_ini$mu_y <- rep(0, 17)


###########################
# 4. Estimation patterns based on Par_ini
#    (which elements are free vs fixed)
###########################

# Gamma pattern: NA = free, 0 = fixed
Estpat$gamma <- matrix(
  c(NA, NA, NA,
    NA, NA, NA,
    0,  0,  0),
  nrow = 3
)

# Initialize lambda_y pattern: 17×3 zeros
Estpat$lambda_y <- matrix(0, nrow = 17, ncol = 3)

# For any initial loading equal to 1, mark that position for further handling
for (r in 1:length(Estpat$lambda_y)) {
  if (Par_ini$lambda_y[r] == 1) {
    Estpat$lambda_y[r] <- 1
  }
}

# Specify which loadings are free to estimate (NA):
# - Items 2–11 on factor 1
# - Items 13–15 on factor 2
# - Item 17 on factor 3
Estpat$lambda_y[c(2:11), 1] <- NA
Estpat$lambda_y[c(13:15), 2] <- NA
Estpat$lambda_y[17, 3]       <- NA

# Pattern for phi (3×3, 6 free parameters: 3 variances + 3 covariances)
Estpat$phi <- array(0, dim = c(3, 3, 6))
for (r in 1:3) {
  Estpat$phi[r, r, r] <- 1
}
m <- r + 1
for (r in 1:2) {
  for (k in (r + 1):3) {
    Estpat$phi[r, k, m] <- Estpat$phi[k, r, m] <- 1
    m <- m + 1
  }
}

# Pattern for psi (only diagonal variances free)
Estpat$psi <- array(0, dim = c(3, 3, 3))
for (r in 1:3) {
  Estpat$psi[r, r, r] <- 1
}

# Pattern for omega_y (only variances free)
Estpat$omega_y <- array(0, dim = c(17, 17, 17))
for (r in 1:17) {
  Estpat$omega_y[r, r, r] <- 1
}

# Control settings for optimization
Control$maxit       <- 1000    # max outer iterations
Control$tol         <- 1e-4    # convergence tolerance
Control$maxit_optim <- 50      # max iterations per optim() call
Control$var_g       <- 500     # variance for Gaussian regularization


###########################
# 5. Estimation patterns based on Truepar (for simulations)
#    (overwrites Estpat with a pattern aligned to the true λ)
###########################

Estpat <- list()

# Same gamma pattern as above
Estpat$gamma <- matrix(
  c(NA, NA, NA,
    NA, NA, NA,
    0,  0,  0),
  nrow = 3
)

# lambda_y pattern determined directly from Truepar$lambda_y:
#  - 1  → constrained to 1
#  - 0  → constrained to 0
#  - other values → free (NA)
Estpat$lambda_y <- matrix(0, nrow = 17, ncol = 3)
for (r in 1:length(Estpat$lambda_y)) {
  if (Truepar$lambda_y[r] == 1) {
    Estpat$lambda_y[r] <- 1
  } else if (Truepar$lambda_y[r] == 0) {
    Estpat$lambda_y[r] <- 0
  } else {
    Estpat$lambda_y[r] <- NA
  }
}

# Same covariance patterns as before
Estpat$phi <- array(0, dim = c(3, 3, 6))
for (r in 1:3) {
  Estpat$phi[r, r, r] <- 1
}
m <- r + 1
for (r in 1:2) {
  for (k in (r + 1):3) {
    Estpat$phi[r, k, m] <- Estpat$phi[k, r, m] <- 1
    m <- m + 1
  }
}

Estpat$psi <- array(0, dim = c(3, 3, 3))
for (r in 1:3) {
  Estpat$psi[r, r, r] <- 1
}

Estpat$omega_y <- array(0, dim = c(17, 17, 17))
for (r in 1:17) {
  Estpat$omega_y[r, r, r] <- 1
}

# Same control settings (duplicated but harmless)
Control$maxit       <- 1000
Control$tol         <- 1e-4
Control$maxit_optim <- 50
Control$var_g       <- 500


#########################################
# 6. Generate data and estimate the model (simulation loop)
#########################################

Elapsed_time_vec      <- c()      # store elapsed time per replication
Est_obj_sim_list      <- list()   # store HM estimation results
Multivariate_est_list <- list()   # store IE regression results
N_samples             <- 100      # number of Monte Carlo replications

Gamma_est_mat <- NULL            # collect estimated gamma (first 6 entries)

for (r in 1:N_samples) {
  # Different seed per replication for reproducibility
  set.seed(20231230 + 100 * r)
  
  # Generate one dataset according to Truepar
  Data_obj <- gen_data_fun(Info, Truepar)
  Gambling_data_sim_list <- Data_obj$gambling_data_sim
  Y_mat                  <- Data_obj$y_mat
  
  # Fit hierarchical model using est_fun()
  # We use Truepar as starting values here (can also use Par_ini)
  try(
    Elapsed_time_vec[r] <-
      system.time(
        Est_obj_sim_list[[r]] <-
          est_fun(Gambling_data_sim_list, Y_mat, Truepar, Estpat, Control)
      )[3]
  )
  
  cat("Simulation ", r, " finished after ",
      round(Elapsed_time_vec[r] / 60, 2), " minutes.\n")
  
  # Store gamma estimates (first 6 elements in vectorized form)
  Gamma_est_mat <- rbind(
    Gamma_est_mat,
    as.vector(Est_obj_sim_list[[r]]$gamma)[1:6]
  )
  
  # Prepare y data for factor analysis (IE approach)
  Df_Y <- as.data.frame(as.matrix(Y_mat))
  
  # CFA model for two latent factors f1 and f2 (using lavaan)
  mla <- 'f1 =~ 1*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11
          f2 =~ 1*V12 + V13 + V14 + V15'
  facanal       <- cfa(mla, data = Df_Y)
  sumfacanal    <- summary(facanal)
  factor_scores <- predict(facanal)
  
  # Multivariate regression of factor scores on xi (IE approach)
  if (!is.null(Est_obj_sim_list[[r]])) {
    mlm <- lm(factor_scores ~ Est_obj_sim_list[[r]]$xi_tilde_mat[ , -3])
    summary(mlm)
    Multivariate_est_list[[r]] <- mlm
  }
}


########################################
# 7. Analyze simulation results: parameter accuracy
##########################################

N_gamma_hm <- length(Est_obj_sim_list[[1]]$par_gamma)  # # of gamma params in HM
N_gamma_fs <- 4                                        # # of gamma params in FS

N_samples <- length(Est_obj_sim_list)

Par_gamma_hm_mat           <- matrix(NA, nrow = N_samples, ncol = N_gamma_hm)
Par_gamma_fs_mat           <- matrix(NA, nrow = N_samples, ncol = N_gamma_fs)
SE_location_gamma_sandwich <- matrix(NA, nrow = N_samples, ncol = N_gamma_hm)
SE_location_gamma_score    <- matrix(NA, nrow = N_samples, ncol = N_gamma_hm)
SE_gamma_sandwich          <- matrix(NA, nrow = N_samples, ncol = N_gamma_hm)
SE_gamma_score             <- matrix(NA, nrow = N_samples, ncol = N_gamma_hm)
SE_gamma_fs                <- matrix(NA, nrow = N_samples, ncol = N_gamma_fs)

N_na_hm  <- 0  # counting failed HM fits
N_na_fs  <- 0  # counting failed IE fits
NA_vec_hm <- NULL
NA_vec_fs <- NULL

# fs_sequence reorders coefficients from the FS regression
fs_sequence <- c(1, 3, 2, 4)

for (r in 1:N_samples) {
  # Collect HM estimates if available
  if (!is.null(Est_obj_sim_list[[r]])) {
    Par_gamma_hm_mat[r, ]            <- Est_obj_sim_list[[r]]$par_gamma
    SE_location_gamma_sandwich[r, ]  <- Est_obj_sim_list[[r]]$se_location_gamma_sandwich
    SE_location_gamma_score[r, ]     <- Est_obj_sim_list[[r]]$se_location_gamma_score
    SE_gamma_sandwich[r, ]           <- Est_obj_sim_list[[r]]$se_gamma_sandwich
    SE_gamma_score[r, ]              <- Est_obj_sim_list[[r]]$se_gamma_score
  } else {
    N_na_hm  <- N_na_hm + 1
    NA_vec_hm <- c(NA_vec_hm, r)
  }
  
  # Collect IE estimates if regression succeeded
  if (!is.null(Multivariate_est_list[[r]])) {
    sumMultivariate_est <- summary(Multivariate_est_list[[r]])
    
    # Extract regression coefficients (excluding intercepts)
    Par_gamma_fs_mat[r, ] <-
      as.vector(Multivariate_est_list[[r]]$coefficients[c(-1, -4), ])[fs_sequence]
    
    # Extract standard errors for same coefficients
    SE_gamma_fs[r, ] <-
      c(sumMultivariate_est$`Response f1`$coefficients[2:3, 2],
        sumMultivariate_est$`Response f2`$coefficients[2:3, 2])[fs_sequence]
  } else {
    N_na_fs  <- N_na_fs + 1
    NA_vec_fs <- c(NA_vec_fs, r)
  }
}

# Number of converged HM samples (ignoring failed runs)
N_conv_samples <- N_samples - N_na_hm

# True gamma parameters (excluding last column which is all zero)
Par_gamma_true     <- as.vector(Truepar$gamma)[-c(7:9)]
Par_gamma_mean_hm  <- colMeans(Par_gamma_hm_mat)
Par_gamma_mean_fs  <- colMeans(Par_gamma_fs_mat)

# Bias of HM and IE estimators
Par_gamma_mean_hm - Par_gamma_true
Par_gamma_mean_fs - Par_gamma_true[1:N_gamma_fs]


# Investigate approximate confidence intervals via Z-scores
Z_gamma_hm_score <- (Par_gamma_hm_mat -
                       t(matrix(rep(Par_gamma_true, N_conv_samples),
                                nrow = N_gamma_hm))) / SE_gamma_score

Z_gamma_fs <- (Par_gamma_fs_mat[, 1:4] -
                 t(matrix(rep(Par_gamma_true[c(1:2, 4:5)], N_conv_samples),
                          nrow = N_gamma_fs))) / SE_gamma_fs

# Compare empirical Z-distributions to standard normal for one parameter
k        <- 4
xi_cont  <- seq(-5, 5, by = 0.01)
plot(xi_cont, dnorm(xi_cont), type = "l", lwd = 3,
     ylim = c(0, 0.5), xlim = c(-9, 9))
lines(density(Z_gamma_hm_score[, k]),  lwd = 2)
lines(density(Z_gamma_fs[, k]),       lwd = 2)


# Summaries: mean, bias, SD, RMSE for gamma (HM vs IE)
k  <- 1
qu <- 0.9
ql <- 0.1

Mean_stimate <- Par_gamma_mean_hm
Bias         <- Par_gamma_mean_hm - Par_gamma_true
SD           <- sqrt(diag(cov(Par_gamma_hm_mat)))
Parameter    <- c("Gamma_11", "Gamma_21", "Gamma_31",
                  "Gamma_12", "Gamma_22", "Gamma_32")
DF_accuracy_results_hm <- data.frame(
  Parameter    = Parameter,
  Trueval      = Par_gamma_true,
  Method       = "HM",
  Mean_stimate = Mean_stimate,
  Bias         = Bias,
  RelBias      = Bias / Par_gamma_true,
  SD           = SD,
  RMSE         = sqrt(SD^2 + Bias^2)
)

Mean_stimate       <- Par_gamma_mean_fs
Par_gamma_true_fs  <- Par_gamma_true[c(1:2, 4:5)]
Bias               <- Par_gamma_mean_fs - Par_gamma_true_fs
SD                 <- sqrt(diag(cov(Par_gamma_fs_mat)))
Parameter          <- c("Gamma_11", "Gamma_21", "Gamma_12", "Gamma_22")
DF_accuracy_results_fs <- data.frame(
  Parameter    = Parameter,
  Trueval      = Par_gamma_true_fs,
  Method       = "IE",
  Mean_stimate = Mean_stimate,
  Bias         = Bias,
  RelBias      = Bias / Par_gamma_true_fs,
  SD           = SD,
  RMSE         = sqrt(SD^2 + Bias^2)
)

DF_accuracy_results <- rbind(DF_accuracy_results_hm, DF_accuracy_results_fs)
DF_accuracy_results[, 4:8] <- round(DF_accuracy_results[, 4:8], digits = 3)
DF_accuracy_results[c(1, 7, 2, 8, 4, 9, 5, 10), ]

xtable(DF_accuracy_results, digits = 3)
xtable(DF_accuracy_results[c(1, 7, 2, 8, 4, 9, 5, 10), ], digits = 3)


################################################
# 8. Correlation structure and partial correlations (true values)
################################################

Gamma_true <- Truepar$gamma[, 1:2]
Phi_true   <- Truepar$phi[1:2, 1:2]
Psi_true   <- Truepar$psi

# Construct true covariance matrix among (eta, xi)
Sigma_true            <- matrix(NA, nrow = 5, ncol = 5)
Sigma_true[1:3, 1:3]  <- Gamma_true %*% Phi_true %*% t(Gamma_true) + Psi_true
Sigma_true[4:5, 1:3]  <- Phi_true %*% t(Gamma_true)
Sigma_true[1:3, 4:5]  <- Gamma_true %*% Phi_true
Sigma_true[4:5, 4:5]  <- Phi_true

SD_theta <- sqrt(diag(Sigma_true))
SD_eta   <- SD_theta[1:3]
SD_xi    <- SD_theta[4:5]

# True correlation matrix
R_true          <- diag(1 / SD_theta) %*% Sigma_true %*% diag(1 / SD_theta)
R_eta_true      <- R_true[1:3, 1:3]
R_xi_true       <- R_true[4:5, 4:5]
R_etaxi_true    <- R_true[1:3, 4:5]

# True partial correlations between eta and xi, controlling for the other xi
R_partial_true <- matrix(NA, nrow = 3, ncol = 2)
for (xi_index in 1:2) {
  for (eta_index in 1:3) {
    rho_xy <- R_etaxi_true[eta_index, xi_index]
    rho_xz <- R_xi_true[xi_index, 3 - xi_index]
    rho_zy <- R_etaxi_true[eta_index, 3 - xi_index]
    
    R_partial_true[eta_index, xi_index] <-
      (rho_xy - rho_xz * rho_zy) /
      (sqrt(1 - rho_xz^2) * sqrt(1 - rho_zy^2))
  }
}

R_partial_true^2


################################################
# 9. Estimated partial correlations from HM and IE
################################################

N_reps <- length(Est_obj_sim_list)

R_partial_hm_mat <- matrix(NA, nrow = N_reps, ncol = 6)
R_partial_fs_mat <- matrix(NA, nrow = N_reps, ncol = 4)

for (r in 1:N_reps) {
  
  ## ---- HM-based partial correlations ----
  Sigma_hm_r <- matrix(NA, nrow = 5, ncol = 5)
  Gamma_hm_r <- Est_obj_sim_list[[r]]$gamma[ , -3]
  Phi_hm_r   <- Est_obj_sim_list[[r]]$phi[1:2, 1:2]
  Psi_hm_r   <- Est_obj_sim_list[[r]]$psi
  
  Sigma_hm_r[1:3, 1:3] <- Gamma_hm_r %*% Phi_hm_r %*% t(Gamma_hm_r) + Psi_hm_r
  Sigma_hm_r[4:5, 1:3] <- Phi_hm_r %*% t(Gamma_hm_r)
  Sigma_hm_r[1:3, 4:5] <- Gamma_hm_r %*% Phi_hm_r
  Sigma_hm_r[4:5, 4:5] <- Phi_hm_r
  
  SD_theta_hm_r <- sqrt(diag(Sigma_hm_r))
  R_hm_r        <- diag(1 / SD_theta_hm_r) %*% Sigma_hm_r %*% diag(1 / SD_theta_hm_r)
  
  R_eta_hm_r   <- R_hm_r[1:3, 1:3]
  R_xi_hm_r    <- R_hm_r[4:5, 4:5]
  R_etaxi_hm_r <- R_hm_r[1:3, 4:5]
  
  R_partial_hm_r <- matrix(NA, nrow = 3, ncol = 2)
  for (xi_index in 1:2) {
    for (eta_index in 1:3) {
      rho_xy <- R_etaxi_hm_r[eta_index, xi_index]
      rho_xz <- R_xi_hm_r[xi_index, 3 - xi_index]
      rho_zy <- R_etaxi_hm_r[eta_index, 3 - xi_index]
      
      R_partial_hm_r[eta_index, xi_index] <-
        (rho_xy - rho_xz * rho_zy) /
        (sqrt(1 - rho_xz^2) * sqrt(1 - rho_zy^2))
    }
  }
  
  ## ---- IE-based partial correlations ----
  Df_Y <- as.data.frame(as.matrix(Est_obj_sim_list[[r]]$y_mat))
  
  mla <- 'f1 =~ 1*V1 + V2 + V3 + V4 +  V5+ V6 + V7 + V8 + V9 + V10 + V11
          f2 =~ 1*V12 + V13 + V14 + V15'
  facanal       <- cfa(mla, data = Df_Y)
  sumfacanal    <- summary(facanal)
  factor_scores <- predict(facanal)
  
  Sigma_fs_r <- cov(cbind(factor_scores,
                          Est_obj_sim_list[[r]]$xi_tilde_mat[ , -3]))
  
  SD_theta_fs_r <- sqrt(diag(Sigma_fs_r))
  R_fs_r        <- diag(1 / SD_theta_fs_r) %*% Sigma_fs_r %*% diag(1 / SD_theta_fs_r)
  
  R_eta_fs_r   <- R_fs_r[1:2, 1:2]
  R_xi_fs_r    <- R_fs_r[3:4, 3:4]
  R_etaxi_fs_r <- R_fs_r[1:2, 3:4]
  
  R_partial_fs_r <- matrix(NA, nrow = 2, ncol = 2)
  for (xi_index in 1:2) {
    for (eta_index in 1:2) {
      rho_xy <- R_etaxi_fs_r[eta_index, xi_index]
      rho_xz <- R_xi_fs_r[xi_index, 3 - xi_index]
      rho_zy <- R_etaxi_fs_r[eta_index, 3 - xi_index]
      
      R_partial_fs_r[eta_index, xi_index] <-
        (rho_xy - rho_xz * rho_zy) /
        (sqrt(1 - rho_xz^2) * sqrt(1 - rho_zy^2))
    }
  }
  
  # Store vectorized partial correlations
  R_partial_hm_mat[r, ] <- as.vector(t(R_partial_hm_r))
  R_partial_fs_mat[r, ] <- as.vector(t(R_partial_fs_r))
  
  print(r)
}

R_partial_true     <- as.vector(t(R_partial_true))
colMeans(R_partial_hm_mat)
colMeans(R_partial_fs_mat)

R_partial_mean_hm <- colMeans(R_partial_hm_mat)
R_partial_mean_fs <- colMeans(R_partial_fs_mat)

# Accuracy summary for partial correlations (HM)
Mean_stimate <- R_partial_mean_hm
Bias         <- R_partial_mean_hm - R_partial_true
SD           <- sqrt(diag(cov(R_partial_hm_mat)))
Parameter    <- c("Gamma_11", "Gamma_12", "Gamma_21",
                  "Gamma_22", "Gamma_31", "Gamma_32")
DF_Rpartial_results_hm <- data.frame(
  Parameter    = Parameter,
  Trueval      = R_partial_true,
  Method       = "HM",
  Mean_stimate = Mean_stimate,
  Bias         = Bias,
  RelBias      = Bias / R_partial_true,
  SD           = SD,
  RMSE         = sqrt(SD^2 + Bias^2)
)

# Accuracy summary for partial correlations (FS)
Mean_stimate        <- c(R_partial_mean_fs, NA, NA)
R_partial_true_fs   <- R_partial_true[c(1:4, NA, NA)]
Bias                <- c(R_partial_mean_fs, NA, NA) - R_partial_true
SD                  <- c(sqrt(diag(cov(R_partial_fs_mat))), NA, NA)
DF_Rpartial_results_fs <- data.frame(
  Parameter    = Parameter,
  Trueval      = R_partial_true_fs,
  Method       = "IE",
  Mean_stimate = Mean_stimate,
  Bias         = Bias,
  RelBias      = Bias / R_partial_true_fs,
  SD           = SD,
  RMSE         = sqrt(SD^2 + Bias^2)
)

# Combine and format results
DF_Rpartial_results <- rbind(DF_Rpartial_results_hm, DF_Rpartial_results_fs)
DF_Rpartial_results[ , c(2, 4:8)] <-
  round(DF_Rpartial_results[ , c(2, 4:8)], digits = 3)

DF_Rpartial_results[c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12), ]

DF_Rpartial_results_N100 <-
  DF_Rpartial_results[c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12), ]

xtable(DF_Rpartial_results_N100[ , -6], digits = 3)

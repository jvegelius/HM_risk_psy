############################################################
# step_by_step_guide.R
#
# Purpose:
#   Example script showing how to:
#     1) Load the example gambling (WOF) data and IDAS data
#     2) Specify starting values and estimation patterns
#     3) Call est_fun() to fit the hierarchical model
#
# Requirements:
#   - Functions_for_hierarchical_WOF_psych.R  (in the same folder)
#   - example_wof_data.csv
#   - example_idas_data.csv
#
# Output:
#   - Est_obj: a list containing parameter estimates and diagnostics
#   - Elapsed_time: run time of the estimation
############################################################

#######
# 1. Load model functions and data
#######

# Load all functions used for hierarchical estimation
source("Functions_for_hierarchical_WOF_psych.R")

# Load the example Wheel-of-Fortune (gambling) dataset.
# Columns are: id, trial, p (probability), amb (ambiguity), gain, y (choice)
DF_wof  <- read.csv2("example_wof_data.csv")

# Load the example IDAS dataset (psychological symptom measures).
# First column is subject id; remaining columns are symptom scores.
DF_idas <- read.csv2("example_idas_data.csv")


#######
# 2. Put trial-wise gambling data into a per-subject list
#######

# Extract unique subject IDs and total number of subjects
IDs <- unique(DF_wof$id)
N   <- length(IDs)

# Gambling_data_list will be a list of length N.
# For each subject j, Gambling_data_list[[j]] is a data.frame with columns:
#   p   – probability,
#   amb – ambiguity level,
#   gain – outcome amount,
#   y   – binary choice (1 = gamble, 0 = safe option)
Gambling_data_list <- list()
for (j in 1:N) {
  DF_wof_j <- DF_wof[DF_wof$id == IDs[j], ]
  Gambling_data_list[[j]] <- data.frame(
    p    = DF_wof_j$p,
    amb  = DF_wof_j$amb,
    gain = DF_wof_j$gain,
    y    = DF_wof_j$y
  )
}


#######
# 3. Specify initial parameter values (Par_ini)
#######

Par_ini <- list()

# phi: covariance matrix of subject-level decision parameters (xi)
# Here: 3×3 diagonal with small variance (0.01) for each dimension
Par_ini$phi      <- diag(0.01, 3)

# psi: covariance matrix of latent factors (eta)
# Also 3×3 diagonal with small variance as starting point
Par_ini$psi      <- diag(0.01, 3)

# omega_y: covariance matrix of observed symptom variables y
# Here: 17×17 diagonal with variance 0.2 for each IDAS subscale
Par_ini$omega_y  <- diag(rep(0.2, 17))

# kappa: mean vector of subject-level decision parameters xi
Par_ini$kappa    <- rep(0, 3)

# gamma: regression matrix linking xi to eta (3×3).
# Initialized to zero (no relationship) and updated by the estimator.
Par_ini$gamma    <- matrix(rep(0, 9), nrow = 3)

# lambda_y: factor loading matrix mapping eta to observed y.
Par_ini$lambda_y <- matrix(
  c(
    1,          rep(0, 10),      # first factor
    rep(0, 6),  rep(0, 11),      # second factor
    c(1, 0, 0, 0),
    c(0, 0),
    rep(0, 15),
    c(1, 0)
  ),
  ncol = 3
)

# mu_y: mean vector for the observed symptom variables y
Par_ini$mu_y     <- rep(0, 17)


#######
# 4. Specify which parameters are free to estimate (Estpat)
#######

Estpat <- list()

# Estpat$gamma: pattern for gamma matrix.
# NA = free parameter; 0 = fixed to zero; Can be fixed to any number.
# (This initial line is overwritten below via explicit setting of lambda_y,
# but is kept here to show how gamma can be constrained if desired.)
Estpat$gamma <- matrix(
  c(NA, NA, NA,
    NA, NA, NA,
    0,  0,  0),
  nrow = 3
)

# Initialize lambda_y pattern as zeros:
# 17 observed variables × 3 latent factors
Estpat$lambda_y <- matrix(0, nrow = 17, ncol = 3)

# For any element of Par_ini$lambda_y that is 1, mark it as a free parameter.
for (r in 1:length(Estpat$lambda_y)) {
  if (Par_ini$lambda_y[r] == 1) {
    Estpat$lambda_y[r] <- 1
  }
}

# Set which loadings are estimated (NA) vs fixed:
#  - Items 2–11 load on factor 1
#  - Items 13–15 load on factor 2
#  - Item 17 loads on factor 3
Estpat$lambda_y[c(2:11), 1] <- NA
Estpat$lambda_y[c(13:15), 2] <- NA
Estpat$lambda_y[17, 3]       <- NA


# Estpat$phi: pattern for the 3×3 covariance matrix phi.
# We represent each free element via a separate 3D slice.
#   - Each slice represents one parameter to be estimated
#   - Ones represent the parameter to be estimated
#   - Several ones per slice means that those are restricted to be the same, e.g., covariances
Estpat$phi <- array(0, dim = c(3, 3, 6))

# Variances (diagonal entries)
for (r in 1:3) {
  Estpat$phi[r, r, r] <- 1
}

# Covariances (off-diagonal entries)
m <- r + 1
for (r in 1:2) {
  for (k in (r + 1):3) {
    Estpat$phi[r, k, m] <- Estpat$phi[k, r, m] <- 1
    m <- m + 1
  }
}

# Estpat$psi: pattern for 3×3 covariance matrix psi.
# Only variances are estimated (diagonal elements).
Estpat$psi <- array(0, dim = c(3, 3, 3))
for (r in 1:3) {
  Estpat$psi[r, r, r] <- 1
}

# Estpat$omega_y: pattern for 17×17 covariance matrix of y.
# Here only diagonal elements (variances) are estimated.
Estpat$omega_y <- array(0, dim = c(17, 17, 17))
for (r in 1:17) {
  Estpat$omega_y[r, r, r] <- 1
}


#######
# 5. Construct Y matrix (symptom data)
#######

# Convert IDAS data (excluding the id column) to a numeric matrix.
# Rows = subjects; columns = symptom scales.
Y_mat <- as.matrix(DF_idas[ , -1])


#######
# 6. Set control options for the optimization
#######

Control <- list()
Control$maxit        <- 1000    # Max number of outer iterations in est_fun()
Control$tol          <- 1e-4    # Convergence tolerance for maximum parameter change
Control$maxit_optim  <- 50      # Max iterations for each inner call to optim()
Control$var_g        <- 500     # Variance for the Gaussian regularization ("super prior")


#######
# 7. Run the hierarchical estimation and record elapsed time
#######

# system.time() returns user, system, and elapsed time.
# We extract the 3rd element (elapsed wall-clock time in seconds).
Elapsed_time <- system.time(
  Est_obj <- est_fun(
    Gambling_data_list,
    Y_mat,
    Par_ini,
    Estpat,
    Control
  )
)[3]

# After running, Est_obj contains:
#  - parameter estimates (gamma, lambda_y, phi, psi, omega_y, kappa, mu_y)
#  - subject-level estimates (xi_tilde_mat)
#  - standard errors and convergence information
# You can inspect it with, e.g.:
#   names(Est_obj)
#   Est_obj$gamma
#   Est_obj$lambda_y
#   Est_obj$converge



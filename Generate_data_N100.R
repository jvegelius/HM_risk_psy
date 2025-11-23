library(tidyverse)
library(xtable)
library(mvtnorm)
library(maxLik)
library(lavaan)

# Load model and data-generation functions (including gen_data_fun())
source("Functions_for_hierarchical_WOF_psych.R")


############################################################
# Generate example data for the hierarchical WOF + IDAS model
#
# Purpose:
#   - Specify true parameters similar to the empirical example
#   - Generate:
#       (1) trial-wise gambling data (Wheel-of-Fortune task)
#       (2) subject-level symptom data (IDAS-like scores)
#   - Save:
#       example_wof_data.csv  (trial-wise choices)
#       example_idas_data.csv (subject-level symptoms)
#
# Usage:
#   Run this script once to create the example CSV files that are
#   used in step_by_step_guide.R.
############################################################


###########################
# 1. Set up model input structures
###########################

Info    <- list()  # Info for generating data (sample size, task design)
Truepar <- list()  # True parameter values for data generation

# Number of subjects
Info$n <- 100

# Number of blocks of trials per subject for the WOF task
Info$n_block_trials <- 8


###########################
# 2. True parameter values (similar to empirical example)
###########################

# These are the parameters that gen_data_fun() will use to generate:
#   - subject-level decision parameters xi
#   - latent symptom factors eta
#   - observed symptom scores y
#   - trial-wise choices in the gambling task

Truepar <- list()

# Covariance matrix of decision parameters xi (3×3)
Truepar$phi <- matrix(
  c(0.74, -0.2,  0.47,
    -0.2, 0.11, -0.11,
    0.47, -0.11, 0.49),
  nrow = 3
)

# Covariance of latent factors eta (3×3, diagonal)
Truepar$psi <- diag(c(0.44, 0.65, 0.66))

# Covariance of observed symptom variables y (17×17, diagonal)
Truepar$omega_y <- diag(c(
  0.22, 0.54, 0.83, 0.73, 1.66, 0.84, 0.71, 0.46, 0.57,
  0.79, 0.78, 0.54, 0.20, 0.48, 0.44, 0.09, 0.70
))

# Mean of decision parameters xi
Truepar$kappa <- c(-1.37, 0.65, -0.86)

# Regression matrix gamma linking xi to eta (3×3)
Truepar$gamma <- matrix(
  c(-0.57, -0.47,  0.33,
    -2.15, -1.50,  1.53,
    0,     0,     0),
  nrow = 3
)

# Factor loadings matrix lambda_y (17 observed variables × 3 factors)
Truepar$lambda_y <- matrix(0, nrow = 17, ncol = 3)
Truepar$lambda_y[1:11, 1]  <- c(1, 0.7, 0.94, 0.86, 0.8, 0.74, 0.82, 0.92, 0.94, 1.11, 1.03)
Truepar$lambda_y[11:14, 2] <- c(1, 0.67, 0.39, 0.49)
Truepar$lambda_y[16:17, 3] <- c(1, 0.51)

# Means of y (IDAS-like symptom scores)
Truepar$mu_y <- c(
  3.03, 2.87, 2.72, 1.85, 2.58, 1.97, 2.43, 2.21, 2.62,
  2.61, 2.53, 2.28, 1.73, 1.60, 1.48, 2.50, 1.72
)

# Optionally override first column of gamma (example values)
Truepar$gamma[1:3, 1] <- c(-1, -0.9, 0.8)


###########################
# 3. Generate a single example dataset (N = 100)
###########################

# Fixed seed for reproducibility of the example files
set.seed(20231231)

# Generate data using the hierarchical model:
#   - gambling_data_sim: list of WOF trials per subject
#   - xi_mat, eta_mat: latent subject-level variables
#   - y_mat: IDAS-like symptom scores (n × 17)
Data_obj <- gen_data_fun(Info, Truepar)

Gambling_data_sim_list <- Data_obj$gambling_data_sim
Y_mat                  <- Data_obj$y_mat

# Create subject IDs 1,…,N
IDs <- 1:Info$n


###########################
# 4. Build example IDAS data frame (subject-level symptoms)
###########################

# DF_idas has:
#   - column "id"   : subject identifier
#   - columns "idas.1", "idas.2", ... "idas.17" (created from the matrix Y_mat)
#     representing 17 symptom scales
DF_idas <- data.frame(id = IDs, idas = Y_mat)


###########################
# 5. Build example WOF data frame (trial-wise choices)
###########################

# DF_wof will contain one row per trial and subject, with:
#   id    – subject ID
#   trial – trial index (1,…, Info$n_block_trials * 10)
#   p     – objective win probability
#   amb   – ambiguity level
#   gain  – gain amount
#   y     – binary choice (1 = gamble, 0 = safe)

DF_wof <- NULL
for (j in 1:Info$n) {
  DF_wof_j <- data.frame(
    id    = IDs[j],
    trial = 1:(Info$n_block_trials * 10),
    p     = Gambling_data_sim_list[[j]]$p,
    amb   = Gambling_data_sim_list[[j]]$amb,
    gain  = Gambling_data_sim_list[[j]]$gain,
    y     = Gambling_data_sim_list[[j]]$y
  )
  DF_wof <- rbind(DF_wof, DF_wof_j)
}

dim(DF_wof)  # Check dimensions: should be (N × n_trials_per_subject) rows

###########################
# 6. Remove all objects except the generated data frames
###########################
rm(list = setdiff(ls(), c("DF_idas", "DF_wof")))

###########################
# 7. Save example data to CSV files
###########################

# Trial-wise WOF data
write.csv2(DF_wof,
           file      = "example_wof_data.csv",
           row.names = FALSE)

# Subject-level IDAS-like data
write.csv2(DF_idas,
           file      = "example_idas_data.csv",
           row.names = FALSE)

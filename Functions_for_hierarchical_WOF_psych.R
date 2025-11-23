library(tidyverse)
library(xtable)
library(mvtnorm)
library(maxLik)
library(lavaan)

#######################################
####################Functions
#######################################

# ------------------------------------------------------------------------------
# prob_gamble_fun()
#
# Computes the probability of choosing to gamble on a risky/ambiguous option.
#
# Arguments:
#   prob  – Objective probability of winning.
#   amb   – Ambiguity level (typically between 0 and 1).
#   gain  – Potential monetary gain.
#   alpha – Subject-specific risk sensitivity parameter.
#   beta  – Subject-specific ambiguity aversion parameter.
#   mu    – Choice noise (inverse temperature) parameter.
#
# Returns:
#   A numeric value between 0 and 1 giving the probability of choosing to gamble.
#
# Notes:
#   - The function computes subjective expected utility for the risky option,
#     compares it to the utility of a safe option, and applies a logistic
#     choice rule.
#   - Extremely small or large probabilities are truncated to avoid numerical
#     underflow or overflow in downstream likelihood calculations.
# ------------------------------------------------------------------------------
prob_gamble_fun <- function(prob,amb,gain,alpha,beta,mu){
  #eu_risky <- (prob-beta*amb/2)*gain^alpha
  eu_risky <- (prob-beta*amb/2)*gain^alpha
  eu_safe <- (0.5)*5^alpha
  prob_gamble <- 1/(1+exp(-(eu_risky-eu_safe)/mu))
  
  if(prob_gamble<1e-6){
    prob_gamble <- 1e-6
  }
  if(prob_gamble>(1-1e-6)){
    prob_gamble <- 1-1e-6
  }
  return(prob_gamble)
}

# ------------------------------------------------------------------------------
# est_fun()
#
# Estimates the hierarchical latent-variable model using a variational-Laplace
# (quadratic) approximation of the marginal log-likelihood.
#
# Arguments:
#   gambling_data_list – A list where each element contains the behavioral
#                        choice data for one subject (probabilities, ambiguity,
#                        gains, and choices).
#
#   y_mat              – An n × k_y matrix of continuous observed variables
#                        (one row per subject).
#
#   par_ini            – A list containing initial values for all model
#                        parameters (phi, psi, gamma, kappa, lambda_y, mu_y,
#                        omega_y).
#
#   estpat             – A list of indicator matrices specifying which elements
#                        of phi, psi, gamma, lambda_y, and omega_y are free
#                        parameters to estimate.
#
#   control            – A list of control settings:
#                           maxit        – Maximum number of outer iterations.
#                           tol          – Convergence tolerance.
#                           maxit_optim  – Max iterations for inner optim()
#                           var_g        – Variance parameter for the Gaussian
#                                          regularization (“super prior”)
#
# Returns:
#   A list containing:
#     par_all         – Final parameter vector.
#     mu_y, lambda_y, gamma, kappa, phi, psi, omega_y – Estimated parameters.
#     xi_tilde_mat    – Mode estimates of subject-specific latent variables.
#     xi_neg_hess_array – Negative Hessian matrices per subject.
#     neg_hess_pos_def– Indicator whether Hessians were positive definite.
#     se_*            – Various standard error estimates (sandwich & score).
#     converge        – TRUE if convergence criteria were met.
#
# Notes:
#   - Alternates between updating gamma, lambda_y, kappa, and covariance
#     matrices using optim().
#   - Computes the Laplace-approximated marginal log-likelihood at each step.
#   - Includes optional numerical gradient and Hessian calculations to compute
#     standard errors.
# ------------------------------------------------------------------------------
est_fun <- function(gambling_data_list,y_mat,par_ini,estpat,control){
  n <- length(gambling_data_list)
  
  #Introduce initial parameters 
  phi_curr <- par_ini$phi
  kappa_curr <- par_ini$kappa
  psi_curr <- par_ini$psi
  omega_y_curr <- par_ini$omega_y
  kappa_curr <- par_ini$kappa
  n_kappa <- length(kappa_curr)
  gamma_curr <- par_ini$gamma
  lambda_y_curr <- par_ini$lambda_y
  mu_y_curr <- par_ini$mu_y
  
  k_xi <- dim(phi_curr)[1]
  k_eta <- dim(psi_curr)[1]
  k_y <- dim(omega_y_curr)[1]
  
  j <- 1
  a_j <- c(0,0,0)
  gambling_data_j <- gambling_data_list[[j]]
  
  negloglik_y_cond_xi_fun <- function(a_j,gambling_data_j){
    value <- 0
    n_trials <- dim(gambling_data_j)[1]
    for(r in 1:n_trials){
      #prob_gamble <-  prob_gamble_fun(gambling_data_j$p[r],gambling_data_j$amb[r],gambling_data_j$gain[r],a_j[1],a_j[2],exp(a_j[3]))
      prob_gamble <-  prob_gamble_fun(gambling_data_j$p[r],gambling_data_j$amb[r],gambling_data_j$gain[r],exp(a_j[1]),a_j[2],exp(a_j[3]))
      value <- value + gambling_data_j$y[r]*log(prob_gamble) + (1-gambling_data_j$y[r])*log(1-prob_gamble)
    }
    value <- value - (1/control$var_g)*n_trials*t(a_j)%*%a_j  #The last term is a regularizing term or a "super prior"
    return(-value)
  }
  
  xi_tilde_mat <- matrix(NA,nrow=n,ncol=k_xi)
  xi_neg_hess_array <- array(NA,dim=c(k_xi,k_xi,n))
  neg_hess_pos_def <- c()
  
  #par(mfrow=c(3,1)) 
  for(j in 1:n){
    obj_fun <- function(xi_j_tmp){
      value <- negloglik_y_cond_xi_fun(xi_j_tmp,gambling_data_list[[j]])
      return(value)
    }
    obj_fun(rep(0,k_xi))
    
    optobj <- optim(rep(0,k_xi),obj_fun)
    #optobj <- optim(c(exp(1),0,0),obj_fun)
    
    
    xi_tilde_mat[j,] <- optobj$par
    
    xi_neg_hess_array[,,j] <- numericHessian(obj_fun,t0=optobj$par)
    
    eigobj <- eigen(xi_neg_hess_array[,,j])
    if(sum(eigobj$values<0)==0){
      neg_hess_pos_def[j] <- 1
    }else{
      neg_hess_pos_def[j] <- 0
      print(cbind(j,"no"))
      xi_neg_hess_array[,,j] <- solve(phi_curr)
      xi_tilde_mat[j,] <- rep(0,k_xi)
    }
    
    # par(mfrow=c(3,1))   
    # if(j>3){
    #   #plot(Data_obj$xi_mat[1:j,1],a_tilde_mat[1:j,1])
    #   plot(density(xi_tilde_mat[1:j,1]))
    #   #abline(0,1)
    #   plot(density(xi_tilde_mat[1:j,2]))
    #   #abline(0,1)
    #   plot(density(xi_tilde_mat[1:j,3]))
    #   #abline(0,1)
    # }
    
  }
  
  
  #colMeans(xi_tilde_mat)
  #cov(xi_tilde_mat)
  #phi_curr
  
  # phi <- phi_curr
  # kappa <- kappa_curr
  log_ev_factor_fun <- function(kappa,phi,xi_tilde_mat,xi_neg_hess_array){
    n <- dim(xi_tilde_mat)[1]
    k_xi <- dim(phi)[1]
    phi_inv <- solve(phi)
    b1 <- matrix(0,nrow=k_xi,ncol=k_xi)
    v1 <- rep(0,k_xi)
    log_ev <- -0.5*n*log(det(phi))
    for(j in 1:n){
      log_ev <- log_ev - 0.5*log(det(phi_inv + xi_neg_hess_array[,,j]))
      log_ev <- log_ev - 0.5*t(xi_tilde_mat[j,])%*%(xi_neg_hess_array[,,j] - xi_neg_hess_array[,,j]%*%solve(phi_inv + xi_neg_hess_array[,,j])%*%xi_neg_hess_array[,,j]    )%*%xi_tilde_mat[j,]
      b1 <- b1 + phi_inv - phi_inv%*%solve(phi_inv + xi_neg_hess_array[,,j])%*%phi_inv
      v1 <- v1 + solve(phi_inv + xi_neg_hess_array[,,j])%*%xi_neg_hess_array[,,j]%*%xi_tilde_mat[j,]
    }
    
    log_ev <- log_ev - 0.5*t(kappa)%*%b1%*%kappa 
    log_ev <- log_ev + t(kappa)%*%phi_inv%*%v1
    neglog_ev <- -log_ev
    return(neglog_ev)
  }
  
  neglog_ev_sem_fun <- function(phi,psi,omega_y,gamma,lambda_y,kappa,mu_y,xi_tilde_mat,xi_neg_hess_array,y_mat){
    n <- dim(y_mat)[1]
    k_xi <- dim(phi)[1]
    k_eta <- dim(psi)[1]
    k_theta <- k_xi + k_eta
    k_y <- dim(y_mat)[2]
    
    phi_inv <- solve(phi)
    psi_inv <- solve(psi)
    omega_y_inv <- solve(omega_y)
    
    lT_om_inv_l <- t(lambda_y)%*%omega_y_inv%*%lambda_y
    lT_om_inv_l__inv <- solve(lT_om_inv_l)
    lT_om_inv <- t(lambda_y)%*%omega_y_inv
    gT_psi_inv_g <- t(gamma)%*%psi_inv%*%gamma
    gT_psi_inv <- t(gamma)%*%psi_inv
    
    lT_om_inv_l_plus_psi_inv <- lT_om_inv_l + psi_inv
    lT_om_inv_l_plus_psi_inv__inv <- solve(lT_om_inv_l_plus_psi_inv)
    
    #Start calculating the log evidence
    log_ev <- -0.5*n*(log(det(omega_y)) + log(det(psi)) + log(det(phi))) 
    
    
    for(j in 1:n){
      
      m_j <- matrix(0,nrow=k_theta,ncol=k_theta)
      m_j[1:k_eta,1:k_eta] <-lT_om_inv_l + psi_inv 
      m_j[(k_eta+1):k_theta,(k_eta+1):k_theta] <- phi_inv + xi_neg_hess_array[,,j] + gT_psi_inv_g
      m_j[(k_eta+1):k_theta,1:k_eta] <- -gT_psi_inv
      m_j[1:k_eta,(k_eta+1):k_theta] <- -t(gT_psi_inv)
      
      m_j_inv <- solve(m_j)
      m_j_inv_11 <- m_j_inv[1:k_eta,1:k_eta]
      m_j_inv_12 <- m_j_inv[(k_eta+1):k_theta,1:k_eta]
      m_j_inv_22 <- m_j_inv[(k_eta+1):k_theta,(k_eta+1):k_theta]
      
      log_ev <- log_ev -0.5*log(det(lT_om_inv_l_plus_psi_inv)*det(-gT_psi_inv%*%lT_om_inv_l_plus_psi_inv__inv%*%t(gT_psi_inv) + phi_inv + xi_neg_hess_array[,,j] + gT_psi_inv_g))
      log_ev <- log_ev - 0.5*t(y_mat[j,]-mu_y)%*%(omega_y_inv - omega_y_inv%*%lambda_y%*%m_j_inv_11%*%t(lambda_y)%*%omega_y_inv)%*%(y_mat[j,] - mu_y)
      log_ev <- log_ev - 0.5*t(xi_tilde_mat[j,] - kappa)%*%(xi_neg_hess_array[,,j] - xi_neg_hess_array[,,j]%*%m_j_inv_22%*%xi_neg_hess_array[,,j])%*%(xi_tilde_mat[j,] - kappa)
      log_ev <- log_ev + t(y_mat[j,]-mu_y)%*%omega_y_inv%*%lambda_y%*%t(m_j_inv_12)%*%xi_neg_hess_array[,,j]%*%(xi_tilde_mat[j,] - kappa)
      #Regularize with respect to alpha_j
      #log_ev <- log_ev - 0.01*t(xi_tilde_mat[j,] - kappa)%*%(xi_tilde_mat[j,] - kappa)
    }
    
    #log_ev <- log_ev - 0.005*n*t(as.vector(gamma))%*%as.vector(gamma)
    log_ev <- log_ev - 0.0001*n*t(as.vector(gamma))%*%as.vector(gamma)
    
    
    neg_log_ev <- -log_ev
    return(neg_log_ev)
  }
  
  neglog_ev_j_sem_fun <- function(phi,psi,omega_y,gamma,lambda_y,kappa,mu_y,xi_tilde_j,xi_neg_hess_j,y_j){
    k_xi <- dim(phi)[1]
    k_eta <- dim(psi)[1]
    k_theta <- k_xi + k_eta
    k_y <- length(y_j)
    
    phi_inv <- solve(phi)
    psi_inv <- solve(psi)
    omega_y_inv <- solve(omega_y)
    
    lT_om_inv_l <- t(lambda_y)%*%omega_y_inv%*%lambda_y
    lT_om_inv_l__inv <- solve(lT_om_inv_l)
    lT_om_inv <- t(lambda_y)%*%omega_y_inv
    gT_psi_inv_g <- t(gamma)%*%psi_inv%*%gamma
    gT_psi_inv <- t(gamma)%*%psi_inv
    
    lT_om_inv_l_plus_psi_inv <- lT_om_inv_l + psi_inv
    lT_om_inv_l_plus_psi_inv__inv <- solve(lT_om_inv_l_plus_psi_inv)
    
    #Start calculating the log evidence
    log_ev <- -0.5*(log(det(omega_y)) + log(det(psi)) + log(det(phi))) 
    
    m_j <- matrix(0,nrow=k_theta,ncol=k_theta)
    m_j[1:k_eta,1:k_eta] <-lT_om_inv_l + psi_inv 
    m_j[(k_eta+1):k_theta,(k_eta+1):k_theta] <- phi_inv + xi_neg_hess_j + gT_psi_inv_g
    m_j[(k_eta+1):k_theta,1:k_eta] <- -gT_psi_inv
    m_j[1:k_eta,(k_eta+1):k_theta] <- -t(gT_psi_inv)
    
    m_j_inv <- solve(m_j)
    m_j_inv_11 <- m_j_inv[1:k_eta,1:k_eta]
    m_j_inv_12 <- m_j_inv[(k_eta+1):k_theta,1:k_eta]
    m_j_inv_22 <- m_j_inv[(k_eta+1):k_theta,(k_eta+1):k_theta]
    
    log_ev <- log_ev -0.5*log(det(lT_om_inv_l_plus_psi_inv)*det(-gT_psi_inv%*%lT_om_inv_l_plus_psi_inv__inv%*%t(gT_psi_inv) + phi_inv + xi_neg_hess_j + gT_psi_inv_g))
    log_ev <- log_ev - 0.5*t(y_j-mu_y)%*%(omega_y_inv - omega_y_inv%*%lambda_y%*%m_j_inv_11%*%t(lambda_y)%*%omega_y_inv)%*%(y_j - mu_y)
    log_ev <- log_ev - 0.5*t(xi_tilde_mat[j,] - kappa)%*%(xi_neg_hess_j - xi_neg_hess_j%*%m_j_inv_22%*%xi_neg_hess_j)%*%(xi_tilde_j - kappa)
    log_ev <- log_ev + t(y_j-mu_y)%*%omega_y_inv%*%lambda_y%*%t(m_j_inv_12)%*%xi_neg_hess_j%*%(xi_tilde_j - kappa)
    
    neg_log_ev <- -log_ev
    return(neg_log_ev)
  }
  
  #Initialize paramter vectors
  
  #mu is calculated directly
  mu_y_curr <- colMeans(y_mat)
  
  #gamma
  gamma_index <- which(is.na(estpat$gamma))
  par_gamma_curr <- gamma_curr[gamma_index]
  n_gamma <- length(gamma_index)
  
  #lambda_y
  lambda_y_index <- which(is.na(estpat$lambda_y))
  par_lambda_y_curr <- lambda_y_curr[is.na(estpat$lambda_y)]
  n_lambda_y <- length(lambda_y_index)
  
  #phi
  n_phi <- dim(estpat$phi)[3]  #Number of free parameters to estimate in Phi
  phi_index <- c()
  for(k in 1:n_phi){
    phi_index[k] <- which(estpat$phi[,,k]==1)[1]
  }
  par_phi_curr <- phi_curr[phi_index]
  
  #psi
  n_psi <- dim(estpat$psi)[3]  #Number of free parameters to estimate in Phi
  psi_index <- c()
  for(k in 1:n_psi){
    psi_index[k] <- which(estpat$psi[,,k]==1)[1]
  }
  par_psi_curr <- phi_curr[psi_index]
  
  #omega_y
  n_omega_y <- dim(estpat$omega_y)[3]  #Number of free parameters to estimate in Phi
  omega_y_index <- c()
  for(k in 1:n_omega_y){
    omega_y_index[k] <- which(estpat$omega_y[,,k]==1)[1]
  }
  par_omega_y_curr <- omega_y_curr[omega_y_index]
  
  par_curr_all <- par_new_all <- c(kappa_curr,par_gamma_curr,par_lambda_y_curr,par_phi_curr,par_psi_curr,par_omega_y_curr)
  n_par_all <- length(par_curr_all)
  
  
  nl_index <- 1
  neglog_ev_vec <- c()
  
  neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_curr,lambda_y_curr,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
  #cat("The log evidence is ",neglog_ev_vec[nl_index]," after initialization.", "\n")
  nl_index <- nl_index + 1
  
  maxparupdate <- 100
  iter <- 1
  #Start optimizing model parameters
  #while loop later...
  
  while(iter<control$maxit & maxparupdate>control$tol){
    
    #Gamma
    obj_fun <- function(par_gamma){
      gamma_tmp <- gamma_curr
      gamma_tmp[gamma_index] <- par_gamma
      value <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_tmp,lambda_y_curr,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      return(value)
    }
    
    optobj_gamma <- optim(par_gamma_curr,obj_fun,control=list(maxit=control$maxit_optim))
    par_gamma_new <- optobj_gamma$par
    gamma_new <- gamma_curr
    gamma_new[gamma_index] <- par_gamma_new
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_new,lambda_y_curr,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating gamma: ",round(par_gamma_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    #Lambda_y
    obj_fun <- function(par_lambda_y){
      lambda_y_tmp <- lambda_y_curr
      lambda_y_tmp[lambda_y_index] <- par_lambda_y
      value <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_new,lambda_y_tmp,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      return(value)
    }
    
    optobj_lambda_y <- optim(par_lambda_y_curr,obj_fun,control=list(maxit=control$maxit_optim))
    par_lambda_y_new <- optobj_lambda_y$par
    lambda_y_new <- lambda_y_curr
    lambda_y_new[lambda_y_index] <- par_lambda_y_new
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_new,lambda_y_new,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating lambda_y :",round(par_lambda_y_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    
    #kappa
    obj_fun <- function(kappa){
      value <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_new,lambda_y_new,kappa,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      return(value)
    }
    
    optobj_kappa <- optim(kappa_curr,obj_fun,control=list(maxit=control$maxit_optim))
    kappa_new <- optobj_kappa$par
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating kappa: ",round(kappa_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    
    #Phi
    obj_fun <- function(par_phi){
      phi_tmp <- phi_curr
      for(k in 1:n_phi){
        index_tmp <- which(estpat$phi[,,k]==1)
        phi_tmp[index_tmp] <- par_phi[k]
      }
      value <- neglog_ev_sem_fun(phi_tmp,psi_curr,omega_y_curr,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      return(value)
    }
    
    optobj_phi <- optim(par_phi_curr,obj_fun,control=list(maxit=control$maxit_optim))
    par_phi_new <- optobj_phi$par
    phi_new <- phi_curr
    for(k in 1:n_phi){
      phi_new[which(estpat$phi[,,k]==1)] <- par_phi_new[k]
    }
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_new,psi_curr,omega_y_curr,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating phi: ",round(par_phi_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    
    #Psi
    #par_psi <- par_psi_curr
    obj_fun <- function(par_psi){
      psi_tmp <- psi_curr
      for(k in 1:n_psi){
        index_tmp <- which(estpat$psi[,,k]==1)
        psi_tmp[index_tmp] <- par_psi[k]
        if(psi_tmp[k]<0.001){
          psi_tmp[k] <- 0.001
        }
      }
      value <- neglog_ev_sem_fun(phi_new,psi_tmp,omega_y_curr,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      
      return(value)
    }
    
    optobj_psi <- optim(par_psi_curr,obj_fun,control=list(maxit=control$maxit_optim))
    par_psi_new <- optobj_psi$par
    
    for(k in 1:n_psi){
      if(par_psi_new[k]<=0){
        par_psi_new[k] <- 0.001
      }
    }
    
    psi_new <- psi_curr
    for(k in 1:n_psi){
      psi_new[which(estpat$psi[,,k]==1)] <- par_psi_new[k]
    }
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_new,psi_new,omega_y_curr,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating psi: ",round(par_psi_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    
    #Omega_y
    obj_fun <- function(par_omega_y){
      omega_y_tmp <- omega_y_curr
      for(k in 1:n_omega_y){
        index_tmp <- which(estpat$omega_y[,,k]==1)
        omega_y_tmp[index_tmp] <- par_omega_y[k]
      }
      value <- neglog_ev_sem_fun(phi_new,psi_new,omega_y_tmp,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
      return(value)
    }
    
    optobj_omega_y <- optim(par_omega_y_curr,obj_fun,control=list(maxit=control$maxit_optim))
    par_omega_y_new <- optobj_omega_y$par
    omega_y_new <- omega_y_curr
    for(k in 1:n_omega_y){
      omega_y_new[which(estpat$omega_y[,,k]==1)] <- par_omega_y_new[k]
    }
    
    neglog_ev_vec[nl_index] <- neglog_ev_sem_fun(phi_new,psi_new,omega_y_new,gamma_new,lambda_y_new,kappa_new,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    #cat("Neg log evidence: ",neglog_ev_vec[nl_index]," Updating omega_y: ",round(par_omega_y_new,digits=2), "\n")
    nl_index <- nl_index + 1
    
    par_new_all <- c(kappa_new,par_gamma_new,par_lambda_y_new,par_phi_new,par_psi_new,par_omega_y_new)
    
    d_par_all <- par_new_all - par_curr_all
    
    maxparupdate <- max(abs(d_par_all))
    
    
    cat("Neg log evidence: ",neglog_ev_vec[nl_index-1]," Iteration:",iter," Max parameter update ",maxparupdate, "\n")
    
    iter <- iter + 1
    kappa_curr <- kappa_new
    gamma_curr <- gamma_new
    lambda_y_curr <- lambda_y_new
    phi_curr <- phi_new
    psi_curr <- psi_new
    omega_y_curr <- omega_y_new
    
    par_gamma_curr <- par_gamma_new
    par_lambda_y_curr <- par_lambda_y_new
    par_phi_curr <- par_phi_new
    par_psi_curr <- par_psi_new
    par_omega_y_curr <- par_omega_y_new
    
    par_curr_all <- par_new_all
    
    # par(mfrow=c(1,1))
    # plot(neglog_ev_vec,type="l")
    
  }
  
  cat(c("Calculating standard errors..."),"\n")
  #########################
  ##############Calculate standard errors
  obj_fun_gradient <- function(par_gamma){
    xi_tilde_j <- xi_tilde_mat[j,]
    xi_neg_hess_j <- xi_neg_hess_array[,,j]
    y_j <- y_mat[j,]
    gamma_tmp <- gamma_curr
    gamma_tmp[gamma_index] <- par_gamma
    value <- neglog_ev_j_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_tmp,lambda_y_curr,kappa_curr,mu_y_curr,xi_tilde_j,xi_neg_hess_j,y_j)
    return(value)
  }
  
  obj_fun_hessian <- function(par_gamma){
    gamma_tmp <- gamma_curr
    gamma_tmp[gamma_index] <- par_gamma
    value <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_tmp,lambda_y_curr,kappa_curr,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    return(value)
  }
  
  
  
  g_gamma_mat <- matrix(0,nrow=n,ncol=n_gamma)
  for(j in 1:n){
    g_gamma_mat[j,] <- numericGradient(obj_fun_gradient,t0=par_gamma_curr)
  }
  
  h_gamma <- numericHessian(obj_fun_hessian,t0=par_gamma_curr)
  h_gamma_inv <- solve(h_gamma)
  
  #Investigate if the Hessian is positive definite
  h_eigen_gamma <- eigen(solve(h_gamma))
  if(sum(h_eigen_gamma$values<0)==0){
    hessian_posdef <- TRUE
    se_gamma_sandwich <- sqrt(diag(h_gamma_inv%*%t(g_gamma_mat)%*%g_gamma_mat%*%h_gamma_inv))
  }else{
    hessian_posdef <- FALSE
    #cat(c("Hessian with respect to gamma is not positive definite."),"\n")
    se_gamma_sandwich <- NA
  }
  
  if(iter<control$maxit & maxparupdate){
    converge <- TRUE
  }else{
    converge <- FALSE
  }
  
  se_gamma_score <- sqrt(diag(solve(t(g_gamma_mat)%*%g_gamma_mat)))
  
  
  #Calculate standard erreors with respect to all location parameters, i.e., gamma, lambda_y, tau and mu
  
  par_location <- c(par_gamma_curr,par_lambda_y_curr,kappa_curr)
  n_location <- length(par_location)
  
  obj_fun_gradient <- function(par_location){
    xi_tilde_j <- xi_tilde_mat[j,]
    xi_neg_hess_j <- xi_neg_hess_array[,,j]
    y_j <- y_mat[j,]
    
    gamma_tmp <- gamma_curr
    gamma_tmp[gamma_index] <- par_location[1:n_gamma]
    par_location <- par_location[-c(1:n_gamma)]
    
    lambda_y_tmp <- lambda_y_curr
    lambda_y_tmp[lambda_y_index] <- par_location[1:n_lambda_y]
    par_location <- par_location[-c(1:n_lambda_y)]
    
    kappa_tmp <- par_location[1:n_kappa]
    par_location <- par_location[-c(1:n_kappa)]
    
    value <- neglog_ev_j_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_tmp,lambda_y_tmp,kappa_tmp,mu_y_curr,xi_tilde_j,xi_neg_hess_j,y_j)
    return(value)
  }
  
  
  g_location_mat <- matrix(0,nrow=n,ncol=n_location)
  for(j in 1:n){
    g_location_mat[j,] <- numericGradient(obj_fun_gradient,t0=par_location)
  }
  
  obj_fun_hessian <- function(par_location){
    
    gamma_tmp <- gamma_curr
    gamma_tmp[gamma_index] <- par_location[1:n_gamma]
    par_location <- par_location[-c(1:n_gamma)]
    
    lambda_y_tmp <- lambda_y_curr
    lambda_y_tmp[lambda_y_index] <- par_location[1:n_lambda_y]
    par_location <- par_location[-c(1:n_lambda_y)]
    
    kappa_tmp <- par_location[1:n_kappa]
    par_location <- par_location[-c(1:n_kappa)]
    
    value <- neglog_ev_sem_fun(phi_curr,psi_curr,omega_y_curr,gamma_tmp,lambda_y_tmp,kappa_tmp,mu_y_curr,xi_tilde_mat,xi_neg_hess_array,y_mat)
    return(value)
  }
  
  h_location <- numericHessian(obj_fun_hessian,t0=par_location)
  h_location_inv <- solve(h_location)
  
  #Investigate if the Hessian is positive definite
  h_eigen_location <- eigen(solve(h_location))
  if(sum(h_eigen_location$values<0)==0){
    hessian_location_posdef <- TRUE
    se_location_sandwich <- sqrt(diag(h_location_inv%*%t(g_location_mat)%*%g_location_mat%*%h_location_inv))
  }else{
    hessian_posdef <- FALSE
    #cat(c("Hessian with respect to gamma is not positive definite."),"\n")
    se_location_sandwich <- NA
  }
  
  
  se_location_score <- sqrt(diag(solve(t(g_location_mat)%*%g_location_mat)))
  
  
  #Score based se
  se_location_gamma_score <- se_location_score[1:n_gamma]
  se_location_score <- se_location_score[-c(1:n_gamma)]
  
  se_location_lambda_y_score <- se_location_score[1:n_lambda_y]
  se_location_score <- se_location_score[-c(1:n_lambda_y)]
  
  se_location_kappa_score <- se_location_score[1:n_kappa]
  se_location_score <- se_location_score[-c(1:n_kappa)]
  
  
  #Location based se
  se_location_gamma_sandwich <- se_location_sandwich[1:n_gamma]
  se_location_sandwich <- se_location_sandwich[-c(1:n_gamma)]
  
  se_location_lambda_y_sandwich <- se_location_sandwich[1:n_lambda_y]
  se_location_sandwich <- se_location_sandwich[-c(1:n_lambda_y)]
  
  se_location_kappa_sandwich <- se_location_sandwich[1:n_kappa]
  se_location_sandwich <- se_location_sandwich[-c(1:n_kappa)]
  
  
  #Gather results
  res <- list()
  res$par_all <- par_curr_all
  res$mu_y <- mu_y_curr
  res$lambda_y <- lambda_y_curr
  res$kappa <- kappa_curr
  res$par_lambda_y <- par_lambda_y_curr
  res$gamma <- gamma_curr
  res$phi <- phi_curr
  res$psi <- psi_curr
  res$omega_y <- omega_y_curr
  res$par_gamma <- par_gamma_curr
  res$par_phi <- par_phi_curr
  res$par_psi <- par_psi_curr
  res$par_omega_y <- par_omega_y_curr
  res$iter <- iter
  res$y_mat <- y_mat
  res$gambling_data_list <- gambling_data_list
  res$xi_tilde_mat <- xi_tilde_mat
  res$xi_neg_hess_array  <- xi_neg_hess_array
  res$neg_hess_pos_def <- neg_hess_pos_def  #1 if the negative Hessian at the mode was positive definite at the optimum
  res$se_gamma_sandwich <- se_gamma_sandwich
  res$se_gamma_score <- se_gamma_score
  res$se_location_gamma_score <- se_location_gamma_score
  res$se_location_lambda_y_score <- se_location_lambda_y_score
  res$se_location_kappa_score <- se_location_kappa_score
  res$se_location_gamma_sandwich <- se_location_gamma_sandwich
  res$se_location_lambda_y_sandwich <- se_location_lambda_y_sandwich
  res$se_location_kappa_sandwich <- se_location_kappa_sandwich
  res$converge <- converge
  
  
  return(res)
  
  
}

# ------------------------------------------------------------------------------
# gen_data_fun()
#
# Generates simulated data (behavioral choices and continuous outcomes) from
# the hierarchical model given true parameters.
#
# Arguments:
#   info     – A list containing simulation specifications:
#                n              – Number of subjects.
#                n_block_trials – Number of repetitions of each probability/
#                                ambiguity block.
#
#   truepar  – A list of true parameters used for data generation:
#                phi, psi, omega_y  – Covariance matrices.
#                kappa              – Mean of xi (subject-level parameters).
#                gamma              – Regression matrix relating xi to eta.
#                lambda_y           – Loading matrix for y.
#                mu_y               – Mean vector for y.
#
# Returns:
#   A list with:
#     gambling_data_sim – A list of simulated choice datasets per subject.
#     xi_mat            – True xi values for all subjects.
#     eta_mat           – True eta values for all subjects.
#     y_mat             – Generated continuous outcomes.
#
# Notes:
#   - Behavioral data are generated trial-by-trial using prob_gamble_fun().
#   - Continuous outcomes y are generated via linear model:
#       y_j = lambda_y * eta_j + epsilon_y.
#   - Useful for testing estimation performance and validating algorithms.
# ------------------------------------------------------------------------------
gen_data_fun <- function(info,truepar){
  
  #Introduce parameters for generating data
  phi <- truepar$phi
  psi <- truepar$psi
  omega_y <- truepar$omega_y
  kappa <- truepar$kappa
  gamma <- truepar$gamma
  lambda_y <- truepar$lambda_y
  mu_y <- truepar$mu_y
  
  
  #Collect some important information
  n <- info$n
  n_block_trials <- info$n_block_trials
  
  k_xi <- dim(phi)[1]
  k_eta <- dim(psi)[1]
  k_y <- length(mu_y)
  xi_mat <- matrix(NA,nrow=n,ncol=k_xi)
  eta_mat <- matrix(NA,nrow=n,ncol=k_eta)
  epsilon_xi <- rmvnorm(n,mean=rep(0,k_xi),sigma=phi)
  epsilon_eta <- rmvnorm(n,mean=rep(0,k_eta),sigma=psi)
  for(j in 1:n){
    xi_mat[j,] <- epsilon_xi[j,] + kappa
    eta_mat[j,] <- gamma%*%(xi_mat[j,] - kappa) + epsilon_eta[j,]
    
  }
  
  
  data_sim <- list()
  y_mat <- matrix(0,nrow=n,ncol=k_y)
  epsilon_y <- rmvnorm(n,mean=rep(0,k_y),sigma=omega_y)
  for(j in 1:n){
    P_tmp <- c(0.125,0.25,0.375,0.5,0.625,0.75,rep(0.5,4))
    A_tmp <- c(rep(0,6),0.25,0.5,0.75,1)
    P_j <- rep(P_tmp,n_block_trials)
    A_j <- rep(A_tmp,n_block_trials)
    
    n_trials <- length(P_j)
    
    data_sim_j <- cbind(P_j,A_j,sample(c(5,8,20,50),n_trials,replace=TRUE))
    data_sim_j <- data.frame(data_sim_j)
    names(data_sim_j) <- c("p","amb","gain")
    #alpha_j <- xi_mat[j,1]
    alpha_j <- exp(xi_mat[j,1])
    beta_j <- xi_mat[j,2]
    mu_j <- exp(xi_mat[j,3])
    prob_gamble <- c()
    gamble <- c()
    for(k in 1:n_trials){
      prob_gamble[k] <- prob_gamble_fun(data_sim_j$p[k],data_sim_j$amb[k],data_sim_j$gain[k],alpha_j,beta_j,mu_j)
      gamble[k] <- as.numeric(runif(1)<prob_gamble[k])
      
    }
    #data_sim_j <- cbind(data_sim_j,prob_gamble,gamble)
    data_sim_j <- data.frame(data_sim_j,prob_gamble,gamble)  #y is the choice to gamble or not (1 or 0)
    
    names(data_sim_j) <- c("p","amb","gain","prob_gamble","y")
    
    y_mat[j,] <- lambda_y%*%eta_mat[j,] + epsilon_y[j,]
    
    data_sim[[j]] <- data_sim_j
    
  }
  
  res <- list()
  res$gambling_data_sim <- data_sim
  res$xi_mat <- xi_mat
  res$eta_mat <- eta_mat
  res$y_mat <- y_mat
  
  return(res)
}
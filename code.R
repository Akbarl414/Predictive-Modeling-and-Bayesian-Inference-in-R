#' Akbar Latif
#' neg_log_lik
#
#' @description Evaluate the negated log-likelihood for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model

neg_log_lik <- function(beta, data, model){
  
  mu <- beta[1] + beta[2]*data[["CAD_Weight"]]
  
  # distinguish between the two models to find the particular standard deviation for the betas
  if(model == "A") {
    sigma <- sqrt(exp(beta[3] + beta[4]*data[["CAD_Weight"]]))
  }else{
    sigma <- sqrt(exp(beta[3])+exp(beta[4]) * (data[["CAD_Weight"]]^2))
  }
  - sum(dnorm(data[["Actual_Weight"]],
              mean = mu,
              sd=sigma,
              log = TRUE))
  
}

#' filament_estimate
#
#' @description Estimate filament models with different variance structure
#' @param data A `data.frame` with the same variables as the `filament1` data set.
#' Must have columns `CAD_Weight` and `Actual_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @return An estimation object suitable for use with [filament1_predict()]

filament1_estimate <- function(data, model) {
  model <- match.arg(model, c("A", "B"))
  if (model == "A") {
    beta_start <- c(-0.1, 1.07, -2, 0.05)
  } else {
    beta_start <- c(-0.15, 1.07, -13.5, -6.5)
  }
  opt <- optim(beta_start,
               neg_log_lik,
               data = data,
               model = model,
               hessian = TRUE,
               method = "Nelder-Mead",
               control = list(maxit = 5000)
  )
  fit <- list(
    model = model,
    par = opt$par,
    hessian = opt$hessian
  )
  class(fit) <- c("filament1_estimate", "list")
  fit
}

#' filament1_aux_EV
#' 
#' @description Evaluate the expectation and variance for model A and B
#' @param beta A vector with the beta parameters
#' @param data A `data.frame` containing the required predictors, including `CAD_Weight`
#' @param model Either "A" for a log-linear variance model, or "B" for a proportional
#' scaling error model
#' @param Sigma_beta : If not NULL, an estimate of the covariance matrix for
#                 the uncertainty of estimated betas
#' @return A list with four elements:
#     E : E(y|beta,x)
#     V : Var(y|beta,x)
#     VE : Var(E(y|beta,x)|x) or NULL
#     EV : E(Var(y|beta,x)|x) or NULL

filament1_aux_EV <- function(beta, data, model = c("A", "B"),
                             Sigma_beta = NULL) {
  
  model <- match.arg(model)
  if (model == "A") {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      EV <- exp(ZV %*% beta + rowSums(ZV * (ZV %*% Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = exp(ZV %*% beta),
      VE = VE,
      EV = EV
    )
  } else {
    
    ZE.0 <- model.matrix( ~ 1 + CAD_Weight, data = data)
    ZV.0 <- model.matrix( ~ 1 + I(CAD_Weight^2), data = data)
    ZE = cbind(ZE.0, ZV.0 * 0) 
    ZV = cbind(ZE.0 * 0, ZV.0)
    
    VE <- EV <- NULL
    if (!is.null(Sigma_beta)) {
      # E(Var(y|beta,x)|x)
      # (pmin: Ignore large Sigma_beta values)
      EV <- ZV %*% exp(beta + pmin(0.5^2, diag(Sigma_beta)) / 2)
      # Var(E(y|beta,x)|x)
      VE <- rowSums(ZE * (ZE %*% Sigma_beta))
    }
    out <- list(
      E = ZE %*% beta,
      V = ZV %*% exp(beta),
      VE = VE,
      EV = EV
    )
  }
  out
}

#' filament1_predict
#'
#' @description in this function will create a predictive data of the data based on newdata
#' @param data the old data to be predicted upon
#' @param model which model is being used
#' @param newdata the newdata which is being used to predict against the old data
#' 
#' @return 

filament1_predict <- function(data, model = c("A", "B"), newdata){
  
  beta <- filament1_estimate(data, model)
  
  Sigma_beta <- solve(beta$hessian)
  
  EV <- filament1_aux_EV(beta$par, newdata, model, Sigma_beta)
  
  sd = sqrt(EV$EV + EV$VE)
  
  # Compute 95% prediction intervals
  alpha = 0.05
  lwr <- EV$E - qnorm(1 - (alpha/2))*sd
  upr <- EV$E + qnorm(1 - (alpha/2))*sd
  
  # Create data frame to store results
  results <- data.frame(
    mean = EV$E,
    sd = sd,
    lwr = c(lwr),
    upr = c(upr)
  )
  
  return(results)
}

#' square_error_score
#'
#' @description
#' @param prediction a data frame with Actual_Weght and mean
#' @return the prediction data frame with a new column
square_error_score <- function(prediction){
  score <- prediction %>% 
    mutate(
      se = ((Actual_Weight - mean))^2)
}


#' ds_score
#'
#' @description
#' @param prediction a data frame with Actual_Weght and mean and sd
#' @return the prediction data frame with a new column
ds_score <- function (prediction){
  score <- prediction %>% 
    mutate(
  ds = (Actual_Weight - mean)^2/sd^2 + 2 * log(sd))

}

#' leave1out
#'
#' @description
#' @param data the data 
#' @param model which model is being used
#' @return the data frame with the 

leave1out <- function(data, model = c("A","B")){
  data_new <- data %>% mutate(mean = NA_real_, sd = NA_real_, se = NA_real_, ds = NA_real_)
 
   for (i in seq_len(nrow(data_new))) {
    pred <- filament1_predict(data_new[-i, , drop = FALSE], model, data_new[i,])
    data_new[i, "mean"] = pred$mean
    data_new[i, "sd"] <- pred$sd
  }
  
  data_final <- data_new %>% mutate(
    se = ((Actual_Weight - mean))^2,
    ds = (Actual_Weight - mean)^2/sd^2 + 2 * log(sd))

  return(data_final)
}

#' monte_p_value
#'
#' @description Compute p values using monte carlo methods
#' @param dataA the first bit of data to use
#' @param dataB the second bit of data to use
#' @param N the number of samples to be used
#' @return a data frame with the computed se and ds p values 

Monte_p <- function(dataA, dataB, N){

  score_diff <- data.frame(se = dataA$se - dataB$se, 
                           ds = dataA$ds - dataB$ds)
  statistic0 <- score_diff %>% summarise(se = mean(se), ds = mean(ds)) 
  statistic <- data.frame(se = numeric(N),
                          ds = numeric(N)) 
  for (loop in seq_len(N)) {
    random_sign <- sample(c(-1, 1), size = nrow(score_diff), replace = TRUE) 
    statistic[loop, ] <- score_diff %>% summarise(se = mean(random_sign * se), 
                                                  ds = mean(random_sign * ds))
                          }
  p_values <-
    statistic %>%summarise(se = mean(se > statistic0$se), 
                           ds = mean(ds > statistic0$ds))
  # Estimates:
  return(p_values)
}



#' arch_loglike
#'
#' @description compute the combined likelihood
#'
#' @param data is a data frame with columns N, phi
#' @param y is the vector of observations
#' @return

arch_loglike <- function(data, y){
  N <- data$N
  phi <- data$phi
 
 
  log_likelihood <- -lgamma(y[1] + 1) - lgamma(y[2] + 1) - lgamma(N - y[1] + 1) - lgamma(N - y[2] + 1) + 
    2 * lgamma(N + 1) + (y[1] + y[2]) * log(phi) + (2 * N - y[1] - y[2]) * log(1 - phi)
  
  
  return(log_likelihood)
}



#' estimate
#'
#' @description
#' @param y the observations 
#' @param xi the prior for n 
#' @param a prior for phi
#' @param b prior for phi
#' @param K the number of samples
#' @return

estimate <- function(y, xi, a, b, K){

  # Sample from the prior distributions
  N_samples <- rgeom(K, xi)
  phi_samples <- rbeta(K, a, b)

  # Calculate the unnormalized posterior probabilities
  data <- data.frame(N = N_samples, phi_samples)
  log_posterior <- arch_loglike(data, y)
  posterior_unnormalized <- exp(log_posterior)

  # Calculate the normalization constant
  normalization_constant <- sum(posterior_unnormalized)

  # Calculate posterior probability p(N|y)
  p_N_given_y <- posterior_unnormalized / normalization_constant

  # Calculate expected value of N given y
  E_N_given_y <- sum(N_samples * p_N_given_y)

  # Calculate expected value of phi given y
  E_phi_given_y <- sum(phi_samples * p_N_given_y)

  # Calculate Monte Carlo estimate of p(y)
  p_y <- normalization_constant / K

  # Return results
  results <- list(
    p_y = p_y,
    E_N_given_y = E_N_given_y,
    E_phi_given_y = E_phi_given_y
  )

  return(results)
}


############################################ Functions created for work in relation to master thesis on NMF ############################################


################################## Regular Poisson NMF updates #######################################

poisson_updates <- function(V, K, Iter = 1000, alpha = 0, epsilon = 1e-4) {
  # If V is 0, make it epsilon
  V[V == 0] <- epsilon
  
  # Define N and M
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Initialize matrices
  # set.seed(123)
  W <- matrix(runif(N*K,0.3,0.5), nrow = N)
  H <- matrix(runif(K*M,0.3,0.5), nrow = K)
  
  # Counter to see convergence
  counter <- 0
  
  # Make empty matrix for Kullback-Leibler divergence
  D_kl <- numeric(Iter)
  
  # Find initial Kullback-Leibler value
  WH <- W %*% H
  WH[WH == 0] <- epsilon  
  D_kl[1] <- sum(V * log(V / WH) - V + WH)
  
  # Find updates iteratively until convergence
  for (i in 2:Iter) {
    # Updating WH
    WH <- W %*% H
    # If zero, make it epsilon
    WH[WH == 0] <- epsilon
    
    # Updating signatures
    H <- H * (t(W) %*% (V / WH)) / colSums(W)
    
    # Normalizing H 
    H <- H / rowSums(H)
    
    # Updating WH with new H
    WH <- W %*% H
    WH[WH == 0] <- epsilon  
    
    # Updating exposures
    W <- W * ((V / WH) %*% t(H)) / rowSums(H)
    
    # Calculating Kullback Leibler divergence
    D_kl[i] <- sum(V * log(V / WH) - V + WH)
    
    # Update counter
    counter <- counter + 1
    
    # Check for convergence
    if (abs(D_kl[i-1] - D_kl[i]) < epsilon) {
      break
    }
  }
  # If V was changed, change it back to zero for use in dpois function
  V[V < 1] <- 0
  
  # Find log-likelihood value of final parameters
  loglik <- sum(dpois(as.vector(V), lambda = WH, log = TRUE))
  
  # Return D_kl, H, W, log-likelihood value and counter as a list
  return(list(D_kl = D_kl, H = H, W = W, loglik = loglik, counter = counter))
}


################################## Regular negative binomial updates #######################################
# Note: These updates where implemented, but not used in the thesis.
NB_updates <- function(V, K, Iter = 1000, epsilon = 1e-10, alpha = 0.1) {
  # If V is 0, make it epsilon
  V[V == 0] <- epsilon
  
  # Define N and M
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Initialize matrices
  W <- matrix(runif(N*K,0.3,0.5), nrow = N)
  H <- matrix(runif(K*M,0.3,0.5), nrow = K)
  
  # Make empty vector for log-likelihood values
  loglik <- numeric(Iter)
  
  # Make empty matrix for Kullback-Leibler divergence
  D_kl <- numeric(Iter)
  
  # Find updates iteratively until convergence 
  for (i in 1:Iter) {
    # Updating WH
    WH <- W %*% H
    # If zero, make it epsilon
    WH[WH == 0] <- epsilon
    
    # Updating signatures
    H <- H * (t(W) %*% (V / WH)) / (t(W) %*% ((V + alpha) / (WH + alpha)))
    
    # Normalizing H 
    H <- H / rowSums(H)
    
    # Updating WH with new H
    WH <- W %*% H
    WH[WH == 0] <- epsilon  
    
    # Updating exposures
    W <- W * ((V / WH) %*% t(H)) / (((V + alpha) / (WH + alpha)) %*% t(H))
    
    # Compute log-likelihood
    loglik[i] <- sum(dnbinom(as.vector(V), size = alpha, mu = WH, log = TRUE))
    
    # Kullback Leibler divergence
    D_kl[i] <- sum(V * log(V / WH) - V + WH)
  }
  # Return D_kl, H, W and log-likelihood value as a list
  return(list(D_kl = D_kl, H = H, W = W, loglik = tail(loglik, n=1)))
}


################################## BIC for choosing number of signatures for Poisson updates #######################################
BIC_val <- function(V,K, Iter, update_function, alpha = 0) {
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Create empty vector for BIC values
  BIC_values <- matrix(0, nrow = length(K))
  parameter_values <- numeric(length(K))
  loglik_values <- numeric(length(K))
  # Number of observations: nxM
  
  # Find BIC value for all given choices of K
  for (i in 1:length(K)) {
    k <- K[i]
    # Print to see how far function has run
    cat("Number of signatures:", k, "\n")
    results <- update_function(V, k, Iter, alpha = alpha) 
    # Find the log-likelihood value
    loglik <- results$loglik
    loglik_values[i] <- loglik
    # Find number of parameters
    parameters <- (N*k + k*M) * log(N*M)
    parameter_values[i] <- parameters
    # Calculate BIC value
    BIC_values[i] <- parameters - 2*loglik
  }
  return(list(BIC_values = BIC_values, loglik_values = loglik_values, parameter_values = parameter_values))
}






################################## BIC for choosing number of signatures for Poisson updates - ONLY NON-ZERO ENTRIES #######################################
BIC_val_new <- function(V,K, Iter, update_function, alpha = 0) {
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Number of non-zero entries in V
  non_zero_entries <- sum(V != 0)
  print(non_zero_entries)
  
  # Create empty vector for BIC values
  BIC_values <- matrix(0, nrow = length(K))
  parameter_values <- numeric(length(K))
  loglik_values <- numeric(length(K))
  # Number of observations: nxM
  
  # Find BIC value for all given choices of K
  for (i in 1:length(K)) {
    k <- K[i]
    # Print to see how far function has run
    cat("Number of signatures:", k, "\n")
    results <- update_function(V, k, Iter, alpha = alpha) 
    # Find the log-likelihood value
    loglik <- results$loglik
    loglik_values[i] <- loglik
    # Find number of parameters
    parameters <- (N*k + k*M) * log(non_zero_entries)
    parameter_values[i] <- parameters
    # Calculate BIC value
    BIC_values[i] <- parameters - 2*loglik
  }
  return(list(BIC_values = BIC_values, loglik_values = loglik_values, parameter_values = parameter_values))
}





################################## Cross validation for choosing number of signatures #######################################
cross_val <- function(V, K, Iter, J, update_function, alpha = 0) {
  # Create empty matrix for the cost
  Cost <- matrix(0, nrow = J, ncol = length(K))
  
  # Calculate the cost for each choice of K
  for (k in K) {
    print("Number of signatures:")
    print(k)
    
    # Find the updates for the full dataset
    results <- update_function(V, k, Iter) 
    H_full <- results$H
    W_full <- results$W
    
    # Repeat calculations J times for stability
    for (j in 1:J) {
      print("Cross validation fold:")
      print(j)
      
      # Determine the number of trios to select (10% of total trios)
      num_rows <- round(0.1 * nrow(V))
      
      # Randomly select indices of trios for cross-validation test set
      indices_for_cv <- sample(nrow(V), num_rows)
      
      # Create training and cross-validation datasets
      V_train <- V[-indices_for_cv, ]
      V_test <- V[indices_for_cv, ]
      V_test[V_test == 0] <- 1e-10
      W_test <- W_full[indices_for_cv, ]
      
      # Find updates for signatures for training set
      results_train <- update_function(V_train, k, Iter, alpha = alpha)
      H_train <- results_train$H
      
      # Calculation of cost by Kullback-Leibler divergence
      WH <- W_test %*% H_train
      WH[WH == 0] <- 1e-10  
      Cost[j,k-1] <- sum(V_test * log(V_test / WH) - V_test + WH)
    }
  }
  # Find the median cost over the J iterations
  Cost_median <- apply(Cost, 2, median)
  return(list(Cost_median = Cost_median, Cost = Cost))
}


################################## Comparing signatures functions #######################################
# Function to find indices and values of max and second largest values in a vector
find_max_second <- function(vec) {
  max_idx <- which.max(vec)
  max_val <- max(vec)
  sorted_idx <- order(vec, decreasing = TRUE)
  second_idx <- sorted_idx[2]
  second_val <- vec[second_idx]
  return(list(max_idx = max_idx, max_val = max_val, 
              second_idx = second_idx, second_val = second_val))
}

find_min_second <- function(vec) {
  min_idx <- which.min(vec)
  min_val <- min(vec)
  sorted_idx <- order(vec, decreasing = FALSE)
  second_idx <- sorted_idx[2]
  second_val <- vec[second_idx]
  return(list(min_idx = min_idx, min_val = min_val, 
              second_idx = second_idx, second_val = second_val))
}


############################################################### Parametric Poisson NMF updates ###############################################################


######################################## Function to find the betas by GLM update ########################################
### To find the betas
glm.update = function(y,X, epsilon) {
  # mu corresponds to exp(X %*% beta) in the thesis
  mu <- y + runif(length(y),0,0.1)
  # v corresponds to v in the thesis
  v <- log(mu) + (y-mu)/mu
  v <- v[v == 0] <- epsilon
  # By syntas of glm.fit, you need to multiply with the squareroot of the weight matrix which is A = mu
  w = sqrt(mu)
  # Use glm.fit to find estimates for the betas
  fit <- glm.fit(X, v*w, family = poisson(link = "log"))
  # Extract the betas
  betas <- fit$coefficients
  # This is to determine log(w)
  param <- as.vector(X %*% betas)
  # Determine w 
  w_param <- exp(param)
  return(w_param)
}

######################################## Function to find the updates for parametric model ########################################
### Updates
poisson_updates_param <- function(V, K, design_matrix, Iter = 1000, alpha = 0, epsilon = 1e-10) {
  # If V is 0, make it epsilon
  V[V == 0] <- epsilon
  
  # Define N and M
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Initialize matrices
  W <- matrix(runif(N*K,0.3,0.5), nrow = N)
  H <- matrix(runif(K*M,0.3,0.5), nrow = K)
  
  # Counter to see convergence
  counter <- 0
  
  # Make empty matrix for Kullback-Leibler divergence
  D_kl <- numeric(Iter)
  
  # Find initial Kullback-Leibler value
  WH <- W %*% H
  WH[WH == 0] <- epsilon  
  D_kl[1] <- sum(V * log(V / WH) - V + WH)
  
  for (i in 2:Iter) {
    if (i %% 1000 == 0) {print(i)}
    # Updating WH
    WH <- W %*% H
    # If zero, make it epsilon
    WH[WH == 0] <- epsilon
    
    # Updating signatures
    H <- H * (t(W) %*% (V / WH)) / colSums(W)
    
    # Normalizing H 
    H <- H / rowSums(H)
    
    # Updating WH with new H
    WH <- W %*% H
    WH[WH == 0] <- epsilon  
    
    # Updating exposures as normal
    W <- W * ((V / WH) %*% t(H)) / rowSums(H)
    
    # Finding betas by applying glm.update function. Use this to find new updates
    W <- t(sapply(1:K, function(k) glm.update(W[,k], design_matrix)))
    W <- t(W)
    
    WH <- W %*% H
    WH[WH == 0] <- epsilon
    
    # Kullback Leibler divergence
    D_kl[i] <- sum(V * log(V / WH) - V + WH)
    
    # Update counter
    counter <- counter + 1
    
    # Check for convergence
    if (abs(D_kl[i-1] - D_kl[i]) < 1e-4) {
      break
    }
  }
  loglik <- sum(dpois(as.vector(V), lambda = WH, log = TRUE))
  # Return D_kl, H, and W as a list
  return(list(D_kl = D_kl, H = H, W = W, loglik = loglik, counter = counter))
}


################################## BIC for Parametric Poisson updates #######################################
BIC_val_param <- function(V,K, design_matrix, Iter, update_function, alpha = 0, save_path = "BIC_results.RDS") {
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Create empty vector for BIC values
  BIC_values <- matrix(0, nrow = length(K))
  loglik_values <- numeric(length(K))
  parameter_values <- numeric(length(K))
  # Number of observations: nxM
  
  # Find BIC value for all given choices of K
  for (i in 1:length(K)) {
    k <- K[i]
    # Print to see how far function has run
    cat("Number of signatures:", k, "\n")
    
    results <- update_function(V, k, design_matrix, Iter, alpha = alpha) 
    # Find the log-likelihood value
    loglik <- results$loglik
    loglik_values[i] <- loglik
    # Find the number of parameters
    parameters <- (2*k + k*(M -1)) * log(N*M)
    parameter_values[i] <- parameters
    # Calculate the BIC value
    BIC_values[i] <- parameters - 2*loglik
    
    # Save results as RDS file, overwriting the same file each time
    result_data <- list(BIC_values = BIC_values, parameter_values = parameter_values, loglik_values = loglik_values)
    saveRDS(result_data, file = save_path)
    
    cat("Results saved for k =", k, "\n")
  }
  cat("All results saved to:", save_path, "\n")
  
  return(list(BIC_values = BIC_values, parameter_values = parameter_values, loglik_values = loglik_values))
}





################################## BIC for Parametric Poisson updates ONLY NON-ZERO ENTRIES #######################################
BIC_val_param_new <- function(V,K, design_matrix, Iter, update_function, alpha = 0, save_path = "BIC_results_test.RDS") {
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  # Create empty vector for BIC values
  non_zero_entries <- sum(V != 0)
  BIC_values <- matrix(0, nrow = length(K))
  loglik_values <- numeric(length(K))
  parameter_values <- numeric(length(K))
  # Number of observations: nxM
  
  # Find BIC value for all given choices of K
  for (i in 1:length(K)) {
    k <- K[i]
    # Print to see how far function has run
    cat("Number of signatures:", k, "\n")
    
    results <- update_function(V, k, design_matrix, Iter, alpha = alpha) 
    # Find the log-likelihood value
    loglik <- results$loglik
    loglik_values[i] <- loglik
    # Find the number of parameters
    parameters <- (2*k + k*(M -1)) * log(non_zero_entries)
    parameter_values[i] <- parameters
    # Calculate the BIC value
    BIC_values[i] <- parameters - 2*loglik
    
    # Save results as RDS file, overwriting the same file each time
    result_data <- list(BIC_values = BIC_values, parameter_values = parameter_values, loglik_values = loglik_values)
    saveRDS(result_data, file = save_path)
    
    cat("Results saved for k =", k, "\n")
  }
  cat("All results saved to:", save_path, "\n")
  
  return(list(BIC_values = BIC_values, parameter_values = parameter_values, loglik_values = loglik_values))
}






################################## Cross validation for Parametric Poisson updates #######################################
cross_val_param <- function(V, K, design_matrix, Iter, J, update_function, save_interval = 1, save_path = "cross_val_progress.RDS") {
  # Create empty matrix for the cost
  Cost <- matrix(0, nrow = J, ncol = length(K))
  
  # Calculate the cost for each choice of K
  for (k_index in 1:length(K)) {
    k <- K[k_index]
    # Print to see how far function has run
    cat("Cross validation. Signature value:", k, "\n")
    # Find the updates for the full dataset
    results <- update_function(V, k, design_matrix, Iter)  # Pass k to the update_function
    H_full <- results$H
    W_full <- results$W
    
    # Repeat calculations J times for stability
    for (j in 1:J) {
      cat("CV fold:", j, "\n")
      
      # Determine the number of trios to select (10% of total trios)
      num_rows <- round(0.1 * nrow(V))
      
      # Randomly select indices of trios for cross-validation test set
      indices_for_cv <- sample(nrow(V), num_rows)
      
      
      # Create training and cross-validation datasets
      V_train <- V[-indices_for_cv, ]
      V_test <- V[indices_for_cv, ]
      V_test[V_test == 0] <- 1e-10
      W_test <- W_full[indices_for_cv, ]
      
      # Change design matrix to match chosen indices
      design_matrix_train <- design_matrix[-indices_for_cv,]
      
      # Find updates for signatures for the training set
      results_train <- update_function(V_train, k, design_matrix_train, Iter)  # Pass k to the update_function
      H_train <- results_train$H
      
      # Calculation of cost by Kullback-Leibler divergence
      WH <- W_test %*% H_train
      WH[WH == 0] <- 1e-10  
      Cost[j, k_index] <- sum(V_test * log(V_test / WH) - V_test + WH)
      
      # Save the progress periodically
      if (j %% save_interval == 0) {
        saveRDS(list(Cost = Cost, k_index = k_index, j = j), save_path)
      }
    }
  }
  # Find the median cost over the J iterations
  Cost_median <- apply(Cost, 2, median)
  return(list(Cost_median = Cost_median, Cost = Cost))
}










############################################### For making signature plots ###############################################
# Function to make the signature plots
create_signature_plot <- function(data, x, y, mutationtype) {
  p <- ggplot(data, aes(x = x, y = y, fill = mutationtype)) + 
    geom_bar(stat = "identity", aes(y = y), color = "black") + 
    scale_fill_manual(values = c("skyblue", "salmon", "green", "purple", "orange", "forestgreen"), labels = unique(mutationtype)) +
    labs(title = sprintf("%s", deparse(substitute(y)))) +
    theme_minimal() + 
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))
}






############################################### For checking slopes of exposures ###############################################
# Calculate t-test for the slope of each exposure signature
calculate_slope_ttest <- function(exposure, age) {
  model <- lm(exposure ~ age)
  summary_model <- summary(model)
  # Find the slope
  slope <- summary_model$coefficients[2, 1] 
  # Find the standard error
  std_error <- summary_model$coefficients[2, 2]  
  # Calculate t-value
  t_value <- slope / std_error 
  # Find degrees of freedom
  df <- summary_model$df[2]  
  # Calculate p-value
  p_value <- 2 * pt(-abs(t_value), df) 
  p_value <- format(p_value, digits = 5, nsmall = 5)
  return(list(slope = slope, std_error = std_error, t_value = t_value, p_value = p_value))
}
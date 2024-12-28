# Load required libraries ------------------------------------------------------
my_packages <- c("quantreg","rqpd", "plm", "psych", "tseries", "pco", "dplyr", 
                 "urca", "lmtest", "car", "MASS", "AER", "rrcov", "Matrix", 
                 "ggplot2", "gridExtra", "tidyr", "stringr", "tibble", "purrr", 
                 "MLmetrics", "patchwork")
lapply(my_packages, library, character.only = TRUE)

# Import dataset ---------------------------------------------------------------
rawdata <- read.csv("newdataset.csv")
str(rawdata)

# Log transformation of wages
data <- rawdata %>%
  mutate(
    log_ump = log(ump)
  )

# Create lagged of variables
lag_data <- data %>%
  group_by(Provinsi) %>%
  mutate(
    log_ump_lag = lag(log_ump),
    exp_lag = lag(expenses, 1),
    tpt_lag = lag(tpt, 1),
    gini_lag = lag(gini, 1)
  ) %>%
  na.omit()

# Summary statistics
describe(pdata.frame(data, index = c("Provinsi", "Time")))
describe(pdata.frame(lag_data, index = c("Provinsi", "Time")))

jarque.bera.test(pdata.frame(data, index = c("Provinsi", "Time"))$ump)
jarque.bera.test(pdata.frame(data, index = c("Provinsi", "Time"))$expenses)
jarque.bera.test(pdata.frame(data, index = c("Provinsi", "Time"))$tpt)
jarque.bera.test(pdata.frame(data, index = c("Provinsi", "Time"))$gini)

jarque.bera.test(pdata.frame(lag_data, index = c("Provinsi", "Time"))$log_ump)
jarque.bera.test(pdata.frame(lag_data, index = c("Provinsi", "Time"))$log_ump_lag)
jarque.bera.test(pdata.frame(lag_data, index = c("Provinsi", "Time"))$exp_lag)
jarque.bera.test(pdata.frame(lag_data, index = c("Provinsi", "Time"))$tpt_lag)
jarque.bera.test(pdata.frame(lag_data, index = c("Provinsi", "Time"))$gini_lag)

# Correlation
cor(pdata.frame(lag_data, index = c("Provinsi", "Time")) 
    %>% dplyr::select(log_ump, exp_lag, tpt_lag, gini_lag))

# Check for multicollinearity --------------------------------------------------
vif(lm(log_ump ~ log_ump_lag + exp_lag + tpt_lag + gini_lag, data = lag_data))

# Detect outliers --------------------------------------------------------------
robust_mahalanobis <- function(data) {
  x <- data %>% 
    dplyr::ungroup() %>%
    dplyr::select(exp_lag, tpt_lag, gini_lag) %>% as.matrix()
  
  # Robust covariance matrix using Fast-MCD
  mcd <- covMcd(x)
  mean_vec <- mcd$center
  cov_matrix <- mcd$cov
  
  # Compute robust Mahalanobis distance
  data$robust_md <- apply(x, 1, function(row) {
    diff <- row - mean_vec
    sqrt(t(diff) %*% solve(cov_matrix) %*% diff)
  })
  return(data)
}

# Apply the function
rmd <- lag_data %>%
  group_by(Provinsi) %>%
  group_modify(~ robust_mahalanobis(.x)) %>%
  ungroup()

# Identify outliers
p <- 3 # Number of variables
threshold <- qchisq(0.95, df = p)
outliers <- rmd %>% filter(robust_md > threshold)
print(threshold)
print(outliers)

# Check for stationarity using Maddala-Wu test --------------------------------
log_y <- data.frame(split(lag_data$log_ump, lag_data$Provinsi))
log_Ly <- data.frame(split(lag_data$log_ump_lag, lag_data$Provinsi))
x1 <- data.frame(split(lag_data$exp_lag, lag_data$Provinsi))
x2 <- data.frame(split(lag_data$tpt_lag, lag_data$Provinsi))
x3 <- data.frame(split(lag_data$gini_lag, lag_data$Provinsi))

purtest(log_y, pmax = 1, test = "madwu")
purtest(log_Ly, pmax = 1, test = "madwu")
purtest(x1, pmax = 1, test = "madwu")
purtest(x2, pmax = 1, test = "madwu")
purtest(x3, pmax = 1, test = "madwu")

# Differencing -----------------------------------------------------------------
log_y_diff <- apply(log_y, 2, diff)
log_Ly_diff <- apply(log_Ly, 2, diff)
x1_diff <- apply(x1, 2, diff)
x2_diff <- apply(x2, 2, diff)
x3_diff <- apply(x3, 2, diff)

# Check for stationarity again
purtest(log_y_diff, pmax = 1, test = "madwu")
purtest(log_Ly_diff, pmax = 1, test = "madwu")
purtest(x1_diff, pmax = 1, test = "madwu")
purtest(x2_diff, pmax = 1, test = "madwu")
purtest(x3_diff, pmax = 1, test = "madwu")

# Convert to long panel format
long_panel <- function(data, var_name) {
  df_long <- data.frame(
    Provinsi = rep(colnames(data), each = nrow(data)),
    Time = rep(2009:2024, times = ncol(data)),
    Value = as.vector(data)
  )
  colnames(df_long)[3] <- var_name  # Rename column
  return(df_long)
}

# Create long-format data frames
log_y_long <- long_panel(log_y_diff, "log_ump")
log_Ly_long <- long_panel(log_Ly_diff, "log_ump_lag")
x1_long <- long_panel(x1_diff, "exp_lag")
x2_long <- long_panel(x2_diff, "tpt_lag")
x3_long <- long_panel(x3_diff, "gini_lag")

# Merge all data frame
diff_data <- Reduce(function(x, y) merge(x, y, by = c("Provinsi", "Time")),
                    list(log_y_long, log_Ly_long, x1_long, x2_long, x3_long))
head(diff_data)
dim(diff_data)

# Create lagged variables ------------------------------------------------------
diff_data <- diff_data %>%
  group_by(Provinsi) %>%
  mutate(
    log_ump_lag2 = lag(log_ump_lag),
    exp_lag2 = lag(exp_lag),
    tpt_lag2 = lag(tpt_lag),
    gini_lag2 = lag(gini_lag),
    log_ump_lag3 = lag(log_ump_lag2),
    exp_lag3 = lag(exp_lag2),
    tpt_lag3 = lag(tpt_lag2),
    gini_lag3 = lag(gini_lag2)
  ) %>%
  ungroup() %>%
  na.omit()
str(diff_data)

# Transform to panel data structure --------------------------------------------
df_panel <- pdata.frame(diff_data, index = c("Provinsi", "Time"))
View(df_panel)

# Instrument test
first_stage <- lm(log_ump_lag ~ log_ump_lag2 + exp_lag2 + tpt_lag2 + gini_lag2 +
                      log_ump_lag3 + exp_lag3 + tpt_lag3 + gini_lag3, data = df_panel)
summary(first_stage)

iv_model <- ivreg(log_ump ~ log_ump_lag + exp_lag + tpt_lag + gini_lag |
                   exp_lag2 + tpt_lag2 + log_ump_lag3 + exp_lag3 + tpt_lag3 + gini_lag3 
                  + exp_lag + tpt_lag + gini_lag, 
                  data = df_panel)
summary(iv_model, diagnostics = TRUE)

# QRPIV for dynamic panel with fixed-effect ------------------------------------
## Step 1: Specify grid search for coefficient alpha
alpha_grid <- seq(-0.99, 0.99, by = 0.01)  # Define the range of alpha values
norm_gamma <- c()  # Placeholder for norm of gamma

## For K > 1
# Define the quantile levels
taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)
tauw <- rep((1/5),5)

# Placeholder for results
results <- list()

# Loop over quantiles
for (tau in taus) {
  # Placeholder for norm_gamma
  norm_gamma <- numeric(length(alpha_grid))
  
  # Perform grid search for alpha
  for (j in seq_along(alpha_grid)) {
    alpha <- alpha_grid[j]
    
    # Transform the dependent variable
    df_panel$y_transformed <- df_panel$log_ump - alpha * df_panel$log_ump_lag
    
    # Fit the model for the current quantile and alpha
    model <- rqpd(
      formula = y_transformed ~ exp_lag + tpt_lag + gini_lag + 
        exp_lag2 + tpt_lag2 + log_ump_lag3 + exp_lag3 + tpt_lag3 + gini_lag3 | Provinsi,
      panel = panel(method = "pfe", taus = taus, tauw = tauw), 
      data = df_panel,
      control = quantreg::sfn.control(tmpmax = 1000000)
    )
    
    # Define the instrument names for the current quantile
    instrument_names <- paste0(c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                                 "exp_lag3", "tpt_lag3", "gini_lag3"), "[", tau, "]")
    
    # Extract coefficients for instruments
    gamma <- model$coef[instrument_names]
    
    # Compute the norm of gamma
    norm_gamma[j] <- sqrt(sum(gamma^2))
  }
  
  # Optimal alpha for current quantile
  optimal_alpha <- alpha_grid[which.min(norm_gamma)]
  
  # Final model with the optimal alpha
  df_panel$y_transformed <- df_panel$log_ump - optimal_alpha * df_panel$log_ump_lag
  
  final_model <- rqpd(
    formula = y_transformed ~ exp_lag + tpt_lag + gini_lag + 
      exp_lag2 + tpt_lag2 + log_ump_lag3 + exp_lag3 + tpt_lag3 + gini_lag3  | Provinsi,
    panel = panel(method = "pfe", taus = taus, tauw = tauw), 
    data = df_panel,
    control = quantreg::sfn.control(tmpmax = 1000000)
  )
  
  # Summary with covariance matrix
  model_summary <- summary.rqpd(final_model, cov = TRUE)
  
  # Compute Pseudo-R2
  residuals <- resid(final_model)
  pseudo_r2 <- 1 - sum(residuals^2) / sum((df_panel$y_transformed - median(df_panel$y_transformed))^2)
  
  # Extract fixed effects
  eta_names <- grep("Provinsi", names(final_model$coef), value = TRUE)
  eta <- final_model$coef[eta_names]
  
  # Extract coefficients for the final quantile model
  quantile_coef_names <- grep(paste0("\\[", tau, "\\]"), names(final_model$coef), value = TRUE)
  quantile_coefs <- final_model$coef[quantile_coef_names]
  
  # Store results
  results[[paste0("tau_", tau)]] <- list(
    alpha = optimal_alpha,
    quantile_coefs = quantile_coefs,
    fixed_effects = eta,
    residuals = residuals,
    PseudoR2 = pseudo_r2
  )
}
## View the results
results

## Confidence interval of estimation
confidence_intervals_list <- list()

# Loop through each tau level
for (tau in taus) {
  # Extract coefficients for the current tau level
  tau_label <- paste0("\\[", tau, "\\]")
  coef_data <- summary.rqpd(final_model, covariance = TRUE)$coefficients[grep(tau_label, rownames(summary.rqpd(final_model, covariance = TRUE)$coefficients)), ]
  
  # Compute 95% confidence intervals
  critical_value <- qnorm(0.975)  # 95% confidence interval
  confidence_intervals <- data.frame(
    Coefficient = rownames(coef_data),
    Estimate = coef_data[, "Value"],
    StdError = coef_data[, "Std. Error"],
    Lower = coef_data[, "Value"] - critical_value * coef_data[, "Std. Error"],
    Upper = coef_data[, "Value"] + critical_value * coef_data[, "Std. Error"],
    Quantile = tau  # Add the tau level
  )
  
  # Store the results
  confidence_intervals_list[[paste0("tau_", tau)]] <- confidence_intervals
}

# Combine results for all tau levels into a single data frame
confidence_intervals_all <- do.call(rbind, confidence_intervals_list)
confidence_intervals_all <- confidence_intervals_all %>%
  mutate(
    Coefficient = str_remove(Coefficient, "\\[\\d\\.\\d+|\\[\\d+\\]"),
    Coefficient = str_remove(Coefficient, "tau_\\d\\.\\d+\\."),
    Coefficient = str_replace(Coefficient, "\\(Intercept\\)", "Intercept"),
    Coefficient = str_trim(Coefficient),
    Coefficient = str_remove(Coefficient, "]")
  )
print(confidence_intervals_all)

# Plot the estimates with confidence intervals
beta_vars <- c("Intercept", "exp_lag", "tpt_lag", "gini_lag")
beta_data <- confidence_intervals_all %>% filter(Coefficient %in% beta_vars)

# Beta plots
beta_plots <- list()
for (var in beta_vars) {
  current_beta_data <- beta_data %>% filter(Coefficient == var)
  
  # Create sequence of quantiles including 0 and 1 for axis breaks
  all_quantiles <- sort(unique(c(0, current_beta_data$Quantile, 1)))
  
  # Determine y-axis limits with padding (same as before)
  y_min <- min(current_beta_data$Lower)
  y_max <- max(current_beta_data$Upper)
  y_range <- y_max - y_min
  y_padding <- y_range * 0.1
  y_limits <- c(y_min - y_padding, y_max + y_padding)
  
  beta_plots[[var]] <- ggplot(current_beta_data, aes(x = Quantile, y = Estimate)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey", alpha = 0.5) +
    geom_line(color = "red") +
    geom_point(color = "red", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_x_continuous(breaks = all_quantiles, labels = all_quantiles, limits=c(0.1,0.9))+
    scale_y_continuous(limits = y_limits) +
    labs(
      title = paste("Coefficient Estimates with 95% CI for", var), # More descriptive title
      x = "Quantile Level",
      y = "Estimate"
    ) +
    theme_bw()
  
  ggsave(filename = file.path("beta_plots", paste0("beta_plot_", var, ".jpg")), 
         plot = beta_plots[[var]], width = 8, height = 4)
}

1# Inference on QR coefficients -------------------------------------------------
p <- 4 # Number of covariates ; correspond to beta including intercept
q <- 6  # Number of instruments ; correspond to gamma

# Compute asymptotic variance-covariance matrix
compute_omega <- function(results, data, tau_levels, tau_weights) {
  NT <- nrow(data)
  k <- length(tau_levels)
  
  # Compute Xi_it for each tau
  compute_xi <- function(data, tau_levels, results) {
    # Extract necessary components from the result
    fixed_effects <- list()
    alpha <- list()
    quantile_coefs <- list()
    beta <- list()
    gamma <- list()
    for (k in seq_along(tau_levels)) {
      tau_k <- tau_levels[k]
      fixed_effects[[k]] <- results[[k]]$fixed_effects
      alpha[[k]] <- results[[k]]$alpha
      coefs <- results[[k]]$quantile_coefs
      beta[[k]] <- coefs[1:4]
      gamma[[k]] <- coefs[5:10]
    }
    
    # Prepare covariates and instruments
    covariate_names <- c("exp_lag", "tpt_lag", "gini_lag")
    instrument_names <- c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                          "exp_lag3", "tpt_lag3", "gini_lag3")
    covariates <- as.matrix(data[, covariate_names, drop = FALSE])
    covariates_with_beta0 <- cbind(1, covariates)
    instruments <- as.matrix(data[, instrument_names, drop = FALSE])
    
    # Initialize a matrix to store xi values for each tau
    xi <- matrix(NA, nrow = NT, ncol = k)
    
    # Loop over taus and compute xi
    for (k in seq_along(tau_levels)) {
      # Extract tau-specific coefficients
      fixed_effects <- list()
      alpha <- list()
      quantile_coefs <- list()
      beta <- list()
      gamma <- list()
      for (k in seq_along(tau_levels)) {
        tau_k <- tau_levels[k]
        fixed_effects[[k]] <- results[[k]]$fixed_effects
        alpha[[k]] <- results[[k]]$alpha
        coefs <- results[[k]]$quantile_coefs
        beta[[k]] <- coefs[1:4]
        gamma[[k]] <- coefs[5:10]
      }
      
      # Compute xi
      for (k in seq_along(tau_levels)) {
        xi[, k] <- results[[k]]$fixed_effects + 
          results[[k]]$alpha  * data$log_ump_lag +
          covariates_with_beta0 %*% beta[[k]] +
          instruments %*% gamma[[k]]
      }
    }
    colnames(xi) <- paste0("tau_", tau_levels)
    return(xi)
  }
  xi <- compute_xi(data, tau_levels, results)
  
  # Compute phi(tau_k)
  compute_phi <- function(xi, tau_levels, data) {
    phi <- lapply(tau_levels, function(tau_k) {
      index <- which(tau_levels == tau_k)
      xi_tau_k <- xi[, index]
      f_hat <- density(xi, bw = "nrd0", n = nrow(data))$y
      phi_matrix <- diag(f_hat)
      return(phi_matrix)
    })
    names(phi) <- paste0("tau_", tau_levels)
    return(phi)
  }
  phi <- compute_phi(xi, tau_levels, data)
  
  # Compute M_Zk
  N <- length(unique(data$Provinsi))
  T <- length(unique(data$Time))
  Z <- kronecker(diag(N), rep(1, T))
  
  compute_MZ <- function(Z, phi, tau_levels) {
    MZ <- list()
    
    for (k in seq_along(phi)) {
      phi_tau_k <- phi[[k]]
      
      # Compute P_Zk = Z (Z' * phi * Z)^(-1) Z' * phi
      Z_phi_Z <- t(Z) %*% phi_tau_k %*% Z
      P_Zk <- Z %*% diag(1/(diag(Z_phi_Z))) %*% t(Z) %*% phi_tau_k
      
      # Compute M_Zk = I - P_Zk
      I <- diag(nrow(P_Zk))
      M_Zk <- I - P_Zk
      
      # Store the result with a descriptive name
      MZ[[paste0("MZ_tau_", tau_levels[k])]] <- M_Zk
    }
    return(MZ)
  }
  M_Zk <- compute_MZ(Z, phi, tau_levels)
  
  # Compute J matrices
  compute_J_matrices <- function(data, results, tau_levels, tau_weights) {
    # Prepare covariates and instruments
    covariate_names <- c("exp_lag", "tpt_lag", "gini_lag")
    instrument_names <- c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                          "exp_lag3", "tpt_lag3", "gini_lag3")
    covariates <- as.matrix(data[, covariate_names, drop = FALSE])
    covariates_with_beta0 <- cbind(1, covariates)
    instruments <- as.matrix(data[, instrument_names, drop = FALSE])
    X_tilde <- cbind(instruments, covariates_with_beta0)
    
    # Initialize matrices to store J_alpha and J_theta for each tau level
    J_alpha_all <- list()
    J_theta_all <- list()
    
    # Loop over each tau_k to retrieve fixed effects, alpha, and quantile coefficients
    fixed_effects <- list()
    alpha <- list()
    quantile_coefs <- list()
    beta <- list()
    gamma <- list()
    for (k in seq_along(tau_levels)) {
      tau_k <- tau_levels[k]
      fixed_effects[[k]] <- results[[k]]$fixed_effects
      alpha[[k]] <- results[[k]]$alpha
      coefs <- results[[k]]$quantile_coefs
      beta[[k]] <- coefs[1:4]
      gamma[[k]] <- coefs[5:10]
    }
    
    # Compute residuals for all tau levels
    residuals <- matrix(0, nrow = NT, ncol = k)
    for (k in seq_along(tau_levels)) {
      residuals[, k] <- data$log_ump - Z %*% results[[k]]$fixed_effects - 
        results[[k]]$alpha * data$log_ump_lag - covariates_with_beta0 %*%beta[[k]]
    }
    
    # Compute bandwidth h dynamically using Silverman's rule of thumb
    bandwidths <- numeric(k)
    for (k in seq_along(tau_levels)) {
      u <- residuals[, k]
      sigma_u <- sd(u)
      bandwidths[k] <- 1.06 * sigma_u * NT^(-0.2)
    }
    
    # Loop over each tau_k to compute J_alpha and J_theta
    for (k in seq_along(tau_levels)) {
      tau_k <- tau_levels[k]
      weight_k <- tau_weights[k]
      M_Z_k <- M_Zk[[k]]
      
      # Extract residuals and bandwidth for the current tau level
      u <- residuals[, k]
      h <- bandwidths[k]
      
      # Indicator function for |u| ≤ h
      indicator <- ifelse(abs(u) <= h, 1, 0)
      
      # Compute J_alpha and J_theta
      y_lag <- as.vector(data$log_ump_lag)
      
      J_alpha_k <- weight_k * t(X_tilde) %*% t(M_Z_k) %*% diag(as.vector(indicator)) %*% M_Z_k %*% y_lag / (2*NT*h)
      J_theta_k <- weight_k * t(X_tilde) %*% t(M_Z_k) %*% diag(as.vector(indicator)) %*% M_Z_k %*% X_tilde / (2*NT*h)
      
      # Store J_alpha and J_theta for the current tau
      colnames(J_alpha_k) <- c("alpha")
      J_alpha_all[[k]] <- J_alpha_k
      J_theta_all[[k]] <- J_theta_k
    }
    
    # Partition J_theta into J_beta and J_gamma for each tau
    p <- ncol(covariates_with_beta0)
    q <- ncol(instruments)
    
    J_beta_bar_all <- list()
    J_gamma_bar_all <- list()
    
    for (k in seq_along(tau_levels)) {
      J_theta_k <- J_theta_all[[k]]
      
      J_gamma_bar <- J_theta_k[1:q, 1:(q + p)]
      J_beta_bar <- J_theta_k[(q + 1):(q + p), 1:(q + p)]
      
      J_beta_bar_all[[k]] <- J_beta_bar
      J_gamma_bar_all[[k]] <- J_gamma_bar
    }
    
    # Return results for all tau levels
    return(list(
      J_alpha = J_alpha_all,
      J_theta = J_theta_all,
      J_beta_bar = J_beta_bar_all,
      J_gamma_bar = J_gamma_bar_all,
      bandwidth = bandwidths
    ))
  }
  
  # Compute K, L, and S matrices
  compute_k_l_s <- function(J_alpha, J_beta_bar, J_gamma_bar, tau_levels, tau_weights, data) {
    NT <- nrow(data)
    k <- length(tau_levels)
    
    compute_V <- function(data, tau_levels, tau_weights) {
      # Initialize V
      V_all <- vector("list", k)
      
      # Prepare covariates and instruments
      covariate_names <- c("exp_lag", "tpt_lag", "gini_lag")
      instrument_names <- c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                            "exp_lag3", "tpt_lag3", "gini_lag3")
      covariates <- as.matrix(data[, covariate_names, drop = FALSE])
      covariates_with_beta0 <- cbind(1, covariates)
      instruments <- as.matrix(data[, instrument_names, drop = FALSE])
      X_tilde <- cbind(instruments, covariates_with_beta0)
      
      # Loop over each tau
      for (k in seq_along(tau_levels)) {
        tau_k <- tau_levels[k]
        weight_k <- tau_weights[k]
        M_Z_k <- M_Zk[[k]]
        
        # Compute the k-th row of V
        V <- weight_k * t(X_tilde) %*% t(M_Z_k)
        
        # Store J_alpha and J_theta for the current tau
        V_all[[k]] <- V
      }
      names(V_all) <-  paste0("tau_", as.character(tau_levels))
      return(V_all)
    }
    V <- compute_V(data, tau_levels, tau_weights)
    
    compute_S <- function(tau_levels, tau_weights, V) {
      k <- length(tau_levels)
      S_all <- vector("list", k)
      
      for (k in seq_along(tau_levels)) {
        tau_k <- tau_levels[k]
        expectation <- (V[[k]] %*% t(V[[k]])) / NT
        S <- (tau_k - tau_k^2) * expectation
        
        # Store the result
        S_all[[k]] <- S
        names(S_all)[k] <- paste0("tau_", tau_k)
      }
      return(S_all)
    }
    S <- compute_S(tau_levels, tau_weights, V)
    
    # Initialize lists for K and L matrices
    K_all <- vector("list", k)
    L_all <- vector("list", k)
    
    # Compute K and L for each tau
    for (k in seq_along(tau_levels)) {
      weight_k <- tau_weights[k]
      
      J_alpha_k <- J_alpha[[k]]
      J_beta_k <- J_beta_bar[[k]]
      J_gamma_k <- J_gamma_bar[[k]]
      S_k <- S[[k]]
      
      A_alpha <- solve(J_gamma_k %*% S_k %*% t(J_gamma_k))
      H <- t(J_gamma_k) %*% A_alpha %*% J_gamma_k
      
      # Compute K
      K_k <- solve(t(J_alpha_k) %*% H %*% J_alpha_k) %*% t(J_alpha_k) %*% H
      
      # Compute L
      M <- diag(J_alpha_k %*% K_k) - (J_alpha_k %*% K_k)
      L_k <- J_beta_k %*% M
      
      # Store results
      K_all[[k]] <- K_k
      L_all[[k]] <- L_k
      S[[k]] <- S_k
    }
    # Return results for all tau levels
    names(K_all) <-  paste0("tau_", as.character(tau_levels))
    names(L_all) <-  paste0("tau_", as.character(tau_levels))
    
    return(list(
      K = K_all,
      L = L_all,
      S = S
    ))
  }
  
  compute_omega_tau <- function(K, L, S, tau_levels) {
    # Initialize list to store omega for each tau
    Omega_all <- list()
    
    # Loop over each tau level
    for (k in seq_along(tau_levels)) {
      K_k <- as.matrix(K[[k]])
      L_k <- as.matrix(L[[k]])
      S_k <- as.matrix(S[[k]])
      
      # Compute Omega(tau_k)
      KL_k <- cbind(t(K_k), t(L_k))
      Omega_k <- t(KL_k) %*% S_k %*% KL_k 
      
      # Store in the list
      Omega_all[[k]] <- Omega_k
      names(Omega_all)[k] <- paste0("tau_", tau_levels[k])
    }
    return(Omega_all)
  }
  
  # Combine all computation
  J <- compute_J_matrices(data, results, tau_levels, tau_weights)
  
  KLS <- compute_k_l_s(J_alpha = J$J_alpha, 
                       J_beta_bar = J$J_beta_bar, 
                       J_gamma_bar = J$J_gamma_bar, 
                       tau_levels, tau_weights, data)
  
  # Compute omega for each tau
  Omega <- compute_omega_tau(KLS$K, KLS$L, KLS$S, tau_levels)
  return(Omega)
}

Omega <- compute_omega(results = results, 
                       data = df_panel, 
                       tau_levels = taus,
                       tau_weights = tauw)
print(Omega)

## Wald test -------------------------------------------------------------------
wald_test <- function(results, omega, tau_levels, data) {
  NT <- nrow(data)
  
  # Initialize a data frame to store results
  wald_results <- data.frame(
    tau = numeric(),
    parameter = character(),
    p_value = numeric(),
    significance = character()
  )
  
  # Loop over each tau level to perform the test
  for (k in seq_along(tau_levels)) {
    tau_k <- tau_levels[k]
    
    # Extract the covariance matrix
    omega_k <- omega[[k]]
    
    # Extract theta (alpha and beta coefficients)
    theta <- c(results[[k]]$alpha, 
               results[[k]]$quantile_coefs[1:4])
    
    # Define parameter names
    param_names <- c("alpha", "intercept", "L1_expenses", "L1_tpt", "L1_gini")
    
    # Loop over each parameter in theta and test its significance
    for (i in seq_along(theta)) {
      param_name <- param_names[i]
      R <- matrix(0, nrow = 1, ncol = length(theta))
      R[1, i] <- 1
      r <- 0
      
      # Calculate the Wald statistic
      R_omega_R <- R %*% omega_k %*% t(R)
      R_theta_r <- (R %*% theta) - r
      wald_statistic <- NT * t(R_theta_r) %*% solve(R_omega_R) %*% R_theta_r
      
      # Calculate the p-value
      p_value <- 1 - pchisq(wald_statistic, df = 1)
      
      # Determine significance
      significant <- ifelse(p_value < 0.05, "Significant", "Not Significant")
      
      # Add results to the data frame
      wald_results <- rbind(wald_results, data.frame(
        tau = tau_k,
        parameter = param_name,
        statistics = wald_statistic,
        p_value = p_value,
        significance = significant))
    }
  }
  return(wald_results)
}
wald_test(results = results, 
          omega = Omega, 
          tau_levels = taus,
          data = df_panel)

# Predictions ------------------------------------------------------------------
generate_predictions <- function(results, data, actual_col = "actual_y", 
                                 initial_values, actual_is_differenced = FALSE) {
  # Placeholder for predictions and metrics
  predictions <- list()
  metrics <- list()
  
  # Loop over each quantile in results
  for (tau_level in names(results)) {
    # Extract quantile-specific components
    tau_result <- results[[tau_level]]
    alpha <- tau_result$alpha
    quantile_coefs <- tau_result$quantile_coefs
    fixed_effects <- tau_result$fixed_effects
    
    # Separate intercept and beta coefficients
    intercept <- quantile_coefs[1]
    beta <- quantile_coefs[2:4]
    
    # Map fixed effects for groups in data
    fixed_effect_names <- sapply(names(fixed_effects), function(name) {
      name <- gsub("Provinsi", "", name)      # Remove "Provinsi"
      name <- gsub("\\.", " ", name)         # Replace remaining dots with spaces
      name <- gsub("KEP\\s{2}", "KEP. ", name) #Replace KEP with KEP.
      name <- gsub("^\\s+|\\s+$", "", name)   # Remove leading/trailing spaces
      return(name)
    })
    names(fixed_effects) <- fixed_effect_names  # Update names to match
    
    data$fixed_effect <- fixed_effects[data$Provinsi]
    
    # Compute differenced predictions
    data$predicted_diff <- intercept + 
      alpha * data$log_ump_lag + 
      as.matrix(data[, c("exp_lag", "tpt_lag", "gini_lag")]) %*% beta + 
      data$fixed_effect
    
    # Undifference the predictions
    data$predicted_y <- initial_values$log_ump_lag + data$predicted_diff
    
    # Store predictions for this quantile
    predictions[[tau_level]] <- data$predicted_y
    
    # Calculate metrics
    if (!is.null(actual_col) && actual_col %in% names(data)) {
      # If actual values are differenced, undifference them
      if (actual_is_differenced) {
        actual_diff <- data[[actual_col]]  # Differenced actual values
        data$actual_y <- initial_values$log_ump_lag + actual_diff  # Undifferencing actual values
      } else {
        data$actual_y <- data[[actual_col]] # Use the values directly
      }
      
      # Compute metrics
      actual <- data$actual_y
      predicted <- data$predicted_y
      
      MAE <- MLmetrics::MAE(exp(predicted), exp(actual))
      RMSE <- MLmetrics::RMSE(exp(predicted), exp(actual))
      MAPE <- MLmetrics::MAPE(exp(predicted), exp(actual)) * 100  # As a percentage
      
      # Store metrics for this quantile
      metrics[[tau_level]] <- list(MAE = MAE, RMSE = RMSE, MAPE = MAPE)
    }
  }
  # Combine predictions into a data frame
  prediction_df <- cbind(data, predictions)
  prediction_df <- prediction_df[, c(names(results))]
  
  # Return the predictions and metrics
  return(list(predictions = exp(prediction_df), metrics = metrics))
}

# Create lagged of wages
new_lag_data <- data %>%
  mutate(
    log_ump_lag = lag(log_ump)
  ) %>%
  na.omit() %>%
  filter(Time >= 2013 & Time <= 2023) 

# Prepare the initial values
initial_values <- new_lag_data %>% 
  filter(Time >= 2013 & Time <= 2023) %>%
  dplyr::select(Provinsi, log_ump_lag)

# Prepare df_panel without 2024
df_panel$Time <- as.numeric(as.character(df_panel$Time))
df_panel2 <- df_panel %>%
  filter(Time >= 2013 & Time <= 2023) %>%
  dplyr::select(Provinsi, Time, log_ump,log_ump_lag, exp_lag, tpt_lag, gini_lag)

## Generate predictions
output <- generate_predictions(results, df_panel2, actual = "log_ump", initial_values = initial_values, actual_is_differenced = TRUE)

predictions <- output$predictions
print(predictions)
print(output$metrics)

# Plot the predictions
predictions$Provinsi <- rownames(predictions)  # Add Province column
predictions <- transform(predictions, 
                         Provinsi = rep(unique(data$Provinsi), each = 11),  # Extract Provinsi
                         Time = as.numeric(sub(".*-", "", rownames(predictions))) + 1)  # Increment year for t+1
predictions_long <- melt(predictions, 
                         id.vars = c("Provinsi", "Time"), 
                         variable.name = "Quantile", 
                         value.name = "Prediction")
actual <- data %>%
  filter(Time >= 2014 & Time <= 2024) %>%
  dplyr::select(Provinsi, Time, ump)

# Merge actual values with predictions
comparison_data <- merge(predictions_long, actual, by = c("Provinsi", "Time"))

# Plot
create_province_plot <- function(province_data) {
  ggplot(province_data, aes(x = Time, y = Prediction, color = Quantile, group = Quantile)) +
    geom_line(size = 1) +  # Quantile prediction lines
    geom_line(aes(y = ump, color = "actual"), size = 1, linetype = "dashed") +  # Actual values
    labs(
      title = paste("Quantile Predictions vs Actual Values for", province_data$Provinsi[1]),
      x = "Year",
      y = "Wages"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      strip.text = element_text(margin = margin(b = 5)), # Keep facet labels with margin
      panel.spacing = unit(1, "lines") # Add spacing between facets
    )
}
province_plots <- split(comparison_data, comparison_data$Provinsi)

# Create a folder to store the plots
dir.create("province_plots", showWarnings = FALSE)
for (province in names(province_plots)) {
  plot <- create_province_plot(province_plots[[province]])
  ggsave(paste0("province_plots/", province, ".jpg"), plot, width = 6, height = 4)
}

## Forecast 2025 ---------------------------------------------------------------
forecast <- function(results, data, initial_values) {
  # Placeholder for predictions and metrics
  predictions <- list()
  
  # Loop over each quantile in results
  for (tau_level in names(results)) {
    # Extract quantile-specific components
    tau_result <- results[[tau_level]]
    alpha <- tau_result$alpha
    quantile_coefs <- tau_result$quantile_coefs
    fixed_effects <- tau_result$fixed_effects
    
    # Separate intercept and beta coefficients
    intercept <- quantile_coefs[1]
    beta <- quantile_coefs[2:4]
    
    # Map fixed effects for groups in data
    fixed_effect_names <- sapply(names(fixed_effects), function(name) {
      name <- gsub("Provinsi", "", name)      # Remove "Provinsi"
      name <- gsub("\\.", " ", name)         # Replace remaining dots with spaces
      name <- gsub("KEP\\s{2}", "KEP. ", name) #Replace KEP with KEP.
      name <- gsub("^\\s+|\\s+$", "", name)   # Remove leading/trailing spaces
      return(name)
    })
    names(fixed_effects) <- fixed_effect_names  # Update names to match
    
    data$fixed_effect <- fixed_effects[data$Provinsi]
    
    # Compute differenced predictions
    data$predicted_diff <- intercept + 
      alpha * data$log_ump_lag + 
      as.matrix(data[, c("exp_lag", "tpt_lag", "gini_lag")]) %*% beta + 
      data$fixed_effect
    
    # Undifference the predictions
    data$predicted_y <- initial_values$log_ump_lag + data$predicted_diff
    
    # Store predictions for this quantile
    predictions[[tau_level]] <- data$predicted_y
    }
  # Combine predictions into a data frame
  prediction_df <- cbind(data, predictions)
  prediction_df <- prediction_df[, c(names(results))]
  
  # Return the predictions and metrics
  return(predictions = exp(prediction_df))
}

# Import actual wages 2024
actual <- data.frame(
  Provinsi = unique(rawdata$Provinsi),
  ump = rawdata[rawdata$Time == 2024, ]$ump)

# Log transformation of wages
actual <- actual %>%
  mutate(
    log_ump = log(ump)
  )

# Create a new dataset to provide the covariates
covariates <- df_panel[df_panel$Time == 2024, ] %>%
  dplyr::select(Provinsi, Time, log_ump,log_ump_lag, exp_lag, tpt_lag, gini_lag)

# Add actual log wages to differenced dataset
covariates$log_ump <- actual$log_ump
head(covariates)

# Prepare the initial values of 2023
initial_2024 <- data  %>%
  mutate(
    log_ump_lag = lag(log_ump)
  ) %>%
  na.omit() %>%
  filter(Time == 2024) %>%
  dplyr::select(Provinsi, log_ump_lag)

# Generate predictions for 2025
predicted_2025 <- forecast(
  results = results, 
  data = covariates, 
  initial_values = initial_2024)
print(predicted_2025)

## Interval --------------------------------------------------------------------
interval_prediction <- function(predictions, results, data, tau_levels, confidence_level = 0.95) {
  # Prepare the components needed
  N <- length(unique(data$Provinsi))
  T <- length(unique(data$Time))
  NT <- nrow(data)
  
  # Define X_check and X_dot
  covariate_names <- c("exp_lag", "tpt_lag", "gini_lag")
  instrument_names <- c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                        "exp_lag3", "tpt_lag3", "gini_lag3")
  covariates <- as.matrix(data[, covariate_names, drop = FALSE])
  X <- cbind(1, covariates)
  W <- as.matrix(data[, instrument_names, drop = FALSE])
  Z <- kronecker(diag(N), rep(1, T))
  X_check <- cbind(Z, W, X)
  
  lastest_data <- data[data$Time == 2024, ]
  z <- diag(1, N)
  y_t <- as.vector(lastest_data$log_ump_lag)
  X_t1 <- cbind(1, as.matrix(lastest_data[, covariate_names, drop = FALSE]))
  X_dot <- cbind(z, y_t, X_t1)
  
  # Define D_NT matrix
  D <- diag(c(
    rep(sqrt(T), N), 
    sqrt(NT), 
    rep(sqrt(NT), ncol(X)))
  )
  
  # Compute Xi_it for each tau
  compute_xi <- function(data, tau_levels, results) {
    # Extract necessary components from the result
    fixed_effects <- list()
    alpha <- list()
    quantile_coefs <- list()
    beta <- list()
    gamma <- list()
    for (k in seq_along(tau_levels)) {
      tau_k <- tau_levels[k]
      fixed_effects[[k]] <- results[[k]]$fixed_effects
      alpha[[k]] <- results[[k]]$alpha
      coefs <- results[[k]]$quantile_coefs
      beta[[k]] <- coefs[1:4]
      gamma[[k]] <- coefs[5:10]
    }
    
    # Prepare covariates and instruments
    covariate_names <- c("exp_lag", "tpt_lag", "gini_lag")
    instrument_names <- c("exp_lag2", "tpt_lag2", "log_ump_lag3",
                          "exp_lag3", "tpt_lag3", "gini_lag3")
    covariates <- as.matrix(data[, covariate_names, drop = FALSE])
    covariates_with_beta0 <- cbind(1, covariates)
    instruments <- as.matrix(data[, instrument_names, drop = FALSE])
    
    # Initialize a matrix to store xi values for each tau
    xi <- matrix(NA, nrow = NT, ncol = k)
    
    # Loop over taus and compute xi
    for (k in seq_along(tau_levels)) {
      # Extract tau-specific coefficients
      fixed_effects <- list()
      alpha <- list()
      quantile_coefs <- list()
      beta <- list()
      gamma <- list()
      for (k in seq_along(tau_levels)) {
        tau_k <- tau_levels[k]
        fixed_effects[[k]] <- results[[k]]$fixed_effects
        alpha[[k]] <- results[[k]]$alpha
        coefs <- results[[k]]$quantile_coefs
        beta[[k]] <- coefs[1:4]
        gamma[[k]] <- coefs[5:10]
      }
      
      # Compute xi
      for (k in seq_along(tau_levels)) {
        xi[, k] <- results[[k]]$fixed_effects + 
          results[[k]]$alpha  * data$log_ump_lag +
          covariates_with_beta0 %*% beta[[k]] +
          instruments %*% gamma[[k]]
      }
    }
    colnames(xi) <- paste0("tau_", tau_levels)
    return(xi)
  }
  xi <- compute_xi(data, tau_levels, results)
  
  # Compute phi(tau_k)
  compute_phi <- function(xi, tau_levels, data) {
    phi <- lapply(tau_levels, function(tau_k) {
      index <- which(tau_levels == tau_k)
      xi_tau_k <- xi[, index]
      f_hat <- density(xi, bw = "nrd0", n = nrow(data))$y
      phi_matrix <- diag(f_hat)
      return(phi_matrix)
    })
    names(phi) <- paste0("tau_", tau_levels)
    return(phi)
  }
  phi <- compute_phi(xi, tau_levels, data)
  
  # Initialize lists to store components
  intervals <- list()
  
  # Loop over each quantile level
  for (tau in tau_levels) {
    # Predicted quantiles for this tau
    tau_col <- paste0("tau_", tau)
    predicted_quantile <- predictions[[tau_col]]
    
    # Compute residual density (phi) at tau
    phi_tau <- phi[[tau_col]]
    
    # Compute J(tau) using the expectation
    y_lag <- as.vector(data$log_ump_lag)
    J_tau <-  t(X_check) %*% phi_tau %*% cbind(Z, y_lag, X) / NT
    
    # Compute S(tau)
    expectation <- (t(X_check) %*% X_check) / NT
    S_tau <- tau * (1 - tau) * expectation
    
    # Compute Sigma(tau)
    if (nrow(J_tau) != ncol(J_tau)) {
      library(MASS)
      J_tau_inv <- ginv(J_tau)
    } else {
      J_tau_inv <- solve(J_tau)
    }
    
    sigma_tau <- J_tau_inv %*% S_tau %*% t(J_tau_inv)
    
    # Compute Omega(tau)
    omega_tau <- X_dot %*% solve(D) %*% sigma_tau %*% solve(D) %*% t(X_dot)
    
    # Compute a_NT
    lambda <- (1 - confidence_level) / 2
    z_lambda <- qnorm(1 - lambda)
    a_NT <- z_lambda * sqrt(diag(omega_tau))
    
    # Compute the prediction interval
    lower_bound <- predicted_quantile - a_NT
    upper_bound <- predicted_quantile + a_NT
    
    # Store the results
    intervals[[tau_col]] <- data.frame(
      Lower = lower_bound,
      Upper = upper_bound
    )
  }
  return(intervals)
}

# Compute prediction intervals
interval <- interval_prediction(predictions = output$predictions, 
                                results = results, 
                                data = df_panel, 
                                tau_levels = taus)
print(interval)

interval2025 <- interval_prediction(predictions = predicted_2025, 
                                results = results, 
                                data = df_panel, 
                                tau_levels = taus)
print(interval2025)

# Combine the datasets
combine_intervals <- function(historical, forecast) {
  combined <- list()
  
  # Loop through quantile levels (keys of the lists)
  for (quantile in names(historical)) {
    # Bind historical and forecast data frames for the current quantile
    combined[[quantile]] <- rbind(
      historical[[quantile]] %>% mutate(Quantile = quantile),
      forecast[[quantile]] %>% mutate(Quantile = quantile)
    )
  }
  
  # Combine all quantiles into a single data frame
  combined_df <- do.call(rbind, combined)
  rownames(combined_df) <- NULL
  return(combined_df)
}

combined_intervals <- combine_intervals(interval, interval2025)
print(combined_intervals)

# Plot the forecast with the interval
forecast_data <- read.csv("forecast.csv")
forecast_data <- forecast_data[, -1]
forecast_data$Time <- as.numeric(forecast_data$Time)
str(forecast_data)

# Define the folder to save the plots
output_folder <- "forecast_plots"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Loop through each province
unique_provinces <- unique(forecast_data$Provinsi)

for (province in unique_provinces) {
  province_data <- forecast_data %>% filter(Provinsi == province)
  
  # Get the last historical ump value
  last_historical_data <- province_data %>%
    filter(Time < 2025) %>%
    slice_max(Time)
  
  plot <- ggplot(province_data, aes(x = Time)) +
    # Historical data line
    geom_line(
      data = province_data %>% filter(Time < 2025),
      aes(y = ump),
      color = "blue", size = 1
    ) +
    # Connecting lines
    geom_segment(
      data = province_data %>% filter(Time == 2025),
      aes(
        xend = Time, yend = ump,
        x = last_historical_data$Time, y = last_historical_data$ump,
        color = Quantile # Color the lines
      ),
      linetype = "dashed"
    ) +
    # Forecast points
    geom_point(
      data = province_data %>% filter(Time == 2025),
      aes(y = ump, color = Quantile),
      size = 2
    ) +
    # Labels and themes
    labs(
      title = paste("Forecast of UMP for", province),
      x = "Year",
      y = "UMP",
      color = "Quantile Levels",
      fill = "Confidence Interval"
    ) +
    scale_color_manual(
      values = c(
        "tau_0.1" = "orange",
        "tau_0.25" = "purple",
        "tau_0.5" = "red",
        "tau_0.75" = "darkgreen",
        "tau_0.9" = "brown"
      )
    ) +
    scale_x_continuous(breaks = unique(province_data$Time))+
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave(
    filename = paste0(output_folder, "/", province, ".jpg"),
    plot = plot,
    width = 12,
    height = 6
  )
}

# Summary of the forecast
summary <- list()
for (quantile in unique(forecast_data$Quantile)) {
  # Filter data for the current quantile
  quantile_data <- forecast_data %>% filter(Quantile == quantile)
  
  # Summary statistics
  description_2025 <- describe(quantile_data %>%
                                 filter(Time == 2025) %>%
                                 group_by(Provinsi))
  
  # Calculate UMP for 2024 and 2025
  ump_2024_2025 <- quantile_data %>%
    filter(Time %in% c(2024, 2025)) %>%
    group_by(Provinsi) %>%
    arrange(Time) %>%
    summarize(
      ump_2024 = first(ump),
      ump_2025 = last(ump)
    )
  
  # Calculate growth from 2024 to 2025
  ump_growth <- ump_2024_2025 %>%
    mutate(
      ump_growth = (ump_2025 - ump_2024) / ump_2024 * 100
    )
  
  # Find top 10 provinces with highest and lowest growth
  top_5_highest_growth <- ump_growth %>%
    arrange(desc(ump_growth)) %>%
    slice(1:5)
  
  top_5_lowest_growth <- ump_growth %>%
    arrange(ump_growth) %>%
    slice(1:5)
  
  # Store the results in a list
  summary[[paste0("Quantile_", quantile)]] <- list(
    description_2025 = description_2025,
    growth_data = ump_growth,
    top_5_highest = top_5_highest_growth,
    top_5_lowest = top_5_lowest_growth
  )
}

# Print the results for each quantile
for (quantile_name in names(summary)) {
  cat("\n--- Results for", quantile_name, "---\n")
  cat("\nSummary Statistics for Forecast 2025:\n")
  print(summary[[quantile_name]]$description_2025)
  cat("\nTop 5 Highest Growth:\n")
  print(summary[[quantile_name]]$top_5_highest)
  cat("\nTop 5 Lowest Growth:\n")
  print(summary[[quantile_name]]$top_5_lowest)
}



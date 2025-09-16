# Load necessary R libraries
library(numDeriv)
library(openxlsx)
library(english)
library(SMPracticals)
library(readxl)
library(plotly)
library(RColorBrewer)
library(survival)
library(processx)
library(orca)
library(BoutrosLab.plotting.general)

# Load the dataset
dataset <- read_excel("dataset_ssm.xlsx")
dataset$time = dataset$time

# Separate data by stress stages
stress_data_stage_1 <- subset(dataset, stress == 1)
stress_data_stage_2 <- subset(dataset, stress == 2)
stress_data_stage_3 <- subset(dataset, stress == 3)

# a list of failure times for each stage
t = list(c(stress_data_stage_1$time), c(stress_data_stage_2$time), c(stress_data_stage_3$time))

#  number of censored observations at each stage
r_1 = 0
r_2 = 0
r_3 = 9

#  sample size, number of stages, and stress change times
n = dim(dataset)[1]
m = 3
tau = c(0, 1.25, 1.41, 1.54)

#  cause identifiers
cause_1 = 2
cause_2 = 3
common_cause = 1

# Count failures for each cause at stage 1
n_11 = table(stress_data_stage_1$cause, stress_data_stage_1$stress)[cause_1, ]
n_12 = table(stress_data_stage_1$cause, stress_data_stage_1$stress)[cause_2, ]
n_10 = table(stress_data_stage_1$cause, stress_data_stage_1$stress)[common_cause, ]

# Count failures for each cause at stage 2
n_21 = table(stress_data_stage_2$cause, stress_data_stage_2$stress)[cause_1, ]
n_22 = table(stress_data_stage_2$cause, stress_data_stage_2$stress)[cause_2, ]
n_20 = table(stress_data_stage_2$cause, stress_data_stage_2$stress)[common_cause, ]

# Count failures for stage 3 (only common cause)
n_31 = 0
n_32 = 0
n_30 = table(stress_data_stage_3$cause, stress_data_stage_3$stress)[common_cause, ]

# Combine failure counts into a list
failure_count = list(c(n_10, n_11, n_12), c(n_20, n_21, n_22), c(n_30, n_31, n_32))

# Combine censored counts into a list
r = list(r_1, r_2, r_3)

# Normalize stress levels
S = c(131.5, 142.5, 158)/150


# Define baseline cumulative distribution function
F_not <- function(t) {
  return(1 - exp(-t))
}

# Define  CDF with parameter lambda
D_F <- function(t, lamb) {
  return(1 - exp(-lamb * t))
}

# Log of baseline PDF
ln_f_not <- function(t) {
  return(-t)
}

# Log of baseline survival function
ln_S_not <- function(t) {
  return(-t)
}

# Compute lambda values for each cause or combined
lambdas <- function(C, S, P, i, j = NULL) {
  if (!is.null(j)) {
    # Lambda for specific cause j
    return(C[j] * (S[i]^P))
  } else {
    # Combined lambda across all causes
    return(C[1] * (S[i]^P) + C[2] * (S[i]^P) + C[3] * (S[i]^P))
  }
}


#  C parameter for cause b based on profile likelihood
C_param <- function(b, P, m, t, failure_count, S , tau, r) {
  
  # 
  t1 = 0
  for (i in 1:m) {
    t1 = t1 + failure_count[[i]][b]
  }
  
  # 
  t2 = 0
  for (i in 2:m) {
    k = i - 1
    for (l in 1:k) {
      a = (sum(as.integer(failure_count[[i]])) + r[[i]]) * log(1 - F_not(tau[l+1])) * ((S[l]^P) - (S[l+1]^P))
      t2 = t2 + a
    }
  }
  
  # 
  t3 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t3 = t3 + (S[i]^P) * log(1 - F_not(t[[i]][[j]]))
      }
    }
  }
  
  # 
  t4 = 0
  for (i in 1:m) {
    t4 = t4 + r[[i]] * (S[i]^P) * log(1 - F_not(tau[i+1]))
  }
  
  # 
  return(-t1 / (t2 + t3 + t4))
}


#  the profile log-likelihood 
prof_log_likelihood <- function(P, m, t, failure_count, S, tau, r) {
  
  # 
  c_0 = C_param(1, P, m, t, failure_count, S, tau, r)
  c_1 = C_param(2, P, m, t, failure_count, S, tau, r)
  c_2 = C_param(3, P, m, t, failure_count, S, tau, r)
  
  # 
  C = c_0 + c_1 + c_2
  C_vec = c(c_0, c_1, c_2)
  
  # 
  t1 = 0
  for (j in 1:3) {
    for (i in 1:m) {
      t1 = t1 + log(C_vec[j]) * failure_count[[i]][j]
    }
  }
  
  # 
  t2 = 0
  for (i in 1:m) {
    t2 = t2 + P * log(S[i]) * sum(as.integer(failure_count[[i]]))
  }
  
  # 
  t3 = 0
  for (i in 2:m) {
    k = i - 1
    for (l in 1:k) {
      a = C * (sum(as.integer(failure_count[[i]])) + r[[i]]) * log(1 - F_not(tau[l+1])) * ((S[l]^P) - (S[l+1]^P))
      t3 = t3 + a
    }
  }
  
  # 
  t4 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t4 = t4 + (C * (S[i]^P) - 1) * log(1 - F_not(t[[i]][[j]]))
      }
    }
  }
  
  # 
  t5 = 0
  for (i in 1:m) {
    t5 = t5 + C * r[[i]] * (S[i]^P) * log(1 - F_not(tau[i+1]))
  }
  
  # 
  t6 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t6 = t6 + ln_f_not(t[[i]][[j]])
      }
    }
  }
  
  # 
  return(t1 + t2 + t3 + t4 + t5 + t6)
}


#  a sequence of candidate P values
P_values <- seq(from = 0.001, to = 20, by = 0.05)

# Initialize vector to store profile log-likelihood 
z <- c()

# Compute profile log-likelihood 
for (i in 1:length(P_values)) {
  pl <- prof_log_likelihood(P_values[i], m, t, failure_count, S, tau, r)
  z <- c(z, pl)
}

#  maximum log-likelihood and corresponding P
max_value <- max(z)
max_position <- which(z == max_value, arr.ind = TRUE)
opt_p <- P_values[max_position[1]]

#  optimal C parameters and corresponding lambda vector
opt_C <- C_param(1, opt_p, m, t, failure_count, S, tau, r) +
  C_param(2, opt_p, m, t, failure_count, S, tau, r) +
  C_param(3, opt_p, m, t, failure_count, S, tau, r)
opt_lamb <- opt_C * (S^opt_p)

#  data for contour plot
contour_dataframe <- data.frame(P_values = P_values, log_likelihood = z)

#  profile log-likelihood plot
contour_gg <- ggplot(contour_dataframe, aes(x = P_values)) +
  geom_line(aes(y = z, color = "Profile log-likelihood"), size = 1) +
  labs(x = "P", y = "Profile log-likelihood") +
  scale_color_manual(name = NULL, values = c("Profile log-likelihood" = "black")) +
  theme_minimal() +
  theme(
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    legend.position = c(0.75, 0.25),
    legend.text = element_text(size = 12, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.line = element_line(color = "black")
  )

# Display the plot
contour_gg

# Save the profile log-likelihood plot 

ggsave(
  filename = "contour_frm_exponential.pdf",
  plot = contour_gg,
  width = 6,            
  height = 4,           
  dpi = 600,            
  units = "in",         
  device = cairo_pdf,   
  bg = "white"          
)


# Set initial parameter for optimization
initial_params <- c(opt_p)

# Define objective function for optimization (negative profile log-likelihood)
objective_func <- function(params) {
  P <- params[1]  
  return(-prof_log_likelihood(P, m, t, failure_count, S , tau, r))
}

# Set bounds for optimization
lower_bound <- c(0.001)
upper_bound <- c(20)

# Perform constrained optimization 
result_1 <- optim(par = initial_params, fn = objective_func, method = "L-BFGS-B", 
                  lower = lower_bound, upper = upper_bound, hessian = TRUE)

# Extract estimated parameter
estimated_P <- result_1$par[1]

#  corresponding C parameters for estimated P
estimated_c_0 <- C_param(1, estimated_P, m, t, failure_count, S , tau, r)
estimated_c_1 <- C_param(2, estimated_P, m, t, failure_count, S , tau, r)
estimated_c_2 <- C_param(3, estimated_P, m, t, failure_count, S , tau, r)

#  full log-likelihood function
log_likelihood <- function(P, c_0, c_1, c_2, m, t, failure_count, S , tau, r) {
  C <- c_0 + c_1 + c_2
  C_vec <- c(c_0, c_1, c_2)
  
  t1 <- 0
  for (j in 1:3) {
    for (i in 1:m) {
      t1 <- t1 + log(C_vec[j]) * failure_count[[i]][j]
    }
  }
  
  t2 <- 0
  for (i in 1:m) {
    t2 <- t2 + P * log(S[i]) * sum(as.integer(failure_count[[i]]))
  }
  
  t3 <- 0
  for (i in 2:m) {
    k <- i - 1
    for (l in 1:k) {
      a <- C * (sum(as.integer(failure_count[[i]])) + r[[i]]) *
        log(1 - F_not(tau[l + 1])) * ((S[l]^P) - (S[l + 1]^P))
      t3 <- t3 + a
    }
  }
  
  t4 <- 0
  for (i in 1:m) {
    k <- sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t4 <- t4 + (C * (S[i]^P) - 1) * log(1 - F_not(t[[i]][[j]]))
      }
    }
  }
  
  t5 <- 0
  for (i in 1:m) {
    t5 <- t5 + C * r[[i]] * (S[i]^P) * log(1 - F_not(tau[i + 1]))
  }
  
  t6 <- 0
  for (i in 1:m) {
    k <- sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t6 <- t6 + ln_f_not(t[[i]][[j]])
      }
    }
  }
  
  return(t1 + t2 + t3 + t4 + t5 + t6)
}

#  function for Hessian-based variance estimation
objective_func_for_hessian <- function(params) {
  P <- params[1] 
  c_0 <- params[2]
  c_1 <- params[3]
  c_2 <- params[4]
  return(-log_likelihood(P, c_0, c_1, c_2, m, t, failure_count, S , tau, r))
}

#  estimated parameters into a single vector
param <- c(estimated_P, estimated_c_0, estimated_c_1, estimated_c_2)

#  variance using numerical Hessian
variance <- diag(SMPracticals::ihess(objective_func_for_hessian, param, ep = 1e-05))

#  95% confidence intervals for parameters
indicator1 <- param - qnorm(p = 0.05 / 2, lower.tail = FALSE) * sqrt(variance)
indicator2 <- param + qnorm(p = 0.05 / 2, lower.tail = FALSE) * sqrt(variance)


#  minimum CDF based on estimated parameters
minimum_cdf <- function(t) {
  P = param[1]
  C = c(param[2], param[3], param[4])
  lamb = sum(C) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t))^lamb2
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    return(1 - t1*t2)
  } else {
    t1 = (1 - F_not(t))^lamb3
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3]))^(lamb2 - lamb3)
    return(1 - t1*t2*t3)
  }
}

#  empirical data with right-censoring
x1_emp_min = c(dataset$time, rep(1.54, 9))
state_emp_min = c(rep(1, length(dataset$time)), rep(0, 9))
km_fit_min <- survfit(Surv(x1_emp_min, state_emp_min) ~ 1)

# 
x1_emp_min = c(0, dataset$time, rep(1.54, 9))
y1_emp_min = c(); y1_emp_lower_conf_min = c(); y1_emp_upper_conf_min = c()

#  KM survival estimates and confidence intervals
for(i in x1_emp_min) {
  temp = summary(km_fit_min, times = i)
  y1_emp_min = c(y1_emp_min, temp$surv)
  y1_emp_lower_conf_min = c(y1_emp_lower_conf_min, temp$lower)
  y1_emp_upper_conf_min = c(y1_emp_upper_conf_min, temp$upper)
}

#  theoretical CDF
y1_theo_min = c()
for(i in x1_emp_min) {
  y1_theo_min = c(y1_theo_min, 1 - minimum_cdf(i))
}
y1_theo_min[1] = 1

#  data for plotting
min_dataframe = data.frame(
  x1_emp_min = x1_emp_min,
  y1_emp_min = y1_emp_min,
  y1_emp_lower_conf_min = y1_emp_lower_conf_min,
  y1_emp_upper_conf_min = y1_emp_upper_conf_min,
  y1_theo_min = y1_theo_min
)

# Plot empirical vs theoretical minimum CDF
fig_minimum_gg <- ggplot(min_dataframe, aes(x = x1_emp_min)) +
  geom_line(aes(y = y1_emp_min, color = "Empirical Probability"), size = 1) +
  geom_ribbon(aes(ymin = y1_emp_lower_conf_min, ymax = y1_emp_upper_conf_min, fill = "Empirical Confidence Interval"), alpha = 0.5) +
  geom_point(aes(y = y1_emp_min), color = "blue", size = 2.2, show.legend = FALSE) +
  geom_line(aes(y = y1_theo_min, color = "Theoretical Probability"), size = 1) +
  labs(x = "Time", y = "Cumulative Distribution Function") +
  scale_color_manual(name = NULL, values = c("Empirical Probability" = "blue", "Theoretical Probability" = "red")) +
  scale_fill_manual(name = NULL, values = c("Empirical Confidence Interval" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    panel.grid = element_blank(),
    legend.position = c(0.2, 0.2),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  )

fig_minimum_gg


#  marginal CDF for cause 1
marginal_1_cdf <- function(t) {
  P = param[1]
  C = c(param[2], param[3], param[4])
  lamb = sum(C[1] + C[2]) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t))^lamb2
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    return(1 - t1*t2)
  } else {
    t1 = (1 - F_not(t))^lamb3
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3]))^(lamb2 - lamb3)
    return(1 - t1*t2*t3)
  }
}

#  empirical data with right-censoring
x1_emp_margin_1 = c(subset(dataset$time, dataset$cause %in% c(common_cause, cause_1)), rep(1.54, 9))
state_emp_margin_1 = c(rep(1, length(subset(dataset$time, dataset$cause %in% c(common_cause, cause_1)))), rep(0, 9))
km_fit_margin_1 <- survfit(Surv(x1_emp_margin_1, state_emp_margin_1) ~ 1)

x1_emp_margin_1 = c(0, x1_emp_margin_1)
y1_emp_margin_1 = c(); y1_emp_lower_conf_margin_1 = c(); y1_emp_upper_conf_margin_1 = c()

#  KM estimates
for(i in x1_emp_margin_1){
  temp = summary(km_fit_margin_1, times = i)
  y1_emp_margin_1 = c(y1_emp_margin_1, temp$surv)
  y1_emp_lower_conf_margin_1 = c(y1_emp_lower_conf_margin_1, temp$lower)
  y1_emp_upper_conf_margin_1 = c(y1_emp_upper_conf_margin_1, temp$upper)
}

#  theoretical probabilities
y1_theo_margin_1 = sapply(x1_emp_margin_1, function(i) 1 - marginal_1_cdf(i))
y1_theo_margin_1[1] = 1

#  data for plotting
margin_1_dataframe = data.frame(
  x1_emp_margin_1 = x1_emp_margin_1,
  y1_emp_margin_1 = y1_emp_margin_1,
  y1_emp_lower_conf_margin_1 = y1_emp_lower_conf_margin_1,
  y1_emp_upper_conf_margin_1 = y1_emp_upper_conf_margin_1,
  y1_theo_margin_1 = y1_theo_margin_1
)

# Plot empirical vs theoretical CDF
fig_margin_1_gg <- ggplot(margin_1_dataframe, aes(x = x1_emp_margin_1)) +
  geom_line(aes(y = y1_emp_margin_1, color = "Empirical Probability"), size = 1) +
  geom_ribbon(aes(ymin = y1_emp_lower_conf_margin_1, ymax = y1_emp_upper_conf_margin_1, fill = "Empirical Confidence Interval"), alpha = 0.5) +
  geom_point(aes(y = y1_emp_margin_1), color = "blue", size = 2.2, show.legend = FALSE) +
  geom_line(aes(y = y1_theo_margin_1, color = "Theoretical Probability"), size = 1) +
  labs(x = "Time", y = "Cumulative Distribution Function") +
  scale_color_manual(name = NULL, values = c("Empirical Probability" = "blue", "Theoretical Probability" = "red")) +
  scale_fill_manual(name = NULL, values = c("Empirical Confidence Interval" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    panel.grid = element_blank(),
    legend.position = c(0.2, 0.2),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  )

fig_margin_1_gg


#  marginal CDF for cause 2
marginal_2_cdf <- function(t) {
  P = param[1]
  C = c(param[2], param[3], param[4])
  lamb = sum(C[1] + C[3]) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t))^lamb2
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    return(1 - t1*t2)
  } else {
    t1 = (1 - F_not(t))^lamb3
    t2 = (1 - F_not(tau[2]))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3]))^(lamb2 - lamb3)
    return(1 - t1*t2*t3)
  }
}

#  empirical data
x1_emp_margin_2 = c(subset(dataset$time, dataset$cause %in% c(common_cause, cause_2)), rep(1.54, 9))
state_emp_margin_2 = c(rep(1, length(subset(dataset$time, dataset$cause %in% c(common_cause, cause_2)))), rep(0, 9))
km_fit_margin_2 <- survfit(Surv(x1_emp_margin_2, state_emp_margin_2) ~ 1)

x1_emp_margin_2 = c(0, x1_emp_margin_2)
y1_emp_margin_2 = c(); y1_emp_lower_conf_margin_2 = c(); y1_emp_upper_conf_margin_2 = c()

#  KM estimates
for(i in x1_emp_margin_2){
  temp = summary(km_fit_margin_2, times = i)
  y1_emp_margin_2 = c(y1_emp_margin_2, temp$surv)
  y1_emp_lower_conf_margin_2 = c(y1_emp_lower_conf_margin_2, temp$lower)
  y1_emp_upper_conf_margin_2 = c(y1_emp_upper_conf_margin_2, temp$upper)
}

#  theoretical probabilities
y1_theo_margin_2 = sapply(x1_emp_margin_2, function(i) 1 - marginal_2_cdf(i))
y1_theo_margin_2[1] = 1

#  data for plotting
margin_2_dataframe = data.frame(
  x1_emp_margin_2 = x1_emp_margin_2,
  y1_emp_margin_2 = y1_emp_margin_2,
  y1_emp_lower_conf_margin_2 = y1_emp_lower_conf_margin_2,
  y1_emp_upper_conf_margin_2 = y1_emp_upper_conf_margin_2,
  y1_theo_margin_2 = y1_theo_margin_2
)

# Plot empirical vs theoretical CDF
fig_margin_2_gg <- ggplot(margin_2_dataframe, aes(x = x1_emp_margin_2)) +
  geom_line(aes(y = y1_emp_margin_2, color = "Empirical Probability"), size = 1) +
  geom_ribbon(aes(ymin = y1_emp_lower_conf_margin_2, ymax = y1_emp_upper_conf_margin_2, fill = "Empirical Confidence Interval"), alpha = 0.5) +
  geom_point(aes(y = y1_emp_margin_2), color = "blue", size = 2.2, show.legend = FALSE) +
  geom_line(aes(y = y1_theo_margin_2, color = "Theoretical Probability"), size = 1) +
  labs(x = "Time", y = "Cumulative Distribution Function") +
  scale_color_manual(name = NULL, values = c("Empirical Probability" = "blue", "Theoretical Probability" = "red")) +
  scale_fill_manual(name = NULL, values = c("Empirical Confidence Interval" = "lightblue")) +
  theme_minimal() +
  theme(
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(0, "cm"),
    panel.grid = element_blank(),
    legend.position = c(0.2, 0.2),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12, color = "black")
  )

fig_margin_2_gg


conf_int <- character(length = length(indicator1))  #  confidence interval strings

for (i in seq_along(indicator1)) {
  conf_int[i] <- paste0("(", indicator1[i], ", ", indicator2[i], ")")  #  bounds
}

MLE_data <- data.frame(Parameters= c("P","$c_0$","$c_1$","$c_2$"))  # parameters table
MLE_data$MLE = param  # MLE estimates
MLE_data$Confidence_interval = conf_int  #  confidence intervals

model_llh <- data.frame(llh= c(-result_1$value))  # log-likelihood
model_llh$AIC <- c(2*length(param)+2*result_1$value)  # AIC calculation

model <- "common_"  
distribution <- "exponential_"  

file_name <- paste(model,distribution,"data_analysis.xlsx")  # output file name

write.xlsx(list(MLE_data = MLE_data, model_llh=model_llh), file =  file_name, rowNames = FALSE)  # export Excel

fig_min_path <- paste(model,distribution,"fig_min.pdf")  # figure paths
fig_margin1_path <- paste(model,distribution,"fig_margin_1.pdf")
fig_margin2_path <- paste(model,distribution,"fig_margin_2.pdf")

ggsave(fig_min_path, plot = fig_minimum_gg, dpi = 1200)  # save figures
ggsave(fig_margin1_path, plot = fig_margin_1_gg, dpi = 1200)
ggsave(fig_margin2_path, plot = fig_margin_2_gg, dpi = 1200)






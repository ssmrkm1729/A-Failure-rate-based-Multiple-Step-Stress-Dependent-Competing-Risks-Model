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

# Load dataset

dataset <- read_excel("dataset_ssm.xlsx")
dataset$time <- dataset$time  # Ensure 'time' column is properly set


# Split data by stress stages

stress_data_stage_1 <- subset(dataset, stress == 1)
stress_data_stage_2 <- subset(dataset, stress == 2)
stress_data_stage_3 <- subset(dataset, stress == 3)


#  time lists for each stage

t <- list(
  c(stress_data_stage_1$time),
  c(stress_data_stage_2$time),
  c(stress_data_stage_3$time)
)


# Censoring  for each stage

r_1 <- 0  # Stage 1 right-censored count
r_2 <- 0  # Stage 2 right-censored count
r_3 <- 9  # Stage 3 right-censored count



n <- dim(dataset)[1]  # Total number of observations
m <- 3                # Number of stress stages
tau <- c(0, 1.25, 1.41, 1.54)  # Stage transition times


# Cause indices

cause_1 <- 2
cause_2 <- 3
common_cause <- 1


#  failures by cause and stage

# Stage 1
n_11 <- table(stress_data_stage_1$cause, stress_data_stage_1$stress)[cause_1, ]
n_12 <- table(stress_data_stage_1$cause, stress_data_stage_1$stress)[cause_2, ]
n_10 <- table(stress_data_stage_1$cause, stress_data_stage_1$stress)[common_cause, ]

# Stage 2
n_21 <- table(stress_data_stage_2$cause, stress_data_stage_2$stress)[cause_1, ]
n_22 <- table(stress_data_stage_2$cause, stress_data_stage_2$stress)[cause_2, ]
n_20 <- table(stress_data_stage_2$cause, stress_data_stage_2$stress)[common_cause, ]

# Stage 3
n_31 <- 0  
n_32 <- 0  
n_30 <- table(stress_data_stage_3$cause, stress_data_stage_3$stress)[common_cause, ]


#  failures into a list

failure_count <- list(
  c(n_10, n_11, n_12),
  c(n_20, n_21, n_22),
  c(n_30, n_31, n_32)
)

# Censoring list
r <- list(r_1, r_2, r_3)

# Stress level ratios (scaled by 150)
S <- c(131.5, 142.5, 158) / 150



# Baseline CDF function

F_not <- function(t, alpha) {
  # Custom CDF function (user-defined)
  a = 1 - exp(1 - exp(t^alpha))
  return(a)
}


#  CDF variation with lambda

D_F <- function(t, lamb, alpha) {
  # Adjusted CDF using stress/failure parameter lambda
  a = 1 - exp(lamb * (1 - exp(t^alpha)))
  return(a)
}


# Log-density function

ln_f_not <- function(t, alpha) {
  # Log of probability density function (user-defined)
  a = 1 - exp(t^alpha) + t^alpha + (alpha - 1) * log(t) + log(alpha)
  return(a)
}


# Log-survival function

ln_S_not <- function(t, alpha) {
  return(1 - exp(t^alpha))
}


# Lambda function for each stress stage

lambdas <- function(C, S, P, i, j = NULL) {
  # C: vector of C parameters
  # S: stress levels
  # P: power parameter
  # i: current stage
  # j: optional cause index
  if (!is.null(j)) {
    return(C[j] * (S[i]^P))  # Single cause
  } else {
    return(C[1] * (S[i]^P) + C[2] * (S[i]^P) + C[3] * (S[i]^P))  # Sum over all causes
  }
}


#  C parameter for a specific cause

C_param <- function(b, P, alpha, m, t, failure_count, S, tau, r) {
  
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
      a = (sum(as.integer(failure_count[[i]])) + r[[i]]) * 
        log(1 - F_not(tau[l+1], alpha)) * ((S[l]^P) - (S[l+1]^P))
      t2 = t2 + a
    }
  }
  
  # 
  t3 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t3 = t3 + (S[i]^P) * log(1 - F_not(t[[i]][[j]], alpha))
      }
    }
  }
  
  # 
  t4 = 0
  for (i in 1:m) {
    t4 = t4 + r[[i]] * (S[i]^P) * log(1 - F_not(tau[i+1], alpha))
  }
  
  # 
  return(-t1 / (t2 + t3 + t4))
}



# Profile log-likelihood function

prof_log_likelihood <- function(P, alpha, m, t, failure_count, S, tau, r) {
  
  # 
  c_0 = C_param(1, P, alpha, m, t, failure_count, S, tau, r)
  c_1 = C_param(2, P, alpha, m, t, failure_count, S, tau, r)
  c_2 = C_param(3, P, alpha, m, t, failure_count, S, tau, r)
  
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
      a = C * (sum(as.integer(failure_count[[i]])) + r[[i]]) * 
        log(1 - F_not(tau[l + 1], alpha)) * ((S[l]^P) - (S[l + 1]^P))
      t3 = t3 + a
    }
  }
  
  # 
  t4 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t4 = t4 + (C * (S[i]^P) - 1) * log(1 - F_not(t[[i]][[j]], alpha))
      }
    }
  }
  
  # 
  t5 = 0
  for (i in 1:m) {
    t5 = t5 + C * r[[i]] * (S[i]^P) * log(1 - F_not(tau[i + 1], alpha))
  }
  
  # 
  t6 = 0
  for (i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if (k > 0) {
      for (j in 1:k) {
        t6 = t6 + ln_f_not(t[[i]][[j]], alpha)
      }
    }
  }
  
  # 
  return(t1 + t2 + t3 + t4 + t5 + t6)
}


#  search over P and alpha

P_values <- seq(from = 0.001, to = 50, by = 0.5)      #  P values
alpha_values <- seq(from = 0.001, to = 3, by = 0.01)  #  alpha values

#  matrix to store log-likelihood values
z <- matrix(0, nrow = length(P_values), ncol = length(alpha_values))

#  log-likelihood for each (P, alpha) combination
for (i in 1:length(P_values)) {
  for (j in 1:length(alpha_values)) {
    pl <- prof_log_likelihood(P_values[i], alpha_values[j], m, t, failure_count, S, tau, r)
    if (is.finite(pl)) {
      z[i, j] <- pl
    } else {
      z[i, j] <- -Inf   
    }
  }
}


#  optimal P and alpha

max_value <- max(z)                           # Maximum log-likelihood
max_position <- which(z == max_value, arr.ind = TRUE)  # 

opt_p = P_values[max_position[1]]            # Optimal P
opt_alpha = alpha_values[max_position[2]]    # Optimal alpha

# Compute optimal C and lambda values
opt_C = C_param(1, opt_p, opt_alpha, m, t, failure_count, S, tau, r) +
  C_param(2, opt_p, opt_alpha, m, t, failure_count, S, tau, r) +
  C_param(3, opt_p, opt_alpha, m, t, failure_count, S, tau, r)

opt_lamb = opt_C * (S^opt_p)


# 2D contour plot of log-likelihood
fig_2d <- plot_ly()
fig_2d <- add_trace(
  fig_2d,
  x = alpha_values,  # x-axis: alpha values
  y = P_values,      # y-axis: P values
  z = z,             # z-axis: log-likelihood values
  type = "contour",  # contour plot
  contours = list(
    size = 10, 
    start = max(z), 
    end = min(z), 
    coloring = "fill", 
    showlabels = TRUE, 
    labelfont = list(size = 30) 
  ),
  colorbar = list(
    len = 1,                               
    tickfont = list(size = 35, color = "black"), 
    title = list(font = list(size = 40, color = "black"))
  )
) %>%
  layout(
    xaxis = list(title = list(text = "\u03B1", font = list(size = 48, color = "black")),
                 tickfont = list(size = 40, color = "black")),
    yaxis = list(title = list(text = "P", font = list(size = 48, color = "black")),
                 tickfont = list(size = 40, color = "black")),
    margin = list(b = 80)
  )

fig_2d  # display 2D contour plot

#  contour plot as PDF
orca(fig_2d, file = "contour_frm_chen.pdf", format = "pdf", width = 1600, height = 1200, scale = 3)

# 3D surface plot of log-likelihood
fig_3d <- plot_ly()
fig_3d <- add_trace(fig_3d, x = alpha_values, y = P_values, z = z, type = "surface")
fig_3d <- fig_3d %>% layout(scene = list(xaxis = list(title = "alpha"),
                                         yaxis = list(title = "P"),
                                         zaxis = list(title = "Log Likelihood")))
fig_3d  # display 3D surface plot

# Initial parameters for optimization
initial_params <- c(opt_p, opt_alpha)

# Objective function for maximization (negative log-likelihood)
objective_func <- function(params) {
  P <- params[1]  
  alpha <- params[2]
  return(-prof_log_likelihood(P, alpha, m, t, failure_count, S , tau, r))
}

# Bounds for parameters
lower_bound <- c(0.001, 0.005)
upper_bound <- c(20, 20)

# Run  optimization
result_1 <- optim(par = initial_params, fn = objective_func,
                  method = "L-BFGS-B", lower = lower_bound, upper = upper_bound, hessian = T)

#  estimated parameters
estimated_P <- result_1$par[1]
estimated_alpha <- result_1$par[2]

#  estimated C parameters
estimated_c_0 <- C_param(1, estimated_P, estimated_alpha, m, t, failure_count, S , tau, r)
estimated_c_1 <- C_param(2, estimated_P, estimated_alpha, m, t, failure_count, S , tau, r)
estimated_c_2 <- C_param(3, estimated_P, estimated_alpha, m, t, failure_count, S , tau, r)


# Function for the full log-likelihood
log_likelihood <- function(P, alpha, c_0, c_1, c_2, m, t, failure_count, S , tau, r) {
  
  C = c_0 + c_1 + c_2             
  C_vec = c(c_0, c_1, c_2)        
  
  # 
  t1 = 0
  for(j in 1:3) {
    for(i in 1:m) {
      t1 = t1 + log(C_vec[j]) * failure_count[[i]][j]
    }
  }
  
  # 
  t2 = 0
  for(i in 1:m) {
    t2 = t2 + P * log(S[i]) * sum(as.integer(failure_count[[i]]))
  }
  
  # 
  t3 = 0
  for(i in 2:m) {
    k = i - 1
    for(l in 1:k) {
      a = C * (sum(as.integer(failure_count[[i]])) + r[[i]]) * 
        log(1 - F_not(tau[l+1], alpha)) * ((S[l]^P) - (S[l+1]^P))
      t3 = t3 + a
    }
  }
  
  # 
  t4 = 0
  for(i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if(k > 0) {
      for(j in 1:k) {
        t4 = t4 + (C*(S[i]^P) - 1) * log(1 - F_not(t[[i]][[j]], alpha))
      }
    }
  }
  
  # 
  t5 = 0
  for(i in 1:m) {
    t5 = t5 + C * r[[i]] * (S[i]^P) * log(1 - F_not(tau[i+1], alpha))
  }
  
  # 
  t6 = 0
  for(i in 1:m) {
    k = sum(as.integer(failure_count[[i]]))
    if(k > 0) {
      for(j in 1:k) {
        t6 = t6 + ln_f_not(t[[i]][[j]], alpha)
      }
    }
  }
  
  return(t1 + t2 + t3 + t4 + t5 + t6)
}

# Objective function for  negative log-likelihood
objective_func_for_hessian <- function(params) {
  P <- params[1] 
  alpha <- params[2]
  c_0 <- params[3]
  c_1 <- params[4]
  c_2 <- params[5]
  
  return(-log_likelihood(P, alpha, c_0, c_1, c_2, m, t, failure_count, S , tau, r))
}

#  estimated parameters 
param = c(estimated_P, estimated_alpha, estimated_c_0, estimated_c_1, estimated_c_2)

#  variance 
variance = diag(SMPracticals::ihess(objective_func_for_hessian, param, ep=1e-05))

#  95% confidence intervals for parameters
indicator1 = param - qnorm(p = .05/2, lower.tail = FALSE) * sqrt(variance)
indicator2 = param + qnorm(p = .05/2, lower.tail = FALSE) * sqrt(variance)


# Function for minimum CDF 
minimum_cdf <- function(t) {
  P = param[1]
  alpha = param[2]
  C = c(param[3], param[4], param[5])
  lamb = sum(C) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  # 
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t, alpha))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t, alpha))^lamb2
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    return(1 - t1 * t2)
  } else {
    t1 = (1 - F_not(t, alpha))^lamb3
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3], alpha))^(lamb2 - lamb3)
    return(1 - t1 * t2 * t3)
  }
}

#  empirical survival data for minimum system
x1_emp_min = c(dataset$time, rep(1.54, 9))   # observed times + censoring
state_emp_min = c(rep(1, length(dataset$time)), rep(0, 9)) # 1=failure, 0=censored
km_fit_min <- survfit(Surv(x1_emp_min, state_emp_min) ~ 1)

# 
x1_emp_min = c(0, dataset$time, rep(1.54, 9))

#  KM estimates and confidence intervals
y1_emp_min = y1_emp_lower_conf_min = y1_emp_upper_conf_min = c()
for(i in x1_emp_min){
  temp = summary(km_fit_min, times = i)
  y1_emp_min = c(y1_emp_min, temp$surv)
  y1_emp_lower_conf_min = c(y1_emp_lower_conf_min, temp$lower)
  y1_emp_upper_conf_min = c(y1_emp_upper_conf_min, temp$upper)
}

#  theoretical CDF for minimum 
y1_theo_min = c()
for(i in x1_emp_min){
  y1_theo_min = c(y1_theo_min, 1 - minimum_cdf(i))
}
y1_theo_min[1] = 1  # survival at t=0

#  data frame for plotting
min_dataframe = data.frame(
  x1_emp_min = x1_emp_min, 
  y1_emp_min = y1_emp_min, 
  y1_emp_lower_conf_min = y1_emp_lower_conf_min, 
  y1_emp_upper_conf_min = y1_emp_upper_conf_min, 
  y1_theo_min = y1_theo_min
)

# Plot empirical and theoretical CDFs
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
    axis.text.x  = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, face = "bold", color = "black"),
    legend.text  = element_text(size = 12,  color = "black"),
    legend.title = element_text(size = 12,  color = "black")
  )  

fig_minimum_gg


# Function for marginal CDF of cause 1
marginal_1_cdf <- function(t) {
  P = param[1]
  alpha = param[2]
  C = c(param[3], param[4], param[5])
  lamb = sum(C[1] + C[2]) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  # 
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t, alpha))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t, alpha))^lamb2
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    return(1 - t1 * t2)
  } else {
    t1 = (1 - F_not(t, alpha))^lamb3
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3], alpha))^(lamb2 - lamb3)
    return(1 - t1 * t2 * t3)
  }
}

#  empirical data for marginal CDF (cause 1)
x1_emp_margin_1 = c(subset(dataset$time, dataset$cause %in% c(common_cause, cause_1)), rep(1.54, 9))
state_emp_margin_1 = c(rep(1, length(subset(dataset$time, dataset$cause %in% c(common_cause, cause_1)))), rep(0, 9))
km_fit_margin_1 <- survfit(Surv(x1_emp_margin_1, state_emp_margin_1) ~ 1)

# 
x1_emp_margin_1 = c(0, x1_emp_margin_1)

#  KM estimates and confidence intervals
y1_emp_margin_1 = y1_emp_lower_conf_margin_1 = y1_emp_upper_conf_margin_1 = c()
for(i in x1_emp_margin_1){
  temp = summary(km_fit_margin_1, times = i)
  y1_emp_margin_1 = c(y1_emp_margin_1, temp$surv)
  y1_emp_lower_conf_margin_1 = c(y1_emp_lower_conf_margin_1, temp$lower)
  y1_emp_upper_conf_margin_1 = c(y1_emp_upper_conf_margin_1, temp$upper)
}

#  theoretical marginal CDF for cause 1
y1_theo_margin_1 = c()
for(i in x1_emp_margin_1){
  y1_theo_margin_1 = c(y1_theo_margin_1, 1 - marginal_1_cdf(i))
}
y1_theo_margin_1[1] = 1

#  data frame for plotting
margin_1_dataframe = data.frame(
  x1_emp_margin_1 = x1_emp_margin_1, 
  y1_emp_margin_1 = y1_emp_margin_1, 
  y1_emp_lower_conf_margin_1 = y1_emp_lower_conf_margin_1, 
  y1_emp_upper_conf_margin_1 = y1_emp_upper_conf_margin_1, 
  y1_theo_margin_1 = y1_theo_margin_1
)

# Plot empirical and theoretical marginal CDF
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
    axis.text.x  = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, face = "bold", color = "black"),
    legend.text  = element_text(size = 12,  color = "black"),
    legend.title = element_text(size = 12,  color = "black")
  )  

fig_margin_1_gg


# Function for marginal CDF of cause 2
marginal_2_cdf <- function(t) {
  P = param[1]
  alpha = param[2]
  C = c(param[3], param[4], param[5])
  lamb = sum(C[1] + C[3]) * S^P
  lamb1 = lamb[1]; lamb2 = lamb[2]; lamb3 = lamb[3]
  
  # 
  if(t > tau[1] && t <= tau[2]) {
    return(1 - (1 - F_not(t, alpha))^lamb1)
  } else if(t > tau[2] && t <= tau[3]) {
    t1 = (1 - F_not(t, alpha))^lamb2
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    return(1 - t1 * t2)
  } else {
    t1 = (1 - F_not(t, alpha))^lamb3
    t2 = (1 - F_not(tau[2], alpha))^(lamb1 - lamb2)
    t3 = (1 - F_not(tau[3], alpha))^(lamb2 - lamb3)
    return(1 - t1 * t2 * t3)
  }
}


#  empirical data for marginal CDF (cause 2)
x1_emp_margin_2 = c(subset(dataset$time, dataset$cause %in% c(common_cause, cause_2)), rep(1.54, 9))
state_emp_margin_2 = c(rep(1, length(subset(dataset$time, dataset$cause %in% c(common_cause, cause_2)))), rep(0, 9))
km_fit_margin_2 <- survfit(Surv(x1_emp_margin_2, state_emp_margin_2) ~ 1)

# 
x1_emp_margin_2 = c(0, x1_emp_margin_2)

#  KM estimates and confidence intervals
y1_emp_margin_2 = y1_emp_lower_conf_margin_2 = y1_emp_upper_conf_margin_2 = c()
for(i in x1_emp_margin_2){
  temp = summary(km_fit_margin_2, times = i)
  y1_emp_margin_2 = c(y1_emp_margin_2, temp$surv)
  y1_emp_lower_conf_margin_2 = c(y1_emp_lower_conf_margin_2, temp$lower)
  y1_emp_upper_conf_margin_2 = c(y1_emp_upper_conf_margin_2, temp$upper)
}

#  theoretical marginal CDF for cause 2
y1_theo_margin_2 = c()
for(i in x1_emp_margin_2){
  y1_theo_margin_2 = c(y1_theo_margin_2, 1 - marginal_2_cdf(i))
}
y1_theo_margin_2[1] = 1

# 
margin_2_dataframe = data.frame(
  x1_emp_margin_2 = x1_emp_margin_2, 
  y1_emp_margin_2 = y1_emp_margin_2, 
  y1_emp_lower_conf_margin_2 = y1_emp_lower_conf_margin_2, 
  y1_emp_upper_conf_margin_2 = y1_emp_upper_conf_margin_2, 
  y1_theo_margin_2 = y1_theo_margin_2
)

# Plot empirical and theoretical marginal CDF for cause 2
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
    axis.text.x  = element_text(size = 10, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 10, face = "bold", color = "black"),
    legend.text  = element_text(size = 12,  color = "black"),
    legend.title = element_text(size = 12,  color = "black")
  )  

fig_margin_2_gg




#  confidence intervals 
conf_int <- character(length = length(indicator1))
for (i in seq_along(indicator1)) {
  conf_int[i] <- paste0("(", indicator1[i], ", ", indicator2[i], ")")
}

#  data frame for MLE estimates and confidence intervals
MLE_data <- data.frame(Parameters = c("P","alpha","$c_0$","$c_1$","$c_2$"))
MLE_data$MLE <- param
MLE_data$Confidence_interval <- conf_int

#  log-likelihood and AIC for the model
model_llh <- data.frame(llh = c(-result_1$value))
model_llh$AIC <- c(2 * length(param) + 2 * result_1$value)

#  model and distribution names 
model <- "frm_"
distribution <- "chen_"


#   file name and write results
file_name <- paste(model, distribution, "data_analysis.xlsx")
write.xlsx(list(MLE_data = MLE_data, model_llh = model_llh), 
           file = file_name, rowNames = FALSE)

#  paths for plots
fig_min_path <- paste(model, distribution, "fig_min.pdf")
fig_margin1_path <- paste(model, distribution, "fig_margin_1.pdf")
fig_margin2_path <- paste(model, distribution, "fig_margin_2.pdf")

# Save plots 
ggsave(fig_min_path, plot = fig_minimum_gg, dpi = 1200)
ggsave(fig_margin1_path, plot = fig_margin_1_gg, dpi = 1200)
ggsave(fig_margin2_path, plot = fig_margin_2_gg, dpi = 1200)



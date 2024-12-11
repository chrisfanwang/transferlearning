library(genlasso)
library(plotrix)
library(ggplot2)
source("transfer_learning_functions.R")
library(tictoc)
library(latex2exp)
library(psych)
library(changepoints)
library(egg)



# Calculate harmonic means for different subsets of source observations
harmonic_sizes = c(
  1*harmonic.mean(c(2000)),
  2*harmonic.mean(c(2000, 1800)),
  3*harmonic.mean(c(2000, 1800, 1600)),
  4*harmonic.mean(c(2000, 1800, 1600, 1400)),
  5*harmonic.mean(c(2000, 1800, 1600, 1400, 1200)),
  6*harmonic.mean(c(2000, 1800, 1600, 1400, 1200, 1000)),
  7*harmonic.mean(c(2000, 1800, 1600, 1400, 1200, 1000, 800)),
  8*harmonic.mean(c(2000, 1800, 1600, 1400, 1200, 1000, 800, 600)),
  9*harmonic.mean(c(2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400)),
  10*harmonic.mean(c(2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200)))

# Create a data frame for plotting harmonic means
sizes  = seq(1, 10, by = 1) 
d_1 <- data.frame(sizes = sizes, harmonic_sizes = harmonic_sizes)

# Plot harmonic means
plott_1 <- ggplot(data = d_1, mapping = aes(x = sizes, y = harmonic_sizes)) +
  geom_line() + 
  scale_x_continuous(breaks = seq(min(d_1$sizes), max(d_1$sizes), by = 1)) +
  labs(y = expression(atop(K %*% "the harmonic mean", "of source observations")), x = "K") +
  theme_classic() +
  theme(
    axis.line.x = element_line(linewidth = 0.5, color = "black"),
    axis.line.y = element_line(linewidth = 1, color = "black")
  )
print(plott_1)



# Simulation parameters
gamma= 0.5
jump.mean = c(2*gamma, 4*gamma, 1*gamma, 5*gamma, 7*gamma, 8*gamma, 2*gamma, 1*gamma, 3*gamma, 4*gamma)
jump.location = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

n_0 = 200
n_k_size = c(2000, 1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200)

H_k = 0.15
sig.delta_1 = 0.2
h_1 = 20
sigma = 0.5
simu = 100

# Initialize matrices to store results
result_fl =  matrix(NA, nrow = length(n_k_size), ncol = simu)
result_l0 = matrix(NA, nrow = length(n_k_size), ncol = simu)
result_fl_est =  matrix(NA, nrow = length(n_k_size), ncol = simu)
result_l0_est = matrix(NA, nrow = length(n_k_size), ncol = simu)


# Perform simulations
for (simu_test in 1:simu){
  set.seed(10*simu_test)
  
  # Generate coefficients and observations
  coef_f  = f_gen_new(jump.mean, jump.location, n_0, n_k_size, H_k, sig.delta_1, h_1,  exact = T)
  f_0 = coef_f$f_0
  W = coef_f$W
  obs = obs_gen_new(f_0, W, sigma)
  y_0 = obs$y_0
  
  for (k in 1:length(n_k_size)){
    k_set = n_k_size[1:k]
    k_subsets = list()
    for (i in 1:length(k_set)) {
      k_combos =  combn(k_set, i, simplify = FALSE)
      k_subsets <- c(k_subsets, k_combos)
    }
    k_harmonic <- lapply(k_subsets, function(x) length(x)^2 / sum(1/x))
    
    k_est_A = k_subsets[[which.max(unlist(k_harmonic))]]
    
    y_harmonic_project = rep(0, n_0)
    
    for (i in k_est_A){
      index = which(n_k_size == i)
      y_project = P_n0_n(n_0, i) %*% obs$obs_y[[index]]
      
      y_harmonic_project = y_harmonic_project + y_project
    }
  
    y_harmonic_project = y_harmonic_project/length(k_est_A)

    # Compute projections using all sources
    y_all_project = rep(0, n_0)
    for (i in 1:k){
      y_project = P_n0_n(n_0, n_k_size[i]) %*% obs$obs_y[[i]]
      
      y_all_project = y_all_project + y_project
    }
    y_all_project = y_all_project/k

    # all multisource l_1-penalized method  
    result_fl_test = fusedlasso1d(y_all_project)
    cv_fl_test = cv.trendfilter(result_fl_test, k=5)
    fit_fl_test = coef(result_fl_test, lambda = cv_fl_test$lambda.min)$beta
    result_fl[k, simu_test] = sum((fit_fl_test - f_0)^2)/n_0
    
     # all multisource l_0-penalized method 
    fit_l0_test = cv_l_0_est(y_all_project) 
    result_l0[k, simu_test] = sum((fit_l0_test-f_0)^2)/n_0

    # selected multisource l_1-penalized method                      
    result_fl_est_test = fusedlasso1d(y_harmonic_project)
    cv_fl_est_test = cv.trendfilter(result_fl_est_test, k=5)    
    fit_fl_est_test = coef(result_fl_est_test, lambda = cv_fl_est_test$lambda.min)$beta
    result_fl_est[k, simu_test] = sum((fit_fl_est_test - f_0)^2)/n_0
    
    # selected multisource l_0-penalized method 
    fit_l0_est_test = cv_l_0_est(y_harmonic_project) 
    result_l0_est[k, simu_test] = sum((fit_l0_est_test-f_0)^2)/n_0
    
    
    
  }
}



sizes  = rep(seq(1, 10, by = 1), 4) 

value = c(apply(result_fl, 1, mean),
          apply(result_fl_est, 1, mean),
          apply(result_l0, 1, mean), 
          apply(result_l0_est, 1, mean) )

Method <- rep(c('fl', 'fl_est',  'l_0', 'l_0_est'), each = 10)

se = c( apply(result_fl, 1, function(column) sd(column) / sqrt(length(column))),  
        apply(result_fl_est, 1, function(column) sd(column) / sqrt(length(column))),  
        apply(result_l0, 1, function(column) sd(column) / sqrt(length(column))),
        apply(result_l0_est, 1, function(column) sd(column) / sqrt(length(column)))
)

d_2 <- data.frame(sizes = sizes, Method = Method, value = value, se = se)

plott_2 <- ggplot(data = d_2, mapping = aes(x =sizes, y = value, colour = Method)) +
  geom_line() + 
  geom_point() + 
  scale_x_continuous(breaks = seq(min(d_1$sizes), max(d_1$sizes), by = 1)) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.3) +
  labs(y = "Estimation error", x = "K") +
  scale_color_manual(
    values = c("#66c2a5", "#fc8d26", "#484B8A", "#e5c494"),   
    labels = c( TeX("$\u2113_1$-T-[K]"), TeX("$\u2113_1$-T-$\\widetilde{A}$"),  TeX("$\u2113_0$-T-[K]"),  TeX("$\u2113_0$-T-$\\widetilde{A}$"))
  ) +
  theme_classic()+  
  theme(
    axis.line.x = element_line(linewidth  = 0.5, color = "black"),
    axis.line.y = element_line(linewidth  = 1, color = "black")
  )

print(plott_2)


ggarrange(plott_1, plott_2, ncol = 2)

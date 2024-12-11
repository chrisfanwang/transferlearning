library(genlasso)
library(plotrix)
library(ggplot2)
source("transfer_learning_functions.R")
library(tictoc)
library(latex2exp)
library(changepoints)
library(egg)

tic()

#Plot 1: configuration 2 and temporal dependence in target data 

n_0 = 200
n_k = 400
K = 10


temp_1_list  = c(0.1, 0.2, 0.3, 0.4)
sigma = 0.5
h_1 = 20
h_2 = 500
t_k = 50

sig.delta_1 = c(0.2)
sig.delta_2  = 2
H_k_list = c(0.075)
A_list = c(4)
gamma_list = c(0.5)

jump.location = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)
simu = 100

result_fl =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_1 =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_ora =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_est =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0 =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_1 =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_ora =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_est =  array(NA, dim = c(length(temp_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))


for (temp_1_index in 1:length(temp_1_list)){
  temp_1 = temp_1_list[temp_1_index]
  for (H_k_index in 1:length(H_k_list)){
    H_k_ratio = H_k_list[H_k_index]
    H_k = H_k_ratio*n_k
    for (A_index in 1:length(A_list)){
      size_A = A_list[A_index]
      for (gamma_index in 1:length(gamma_list)) {
        gamma = gamma_list[gamma_index]
        jump.mean = c(2*gamma, 4*gamma, 1*gamma, 5*gamma, 7*gamma, 8*gamma, 2*gamma, 1*gamma)
        for (simu_test in 1:simu){
        set.seed(10*temp_1_index*H_k_index*A_index*gamma_index*simu_test)
          
        coef_f  = f_gen(jump.mean, jump.location, n_0, n_k, K, size_A, H_k, sig.delta_1, sig.delta_2, h_1, h_2, exact = F)
        f_0 = coef_f$f_0
        W = coef_f$W
        obs = obs_gen_tempde(f_0, W, sigma, temp_1)
        y_0 = obs$y_0
        y_1 = obs$obs_y[1, ]
 
        
        
        result_fl_test = fusedlasso1d(y_0)
        
        cv_fl_test = cv.trendfilter(result_fl_test, k=5)
        
        fit_fl_test = coef(result_fl_test, lambda = cv_fl_test$lambda.min)$beta
        
        result_fl[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_test-f_0)^2)/n_0
        
        
        y_1_project = P_n0_n(n_0, n_k) %*% y_1
        
        result_fl_T_1_test = fusedlasso1d(y_1_project)
        
        cv_fl_T_1_test = cv.trendfilter(result_fl_T_1_test, k=5)
        
        fit_fl_T_1_test = coef(result_fl_T_1_test, lambda = cv_fl_T_1_test$lambda.min)$beta
        
        result_fl_T_1[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_1_test -f_0)^2)/n_0
        
        
        y_A_project = rep(0, n_0)
  
        for (k in 1:size_A){
        y_project = P_n0_n(n_0, n_k) %*% obs$obs_y[k, ]
        y_A_project = y_A_project + y_project
        }
        y_A_project = y_A_project/size_A
        
        result_fl_T_ora_test = fusedlasso1d(y_A_project)
        
        cv_fl_T_ora_test = cv.trendfilter(result_fl_T_ora_test, k=5)
        
        fit_fl_T_ora_test = coef(result_fl_T_ora_test, lambda = cv_fl_T_ora_test$lambda.min)$beta
        
        result_fl_T_ora[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_ora_test -f_0)^2)/n_0
        
        
        
        fit_l0_test = cv_l_0_est(y_0) 
        
        result_l0[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_test-f_0)^2)/n_0
        
        
        fit_l0_T_1_test =  cv_l_0_est(y_1_project)
        
        result_l0_T_1[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_1_test-f_0)^2)/n_0
        
        
        fit_l0_T_ora_test =  cv_l_0_est(y_A_project)
        
        result_l0_T_ora[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_ora_test-f_0)^2)/n_0
        
        
        
        index_k = A_detect_index(obs$obs_y, y_0)
        
        est_tau = A_detect_threshold(obs$obs_y[index_k,], y_0, t_hat_k = t_k)
        
        print(index_k)
        
        est_A = A_detect_alg(obs$obs_y, y_0, t_hat_k = t_k, est_tau)
        
        print(est_A)
        
        y_A_est_project = rep(0, n_0)
        
        for (k in 1:length(est_A)){
          y_project = P_n0_n(n_0, n_k) %*% obs$obs_y[est_A[k], ]
          y_A_est_project = y_A_est_project + y_project
        }
        y_A_est_project = y_A_est_project/length(est_A)
        
        
        result_fl_T_est_test = fusedlasso1d(y_A_est_project)
        
        cv_fl_T_est_test = cv.trendfilter(result_fl_T_est_test, k=5)
        
        fit_fl_T_est_test = coef(result_fl_T_est_test, lambda = cv_fl_T_est_test$lambda.min)$beta
        
        result_fl_T_est[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_est_test -f_0)^2)/n_0
        
        
        fit_l0_T_est_test =  cv_l_0_est(y_A_est_project)
        
        result_l0_T_est[temp_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_est_test-f_0)^2)/n_0
        
        
        
        }
      }
    }
  }
}

toc()


value_1 = c(mean(result_fl[1,1,1,1,]), mean(result_fl_T_1[1,1,1,1,]), mean(result_fl_T_ora[1,1,1,1,]), mean(result_fl_T_est[1,1,1,1,]), 
    mean(result_l0[1,1,1,1,]), mean(result_l0_T_1[1,1,1,1,]), mean(result_l0_T_ora[1,1,1,1,]), mean(result_l0_T_est[1,1,1,1,]),
    mean(result_fl[2,1,1,1,]), mean(result_fl_T_1[2,1,1,1,]), mean(result_fl_T_ora[2,1,1,1,]), mean(result_fl_T_est[2,1,1,1,]), 
    mean(result_l0[2,1,1,1,]), mean(result_l0_T_1[2,1,1,1,]), mean(result_l0_T_ora[2,1,1,1,]), mean(result_l0_T_est[2,1,1,1,]),
    mean(result_fl[3,1,1,1,]), mean(result_fl_T_1[3,1,1,1,]), mean(result_fl_T_ora[3,1,1,1,]), mean(result_fl_T_est[3,1,1,1,]), 
    mean(result_l0[3,1,1,1,]), mean(result_l0_T_1[3,1,1,1,]), mean(result_l0_T_ora[3,1,1,1,]), mean(result_l0_T_est[3,1,1,1,]),
    mean(result_fl[4,1,1,1,]), mean(result_fl_T_1[4,1,1,1,]), mean(result_fl_T_ora[4,1,1,1,]), mean(result_fl_T_est[4,1,1,1,]),
    mean(result_l0[4,1,1,1,]), mean(result_l0_T_1[4,1,1,1,]), mean(result_l0_T_ora[4,1,1,1,]), mean(result_l0_T_est[4,1,1,1,]))


se_1 = c(std.error(result_fl[1,1,1,1,]), std.error(result_fl_T_1[1,1,1,1,]), std.error(result_fl_T_ora[1,1,1,1,]), std.error(result_fl_T_est[1,1,1,1,]), 
  std.error(result_l0[1,1,1,1,]), std.error(result_l0_T_1[1,1,1,1,]), std.error(result_l0_T_ora[1,1,1,1,]), std.error(result_l0_T_est[1,1,1,1,]),
  std.error(result_fl[2,1,1,1,]), std.error(result_fl_T_1[2,1,1,1,]), std.error(result_fl_T_ora[2,1,1,1,]), std.error(result_fl_T_est[2,1,1,1,]), 
  std.error(result_l0[2,1,1,1,]), std.error(result_l0_T_1[2,1,1,1,]), std.error(result_l0_T_ora[2,1,1,1,]), std.error(result_l0_T_est[2,1,1,1,]),
  std.error(result_fl[3,1,1,1,]), std.error(result_fl_T_1[3,1,1,1,]), std.error(result_fl_T_ora[3,1,1,1,]), std.error(result_fl_T_est[3,1,1,1,]), 
  std.error(result_l0[3,1,1,1,]), std.error(result_l0_T_1[3,1,1,1,]), std.error(result_l0_T_ora[3,1,1,1,]), std.error(result_l0_T_est[3,1,1,1,]),
  std.error(result_fl[4,1,1,1,]), std.error(result_fl_T_1[4,1,1,1,]), std.error(result_fl_T_ora[4,1,1,1,]), std.error(result_fl_T_est[4,1,1,1,]),
  std.error(result_l0[4,1,1,1,]), std.error(result_l0_T_1[4,1,1,1,]), std.error(result_l0_T_ora[4,1,1,1,]), std.error(result_l0_T_est[4,1,1,1,]))





Method <- rep(c('FL', 'FL-T-1', 'FL-T-ora', 'FL-T-est', 'l_0', 'l_0-T-1', 'l_0-T-ora', 'l_0-T-est'), times = 4)
rho <- rep(c(0.1, 0.2, 0.3, 0.4), each = 8)

d_1 <- data.frame(rho = rho, Method = Method, value = value_1, se = se_1)


d_1$Method <- factor(d_1$Method, levels = c('FL', 'FL-T-1', 'FL-T-ora', 'FL-T-est', 'l_0', 'l_0-T-1', 'l_0-T-ora', 'l_0-T-est'))

# Define the colors for each method in the same order as the levels
color_list <- setNames(c("#66c2a5", "#fc8d26", "#8da0cb", "#e78ac3", "#484B8A", "#e5c494", "#ffd92f", "#b3b3bc"), levels(d_1$Method))



plott_1 <- ggplot(data = d_1, mapping = aes(x = rho, y = value, colour = Method)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.01) +
  labs(y = "Estimation error", x = TeX("$\\rho_1$")) +
  scale_color_manual(
    values = color_list,
    labels = c(TeX("$\u2113_1$"), TeX("$\u2113_1$-T-1"),TeX("$\u2113_1$-T-$A$"), TeX("$\u2113_1$-T-$\\widehat{A}$"), 
               TeX("$\u2113_0$"), TeX("$\u2113_0$-T-1"),TeX("$\u2113_0$-T-$A$"), TeX("$\u2113_0$-T-$\\widehat{A}$"))
  )  +
  theme_classic()

print(plott_1)



tic()
#Plot 2: configuration 2 and temporal dependence in difference 

n_0 = 200
n_k = 400
K = 10



sigma = 0.5
h_1 = 20
temp_h_list  = c(0.1, 0.2, 0.3, 0.4)
h_2 = 500
t_k = 50

sig.delta_1 = c(0.2)
sig.delta_2  = 2
H_k_list = c(0.075)
A_list = c(4)
gamma_list = c(0.5)

jump.location = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)
simu = 100

result_fl =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_1 =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_ora =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_est =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0 =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_1 =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_ora =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_est =  array(NA, dim = c(length(temp_h_list), length(H_k_list), length(A_list),  length(gamma_list), simu))


for (temp_h_index in 1:length(temp_h_list)){
  temp_h = temp_h_list[temp_h_index]
  for (H_k_index in 1:length(H_k_list)){
    H_k_ratio = H_k_list[H_k_index]
    H_k = H_k_ratio*n_k
    for (A_index in 1:length(A_list)){
      size_A = A_list[A_index]
      for (gamma_index in 1:length(gamma_list)) {
        gamma = gamma_list[gamma_index]
        jump.mean = c(2*gamma, 4*gamma, 1*gamma, 5*gamma, 7*gamma, 8*gamma, 2*gamma, 1*gamma)
        for (simu_test in 1:simu){
          set.seed(10*temp_h_index*H_k_index*A_index*gamma_index*simu_test)
          coef_f  = f_gen_tempde(jump.mean, jump.location, n_0, n_k, K, size_A, H_k, h_1, h_2, temp_h)
          f_0 = coef_f$f_0
          W = coef_f$W
          obs = obs_gen(f_0, W, sigma)
          y_0 = obs$y_0
          y_1 = obs$obs_y[1, ]
          
          
          
          result_fl_test = fusedlasso1d(y_0)
          
          cv_fl_test = cv.trendfilter(result_fl_test, k=5)
          
          fit_fl_test = coef(result_fl_test, lambda = cv_fl_test$lambda.min)$beta
          
          result_fl[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_test-f_0)^2)/n_0
          
          
          y_1_project = P_n0_n(n_0, n_k) %*% y_1
          
          result_fl_T_1_test = fusedlasso1d(y_1_project)
          
          cv_fl_T_1_test = cv.trendfilter(result_fl_T_1_test, k=5)
          
          fit_fl_T_1_test = coef(result_fl_T_1_test, lambda = cv_fl_T_1_test$lambda.min)$beta
          
          result_fl_T_1[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_1_test -f_0)^2)/n_0
          
          
          y_A_project = rep(0, n_0)
          
          for (k in 1:size_A){
            y_project = P_n0_n(n_0, n_k) %*% obs$obs_y[k, ]
            y_A_project = y_A_project + y_project
          }
          y_A_project = y_A_project/size_A
          
          result_fl_T_ora_test = fusedlasso1d(y_A_project)
          
          cv_fl_T_ora_test = cv.trendfilter(result_fl_T_ora_test, k=5)
          
          fit_fl_T_ora_test = coef(result_fl_T_ora_test, lambda = cv_fl_T_ora_test$lambda.min)$beta
          
          result_fl_T_ora[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_ora_test -f_0)^2)/n_0
          
          
          
          fit_l0_test = cv_l_0_est(y_0) 
          
          result_l0[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_test-f_0)^2)/n_0
          
          
          fit_l0_T_1_test =  cv_l_0_est(y_1_project)
          
          result_l0_T_1[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_1_test-f_0)^2)/n_0
          
          
          fit_l0_T_ora_test =  cv_l_0_est(y_A_project)
          
          result_l0_T_ora[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_ora_test-f_0)^2)/n_0
          
          
          
          index_k = A_detect_index(obs$obs_y, y_0)
          
          est_tau = A_detect_threshold(obs$obs_y[index_k,], y_0, t_hat_k = t_k)
          
          print(index_k)
          
          est_A = A_detect_alg(obs$obs_y, y_0, t_hat_k = t_k, est_tau)
          
          print(est_A)
          
          y_A_est_project = rep(0, n_0)
          
          for (k in 1:length(est_A)){
            y_project = P_n0_n(n_0, n_k) %*% obs$obs_y[est_A[k], ]
            y_A_est_project = y_A_est_project + y_project
          }
          y_A_est_project = y_A_est_project/length(est_A)
          
          
          result_fl_T_est_test = fusedlasso1d(y_A_est_project)
          
          cv_fl_T_est_test = cv.trendfilter(result_fl_T_est_test, k=5)
          
          fit_fl_T_est_test = coef(result_fl_T_est_test, lambda = cv_fl_T_est_test$lambda.min)$beta
          
          result_fl_T_est[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_est_test -f_0)^2)/n_0
          
          
          fit_l0_T_est_test =  cv_l_0_est(y_A_est_project)
          
          result_l0_T_est[temp_h_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_est_test-f_0)^2)/n_0
          
          
          
        }
      }
    }
  }
}

toc()


value_2 = c(mean(result_fl[1,1,1,1,]), mean(result_fl_T_1[1,1,1,1,]), mean(result_fl_T_ora[1,1,1,1,]), mean(result_fl_T_est[1,1,1,1,]), 
            mean(result_l0[1,1,1,1,]), mean(result_l0_T_1[1,1,1,1,]), mean(result_l0_T_ora[1,1,1,1,]), mean(result_l0_T_est[1,1,1,1,]),
            mean(result_fl[2,1,1,1,]), mean(result_fl_T_1[2,1,1,1,]), mean(result_fl_T_ora[2,1,1,1,]), mean(result_fl_T_est[2,1,1,1,]), 
            mean(result_l0[2,1,1,1,]), mean(result_l0_T_1[2,1,1,1,]), mean(result_l0_T_ora[2,1,1,1,]), mean(result_l0_T_est[2,1,1,1,]),
            mean(result_fl[3,1,1,1,]), mean(result_fl_T_1[3,1,1,1,]), mean(result_fl_T_ora[3,1,1,1,]), mean(result_fl_T_est[3,1,1,1,]), 
            mean(result_l0[3,1,1,1,]), mean(result_l0_T_1[3,1,1,1,]), mean(result_l0_T_ora[3,1,1,1,]), mean(result_l0_T_est[3,1,1,1,]),
            mean(result_fl[4,1,1,1,]), mean(result_fl_T_1[4,1,1,1,]), mean(result_fl_T_ora[4,1,1,1,]), mean(result_fl_T_est[4,1,1,1,]),
            mean(result_l0[4,1,1,1,]), mean(result_l0_T_1[4,1,1,1,]), mean(result_l0_T_ora[4,1,1,1,]), mean(result_l0_T_est[4,1,1,1,]))


se_2 = c(std.error(result_fl[1,1,1,1,]), std.error(result_fl_T_1[1,1,1,1,]), std.error(result_fl_T_ora[1,1,1,1,]), std.error(result_fl_T_est[1,1,1,1,]), 
         std.error(result_l0[1,1,1,1,]), std.error(result_l0_T_1[1,1,1,1,]), std.error(result_l0_T_ora[1,1,1,1,]), std.error(result_l0_T_est[1,1,1,1,]),
         std.error(result_fl[2,1,1,1,]), std.error(result_fl_T_1[2,1,1,1,]), std.error(result_fl_T_ora[2,1,1,1,]), std.error(result_fl_T_est[2,1,1,1,]), 
         std.error(result_l0[2,1,1,1,]), std.error(result_l0_T_1[2,1,1,1,]), std.error(result_l0_T_ora[2,1,1,1,]), std.error(result_l0_T_est[2,1,1,1,]),
         std.error(result_fl[3,1,1,1,]), std.error(result_fl_T_1[3,1,1,1,]), std.error(result_fl_T_ora[3,1,1,1,]), std.error(result_fl_T_est[3,1,1,1,]), 
         std.error(result_l0[3,1,1,1,]), std.error(result_l0_T_1[3,1,1,1,]), std.error(result_l0_T_ora[3,1,1,1,]), std.error(result_l0_T_est[3,1,1,1,]),
         std.error(result_fl[4,1,1,1,]), std.error(result_fl_T_1[4,1,1,1,]), std.error(result_fl_T_ora[4,1,1,1,]), std.error(result_fl_T_est[4,1,1,1,]),
         std.error(result_l0[4,1,1,1,]), std.error(result_l0_T_1[4,1,1,1,]), std.error(result_l0_T_ora[4,1,1,1,]), std.error(result_l0_T_est[4,1,1,1,]))





Method <- rep(c('FL', 'FL-T-1', 'FL-T-ora', 'FL-T-est', 'l_0', 'l_0-T-1', 'l_0-T-ora', 'l_0-T-est'), times = 4)
rho <- rep(c(0.1, 0.2, 0.3, 0.4), each = 8)

d_2 <- data.frame(rho = rho, Method = Method, value = value_2, se = se_2)


d_2$Method <- factor(d_2$Method, levels = c('FL', 'FL-T-1', 'FL-T-ora', 'FL-T-est', 'l_0', 'l_0-T-1', 'l_0-T-ora', 'l_0-T-est'))

# Define the colors for each method in the same order as the levels
color_list <- setNames(c("#66c2a5", "#fc8d26", "#8da0cb", "#e78ac3", "#484B8A", "#e5c494", "#ffd92f", "#b3b3bc"), levels(d_2$Method))



plott_2 <- ggplot(data = d_2, mapping = aes(x = rho, y = value, colour = Method)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.01) +
  labs(y = "Estimation error", x = TeX("$\\rho_2$")) +
  scale_color_manual(
    values = color_list,
    labels = c(TeX("$\u2113_1$"), TeX("$\u2113_1$-T-1"),TeX("$\u2113_1$-T-$A$"), TeX("$\u2113_1$-T-$\\widehat{A}$"), 
               TeX("$\u2113_0$"), TeX("$\u2113_0$-T-1"),TeX("$\u2113_0$-T-$A$"), TeX("$\u2113_0$-T-$\\widehat{A}$"))
  )  +
  theme_classic()

print(plott_2)



ggarrange(plott_1,  plott_2, ncol = 2)


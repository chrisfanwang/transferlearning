library(genlasso)
library(plotrix)
library(ggplot2)
source("transfer_learning_functions.R")
library(tictoc)
library(latex2exp)
library(changepoints)
library(egg)


#Plot 1: configuration 1 and vary a  

n_0 = 200
n_k = 400
K = 10


sigma  = 0.5
h_1 = 20
h_2 = 500
t_k = 50

sig.delta_1_list = c(0.2)
sig.delta_2  = 2
H_k_list = c(0.15)
A_list = c(2, 4, 6, 8)
gamma_list = c(0.5)

jump.location = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875)
simu = 100

result_fl =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_1 =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_1_target =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_ora =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_fl_T_ora_target =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))

result_l0 =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_1 =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_1_target =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_ora =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))
result_l0_T_ora_target =  array(NA, dim = c(length(sig.delta_1_list), length(H_k_list), length(A_list),  length(gamma_list), simu))

for (sig.delta_1_index in 1:length(sig.delta_1_list)){
  sig.delta_1 = sig.delta_1_list[sig.delta_1_index]
  for (H_k_index in 1:length(H_k_list)){
    H_k_ratio = H_k_list[H_k_index]
    H_k = H_k_ratio*n_k
    for (A_index in 1:length(A_list)){
      size_A = A_list[A_index]
      for (gamma_index in 1:length(gamma_list)) {
        gamma = gamma_list[gamma_index]
        jump.mean = c(2*gamma, 4*gamma, 1*gamma, 5*gamma, 7*gamma, 8*gamma, 2*gamma, 1*gamma)
        for (simu_test in 1:simu){
          set.seed(10*sig.delta_1_index*H_k_index*A_index*gamma_index*simu_test)
          
          coef_f  = f_gen(jump.mean, jump.location, n_0, n_k, K, size_A, H_k, sig.delta_1, sig.delta_2, h_1, h_2, exact = T)
          f_0 = coef_f$f_0
          W = coef_f$W
          obs = obs_gen(f_0, W, sigma)
          y_0 = obs$y_0
          y_1 = obs$obs_y[1, ]
          
          
          result_fl_test = fusedlasso1d(y_0)
          
          cv_fl_test = cv.trendfilter(result_fl_test, k=5)
          
          fit_fl_test = coef(result_fl_test, lambda = cv_fl_test$lambda.min)$beta
          
          result_fl[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_test-f_0)^2)/n_0
          
          
          y_1_project = P_n0_n(n_0, n_k) %*% y_1
          
          result_fl_T_1_test = fusedlasso1d(y_1_project)
          
          cv_fl_T_1_test = cv.trendfilter(result_fl_T_1_test, k=5)
          
          fit_fl_T_1_test = coef(result_fl_T_1_test, lambda = cv_fl_T_1_test$lambda.min)$beta
          
          result_fl_T_1[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_1_test -f_0)^2)/n_0
          
    
          y_1_target = insert_values (y_0, y_1)
          
          y_1_target_project = P_n0_n(n_0, (n_k+n_0)) %*%  y_1_target
          
          result_fl_T_1_target_test = fusedlasso1d(y_1_target_project)
          
          cv_fl_T_1_target_test = cv.trendfilter(result_fl_T_1_target_test, k=5)
          
          fit_fl_T_1_target_test = coef(result_fl_T_1_target_test, lambda = cv_fl_T_1_target_test$lambda.min)$beta
          
          result_fl_T_1_target[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_1_target_test -f_0)^2)/n_0
          
          
          y_A_project = rep(0, n_0)
          
          for (k in 1:size_A){
            y_project = P_n0_n(n_0, n_k) %*% obs$obs_y[k, ]
            y_A_project = y_A_project + y_project
          }
          y_A_project = y_A_project/size_A
          
          result_fl_T_ora_test = fusedlasso1d(y_A_project)
          
          cv_fl_T_ora_test = cv.trendfilter(result_fl_T_ora_test, k=5)
          
          fit_fl_T_ora_test = coef(result_fl_T_ora_test, lambda = cv_fl_T_ora_test$lambda.min)$beta
          
          result_fl_T_ora[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_ora_test -f_0)^2)/n_0
          
          
          list_y_A <- list()
          
          for (k in 1:size_A){
            list_y_A[[k]] <- obs$obs_y[k, ]
          }
          
          y_ora_target = construct_vector(y_0, list_y_A)
          
          P_alt = alignment_matrix(n_0, size_A, c(n_0, rep(n_k, size_A)))
         
          y_ora_target_project = P_alt %*%  y_ora_target
          
          result_fl_T_ora_target_test = fusedlasso1d(y_ora_target_project)
          
          cv_fl_T_ora_target_test = cv.trendfilter(result_fl_T_ora_target_test, k=5)
          
          fit_fl_T_ora_target_test = coef(result_fl_T_ora_target_test, lambda =  cv_fl_T_ora_target_test$lambda.min)$beta
          
          result_fl_T_ora_target[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_fl_T_ora_target_test -f_0)^2)/n_0
          
          
          
          fit_l0_test =  cv_l_0_est(y_0)
          
          result_l0[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_test-f_0)^2)/n_0
          
          
          fit_l0_T_1_test =  cv_l_0_est(y_1_project) 
          
          result_l0_T_1[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_1_test-f_0)^2)/n_0
          
          
          fit_l0_T_1_target_test =  cv_l_0_est(y_1_target_project) 
          
          result_l0_T_1_target[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_1_target_test-f_0)^2)/n_0
          
          
          fit_l0_T_ora_test =  cv_l_0_est(y_A_project) 
          
          result_l0_T_ora[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_ora_test-f_0)^2)/n_0
          
          
          fit_l0_T_ora_target_test =  cv_l_0_est(y_ora_target_project) 
          
          result_l0_T_ora_target[sig.delta_1_index, H_k_index, A_index,gamma_index, simu_test] = sum((fit_l0_T_ora_target_test-f_0)^2)/n_0
          
          
          

          
          
        }
      }
    }
  }
}




value_1 = c(mean(result_fl[1,1,1,1,]), mean(result_fl_T_1[1,1,1,1,]), mean(result_fl_T_1_target[1,1,1,1,]), mean(result_fl_T_ora[1,1,1,1,]),  mean(result_fl_T_ora_target[1,1,1,1,]), 
            mean(result_l0[1,1,1,1,]), mean(result_l0_T_1[1,1,1,1,]), mean(result_l0_T_1_target[1,1,1,1,]), mean(result_l0_T_ora[1,1,1,1,]),  mean(result_l0_T_ora_target[1,1,1,1,]),
            mean(result_fl[1,1,2,1,]), mean(result_fl_T_1[1,1,2,1,]), mean(result_fl_T_1_target[1,1,2,1,]), mean(result_fl_T_ora[1,1,2,1,]),  mean(result_fl_T_ora_target[1,1,2,1,]),
            mean(result_l0[1,1,2,1,]), mean(result_l0_T_1[1,1,2,1,]), mean(result_l0_T_1_target[1,1,2,1,]), mean(result_l0_T_ora[1,1,2,1,]),  mean(result_l0_T_ora_target[1,1,2,1,]),
            mean(result_fl[1,1,3,1,]), mean(result_fl_T_1[1,1,3,1,]), mean(result_fl_T_1_target[1,1,3,1,]), mean(result_fl_T_ora[1,1,3,1,]),  mean(result_fl_T_ora_target[1,1,3,1,]), 
            mean(result_l0[1,1,3,1,]), mean(result_l0_T_1[1,1,3,1,]), mean(result_l0_T_1_target[1,1,3,1,]), mean(result_l0_T_ora[1,1,3,1,]),  mean(result_l0_T_ora_target[1,1,3,1,]),
            mean(result_fl[1,1,4,1,]), mean(result_fl_T_1[1,1,4,1,]), mean(result_fl_T_1_target[1,1,4,1,]), mean(result_fl_T_ora[1,1,4,1,]),  mean(result_fl_T_ora_target[1,1,4,1,]),
            mean(result_l0[1,1,4,1,]), mean(result_l0_T_1[1,1,4,1,]), mean(result_l0_T_1_target[1,1,4,1,]), mean(result_l0_T_ora[1,1,4,1,]),  mean(result_l0_T_ora_target[1,1,4,1,]))


se_1 = c(std.error(result_fl[1,1,1,1,]), std.error(result_fl_T_1[1,1,1,1,]), std.error(result_fl_T_1_target[1,1,1,1,]), std.error(result_fl_T_ora[1,1,1,1,]),  std.error(result_fl_T_ora_target[1,1,1,1,]), 
         std.error(result_l0[1,1,1,1,]), std.error(result_l0_T_1[1,1,1,1,]), std.error(result_l0_T_1_target[1,1,1,1,]), std.error(result_l0_T_ora[1,1,1,1,]),  std.error(result_l0_T_ora_target[1,1,1,1,]),
         std.error(result_fl[1,1,2,1,]), std.error(result_fl_T_1[1,1,2,1,]), std.error(result_fl_T_1_target[1,1,2,1,]), std.error(result_fl_T_ora[1,1,2,1,]),  std.error(result_fl_T_ora_target[1,1,2,1,]), 
         std.error(result_l0[1,1,2,1,]), std.error(result_l0_T_1[1,1,2,1,]), std.error(result_l0_T_1_target[1,1,2,1,]), std.error(result_l0_T_ora[1,1,2,1,]),  std.error(result_l0_T_ora_target[1,1,2,1,]),
         std.error(result_fl[1,1,3,1,]), std.error(result_fl_T_1[1,1,3,1,]), std.error(result_fl_T_1_target[1,1,3,1,]), std.error(result_fl_T_ora[1,1,3,1,]),  std.error(result_fl_T_ora_target[1,1,3,1,]), 
         std.error(result_l0[1,1,3,1,]), std.error(result_l0_T_1[1,1,3,1,]), std.error(result_l0_T_1_target[1,1,3,1,]), std.error(result_l0_T_ora[1,1,3,1,]),  std.error(result_l0_T_ora_target[1,1,3,1,]),
         std.error(result_fl[1,1,4,1,]), std.error(result_fl_T_1[1,1,4,1,]), std.error(result_fl_T_1_target[1,1,4,1,]), std.error(result_fl_T_ora[1,1,4,1,]),  std.error(result_fl_T_ora_target[1,1,4,1,]),
         std.error(result_l0[1,1,4,1,]), std.error(result_l0_T_1[1,1,4,1,]), std.error(result_l0_T_1_target[1,1,4,1,]), std.error(result_l0_T_ora[1,1,4,1,]), std.error(result_l0_T_ora_target[1,1,4,1,]))

Method = rep(c('FL', 'FL-T-1', 'FL-T-1-target', 'FL-T-ora', 'FL-T-ora-target', 
               'L0','L0-T-1', 'L0-T-1-target', 'L0-T-ora',  'L0-T-ora-target'), times = 4)



A = rep(c(2, 4, 6, 8), each = 10)

d_1 = data.frame(A = A, Method = Method, value = value_1, se = se_1)

d_1$Method <- factor(d_1$Method, levels = c('FL', 'FL-T-1', 'FL-T-1-target', 'FL-T-ora', 'FL-T-ora-target', 
                                            'L0','L0-T-1', 'L0-T-1-target', 'L0-T-ora',  'L0-T-ora-target'))

color_list <- setNames(c("#66c2a5", "#fc8d26", "#e78ac3", "#8da0cb", "#a6d854", "#484B8A", "#e5c494", "#b3b3bc", "#ffd92f", "#C35817"),
                       levels(d_1$Method))


plott_1 = ggplot(data = d_1, mapping = aes(x = A, y = value, colour = Method)) +
  geom_line() + 
  geom_point()+ 
  geom_errorbar(aes(ymin=value-se, ymax= value+se), width= 0.2)+
    labs( y ="Estimation error", x = TeX("$a$"))+
  scale_color_manual(
    values = color_list,
    labels = c(TeX("$\u2113_1$"), TeX("$\u2113_1$-T-1"), TeX("$\u2113_1$-T-$1$-target"), TeX("$\u2113_1$-T-$A$"), TeX("$\u2113_1$-T-$A$-target"),  
               TeX("$\u2113_0$"), TeX("$\u2113_0$-T-1"), TeX("$\u2113_0$-T-$1$-target"), TeX("$\u2113_0$-T-$A$"), TeX("$\u2113_0$-T-$A$-target"))
  ) +
  theme_classic()

print(plott_1)


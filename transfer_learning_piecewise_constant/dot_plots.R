library(genlasso)
library(plotrix)
library(ggplot2)
source("transfer_learning_functions.R")
library(tictoc)
library(latex2exp)
library(changepoints)
library(egg)


#Plot 1: configuration 1 and vary alpha  

n_0 = 400
n_k = 800
K = 20


sigma  = 0.5
h_1 = 100
h_2 = 400
t_k = 150

sig.delta_1_list = c(0.25, 0.5, 0.75, 1)
sig.delta_2  = 3
H_k_list = c(0.2)
A_list = c(10)
gamma_list = c(0.5)

jump.location_equal = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
jump.location_unequal = c(0, 0.2,  0.22,  0.24,   0.4,  0.42,  0.6,  0.62,   0.8,   0.82)


gamma = 0.5
jump.mean = c(2*gamma, 4*gamma, 1*gamma, 5*gamma, 7*gamma, 8*gamma, 2*gamma, 1*gamma, 3*gamma, 4*gamma)

set.seed(10)


coef_f_eq  = f_gen(jump.mean, jump.location_equal, n_0, n_k, K, size_A = 10, H_k = 0.2*n_k, sig.delta_1 = 0.5, sig.delta_2 = 3, h_1 = 100, h_2 = 400, exact = T)
f_0_eq = coef_f_eq$f_0
W_eq = coef_f_eq$W
obs_eq = obs_gen(f_0_eq, W_eq, sigma)$y_0

coef_f_un  = f_gen(jump.mean, jump.location_unequal, n_0, n_k, K, 10, 0.2*n_k, 0.5, 3, 100, 400, exact = T)
f_0_un = coef_f_un$f_0
W_un = coef_f_un$W
obs_un = obs_gen(f_0_un, W_un, sigma)$y_0


data_equal <- data.frame(
  x = c(1:n_0),
  obs = obs_eq,
  f_0 = f_0_eq
)

data_unequal <- data.frame(
  x = c(1:n_0),
  obs = obs_un,
  f_0 = f_0_un
)

plot_uneq = ggplot(data_unequal, aes(x = x)) +
  geom_point(aes(y = obs), size = 0.5, color = 'black') +
  geom_line(aes(y = f_0), color = 'black', linewidth = 1) +
  labs(title = "", x = "i", y = expression(y[i])) +
  theme_classic() 

plot_uneq 

plot_eq = ggplot(data_equal, aes(x = x)) +
  geom_point(aes(y = obs), size = 0.5, color = 'black') +
  geom_line(aes(y = f_0), color = 'black', linewidth = 1) +
  labs(title = "", x = "i", y = expression(y[i])) +
    theme_classic() 
  
plot_eq 


ggarrange(plot_eq, plot_uneq, ncol = 2)

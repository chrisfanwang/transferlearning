form.truth <- function(jump.mean_new, jump.location_new, n_new){
  tmp = seq(0,1,length.out=n_new)
  jump.location_new_2 = sapply(jump.location_new,function(x){max(min(which(tmp>=x)),1)-1})
  jump.location_new_2 = sort(jump.location_new_2)
  jump.location_new_3 = c(jump.location_new_2,n_new)
  jump.location_new_3[1] = 0
  return(rep(jump.mean_new, times=diff(jump.location_new_3)))
}

P_n0_n<-function(n0_new, n_new){
  P = matrix(0,nrow = n0_new, ncol = n_new)
  m_ratio = n_new/n0_new
  for (i in 1:n0_new){
    for (j in 1:n_new){
      if (j %in% c((ceiling((i-1)*m_ratio)+1): ceiling(i*m_ratio)) ){
        P[i,j] = 1/( ceiling(i*m_ratio)- ceiling((i-1)*m_ratio))
      }
    }
  }
  return(P)
}

P_n_n0<-function(n0_new, n_new){
  P = matrix(0,ncol = n0_new, nrow = n_new)
  m = n_new/n0_new
  for (j in 1:n0_new){
    for (i in 1:n_new){
      if (i %in% c((ceiling((j-1)*m)+1): ceiling(j*m)) ){
        P[i,j] = 1
      }
    }
  }
  return(P)
}

f_gen <- function(jump.mean_new, jump.location_new, n0_new, n_new, size_K, size_A, H_k, sig.delta_1, sig.delta_2, h_1, h_2, exact = T){
  f_0 = form.truth(jump.mean_new, jump.location_new, n0_new)
  W = matrix(NA, nrow = size_K, ncol = n_new)
  for (k in 1:size_K){
    sample_k = sample(1:n_new, H_k, replace = F) 
    if(k <= size_A){
      if(exact){
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + rep(sig.delta_1, H_k)
      }else{
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + rnorm(H_k, 0, h_1/100)
      }
    }else{
      if(exact){
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + rep(sig.delta_2, H_k)
      }else{
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + rnorm(H_k, 0, h_2/100)
      }
    }
  }
  return(list(f_0 = f_0, W=W))
}

obs_gen <- function(f_0, W, sigma){
  y_0 = f_0 + rnorm(length(f_0), 0, sigma^2)
  obs_y = matrix(NA, nrow = dim(W)[1], ncol = dim(W)[2])
  K.length = dim(W)[1]
  for (k in 1: K.length){
    obs_y[k, ] =  W[k, ] + rnorm(dim(W)[2], 0, sigma^2)
  }
  return(list(y_0 = y_0, obs_y =obs_y))
}


# Returns a sequence of integers from a to b if a <= b, otherwise nothing.
Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
  
}



cv_l_0_est <-function(y_tar){
  lamb.init = max((y_tar-mean(y_tar))^2)
  max_in = 100
  gamma_ratio = 0.95
  gamma_set = rep(NA, max_in) 
  for(i in (1:max_in)){
    if (i==1){
      gamma_set[i] = lamb.init
    }
    else{
      gamma_set[i] = gamma_set[i-1]*gamma_ratio
    }
  }
  cv_result = CV.search.DP.univar(y_tar, gamma_set, delta = 1)
  min_idx = which.min(cv_result$test_error)
  final_fit = DP.univar(y_tar, gamma_set[min_idx], delta = 1)$yhat 
  
}

A_detect_alg_calculation <- function(y_new, y_old, t_hat_k){
  n_new = length(y_new)
  n0_new = length(y_old)
  Delta_hat =  y_new/sqrt(n_new) -  (P_n_n0(n0_new, n_new) %*% y_old)/sqrt(n_new)
  Delta_hat_subset = sort(Delta_hat, decreasing = TRUE)[1:t_hat_k]
  Delta_check = sum(Delta_hat_subset^2)
  return(Delta_check)
}


A_detect_threshold <- function(y_new, y_old, t_hat_k, trials = 100, quant = 0.99){
  
  result_test =fusedlasso1d(y_new)
  cv_test = cv.trendfilter(result_test, k=5)
  fit_test = coef(result_test, lambda =cv_test$lambda.min)$beta

  
  residual = y_new - fit_test
  
  detect_calculation = rep(NA, trials)
  
  for (trial in 1:trials) {
    set.seed(100*trial)
    y_new_2 =fit_test + residual[sample(length(residual))]
    detect_calculation[trial] = A_detect_alg_calculation(y_new_2, y_old, t_hat_k)
  }
  detect_ini = A_detect_alg_calculation(y_new, y_old, t_hat_k)
  
  return(max(quantile(detect_calculation, probs = quant),  detect_ini ))
}



A_detect_alg <- function(aux_matrix, y_old,t_hat_k, tau){
  check_length = dim(aux_matrix)[1]
  A_est = c()
  for (k in 1:check_length){
    check_Delta =  A_detect_alg_calculation(aux_matrix[k, ], y_old, t_hat_k)
    if (check_Delta  <= tau){
      A_est = c(A_est, k)
    }
    
  }
  return(A_est)
}
  


A_detect_index <- function(y_matrix, y_old){
  source_size = dim(y_matrix)[2]
  source_number = dim(y_matrix)[1]
  target_size = length(y_old)
  result_test = rep(0, source_number)
  for (source_index in 1:source_number){
    y_new_project = P_n0_n(target_size, source_size) %*% y_matrix[source_index,]
    source_result =fusedlasso1d(y_new_project)
    source_cv = cv.trendfilter(source_result, k=5)
    source_fit = coef(source_result, lambda = source_cv$lambda.min)$beta


    result_test[source_index] =  sum((source_fit-y_old)^2)/target_size
  }
  return(which.min( result_test))
  

}


f_gen_new <- function(jump.mean_new, jump.location_new, n0_new, n_new_size, H_k, sig.delta_1,  h_1,  exact = T){
  f_0 = form.truth(jump.mean_new, jump.location_new, n0_new)
  W = list()
  for (k in 1:length( n_new_size)){
    sample_k = sample(1:n_new_size[k], H_k*n_new_size[k], replace = F) 
    new_signal = rep(NA, n_new_size[k])
    if(exact){
      new_signal[-sample_k] = (P_n_n0(n0_new, n_new_size[k]) %*% f_0)[-sample_k, ] 
      new_signal[sample_k] = (P_n_n0(n0_new, n_new_size[k]) %*% f_0)[sample_k, ]  + rep(sig.delta_1, H_k*n_new_size[k])
    }else{
      new_signal[-sample_k]= (P_n_n0(n0_new,  n_new_size[k]) %*% f_0)[-sample_k, ] 
      new_signal[sample_k] = (P_n_n0(n0_new,  n_new_size[k]) %*% f_0)[sample_k, ]  + rnorm(H_k*n_new_size[k], 0, h_1/100)
    }
    W = c(W, list(new_signal))
  }
  return(list(f_0 = f_0, W=W))
}


obs_gen_new <- function(f_0, W, sigma){
  y_0 = f_0 + rnorm(length(f_0), 0, sigma^2)
  obs_y = list()
  K.length = length(W)
  for (k in 1: K.length){
    obs_y_new =  W[[k]] + rnorm(length(W[[k]]), 0, sigma^2)
    obs_y = c(obs_y, list(obs_y_new))
  }
  
  return(list(y_0 = y_0, obs_y =obs_y))
}



  
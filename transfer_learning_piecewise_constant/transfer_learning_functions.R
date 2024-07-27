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

f_gen_tempde <- function(jump.mean_new, jump.location_new, n0_new, n_new, size_K, size_A, H_k,  h_1, h_2, temp_h){
  f_0 = form.truth(jump.mean_new, jump.location_new, n0_new)
  W = matrix(NA, nrow = size_K, ncol = n_new)
  for (k in 1:size_K){
    sample_k = sort(sample(1:n_new, H_k, replace = F)) 
    if(k <= size_A){
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        error = rep(NA, length(sample_k))
        error[1] = rnorm(1, 0, h_1/100)
        for (i in 1:(length(sample_k)-1)) {
          error[i+1] = temp_h*error[i] + (1-temp_h)*rnorm(1, 0, h_1/100)
        }
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + error
      }
    else{
        W[k, -sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[-sample_k, ] 
        error = rep(NA, length(sample_k))
        error[1] = rnorm(1, 0, h_2/100)
        for (i in 1:(length(sample_k)-1)) {
          error[i+1] = temp_h*error[i] + (1-temp_h)*rnorm(1, 0, h_2/100)
        }
        W[k, sample_k] = (P_n_n0(n0_new, n_new) %*% f_0)[sample_k, ]  + error
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


obs_gen_tempde<- function(f_0, W, sigma, coef_temp_1, coef_temp_2 = 0){
  error = rep(NA, length(f_0))
  error[1] = rnorm(1, 0, sigma^2)
  for (i in 1:(length(f_0)-1)) {
    error[i+1] = coef_temp_1*error[i] + (1-coef_temp_1)*rnorm(1, 0, sigma^2)
  }
  y_0 = f_0 + error
  
  error_y = matrix(NA, nrow = dim(W)[1], ncol = dim(W)[2])
  K.length = dim(W)[1]
  
  for (k in 1: K.length){
    error_y[k, 1] =  rnorm(1, 0, sigma^2)
    for (i in 1:( dim(W)[2]-1)) {
      error_y[k, i+1] = coef_temp_2*error_y[k, i] +  (1-coef_temp_2)*rnorm(1, 0, sigma^2)
    } 
  }
  obs_y = W + error_y
  return(list(y_0 = y_0, obs_y =obs_y))
}



# Returns a sequence of integers from a to b if a <= b, otherwise nothing.
Seq <- function(a, b, ...) {
  if (a<=b) return(seq(a,b,...))
  else return(numeric(0))
  
}


CV.l0 = function(y, k = 5, gamma, delta = 1){
    n = length(y)
    foldid = c(0,rep(Seq(1,k),n-2)[Seq(1,n-2)],0)
    cvall = rep(0,k)
    pos = 1:n
    for (i in Seq(1,k)) {
  
    otr = which(foldid!=i)
    ntr = length(otr)
    ytr = y[otr]
    ptr = pos[otr]
    
    b = DP.univar(ytr, gamma, delta = 1)$yhat 
    
    
    ote = which(foldid==i)
    yte = y[ote]
    pte = pos[ote]
    ilo = which((Seq(1,n)%in%(ote-1))[otr])
    ihi = which((Seq(1,n)%in%(ote+1))[otr])
    a = (pte - ptr[ilo])/(ptr[ihi]-ptr[ilo])
    pred = b[ilo]*(1-a) + b[ihi]*a
    cvall[i] = mean((yte-pred)^2)
    }
    
  cverr = mean(cvall)
  
  result = list(test_error = cverr)
  return(result)
  
}  


cv_l_0_est <-function(y_tar){
  lamb.init = max((y_tar-mean(y_tar))^2)
  max_in = 1000
  gamma_ratio = 0.99
  gamma_set = rep(NA, max_in) 
  for(i in (1:max_in)){
    if (i==1){
      gamma_set[i] = lamb.init
    }
    else{
      gamma_set[i] = gamma_set[i-1]*gamma_ratio
    }
  }
  cv_result = sapply(1:length(gamma_set), function(j) CV.l0(y_tar, k = 5, gamma_set[j], delta = 1)$test_error)
  min_idx = which.min(cv_result)
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


insert_values <- function(y, y1) {
  n0 <- length(y)
  n1 <- length(y1)
  n_total <- n0 + n1
  
  yt <- numeric(n_total)
  
  insert_positions <- ceiling((1:n0) * n1 / n0) + (1:n0)
  
  yt[insert_positions] <- y
  
  remaining_positions <- setdiff(1:n_total, insert_positions)
  yt[remaining_positions] <- y1[1:length(remaining_positions)]
  
  return(yt)
}


alignment_matrix <- function(n, h, m) {
  if (length(m) != (h + 1)) {
    stop("Length of m must be h + 1.")
  }
  
  total_m <- sum(m)
  
  P <- matrix(0, nrow = n, ncol = total_m)
  
  for (i in 1:n) {
    start_indices <- sapply(m, function(mk) ceiling((i-1) * mk / n))
    end_indices <- sapply(m, function(mk) ceiling(i * mk / n))
    
    start_sum <- cumsum(start_indices)
    end_sum <- cumsum(end_indices)
    
    start_j <- start_sum[h+1] + 1
    end_j <- end_sum[h+1]
    
    if (start_j <= end_j) {
      denominator <- sum(end_indices - start_indices)
      
      P[i, start_j:end_j] <- 1 / denominator
    }
  }
  
  return(P)
}


construct_vector <- function(y, y_list) {
  n0 <- length(y)
  
  n_k <- sapply(y_list, length)
  
  total_n <- sum(c(n0, n_k))
  
  y_all <- numeric(total_n)
  

  calc_J <- function(j, k, n0, n_k) {
    
    if (k == 1) {
      sum_lk <- j
    } else {
    sum_lk <- sum(sapply(1:(k-1), function(t) ceiling(j * n_k[t] / n0))) +j
    }
    if (k == length(n_k) + 1) {
      sum_rk = 0
    }
    else {
    sum_rk <- sum(sapply(k:length(n_k), function(t) ceiling((j-1) * n_k[t] / n0)))
    }
    return(sum_lk + sum_rk)
  }
  
  for (j in 1:n0) {
    J1 <- calc_J(j, 1, n0, n_k)
    y_all[J1 ] <- y[j]
  }
  
  for (k in 1:length(y_list)) {
    y_k <- y_list[[k]]
    
    nk_len <- length(y_k)
    
    
    
    for (j in 1:n0) {
      Jk <- calc_J(j, k, n0, n_k)
      Jk1 <- calc_J(j, k + 1, n0, n_k)
      indices <- (Jk + 1):Jk1
      y_k_index <- (ceiling((j-1)* nk_len / n0) +1): ceiling(j * nk_len / n0)
      if (length(indices) > 0 && length(y_k_index)> 0) {
        y_all[indices] <- y_k[y_k_index]
        y_k_index <- y_k_index + length(indices)
      }
    }
  }
  
  return(y_all)
}

# Example usage
y <- c(1, 2, 3)
y1 <- c(10, 11, 12, 13, 14)
y2 <- c(20, 21, 22, 23, 24, 25, 26)

y_all <- construct_vector(y, list(y1, y2))
print(y_all)





  
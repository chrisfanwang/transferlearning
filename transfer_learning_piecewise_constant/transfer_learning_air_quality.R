library(lubridate) #to use function ymd
library(genlasso)
library(plotrix)
library(zoo)
source("/Users/transfer_learning_piecewise_constant/transfer_learning_functions.R")
library(changepoints)




path = "/Users/chris/air quality"
fileNames = dir(path)
filepath = sapply(fileNames, function(x){paste(path, x, sep = "/")})
filepath
data_air_quality = lapply(filepath, function(x){read.csv2(x,header = T, sep = ",")})

data_air_quality = lapply(data_air_quality, function(x){
  x$date = as.Date(x$date)
  x$date = ymd(x$date)
  x = x[order(x$date,decreasing = T),]
  x = x[, c(1,2)]
  x
})

date_1 = ymd(as.Date("2023-07-01"))
date_2 = ymd(as.Date("2020-07-02"))

data_air_quality = lapply(data_air_quality, function(x){
  x[which(x$date == date_1): which(x$date == date_2), ]
})

date_interval = data_air_quality$`beijing-air-quality.csv`$date
  

date_miss =  lapply(data_air_quality, function(x){
  which(is.na(match(date_interval, x$date)))
})


data_full =  lapply(data_air_quality, function(x){
  date_miss = which(is.na(match(date_interval, x$date)))
  col_1 = date_interval[date_miss]
  col_2 = rep(NA, length(date_miss))
  data_new = data.frame(date = col_1, pm25 = col_2)
  y = rbind(x, data_new)
})

data_full = lapply(data_full, function(x){
  x = x[order(x$date,decreasing = T),]
  x
})

data_full = lapply(data_full, function(x){
  x[, 2]
})

data_target =    data_full$`paris-air-quality.csv` 

data_auxiliary = rbind(     data_full$`london-air-quality.csv`,
                            data_full$`haerbin-air-quality.csv`,
                            data_full$`bangkok-air-quality.csv`,
                            data_full$`los-angeles-north main street-air-quality.csv`,
                            data_full$`kunming-air-quality.csv`,
                            data_full$`dalian-air-quality.csv`, 
                            data_full$`xian-air-quality.csv`,
                            data_full$`tianjin-air-quality.csv`, 
                            data_full$`central,-singapore-air-quality.csv`,
                            data_full$`amsterdam-air-quality.csv`, 
                            data_full$`beijing-air-quality.csv`,
                            data_full$`chongqing-air-quality.csv`,
                            data_full$`hefei-air-quality.csv`,
                            data_full$`seoul-air-quality.csv`,
                            data_full$`sanya-air-quality.csv`,
                            data_full$`shanghai-air-quality.csv`,
                            data_full$`hongkong-air-quality.csv`, 
                            data_full$`sternschanze,-hamburg, germany-air-quality.csv` 
                            )
data_auxiliary = data.matrix(data_auxiliary)





data_auxiliary = t(apply(data_auxiliary, 1, function(x) na.locf(x, na.rm = FALSE)))

data_target = na.locf(data_target, na.rm = FALSE)


week_length = length(data_target)/7

data_target_weekly = rep(NA, week_length)



for (i in 1:week_length){
  data_target_weekly[i] = data_target[(i-1)*7+1]
}  


train_set  = seq(1, week_length, by = 2)
test_set= seq(2, week_length, by = 2)

data_target_weekly_train = data_target_weekly[train_set]

data_target_weekly_test = data_target_weekly[test_set]

dim_1 = dim(data_auxiliary)[1]

dim_2 = dim(data_auxiliary)[2]

data_auxiliary = t(scale(t(data_auxiliary)))


data_target_weekly_train= scale(data_target_weekly_train)

data_target_weekly_test = scale(data_target_weekly_test)




n_0 = length(data_target_weekly_train)

result_fl_train = fusedlasso1d(data_target_weekly_train)

cv_fl_train = cv.trendfilter(result_fl_train, k=5)

fit_fl_train = coef(result_fl_train, lambda = cv_fl_train$lambda.min)$beta

result_fl = sum((fit_fl_train-data_target_weekly_test)^2)/n_0



index_k = A_detect_index(data_auxiliary, data_target_weekly_train)


n_k = dim(data_auxiliary)[2]

y_1_project = as.vector(P_n0_n(n_0, n_k) %*% data_auxiliary[1, ])

result_fl_T_1_train = fusedlasso1d(y_1_project)

cv_fl_T_1_train = cv.trendfilter(result_fl_T_1_train, k=5)

fit_fl_T_1_train = coef(result_fl_T_1_train, lambda = cv_fl_T_1_train$lambda.min)$beta

result_fl_T_1 = sum((fit_fl_T_1_train -data_target_weekly_test)^2)/n_0





fit_l0_train =  cv_l_0_est(data_target_weekly_train)

result_l0 = sum((fit_l0_train-data_target_weekly_test)^2)/n_0




fit_l0_T_1_train =  cv_l_0_est(y_1_project)

result_l0_T_1 = sum((fit_l0_T_1_train-data_target_weekly_test)^2)/n_0




t_k = 200

est_tau = A_detect_threshold(data_auxiliary[index_k, ], data_target_weekly_train, t_hat_k = t_k)

est_A = A_detect_alg(data_auxiliary, data_target_weekly_train, t_hat_k = t_k, est_tau)

print(est_A)


y_A_est_project = rep(0, n_0)

for (k in 1:length(est_A)){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[est_A[k], ]
  y_A_est_project = y_A_est_project + y_project
}
y_A_est_project = y_A_est_project/length(est_A)


result_fl_T_est_train = fusedlasso1d(y_A_est_project)

cv_fl_T_est_train = cv.trendfilter(result_fl_T_est_train, k=5)

fit_fl_T_est_train = coef(result_fl_T_est_train, lambda = cv_fl_T_est_train$lambda.min)$beta

result_fl_T_est = sum((fit_fl_T_est_train -data_target_weekly_test)^2)/n_0





fit_l0_T_est_train =  cv_l_0_est(y_A_est_project) 

result_l0_T_est = sum((fit_l0_T_est_train-data_target_weekly_test)^2)/n_0


y_all_project = rep(0, n_0)

for (k in 1:(dim(data_auxiliary)[1])){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[k, ]
  y_all_project = y_all_project + y_project
}
y_all_project = y_all_project/dim(data_auxiliary)[1]


result_fl_T_all_train = fusedlasso1d(y_all_project)

cv_fl_T_all_train = cv.trendfilter(result_fl_T_all_train, k=5)

fit_fl_T_all_train = coef(result_fl_T_all_train, lambda = cv_fl_T_all_train$lambda.min)$beta

result_fl_T_all = sum((fit_fl_T_all_train -data_target_weekly_test)^2)/n_0



fit_l0_T_all_train =  cv_l_0_est(y_all_project) 

result_l0_T_all = sum((fit_l0_T_all_train-data_target_weekly_test)^2)/n_0




result_fl
result_fl_T_1
result_fl_T_est
result_fl_T_all
result_l0
result_l0_T_1
result_l0_T_est
result_l0_T_all 
est_A
index_k 

# London


data_target = data_full$`london-air-quality.csv` 

data_auxiliary = rbind(     data_full$`paris-air-quality.csv`,
                            data_full$`amsterdam-air-quality.csv`, 
                            data_full$`dalian-air-quality.csv`, 
                            data_full$`xian-air-quality.csv`,
                            data_full$`tianjin-air-quality.csv`, 
                            data_full$`bangkok-air-quality.csv`, 
                            data_full$`beijing-air-quality.csv`,
                            data_full$`central,-singapore-air-quality.csv`,
                            data_full$`chongqing-air-quality.csv`,
                            data_full$`haerbin-air-quality.csv`,
                            data_full$`hefei-air-quality.csv`,
                            data_full$`hongkong-air-quality.csv`, 
                            data_full$`kunming-air-quality.csv`,
                            data_full$`los-angeles-north main street-air-quality.csv`,
                            data_full$`sanya-air-quality.csv`,
                            data_full$`seoul-air-quality.csv`,
                            data_full$`shanghai-air-quality.csv`, 
                            data_full$`sternschanze,-hamburg, germany-air-quality.csv` 
)
data_auxiliary = data.matrix(data_auxiliary)






data_auxiliary = t(apply(data_auxiliary, 1, function(x) na.locf(x, na.rm = FALSE)))

data_target = na.locf(data_target, na.rm = FALSE)


week_length = length(data_target)/7

data_target_weekly = rep(NA, week_length)



for (i in 1:week_length){
  data_target_weekly[i] = data_target[(i-1)*7+1]
}  


train_set  = seq(1, week_length, by = 2)
test_set= seq(2, week_length, by = 2)

data_target_weekly_train = data_target_weekly[train_set]

data_target_weekly_test = data_target_weekly[test_set]

dim_1 = dim(data_auxiliary)[1]

dim_2 = dim(data_auxiliary)[2]

data_auxiliary = t(scale(t(data_auxiliary)))


data_target_weekly_train= scale(data_target_weekly_train)

data_target_weekly_test = scale(data_target_weekly_test)




n_0 = length(data_target_weekly_train)

result_fl_train_new = fusedlasso1d(data_target_weekly_train)

cv_fl_train_new = cv.trendfilter(result_fl_train_new, k=5)

fit_fl_train_new = coef(result_fl_train_new, lambda = cv_fl_train_new$lambda.min)$beta

result_fl = sum((fit_fl_train_new-data_target_weekly_test)^2)/n_0





n_k = dim(data_auxiliary)[2]

y_1_project = as.vector(P_n0_n(n_0, n_k) %*% data_auxiliary[1, ])

result_fl_T_1_train = fusedlasso1d(y_1_project)

cv_fl_T_1_train = cv.trendfilter(result_fl_T_1_train, k=5)

fit_fl_T_1_train = coef(result_fl_T_1_train, lambda = cv_fl_T_1_train$lambda.min)$beta

result_fl_T_1 = sum((fit_fl_T_1_train -data_target_weekly_test)^2)/n_0





fit_l0_train =  cv_l_0_est(data_target_weekly_train)

result_l0 = sum((fit_l0_train-data_target_weekly_test)^2)/n_0




fit_l0_T_1_train =  cv_l_0_est(y_1_project)

result_l0_T_1 = sum((fit_l0_T_1_train-data_target_weekly_test)^2)/n_0


index_k = A_detect_index(data_auxiliary, data_target_weekly_train)


t_k = 200

est_tau = A_detect_threshold(data_auxiliary[index_k, ], data_target_weekly_train, t_hat_k = t_k)

est_A = A_detect_alg(data_auxiliary, data_target_weekly_train, t_hat_k = t_k, est_tau)

print(est_A)


y_A_est_project = rep(0, n_0)

for (k in 1:length(est_A)){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[est_A[k], ]
  y_A_est_project = y_A_est_project + y_project
}
y_A_est_project = y_A_est_project/length(est_A)


result_fl_T_est_train = fusedlasso1d(y_A_est_project)

cv_fl_T_est_train = cv.trendfilter(result_fl_T_est_train, k=5)

fit_fl_T_est_train = coef(result_fl_T_est_train, lambda = cv_fl_T_est_train$lambda.min)$beta

result_fl_T_est = sum((fit_fl_T_est_train -data_target_weekly_test)^2)/n_0




fit_l0_T_est_train =  cv_l_0_est(y_A_est_project) 

result_l0_T_est = sum((fit_l0_T_est_train-data_target_weekly_test)^2)/n_0


y_all_project = rep(0, n_0)

for (k in 1:(dim(data_auxiliary)[1])){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[k, ]
  y_all_project = y_all_project + y_project
}
y_all_project = y_all_project/dim(data_auxiliary)[1]


result_fl_T_all_train = fusedlasso1d(y_all_project)

cv_fl_T_all_train = cv.trendfilter(result_fl_T_all_train, k=5)

fit_fl_T_all_train = coef(result_fl_T_all_train, lambda = cv_fl_T_all_train$lambda.min)$beta

result_fl_T_all = sum((fit_fl_T_all_train -data_target_weekly_test)^2)/n_0



fit_l0_T_all_train =  cv_l_0_est(y_all_project) 

result_l0_T_all = sum((fit_l0_T_all_train-data_target_weekly_test)^2)/n_0




result_fl
result_fl_T_1
result_fl_T_est
result_fl_T_all
result_l0
result_l0_T_1
result_l0_T_est
result_l0_T_all 
est_A
index_k


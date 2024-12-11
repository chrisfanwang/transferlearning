library(lubridate) #to use function ymd
library(genlasso)
library(plotrix)
library(zoo)
source("transfer_learning_functions.R")
library(changepoints)



# Load GDP data
filepath = "gdp_q.csv"
data_GDP  = read.csv2(filepath,header = T, sep = ",")





# Extract target country data ("Hungary") and preprocess
data_target =  data_GDP[which(data_GDP$Country =="Hungary" )  , ]
data_target = as.numeric(as.vector(unlist(data_target))[-c(1)])




# Auxiliary countries considered for analysis
# "Bulgaria", "Croatia", "Czechia", "Estonia", "Greece", "Hungary", "Latvia", 
# "Lithuania", "Poland (kind of)", "Romania", "Slovak Republic", "Slovenia"

# Load Industrial Production (IP) data


filepath = "ip.csv"
data_ip  = read.csv2(filepath,header = T, sep = ",")
data_ip = data_ip[, -c(278:286)]
dim(data_ip)

# Convert IP data to numeric and matrix formats

data_ip_df <- as.data.frame(data_ip)

# Remove the first column (Country names)
data_ip_df <- data_ip_df[, -1]
# Convert all values to numeric
data_ip_numeric <- as.data.frame(lapply(data_ip_df, as.numeric))
data_ip_numeric <- as.data.frame(lapply(data_ip_df, as.numeric))
# Convert the data frame to a matrix
data_ip_matrix <- as.matrix(data_ip_numeric)
data_auxiliary = data_ip_matrix

# Compute auxiliary data by taking 12-month differences (seasonal adjustment)
data_auxiliary<- t(apply(data_auxiliary, 1, function(x) diff(x, lag = 12)))
dim(data_auxiliary)


# Standardize auxiliary data
year_length = length(data_target)
dim_1 = dim(data_auxiliary)[1]
dim_2 = dim(data_auxiliary)[2]
data_auxiliary = t(scale(t(data_auxiliary)))


# Split data into training and testing sets
train_set  = seq(1,year_length, by = 2)
test_set= seq(2, year_length, by = 2)
data_target_train = data_target[train_set]
data_target_test = data_target[test_set]

# Scale training and testing data
data_target_train = scale(data_target_train)
data_target_test = scale(data_target_test)




n_0 = length(data_target_train)

# Fused lasso on target training data
result_fl_train = fusedlasso1d(data_target_train)
cv_fl_train = cv.trendfilter(result_fl_train, k=5)
fit_fl_train = coef(result_fl_train, lambda = cv_fl_train$lambda.min)$beta
result_fl = sum((fit_fl_train-data_target_test)^2)/n_0


# unisource l_0-penalized method 
n_k = dim(data_auxiliary)[2]
y_1_project = as.vector(P_n0_n(n_0, n_k) %*% data_auxiliary[6, ])
result_fl_T_1_train = fusedlasso1d(y_1_project)
cv_fl_T_1_train = cv.trendfilter(result_fl_T_1_train, k=5)
fit_fl_T_1_train = coef(result_fl_T_1_train, lambda = cv_fl_T_1_train$lambda.min)$beta
result_fl_T_1 = sum((fit_fl_T_1_train -data_target_test)^2)/n_0


# l_0-penalized method  on target training data
fit_l0_train =  cv_l_0_est(data_target_train)
result_l0 = sum((fit_l0_train-data_target_test)^2)/n_0
fit_l0_T_1_train =  cv_l_0_est(y_1_project)
result_l0_T_1 = sum((fit_l0_T_1_train-data_target_test)^2)/n_0
index_k = A_detect_index(data_auxiliary, data_target_train)

# informative source selection
t_k = 50
est_tau = A_detect_threshold(data_auxiliary[index_k, ], data_target_train, t_hat_k = t_k)
est_A = A_detect_alg(data_auxiliary, data_target_train, t_hat_k = t_k, est_tau)
print(est_A)
y_A_est_project = rep(0, n_0)

for (k in 1:length(est_A)){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[est_A[k], ]
  y_A_est_project = y_A_est_project + y_project
}
y_A_est_project = y_A_est_project/length(est_A)


# selected multisource l_1-penalized method 
result_fl_T_est_train = fusedlasso1d(y_A_est_project)
cv_fl_T_est_train = cv.trendfilter(result_fl_T_est_train, k=5)
fit_fl_T_est_train = coef(result_fl_T_est_train, lambda = cv_fl_T_est_train$lambda.min)$beta
result_fl_T_est = sum((fit_fl_T_est_train-data_target_test)^2)/n_0




# selected multisource l_0-penalized method 
fit_l0_T_est_train =  cv_l_0_est(y_A_est_project) 
result_l0_T_est= sum((fit_l0_T_est_train-data_target_test)^2)/n_0




y_all_project = rep(0, n_0)

for (k in 1:(dim(data_auxiliary)[1])){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[k, ]
  y_all_project = y_all_project + y_project
}
y_all_project = y_all_project/dim(data_auxiliary)[1]

# all multisource l_1-penalized method 
result_fl_T_all_train = fusedlasso1d(y_all_project)
cv_fl_T_all_train = cv.trendfilter(result_fl_T_all_train, k=5)
fit_fl_T_all_train = coef(result_fl_T_all_train, lambda = cv_fl_T_all_train$lambda.min)$beta
result_fl_T_all = sum((fit_fl_T_all_train-data_target_test)^2)/n_0


# all multisource l_0-penalized method 

fit_l0_T_all_train =  cv_l_0_est(y_all_project) 
result_l0_T_all = sum((fit_l0_T_all_train-data_target_test)^2)/n_0

                         

# Print results
result_fl 
result_fl_T_1
result_fl_T_est 
result_fl_T_all
result_l0
result_l0_T_1
result_l0_T_est
result_l0_T_all


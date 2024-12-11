library(lubridate) #to use function ymd
library(genlasso)
library(plotrix)
library(zoo)
source("transfer_learning_functions.R")
library(changepoints)
library(httr)
library(jsonlite)

# Get the dataset for each region
url_west <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONA&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_west <- GET(url_west)
if (http_status(response_west)$category == "Success") {
  # Parse the JSON response
  data_west <- fromJSON(content(response_west, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_genesee <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONB&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_genesee <- GET(url_genesee)
if (http_status(response_genesee)$category == "Success") {
  # Parse the JSON response
  data_genesee <- fromJSON(content(response_genesee, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_central <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONC&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_central <- GET(url_central)
if (http_status(response_central)$category == "Success") {
  # Parse the JSON response
  data_central <- fromJSON(content(response_central, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_north <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZOND&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_north <- GET(url_north)
if (http_status(response_north)$category == "Success") {
  # Parse the JSON response
  data_north <- fromJSON(content(response_north, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_MV <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONE&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_MV <- GET(url_MV)
if (http_status(response_MV)$category == "Success") {
  # Parse the JSON response
  data_MV <- fromJSON(content(response_MV, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_capital <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONF&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_capital <- GET(url_capital)
if (http_status(response_capital)$category == "Success") {
  # Parse the JSON response
  data_capital <- fromJSON(content(response_capital, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_HV <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONG&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_HV <- GET(url_HV)
if (http_status(response_HV)$category == "Success") {
  # Parse the JSON response
  data_HV <- fromJSON(content(response_HV, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_millwood <- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONH&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_millwood <- GET(url_millwood)
if (http_status(response_millwood)$category == "Success") {
  # Parse the JSON response
  data_millwood <- fromJSON(content(response_millwood, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_dunwoodie<- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONI&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_dunwoodie<- GET(url_dunwoodie)
if (http_status(response_dunwoodie)$category == "Success") {
  # Parse the JSON response
  data_dunwoodie<- fromJSON(content(response_dunwoodie, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_NYcity<- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONJ&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_NYcity<- GET(url_NYcity)
if (http_status(response_NYcity)$category == "Success") {
  # Parse the JSON response
  data_NYcity<- fromJSON(content(response_NYcity, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}
url_longisland<- "https://api.eia.gov/v2/electricity/rto/daily-region-sub-ba-data/data/?frequency=daily&data[0]=value&facets[subba][]=ZONK&facets[timezone][]=Pacific&start=2020-07-02&end=2023-07-01&sort[0][column]=period&sort[0][direction]=desc&offset=0&length=5000&api_key=utmJaP3d2cheILHEHp314wb246ijsmmFpTb2vQ5f"
response_longisland<- GET(url_longisland)
if (http_status(response_longisland)$category == "Success") {
  # Parse the JSON response
  data_longisland<- fromJSON(content(response_longisland, as = "text", encoding = "UTF-8"), flatten = TRUE)
} else {
  print("API request failed")
}

# Extract the target variable for the central region
data_target = as.numeric(data_central$response$data$value[3:1094])

# Calculate weekly averages for the target data
day_length = length(data_target)/7
data_target_weekly = rep(NA, day_length)

# Aggregate daily data into weekly data
for (i in 1:day_length){
  data_target_weekly[i ] = data_target[(i-1)*7+1]
}  

# Combine auxiliary data from all regions
data_auxiliary = rbind(as.numeric(data_west$response$data$value[3:1094]),
                       as.numeric(data_genesee$response$data$value[3:1094]), 
                       as.numeric(data_NYcity$response$data$value[3:1094]),
                       as.numeric(data_north$response$data$value[3:1094]),
                       as.numeric(data_MV$response$data$value[3:1094]), 
                       as.numeric(data_capital$response$data$value[3:1094]),
                       as.numeric(data_HV$response$data$value[3:1094]), 
                       as.numeric(data_millwood$response$data$value[3:1094]), 
                       as.numeric(data_dunwoodie$response$data$value[3:1094]),
                       as.numeric(data_longisland$response$data$value[3:1094]))

# Standardize the auxiliary data
data_auxiliary = t(scale(t(data_auxiliary)))

# Split the target data into training and testing sets
train_set  = seq(1, day_length, by = 2)
test_set= seq(2, day_length, by = 2)
data_target_weekly_train = data_target_weekly[train_set]
data_target_weekly_test = data_target_weekly[test_set]
data_target_weekly_train= scale(data_target_weekly_train)
data_target_weekly_test = scale(data_target_weekly_test)



n_0 = length(data_target_weekly_train)

# Fused lasso on target training data
result_fl_train = fusedlasso1d(data_target_weekly_train)
cv_fl_train = cv.trendfilter(result_fl_train, k=5)
fit_fl_train = coef(result_fl_train, lambda = cv_fl_train$lambda.min)$beta
result_fl = sum((fit_fl_train-data_target_weekly_test)^2)/n_0

# l_0-penalized method  on target training data
fit_l0_train =  cv_l_0_est(data_target_weekly_train)
result_l0 = sum((fit_l0_train-data_target_weekly_test)^2)/n_0



# unisource l_1-penalized method 
n_k = dim(data_auxiliary)[2]
y_1_project = P_n0_n(n_0, n_k) %*% data_auxiliary[5,]
result_fl_T_1_train = fusedlasso1d(y_1_project)
cv_fl_T_1_train = cv.trendfilter(result_fl_T_1_train, k=5)
fit_fl_T_1_train = coef(result_fl_T_1_train, lambda = cv_fl_T_1_train$lambda.min)$beta
result_fl_T_1 = sum((fit_fl_T_1_train -data_target_weekly_test)^2)/n_0


# l_0-penalized method  on target training data
fit_l0_T_1_train =  cv_l_0_est(y_1_project)
result_l0_T_1 = sum((fit_l0_T_1_train-data_target_weekly_test)^2)/n_0



# informative source selection
index_k = A_detect_index(data_auxiliary, data_target_weekly_train)
index_k
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



# selected multisource l_1-penalized method 
result_fl_T_est_train = fusedlasso1d(y_A_est_project)
cv_fl_T_est_train = cv.trendfilter(result_fl_T_est_train, k=5)
fit_fl_T_est_train = coef(result_fl_T_est_train, lambda = cv_fl_T_est_train$lambda.min)$beta
result_fl_T_est = sum((fit_fl_T_est_train -data_target_weekly_test)^2)/n_0


# selected multisource l_0-penalized method 
fit_l0_T_est_train =  cv_l_0_est(y_A_est_project) 
result_l0_T_est = sum((fit_l0_T_est_train-data_target_weekly_test)^2)/n_0


# all multisource l_1-penalized method 
y_all_project = rep(0, n_0)
for (k in 1:dim(data_auxiliary)[1]){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[k, ]
  y_all_project = y_all_project + y_project
}
y_all_project = y_all_project/dim(data_auxiliary)[1]

result_fl_T_all_train = fusedlasso1d(y_all_project)
cv_fl_T_all_train = cv.trendfilter(result_fl_T_all_train, k=5)
fit_fl_T_all_train = coef(result_fl_T_all_train, lambda = cv_fl_T_all_train$lambda.min)$beta
result_fl_T_all = sum((fit_fl_T_all_train -data_target_weekly_test)^2)/n_0


# all multisource l_0-penalized method 
fit_l0_T_all_train =  cv_l_0_est(y_all_project) 
result_l0_T_all = sum((fit_l0_T_all_train-data_target_weekly_test)^2)/n_0


# Print results
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





















# Extract the target variable for NYcity

data_target = as.numeric(data_NYcity$response$data$value[3:1094])

# Calculate weekly averages for the target data
day_length = length(data_target)/7
data_target_weekly = rep(NA, day_length)
for (i in 1:day_length){
  data_target_weekly[i ] = data_target[(i-1)*7+1]
}  


# Combine auxiliary data from all regions
data_auxiliary = rbind(as.numeric(data_west$response$data$value[3:1094]),
                       as.numeric(data_genesee$response$data$value[3:1094]), 
                       as.numeric(data_central$response$data$value[3:1094]),
                       as.numeric(data_north$response$data$value[3:1094]),
                       as.numeric(data_MV$response$data$value[3:1094]), 
                       as.numeric(data_capital$response$data$value[3:1094]),
                       as.numeric(data_HV$response$data$value[3:1094]), 
                       as.numeric(data_millwood$response$data$value[3:1094]), 
                       as.numeric(data_dunwoodie$response$data$value[3:1094]),
                       as.numeric(data_longisland$response$data$value[3:1094]))
data_auxiliary = t(scale(t(data_auxiliary)))

# Split the target data into training and testing sets
train_set  = seq(1, day_length, by = 2)
test_set= seq(2, day_length, by = 2)
data_target_weekly_train = data_target_weekly[train_set]
data_target_weekly_test = data_target_weekly[test_set]
data_target_weekly_train= scale(data_target_weekly_train)
data_target_weekly_test = scale(data_target_weekly_test)





n_0 = length(data_target_weekly_train)

# Fused lasso on target training data
result_fl_train = fusedlasso1d(data_target_weekly_train)
cv_fl_train = cv.trendfilter(result_fl_train, k=5)
fit_fl_train = coef(result_fl_train, lambda = cv_fl_train$lambda.min)$beta
result_fl = sum((fit_fl_train-data_target_weekly_test)^2)/n_0


# l_0-penalized method  on target training data
fit_l0_train =  cv_l_0_est(data_target_weekly_train)
result_l0 = sum((fit_l0_train-data_target_weekly_test)^2)/n_0



# unisource l_1-penalized method 
n_k = dim(data_auxiliary)[2]
y_1_project = P_n0_n(n_0, n_k) %*% data_auxiliary[9,]
result_fl_T_1_train = fusedlasso1d(y_1_project)
cv_fl_T_1_train = cv.trendfilter(result_fl_T_1_train, k=5)
fit_fl_T_1_train = coef(result_fl_T_1_train, lambda = cv_fl_T_1_train$lambda.min)$beta
result_fl_T_1 = sum((fit_fl_T_1_train -data_target_weekly_test)^2)/n_0



# unisource l_0-penalized method 
fit_l0_T_1_train =  cv_l_0_est(y_1_project)
result_l0_T_1 = sum((fit_l0_T_1_train-data_target_weekly_test)^2)/n_0



# informative source selection
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


# selected multisource l_1-penalized method 
result_fl_T_est_train = fusedlasso1d(y_A_est_project)
cv_fl_T_est_train = cv.trendfilter(result_fl_T_est_train, k=5)
fit_fl_T_est_train = coef(result_fl_T_est_train, lambda = cv_fl_T_est_train$lambda.min)$beta
result_fl_T_est = sum((fit_fl_T_est_train -data_target_weekly_test)^2)/n_0


# selected multisource l_0-penalized method 
fit_l0_T_est_train =  cv_l_0_est(y_A_est_project) 
result_l0_T_est = sum((fit_l0_T_est_train-data_target_weekly_test)^2)/n_0

# all multisource l_1-penalized method 
y_all_project = rep(0, n_0)
for (k in 1:dim(data_auxiliary)[1]){
  y_project = P_n0_n(n_0, n_k) %*% data_auxiliary[k, ]
  y_all_project = y_all_project + y_project
}
y_all_project = y_all_project/dim(data_auxiliary)[1]
result_fl_T_all_train = fusedlasso1d(y_all_project)
cv_fl_T_all_train = cv.trendfilter(result_fl_T_all_train, k=5)
fit_fl_T_all_train = coef(result_fl_T_all_train, lambda = cv_fl_T_all_train$lambda.min)$beta
result_fl_T_all = sum((fit_fl_T_all_train -data_target_weekly_test)^2)/n_0


# all multisource l_0-penalized method 
fit_l0_T_all_train =  cv_l_0_est(y_all_project) 
result_l0_T_all = sum((fit_l0_T_all_train-data_target_weekly_test)^2)/n_0


# Print results
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




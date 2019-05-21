###########################################################################
# Estimating ARMA models
###########################################################################


# the process
arma_proc <- data[,2]


# AR(1) Estimation 
ar1_estimates <- sarima (arma_proc, p=1,d=0,q=0)
ar1_estimates$ttable


# MA(3) Estimation
ma3_estimates <- sarima (arma_proc, p=0,d=0,q=3)
ma3_estimates$ttable



# ARMA(1,1) Estimation combines the two inputs.
arma11_estimates <- sarima (arma_proc, p=1,d=0,q=1)
arma11_estimates$ttable

###########################################################################
# Estimating ARMA models
###########################################################################
P<-4
Q<-4
AIC_MODELS<-matrix(0,P+1,Q+1)
BIC_MODELS<-matrix(0,P+1,Q+1)

# Estimate ARMA models
for (x in 0:P){
  for (y in 0:Q){
    model_fit<- sarima (arma_proc, p=x,d=0,q=y)
    AIC_MODELS[x+1,y+1]<- model_fit$AIC 
    BIC_MODELS[x+1,y+1]<- model_fit$BIC
  }
}

coeffs<-which(BIC_MODELS == min(BIC_MODELS), arr.ind=TRUE)
ar_coeff<-coeffs[[1]]-1
ma_coeff<-coeffs[[2]]-1

final_model <- sarima(arma_proc, p=ar_coeff,d=0,q=ma_coeff)
print (final_model$ttable)

coeffs_A <- which(AIC_MODELS == min(AIC_MODELS),arr.ind = TRUE)
ar_coeffA<-coeffs_A[[1]]-1
ma_coeffA<-coeffs_A[[2]]-1

final_modelA <- sarima(arma_proc, p=ar_coeffA,d=0,q=ma_coeffA)
print (final_modelA$ttable)

# ARMA(1,1) Estimation combines the two inputs.
arma11_estimates <- sarima (arma_proc, p=1,d=0,q=1)
arma11_estimates$ttable


arma_rec_forecast<-vector(mode="double",length=length(arma_proc))

rw_rec_forecasts<-c(arma_proc[1],arma_proc[1:
                                    (length(arma_proc)-1)])

residuals_ARMA <- arma_proc - sim_time_series
residuals_rw <- arma_proc - rw_rec_forecasts


# Ljung-Box Test
# Ljung-Box Test for the random walk model
max_lag <- 20
errors_rw<-arma_proc[2:length(arma_proc)]-
  arma_proc[1:(length(arma_proc)-1)] 

rw_pval<-matrix(0,max_lag,2)
colnames(rw_pval)<-c("pval","lag")
for (i in 1:max_lag){
  rw_pval[i,1]<-Box.test(errors_rw, 
                         i, type = "Ljung-Box")$p.value
  rw_pval[i,2]<-i    
} 
print('Ljung-Box Test (Test stat, p-value)' )
print(rw_pval)

# The LJ test look at multiple correlations simultaneously.  These
# results output the test statistic and the p-value.  
arma_pval<-matrix(0,max_lag,2)
colnames(arma_pval)<-c("pval","lag")
for (i in 1:max_lag){
  arma_pval[i,1]<-Box.test(final_model$fit$residuals, 
                           ar_coeff+ma_coeff+i, 
                           type = "Ljung-Box",
                           fitdf=ar_coeff+ma_coeff)$p.value
  arma_pval[i,2]<-ar_coeff+ma_coeff+i    
} 

print('Lung-Box Test (Test stat, p-value)' )
print(arma_pval)


#Question 6

N <- length(arma_proc); phi <- 0.9494; theta<-0.4243; hold_out<-364
AR_order <- 1; MA_order<- 1; Last_in_sample_obs <- N-hold_out
estimation_sample<- arma_proc[1:Last_in_sample_obs]

# Estimating the model on the observations before the holdout period starts
estimated_model<-arima(estimation_sample,c(1,0,1))

# Predict the time-series over the hold-out period using an ARMA model
prediction<-predict(estimated_model, n.ahead = hold_out)$pred
arma_predictions<-rbind(transpose(estimation_sample),
                        transpose(prediction))

# Predict the time-series over the hold-out period using a random walk model
rw_forecast <- rwf(estimation_sample,h=hold_out)$mean
rw_predictions<-rbind(transpose(estimation_sample),
                      transpose(rw_forecast))

# Plotting the results
plot(arma_proc[(N-hold_out-10+1):N], col="black", type='l',lwd=2)
lines(arma_predictions[(N-hold_out-10+1):N],col="red", type="l",lwd=2)
lines(rw_predictions[(N-hold_out-10+1):N],col="green", type="l",lwd=2)

# Estimating the model on the observations before the holdout period starts
arma_rec_forecast<-vector(mode="double",length=N-Last_in_sample_obs)

for (i in 0:(hold_out-1)){
  temp<-arima(arma_proc[1:(Last_in_sample_obs+i)],c(1,0,1))
  arma_rec_forecast[i+1]<-predict(temp, n.ahead = 1)$pred
  print(i)
}

# Keep realizations and random walk forecasts  
realizations<-arma_proc[(Last_in_sample_obs+1):
                                length(arma_proc)]

rw_rec_forecasts<-arma_proc[(Last_in_sample_obs):
                                    (length(arma_proc)-1)]

# Plot the forecasts against the realizations
plot(realizations,type="l",lwd=2)
lines(arma_rec_forecast,col="red",lwd=2)
lines(rw_rec_forecasts,col="blue",lwd=2)

###########################################################################
###########################################################################
# MZ Regressions
###########################################################################
###########################################################################


###########################################################################
# Computing the tests for the ARMA forecasts
###########################################################################

lm_model <- lm(realizations~arma_rec_forecast)

# Display parameters and p-values
results <- rbind(
  round(summary(lm_model)$coefficients[,1],digits=3),
  round(summary(lm_model)$coefficients[,4],digits=3))
colnames(results)<- c("Alpha","Beta")
rownames(results)<- c("Coeff", "Pval")

beta=
  transpose(as.vector(
    summary(lm_model)$coefficients[,1]))

# Display the F-test
r=2; A=matrix( 
  c(1, 0, 0, 1),
  nrow = 2, 
  ncol = 2)
c =rbind(0,1)


# Here I construct the F-test
Wald_Stat= (transpose(A%*%beta-c) %*%
              solve(A %*% vcov(lm_model) %*% transpose(A)) %*%
              (A%*%beta-c)) 
Pval_Wald_Stat=1-pchisq(Wald_Stat, 2)

print(results)
print(paste('Wald Stat=', round(Wald_Stat,digits=2)))
print(paste('P-value Wald Stat=',round(Pval_Wald_Stat, digits=2)))


###########################################################################
# Computing the tests for the RW forecasts
###########################################################################

lm_model <- lm(realizations~rw_rec_forecasts)

# Display parameters and p-values
results <- rbind(
  round(summary(lm_model)$coefficients[,1],digits=2),
  round(summary(lm_model)$coefficients[,4],digits=2))
colnames(results)<- c("Alpha","Beta")
rownames(results)<- c("Coeff", "Pval")

beta=
  transpose(as.vector(
    summary(lm_model)$coefficients[,1]))

# Display the F-test
r=2; A=matrix( 
  c(1, 0, 0, 1),
  nrow = 2, 
  ncol = 2)
c =rbind(0,1)


# Here I construct the F-test
Wald_Stat= (transpose(A%*%beta-c) %*%
              solve(A %*% vcov(lm_model) %*% transpose(A)) %*%
              (A%*%beta-c)) 
Pval_Wald_Stat=1-pchisq(Wald_Stat, 2)

print(results)
print(paste('Wald Stat=', round(Wald_Stat,digits=2)))
print(paste('P-value Wald Stat=',round(Pval_Wald_Stat, digits=2)))

###########################################################################
# Diebold Mariano tests
# DM tests use the squarred errors and an OLS regression 
# Positive indicate the random walk model is preferred
# Negative indicate the ARMA model is preferred,
###########################################################################

squared_diff=((realizations-arma_rec_forecast)^2 - (realizations-rw_rec_forecasts)^2)
summary(lm(squared_diff~1))
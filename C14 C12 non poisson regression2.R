####################     FUNCTIONS    ########################################
mu_fitting<-function(raw_cts,raw_curr) {
  if (length(raw_cts)<=5) {
    mu_fit<-rep(mean(raw_cts),length(raw_cts))
    ideal_curr_best<-rep(mean(raw_curr),length(raw_curr))
  }  
  if (length(raw_cts)>5) {
    #find secciones,lambda of curr spline, k_scale that minimize AIC of scaled (ideal) cts. Non-Poisson is not included yet.
    #first optimize lambda for spline fit of the current fixed (anchored) to the center (mean cts)
    #ideal cts: normalize the curr. fit and scale the fix shape to cts scale. Find optimum scale for best AIC of cts.
    x<<-c(1:length(raw_cts)) ; xy<<-data.frame(x,raw_cts); xcurr<<-data.frame(x,raw_curr)
    aic_vect<<-c(); ideal_cts_list<<-list(); curr_spline_list<<-list(); lambda_vect<<-c()
    nknots_vect=c(2,3,4)
    
    for (nk in c(1:length(nknots_vect))) {
      nknots_test<<-nknots_vect[nk]; optim_lambda<<-optimise(aic_lambda, lower=5E-4, upper=5E-3, tol=1E-4)
      save_aic_data(nk)
    }

    #best model in function of nknots
    best<-which.min(aic_vect)
    my_aic<-aic_vect[best]; ideal_cts<-ideal_cts_list[[best]]; ideal_curr_best<-curr_spline_list[[best]]
    mu_fit<-ideal_cts

  }
  fit_out<-list(my_aic,mu_fit,ideal_curr_best)
  return(fit_out)
}
save_aic_data<-function(iter){
  aic_vect[iter]<<-optim_lambda$objective; lambda_vect[iter]<<-optim_lambda$minimum
  ideal_cts_list[[iter]]<<-ideal_cts; curr_spline_list[[iter]]<<-as.vector(curr_spline_sm$y)
}
aic_lambda<-function(lambda) {
  #first optimize lambda for spline fit the current
  curr_spline_sm<<-smooth.spline(x,raw_curr,nknots=nknots_test,lambda=lambda) #trend to be gaussian
  curr_spline<-as.vector(curr_spline_sm$y)
  center_cts<<-mean(raw_cts)           ; center_curr<-mean(curr_spline)
  norm_curr_spline<-curr_spline/center_curr
  ideal_cts_unit<<-(1*center_cts*norm_curr_spline)
  #ideal cts: normalize the curr. fit -0.5 to 0.5 and scale the fix shape to cts scale. Find optimum scale for best AIC of cts.
  ideal_cts_model<-optimise(ideal_cts_calc_c,lower=0.02,upper=1.5,tol=0.01)
  k_scale_good<-ideal_cts_model$minimum
  my_aic<-ideal_cts_model$objective
  ideal_cts<<- (k_scale_good*ideal_cts_unit)
  
  return(my_aic)
}

ideal_cts_calc_c<-function(k_scale) {
  ideal_cts<- round(k_scale*ideal_cts_unit)
  my_aic<-my_pois_aic(raw_cts,ideal_cts,4)
  return(my_aic)
}

my_pois_aic<-function(data_y,mu_set,n_par) {
  log_lkh<-my_pois_likelihood(data_y,mu_set)
  aic<-(2*n_par)+(-2*log_lkh)
  return(aic)
}

my_pois_likelihood<-function(data_y,mu_set) {
  log_lkh<-0 
  for (l in c(1:length(data_y))) {
    lkh<-0
    try((lkh<-dpois(data_y[l],lambda=mu_set[l])),silent=T)
    if (lkh>0) {(log_lkh<-log_lkh+log(lkh))}
  }
  return(log_lkh)
}

chi2_transformation<-function(raw_cts,mu_fit) {
  global_avgcts<-mean(raw_cts)
  chi2<-((raw_cts-mu_fit)/sqrt(mu_fit))^2
  E_diff<-sign((raw_cts-mu_fit))*sqrt(chi2)*sqrt(global_avgcts)
  O_transform<-round(global_avgcts + E_diff)
  return(O_transform)
}
##############################################################################################
library(gamlss)
library(bigsplines)
#read data
data=read.csv("C:/Project_results/13C  beam correction/C14_C12_example_2.csv", header=TRUE)
# unique passes labels
label_vec=unique(data$label)
total_passes=length(label_vec)
mu_out=c() ; ideal_curr=c() 
D_qp=0; D_NB=0
cts_unc=data$raw_cts

#iterate thru each pass
for (pass in c(1:length(label_vec))) {
  #read data for corresponding pass
  pass_data=data[which(data$label==label_vec[pass]),]
  raw_curr<<-pass_data$raw_curr
  raw_cts<<-pass_data$raw_cts
  #poisson fitting
  fit_out=mu_fitting(raw_cts, raw_curr)
  #save fitting cts and current data
  mu_out<-c(mu_out,fit_out[[2]]); ideal_curr<-c(ideal_curr,fit_out[[3]]) 
}

D_qp=sqrt(sum(((cts_unc-mu_out)/sqrt(mu_out))^2)/(length(cts_unc)-(total_passes))) # traditional quasi poisson dispersion D
cts_transformation<-chi2_transformation(cts_unc,mu_out) #non stationary to stationary transformation

data_out=data.frame("raw_curr"=data$raw_curr,"raw_cts"=cts_unc,"stationary"=cts_transformation,"mu_out"=mu_out,"ideal_curr"=ideal_curr)

# try to carry out negative binomial fitting
try(NB_fitII<-gamlss(data_out$stationary ~ 1, data=data_out, family = NBII(mu.link = "log", sigma.link = "log"), method=mixed(),  
                     control=gamlss.control(trace=FALSE)),silent=TRUE)
try(D_NB<-sqrt(1+exp(NB_fitII$sigma.coefficients[1])), silent=TRUE) #NB dispersion D

plot(data_out$raw_cts)
lines(data_out$mu_out, type="l", col="red")

plot(data_out$raw_curr)
lines(data_out$ideal_curr, type="l", col="red")

hist(data_out$stationary, main=c("D_qp=",phi_out,"D_NB=",D_NB))

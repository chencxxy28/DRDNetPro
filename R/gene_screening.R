#'@title Estimating equation for ELCIC under GLM
#'@description A specified estimating equation for ELCIC under GLM. This estimating equation is used for marginal mean selection.
#'@usage spearman_screen(index_data,data_vessel,size)
#'@param index_data A matrix containing covariates. The first column should be all ones corresponding to the intercept. See more details in
#'@param data_vessel A vector containing outcomes.
#'@param size A plug-in estimator solved by an external estimating procedure.
#'
#'@return A matrix containing values of calculated estimating equations.
#'@export

spearman_screen<-function(index_data,data_vessel,size)
{
  index<-as.numeric(as.matrix(index_data$rate,ncol=1))
  
  #order the data by the index:
  exp_index<-index[order(index)]
  X<-log(exp_index)
  data_vessel_order<-data_vessel[,order(index)]
  index_data<-index_data[order(index),]
  dim(data_vessel_order)
  Y<-data_vessel_order
  
  ##screening by spearman correlation (seperate screening based on smoking groups)
  #do smoking group 
  x_original<-t(Y)
  data_vessel_order_smoking<-data_vessel_order
  exp_index_smoking<-exp_index
  Y_smoking<-Y
  p<-nrow(data_vessel_order_smoking)
  n<-ncol(data_vessel_order_smoking)
  y<-exp_index
  x<-t(Y_smoking)
  spearman<-rep()
  for (i in 1:ncol(x))
  {
    spearman_i<-abs(cor(x=exp_index_smoking, y=x[,i], method = 'spearman'))
    spearman<-c(spearman,spearman_i)
  }
  size<-size
  sis_smoking<-cbind(spearman,1:p)
  sis_smoking<-sis_smoking[order(sis_smoking[,1],decreasing = TRUE),2][1:size]
  return(final_index=sis_smoking)
}



#'@title Estimating equation for ELCIC under GLM
#'@description A specified estimating equation for ELCIC under GLM. This estimating equation is used for marginal mean selection.
#'@usage test_screen(data_vessel,group,covariate)
#'@param data_vessel A matrix containing covariates. The first column should be all ones corresponding to the intercept. See more details in
#'@param group A vector containing outcomes.
#'@param covariate A plug-in estimator solved by an external estimating procedure.
#'
#'@return A matrix containing values of calculated estimating equations.
#'@export

test_screen<-function(data_vessel,group,covariate){
  Y<-data_vessel
  sig_all<-NULL
  for(i in 1:nrow(Y))
  {
    outcome_i<-as.vector(as.matrix(Y[i,]))
    #gam_fit<-glm(outcome_i~index_data$hyper+index_data$smoke+index_data$sex+as.numeric(index_data$height)+as.numeric(index_data$bmi))
    gam_fit<-glm(outcome_i~group+covariate)
    sum_stat<-summary(gam_fit)
    sig_all_i<-sum_stat$coefficients[2,4]
    sig_all<-c(sig_all,sig_all_i)
  }
  sig_all_adjust<-p.adjust(sig_all,method = "hochberg")
  sig_all_index<-which(sig_all_adjust<=0.05)
  return(sig_all_index)
}

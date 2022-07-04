#'@title Gene screening based on the Spearman correlation.
#'@description This function can be used to screen the genes based on Spearman correlation coefficient.
#'@usage spearman_screen(agent,data_vessel,size)
#'@param agent The imputed disease risk.
#'@param data_vessel The gene expression matrix.
#'@param size The number of selected genes.
#'
#'@return The vector containing the index for the selected genes.
#'@export

spearman_screen<-function(agent,data_vessel,size)
{
  #index<-as.numeric(as.matrix(index_data$rate,ncol=1))
  index<-agent

  #order the data by the index:
  exp_index<-index[order(index)]
  X<-log(exp_index)
  data_vessel_order<-data_vessel[,order(index)]
  #index_data<-index_data[order(index),]
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



#'@title Gene selection based on the test of differential gene expression.
#'@description This function can be used to select the genes based on the test of differential gene expression.
#'@usage test_screen(data_vessel,group,covariate)
#'@param data_vessel The gene expression matrix.
#'@param group The vector including disease status (0 and 1).
#'@param covariate The vector including covariate values (eg., smoking: 0 and 1).
#'@return The vector containing the index for the selected genes.
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

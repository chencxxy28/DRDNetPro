
#'@title Fit a varying coefficient model
#'@description This function fit the gene expression data into a varying coefficient model.
#'@usage vc.fit(agent,data_observe,x_cov)
#'@param data_observe The gene expression matrix.
#'@param agent The imputed disease risk.
#'@param x_cov The vector including covariate values (eg., smoking: 0 and 1).
#'@return A list of two matrices including varying intercept and varying covariate effect, respective.
#'@export
#'@import np splines2 grpreg Matrix pROC
#'@importFrom graphics lines par
#'@importFrom stats coef cor glm p.adjust

vc.fit<-function(agent,data_observe,x_cov)
{
  t<-agent
  #get fitted values of observations using varying coefficient kernel regression
  data_fitted<-rep()
  data_fitted_cov<-rep()
  #bws<-0.01
  for(i in 1:ncol(data_observe))
  {
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=x_cov),ezdat=data.frame(t=t),regtype="ll")
    fitted<-coef(bw)[,1]
    data_fitted<-cbind(data_fitted,fitted)
    fitted_cov<-coef(bw)[,2]
    data_fitted_cov<-cbind(data_fitted_cov,fitted_cov)
  }
  return(list(data_fitted=data_fitted,data_fitted_cov=data_fitted_cov))
}



#'@title The base generation function
#'@description This function generates the base matrices used in the network learning.
#'@usage base.construct(data_observe,data_fitted, degree=3,
#'len.knots=3,data_fitted_cov,agent,x_cov)
#'@param data_observe The gene expression matrix.
#'@param data_fitted The matrix containing varying intercepts.
#'@param degree The degree in the B-spline base. The default is 3
#'@param len.knots The number of knots. The default is 3
#'@param data_fitted_cov The matrix containing varying covariate effects
#'@param agent The imputed disease risk.
#'@param x_cov The vector including covariate values (eg., smoking: 0 and 1).
#'@return The list of four matrices based on varying intercepts, varying covariate effects, varying
#'intercepts with the column containing ones, and varying covariates with the column containing ones.
#'@export
#'@import np splines2 grpreg Matrix pROC
#'@importFrom graphics lines par
#'@importFrom stats coef cor glm p.adjust


base.construct<-function(data_observe=data_observe,
                         data_fitted=data_fitted, degree=3,len.knots=3,
                         data_fitted_cov=data_fitted_cov,
                         agent=agent,
                         x_cov=x_cov)
{
  t=agent
  #some useful fct
  step_fct<-function(x,step)
  {
    t<-length(x)
    x_new<-x[1]
    for (i in 1:(t-1))
    {
      x_int<-seq(x[i],x[i+1],by=step)[-1]
      x_new<-c(x_new,x_int)
    }
    x_new
  }

  #generate B-spline and its integration for baseline
  #generate B-spline and its integration for cov
  X_big<-rep()
  X_big_cov<-rep()
  X_big_int<-rep()
  X_big_int_cov<-rep()
  degree<-degree
  round_num<-attr(regexpr("(?<=\\.)0+", format(min(diff(t,lag=1)),scientific = FALSE), perl = TRUE), "match.length")+1
  t<-round(t,round_num)
  if(min(diff(t,lag=1))>=2*10^(-round_num))
  {
    step<-10^(-round_num)
  }else
  {
    step<-5*10^(-round_num-1)
    t<-round(t,round_num+1)
  }
  cluster<-ncol(data_observe)
  lambda<-seq()
  for (i in 1:cluster)
  {
    knots <- sort(data_fitted[,i])[c(round(nrow(data_fitted)/len.knots),round(nrow(data_fitted)*2/len.knots),round(nrow(data_fitted)*(len.knots-1)/len.knots))]
    knots_cov <-sort(data_fitted_cov[,i])[c(round(nrow(data_fitted)/len.knots),round(nrow(data_fitted)*2/len.knots),round(nrow(data_fitted)*(len.knots-1)/len.knots))]
    bsMat <- bSpline(data_fitted[,i], knots = knots, degree = degree, intercept=F)
    bsMat_cov <- bSpline(data_fitted_cov[,i], knots = knots_cov, degree = degree, intercept=F)
    x_int<-step_fct(as.vector(t),step)
    if(min(diff(t,lag=1))<2*10^(-round_num))
    {x_int<-round(x_int,round_num+1)}else
    {x_int<-round(x_int,round_num)}
    x_cov_int<-rep(0,length(x_int))
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=c(x_cov_int)),ezdat=data.frame(t=c(x_int)),regtype="ll")
    x_hat<-coef(bw)[,1]
    x_hat_cov<-coef(bw)[,2]

    #x_hat<-exp(cbind(1,log(x_int)) %*% beta[,i])
    basis_int<-bSpline(x_hat, knots = knots, degree = degree, intercept=F)
    basis_int_cov<-bSpline(x_hat_cov, knots = knots_cov, degree = degree, intercept=F)
    base_int<-0
    base_int_cov<-0
    for (j in 1:(length(t)-1))
    {
      int_row<-apply(basis_int[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
      base_int<-rbind(base_int,int_row)

      int_row_cov<-apply(basis_int_cov[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
      base_int_cov<-rbind(base_int_cov,int_row_cov)
    }
    base_int<-apply(base_int,2,cumsum)
    base_int_cov<-apply(base_int_cov,2,cumsum)
    #base_int<-ibs(data_fitted[,i], knots = knots, degree = degree, intercept=F)*(max(t)-min(t)) #not correct
    X_big <- cbind(X_big,bsMat)
    X_big_cov <- cbind(X_big_cov,bsMat_cov)

    X_big_int<-cbind(X_big_int,base_int)
    X_big_int_cov<-cbind(X_big_int_cov,base_int_cov)
  }
  return(list(X_big=X_big,
              X_big_cov=X_big_cov,
              X_big_int=X_big_int,
              X_big_int_cov=X_big_int_cov))
}





#'@title Estimating equation for ELCIC under GLM
#'@description A specified estimating equation for ELCIC under GLM. This estimating equation is used for marginal mean selection.
#'@usage network.learn(data_observe,x_cov,X_big_int,X_big_int_cov,
#'agent,degree=3,len.knots=3,cv=TRUE,nfolds=20,alpha=1)
#'@param data_observe A matrix containing covariates. The first column should be all ones corresponding to the intercept. See more details in
#'@param x_cov A plug-in estimator solved by an external estimating procedure.
#'@param X_big_int A vector containing outcomes.
#'@param X_big_int_cov A vector containing outcomes.
#'@param agent A plug-in estimator solved by an external estimating procedure.
#'@param degree A plug-in estimator solved by an external estimating procedure.
#'@param len.knots A plug-in estimator solved by an external estimating procedure.
#'@param cv A plug-in estimator solved by an external estimating procedure.
#'@param nfolds A plug-in estimator solved by an external estimating procedure.
#'@param alpha A plug-in estimator solved by an external estimating procedure.
#'@return A matrix containing values of calculated estimating equations.
#'@export

network.learn<-function(data_observe=data_observe,
                        x_cov=x_cov,
                        X_big_int=X_big_int,
                        X_big_int_cov=X_big_int_cov,
                        agent=agent,
                        degree=3,
                        len.knots=3,
                        cv=TRUE,
                        nfolds=20,
                        alpha=1)
{
  t<-agent
  cluster<-ncol(data_observe)
  #X_big<-apply(X_big,2,function(x) x/sum(x))
  #X_big_int<-apply(X_big_int,2,function(x) x/sum(x))
  num_cov<-degree+len.knots
  #X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov)

  #some setup for final analysis
  X_big_int_exp<-cbind(t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
  X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
  X_big_int_exp_intcep_0<-cbind(1,t,X_big_int,0*X_big_int_cov,t*0,0)
  X_big_int_exp_intcep_1<-cbind(1,t,X_big_int,1*X_big_int_cov,t*1,1)

  x_cov<-matrix(x_cov,ncol=1)
  group<-c(0,rep(1:(2*ncol(data_observe)),each=num_cov),rep(0,ncol(x_cov)),rep(0,ncol(x_cov)))
  self_size<-rep()
  self_size_cov<-rep()
  self_size_all<-rep()
  gene_whole<-rep()
  gene_whole_cov<-rep()
  gene_whole_all<-rep()
  fitted_1_all<-rep()
  fitted_0_all<-rep()
  smoking_effect_all<-rep()
  trend_base_all<-rep()
  trend_smoking_effect_all<-rep()

  for (j in 1:ncol(data_observe))  #
  {
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=5,alpha=alpha1)
    #beta_ini<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha1)$beta
    #group_id<-c(1,2,rep(3:(2+ncol(data_observe)),each=num_cov))
    #beta_each_g<-unlist(tapply(beta_ini,group_id,function(x) x))
    #weight_ini<-rep()
    #for(i in 1:ncol(data_observe))
    #{
    #  m<-sum((beta_each_g[(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)])^2)
    #  weight_ini<-c(weight_ini,m)
    #}
    #weight<-(1/(weight_ini+0.0000001))^0.3
    cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=nfolds,alpha=alpha,max.iter=100000,seed=23556)
    fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda =cvfit$lambda.min,alpha=alpha,max.iter=100000)
    #beta_group<-fit$beta
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda.min = lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda = cvfit$lambda.min)
    #beta_group<-select(fit,"EBIC")$beta
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",nfolds=50)
    #fit<-grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",lambda = cvfit$lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,2], group=group, penalty="grLasso")
    #beta_group<-select(fit,"EBIC")$beta
    #weight

    beta_select<-fit$beta
    index_nonzero<-which(beta_select!=0)
    beta_select[index_nonzero]

    find_index<-index_nonzero[-c(1:2,length(index_nonzero)-1,length(index_nonzero))]-2
    find_index_comb<-unique(ifelse(find_index>num_cov*cluster,find_index-num_cov*cluster,find_index))
    #influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    #m<-which(influ==j)
    #if (length(m)==0)
    #{
    fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_1<-X_big_int_exp_intcep_1[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_0<-X_big_int_exp_intcep_0[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_1_all<-cbind(fitted_1_all,fitted_1)
    fitted_0_all<-cbind(fitted_0_all,fitted_0)
    fitted_self<-X_big_int_exp_intcep[,1:2]%*%beta_select[1:2]
    fitted_self<-t(fitted_self)
    self_size<-rbind(self_size,fitted_self)
    smoking_effect<-(X_big_int_exp_intcep_1[,ncol(X_big_int_exp_intcep_1)]*beta_select[ncol(X_big_int_exp_intcep_1)])[1]
    smoking_effect_all<-c(smoking_effect_all,smoking_effect)
    trend_base<-X_big_int_exp_intcep_0[,2]*beta_select[2]
    trend_base_all<-cbind(trend_base_all,trend_base)
    trend_smoking_effect<-X_big_int_exp_intcep_1[,ncol(X_big_int_exp_intcep_1)-1]*beta_select[ncol(X_big_int_exp_intcep_1)-1]
    trend_smoking_effect_all<-cbind(trend_smoking_effect_all,trend_smoking_effect)

    fitted_self_cov<-X_big_int_exp_intcep_1[,c(length(beta_select)-1,length(beta_select))]%*%beta_select[c(length(beta_select)-1,length(beta_select))]
    fitted_self_cov<-t(fitted_self_cov)
    self_size_cov<-rbind(self_size_cov,fitted_self_cov)

    fitted_self_all<-X_big_int_exp_intcep_1[,c(1,2,length(beta_select)-1,length(beta_select))]%*%beta_select[c(1,2,length(beta_select)-1,length(beta_select))]
    fitted_self_all<-t(fitted_self_all)
    self_size_all<-rbind(self_size_all,fitted_self_all)


    #}
    #else
    #{
    #  fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
    #  fitted_self<-X_big_int_exp_intcep[,c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]%*%beta_select[c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]
    #  fitted_self<-t(fitted_self)
    #  self_size<-rbind(self_size,fitted_self)
    #  find_self<-which(influ==j)
    #  find_index<-find_index[-c((num_cov*(find_self-1)+1):(num_cov*(find_self-1)+num_cov))]
    #}
    par(mfrow = c(1, 1))
    gene_index<-j
    plot(x=t,y=data_observe[,gene_index],ylim=c(-1,max(data_observe[,gene_index])+0.2))
    lines(x=t,y=fitted_self,col="black")
    lines(x=t,y=fitted,col="red")

    num_gene<-length(find_index)/num_cov
    num_gene_comb<-length(find_index_comb)/num_cov
    gene_cluster<-rep()
    gene_cluster_cov<-rep()
    gene_cluster_all<-rep()
    fitted_gene_matrix<-rep()
    fitted_gene_matrix_cov<-rep()
    fitted_gene_matrix_all<-rep()
    if (num_gene>0)
    {
      for (i in 1:num_gene)
      {
        gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
        if (gene_relate<=ncol(data_observe))
        {
          gene_cluster<-c(gene_cluster,gene_relate)
          fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
          fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
          lines(x=t,y=fitted_gene,ylab="y",col=i+2)
        }else
        {
          gene_cluster_cov<-c(gene_cluster_cov,gene_relate)
          fitted_gene_cov<-X_big_int_exp_intcep_1[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
          fitted_gene_matrix_cov<-cbind(fitted_gene_matrix_cov,fitted_gene_cov)
          lines(x=t,y=fitted_gene_cov,ylab="y",col=i+2)
        }
      }
      if(!is.null(fitted_gene_matrix))
      {
        gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
        gene_whole<-rbind(gene_whole,gene_cluster)
      }
      if(!is.null(fitted_gene_matrix_cov))
      {
        gene_cluster_cov<-cbind(gene_cluster_cov,j,t(fitted_gene_matrix_cov))
        gene_whole_cov<-rbind(gene_whole_cov,gene_cluster_cov)
      }
      for (i in 1:num_gene_comb)
      {
        gene_relate_comb<-find_index_comb[(i-1)*num_cov+1]%/%num_cov+1
        gene_cluster_all<-c(gene_cluster_all,gene_relate_comb)
        fitted_gene_all<-X_big_int_exp_intcep_1[,c(find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2,find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2+cluster*num_cov)]%*%beta_select[c(find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2,find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2+cluster*num_cov)]
        fitted_gene_matrix_all<-cbind(fitted_gene_matrix_all,fitted_gene_all)
      }
      gene_cluster_all<-cbind(gene_cluster_all,j,t(fitted_gene_matrix_all))
      gene_whole_all<-rbind(gene_whole_all,gene_cluster_all)
    }
    print(j)

    ###capture the beta_cov and overall

  }

  self_size<-cbind(rep(1:cluster),self_size)
  #dim(self_size)
  self_size_cov<-cbind(rep(1:cluster),self_size_cov)
  #dim(self_size_cov)
  self_size_all<-cbind(rep(1:cluster),self_size_all)
  #dim(self_size_all)
  gene_whole_cov[,1]<-gene_whole_cov[,1]-ncol(data_observe)

  return(list(self_size=self_size,
              self_size_cov=self_size_cov,
              self_size_all=self_size_all,
              gene_whole_cov=gene_whole_cov,
              gene_whole_all=gene_whole_all,
              gene_whole=gene_whole))
}




# fucntions for CDR3 aa probability postanalysis. 

#function to fit Q-offset
get_q<-function(x,total=2e9){
  coef(lm(log10(x$ML)~offset(log10((x$sim_num+1)/(total/3))),singular.ok = F))
}

#accesory function to calculate posterior dist of P_data given occurences vector x_vec and sample sizes vector n_vec
get_P_Data_posterior<-function(x_vec,n_vec,ps=10^-seq(8,0,length.out = 100),quantiles=c(0.025,0.975)){
  likelihood=sapply(ps,function(p)prod(c(1-exp(-n_vec[x_vec]*p)),exp(-n_vec[!x_vec]*p)))
  posteriorP=prop.table(diff(ps)*likelihood[-1])
  cs<-cumsum(posteriorP)
  list(probs=posteriorP,
       cumsum=cs,
       CI=c(ps[which(cs>quantiles[1])[1]],ps[which(cs>quantiles[2])[1]]),
       maxP=ps[which.max(likelihood)])
}

#function for ML estimation of P_data, and confidence interval for P_data
add_p_data_intervals<-function(data,total=2e9,indvec,sizevec){
  left<-numeric(nrow(data))
  right<-numeric(nrow(data))
  ML<-numeric(nrow(data))
  donors<-numeric(nrow(data))
  for (i in 1:nrow(data))
  {
    x_vec<-(data[i,]!=0)[which(indvec)+2] #what??? Not perfect design
    n_vec<-sizevec[indvec]
    posterior=get_P_Data_posterior(x_vec = x_vec,n_vec = n_vec)
    left[i]=posterior$CI[1]
    right[i]=posterior$CI[2]
    ML[i]=posterior$maxP
    donors[i]=sum(x_vec)
  }
  data$left<-left
  data$right<-right
  data$ML<-ML
  data$donors<-donors
  data
}

#function to calculate p-value from P_data posterior distribution.
add_p_value<-function(data,total=2e9,sizevec,indvec){
  data<-data[data$donors>1,]
  q_offset=get_q(data,total=total)
  data$P_post=(data$sim_num*3/total)*10^q_offset
  ps=10^-seq(8,0,length.out = 100)
  pval<-numeric(nrow(data))
  effect_size<-numeric(nrow(data))
  for (i in 1:nrow(data))
  {
    x_vec<-(data[i,]!=0)[which(indvec)+2]
    n_vec<-sizevec[indvec]
    posterior=get_P_Data_posterior(x_vec = x_vec,n_vec = n_vec,ps = ps)
    pval[i]=sum(posterior$cumsum[which(ps>data$P_post[i])[1] ])
    effect_size[i]<-log10(data$ML[i])-log10(data$P_post[i])
  }
  data$pval_post=pval
  data$effect_size<-effect_size
  data
}
  
#main function for automated analysis. Data is dataset with sim_num and CDR3 amino acid sequence in first two columns.
#Next columns indicate presence (or abundance in number of reads) of CDR3 in each sample. 
#sizevec - is vector of samplesize, indvec is logical vector indicating, if sample corresponds to patient cohort. 
do_analysis<-function(data,total,sizevec,indvec){
  tmp<-add_p_data_intervals(data,total=total,sizevec=sizevec, indvec=indvec )
  add_p_value(tmp,total=total,sizevec=sizevec, indvec=indvec)
}
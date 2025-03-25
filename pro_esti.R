library(mclust)
library(cluster)
#library('ggplot2')
#library(scatterplot3d)
library(MASS)
library(network)

source('models.R')

###############################################################################
###############################################################################
###############################################################################

#table1~3

{#init
  K_u=3;K_v=3 #ture communities
  n1=200
  n2=200 #nodes
  T_=15  #times
  #b=seq(0,8,1) #autoterm
  #b=0
  rho_0=seq(5,30,5)  #strength of signal
  
  r=3 #dim of embedding
  
  BB=100
}

pro_est_error<-function(prob_list,u0,v0_list,u0_list,v0){
  result<-0
  for (t in 1:T_) {
    result=result+sum((u0%*%t(v0_list[[t]])+u0_list[[t]]%*%t(v0) - 2*prob_list[[t]])^2)/(4*sum(prob_list[[t]]^2))
  }
  return(result/T_)
}

for (b in c(0,0.5,1)) {
  for (rrho in 1:length(rho_0) ) {
    
    PPmat<-matrix(c(0.12, 0.03, 0.05,
                    0.02, 0.09, 0.05,
                    0.03, 0.02, 0.12), nrow = K_u, ncol = K_v)*rho_0[rrho]*(log(n1)+log(n2))/ (2*sqrt(n1*n2))
    PPmat_last<-matrix(c(0.05, 0.14, 0.12,
                         0.02, 0.09, 0.05,
                         0.14, 0.15, 0.05), nrow = K_u, ncol = K_v)*rho_0[rrho]*(log(n1)+log(n2))/ (2*sqrt(n1*n2))
    
    epsi=0.00008 #converge error
    lambda=(10^-12) #penalty
    
    est_error<-numeric(BB)
    ture_b<-numeric(BB)
    train_b<-numeric(BB)
    
    for (iter in 1:BB) {
      data_result<-get_A_meb_list(n1,n2,T_, b, rho_0[rrho],K_u,K_v)
      A_list=data_result$A_list
      membership_u_list=data_result$membership_u_list
      membership_v_list=data_result$membership_v_list
      
      {
        P_0_list<-list()
        
        prob_list<-list()
        for (t in 1:(T_)) {
          if (t%%2 == 1) {
            P_mat<-matrix(0,n,n)
            for (i in 1:n) {
              for (j in 1:n) {
                P_mat[i,j]<-logit(PPmat[membership_u_list[[t]][i] ,membership_v_list[[t]][j] ])
              }
            }
            prob_list[[t]]<-P_mat
          }
          if (t%%2 == 0) {
            P_mat<-matrix(0,n,n)
            for (i in 1:n) {
              for (j in 1:n) {
                P_mat[i,j]<-logit(PPmat_last[membership_u_list[[t]][i] ,membership_v_list[[t]][j] ])
              }
            }
            prob_list[[t]]<-P_mat
          }
        }
      }
      
      result=stable_embedding(n1,n2,T_,A_list,K_u,K_v,r,epsi,lambda)
      
      b0=result$b0
      u0=result$u0
      v0=result$v0
      v0_list=result$v0_list
      u0_list=result$u0_list
      
      ture_b[iter]<-b
      train_b[iter]<-b0
      est_error[iter]<-pro_est_error(prob_list,u0,v0_list,u0_list,v0)
    }
    final_result<-data.frame(ture_b,train_b,est_error)
    write.table(final_result,paste("final_pro_esti_result_n1=",n1,"_n2=",n2,"_T=",T_,"_b=",b,
                                   "_rho_0=",rho_0[rrho],".csv",sep = ''),row.names=FALSE,col.names=TRUE,sep=",")
    
  }
}














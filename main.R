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



{#init
  K_u=3;K_v=3 #ture communities
  n1=200
  n2=200 #nodes
  
  T_=15 #times
  b=0
  rho_0=seq(8,15,1)  #strength of signal
  
  r=3 #dim of embedding
  
  BB=100
}


for (rrho in 1:length(rho_0) ) {
  epsi=0.00008 #converge error
  lambda=(10^-12)#penalty
  
  
  ture_b<-numeric(BB)
  train_b<-numeric(BB)
  Embedding_SoS<-numeric(BB)
  Unbias_SoS<-numeric(BB)
  Uase<-numeric(BB)
  Sum_spectral<-numeric(BB)
  
  for (iter in 1:BB) {
    data_result<-get_A_meb_list(n1,n2,T_, b, rho_0[rrho],K_u,K_v)
    A_list=data_result$A_list
    membership_u_list=data_result$membership_u_list
    membership_v_list=data_result$membership_v_list
    
    result=stable_embedding(n1,n2,T_,A_list,K_u,K_v,r,epsi,lambda)
    
    b0=result$b0
    u0=result$u0
    v0=result$v0
    v0_list=result$v0_list
    u0_list=result$u0_list
    
    ARI_embedding_SoS<-embedding_SoS(u0,v0_list,u0_list,v0,membership_u_list,membership_v_list,K_u,K_v,r=3)
    ARI_unbias_SoS<-unbias_SoS(A_list,membership_u_list,membership_v_list,K_u,K_v,r=3)
    ARI_UASE<-UASE(A_list,membership_u_list,membership_v_list)
    ARI_sum_spectral<-sum_spectral(A_list,membership_u_list,membership_v_list)
    
    ture_b[iter]<-b
    train_b[iter]<-b0
    Embedding_SoS[iter]<-(ARI_embedding_SoS$acc_Z+ARI_embedding_SoS$acc_Y)/2
    Unbias_SoS[iter]<-(ARI_unbias_SoS$acc_Z+ARI_unbias_SoS$acc_Y)/2
    Uase[iter]<-(ARI_UASE$acc_Z+ARI_UASE$acc_Y)/2
    Sum_spectral[iter]<-(ARI_sum_spectral$acc_Z+ARI_sum_spectral$acc_Y)/2
  }
  final_result<-data.frame(ture_b,train_b,Embedding_SoS,Unbias_SoS,Uase,Sum_spectral)
  write.table(final_result,paste("final_result_n1=",n1,"_n2=",n2,"_T=",T_,"_b=",b,
                                 "_rho_0=",rho_0[rrho],".csv",sep = ''),row.names=FALSE,col.names=TRUE,sep=",")
  
}

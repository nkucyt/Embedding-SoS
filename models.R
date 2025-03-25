
library(mclust)
library(cluster)
library('ggplot2')
library(scatterplot3d)
library(MASS)
library(network)

{
  Probmat_new<-function(u,v,b,A_last){
    
    ee<-10^-50
    
    result<-u%*%t(v)+b*(A_last)
    
    result<-1/((1+exp(-result)))
    result[result>=1]=1-ee
    result[result<=0]=0+ee
    
    return(result)
  }
  
  negloglikelihood <- function(A_list,p_list, u, v_list, b, lambda) {
    TT=length(A_list)
    n1=dim(A_list[[1]])[1]
    n2=dim(A_list[[1]])[2]
    
    term1=matrix(0,n1,n2)
    term2=0
    
    term1<-A_list[[1]]*( u%*%t(v_list[[1]]) )-log(1+exp(u%*%t(v_list[[1]]) ) )
    for (t in 2:TT) {
      temp1<-u%*%t(v_list[[t]])+b*(A_list[[t-1]])
      temp2<-log( 1+exp(u%*%t(v_list[[t]])+b*(A_list[[t-1]]) ) )

      term1<-term1+( A_list[[t]]*temp1-temp2 )
    }
    for (t in 1:TT) {
      term2<-term2+sum(v_list[[t]]^2)
    }
    term3=lambda*sum(u^2)
   
    return(-sum(term1)/(n1*n2*TT) + term3+lambda*term2/TT)
  }
  
  partiallik_u <- function(s, u, v_list, A_list,p_list, lambda) { #r表示截距项embedding的维数
    TT<-length(A_list)
    n2<-dim(A_list[[1]])[2]
    term1=0
    for (t in 1:TT) {
      term1<-term1+(A_list[[t]][s,] -p_list[[t]][s,] )%*%v_list[[t]]
    }
    
    result<-(-term1)/(n2*TT) + 2*lambda*u[s,] 
    return(result)
  }
  
  
  partiallik_v_list <- function(s, u, v=v_list[[t]], A=A_list[[t]], p=p_list[[t]], lambda){
    TT<-length(A_list)
    n1<-dim(A)[1]
    term1<-((A[,s]-p[,s]) )%*%u
    result<-(-term1/(n1) + 2*lambda*v[s,]/TT )
    return(result)
  }
  
  partiallik_b <- function(u, v_list, A_list,p_list){
    TT<-length(A_list)
    n1=dim(A_list[[1]])[1]
    n2=dim(A_list[[1]])[2]
    term1=0
    for (t in 2:TT) {
      term1<-term1+sum( ((A_list[[t]] - p_list[[t]] ))*(A_list[[t-1]]) )
    }
    result <- (-term1)
    return(result)
  }
  
  partiallik_b_2 <- function(A_list,p_list){
    TT<-length(A_list)
    term1=0
    for (t in 2:TT) {
      term1<-term1+sum( (  p_list[[t]]*(1-p_list[[t]]) )*(A_list[[t-1]]) )
    }
    
    result <- (term1)/(2*n1*n2*TT) 
    return(result)
  }

}




{
  sim_e<-function(prob_K,n){
    cum_prob<-cumsum(prob_K)
    
    u<-runif(n)
    elist<-matrix(0, 1, n)
    for (i in 1:n) {
      elist[1,i]<-rank(c(u[i],cum_prob))[1]
    }
    return(elist)
  }
  
  rho<-function(n, lambda){
    return(lambda/n)
  }
  
  get_membership<-function(n,pi_1,k){
    pi_list<-rep((1-pi_1)/(k-1), k)
    pi_list[1]<-pi_1
    
    c<-sim_e(pi_list,n)
    return(c)
  }
  
}


{
  get_adjmat_dynamic<-function(n1,n2,Pmat){
    A<-matrix(0,n1,n2)
    
    for (i in 1:n1 ) {
      for (j in 1:n2) {
        A[i,j]<-rbinom(1,1,min(max(Pmat[i,j],0),1))
      }
    }
    return(A)
  }
  
  # generate_net_dynamic<-function(u, v, n, lambda, b, A_last, K_u,K_v){ #生成网络数据
  #   P_0<-Probmat_new(u,v,b,A_last)*lambda/n
  #   
  #   A<-get_adjmat_dynamic(n,P_0)
  #   A <- as.network(A,matrix.type = 'adjacency',directed = TRUE)
  #   
  #   cl_u<-kmeans(u,K_u)
  #   cl_v<-kmeans(v,K_v)
  #   
  #   #center_u<-t(cl_u$center)
  #   uc<-cl_u$cluster
  #   #center_v<-t(cl_v$center)
  #   vc<-cl_v$cluster
  #   
  #   result<-list()
  #   result$net<-A
  #   result$member_u<- c(uc)
  #   result$member_v<- c(vc)
  #   return(result)
  # }
  
  logit<-function(p){ log(p/(1-p)) }
  antilogit<-function(x) {1/(1+exp(-x))}
  
  generate_net_dynamic_by_mat<-function(P_0, n1,n2, b, A_last, membership_u,membership_v){ #generate networks data
    P_mat<-matrix(0,n1,n2)
    for (i in 1:n1) {
      for (j in 1:n2) {
        P_mat[i,j]<-antilogit( logit(P_0[membership_u[i] ,membership_v[j] ]) + b*(A_last[i,j]) )
      }
    }
    A<-get_adjmat_dynamic(n1,n2, P_mat)
    
    A <- as.network(A , matrix.type = 'adjacency', directed = TRUE ,bipartite = TRUE)
    
    result<-list()
    result$net<-A
    result$member_u<- membership_u
    result$member_v<- membership_v
    return(result)
  }
  
}


################################################################################
################################################################################
################################################################################

get_A_meb_list<-function(n1,n2,T_, b, rho_0,K_u,K_v){ #parematers setting
  r=min(K_u,K_v)
  
  membership_u_list<-list()
  membership_v_list<-list()
  membership_u<-get_membership(n1,1/K_u,K_u)
  membership_v<-get_membership(n2,1/K_v,K_v)
  
  membership_u_list[[1]]<-membership_u
  membership_v_list[[1]]<-membership_v 
  
  gam=0.0
  transition_u<-matrix(gam/(K_u-1),K_u,K_u )
  diag(transition_u)<-1-gam
  transition_v<-matrix(gam/(K_v-1),K_v,K_v )
  diag(transition_v)<-1-gam
  
  for (t in 1:(T_-1)) {
    for (i in 1:n1) {
      membership_u[i]=sample(1:K_u,size=1,prob=transition_u[membership_u[i],])
    }
    for (j in 1:n2) {
      membership_v[j]=sample(1:K_v,size=1,prob=transition_v[membership_v[j],]) 
    }
    
    membership_u_list[[t+1]]<-membership_u
    membership_v_list[[t+1]]<-membership_v
    
  }

  A_list<-list()
  P_0_list<-list()

  
  PPmat<-matrix(c(0.12, 0.03, 0.05,
                  0.02, 0.09, 0.05,
                  0.03, 0.02, 0.12), nrow = K_u, ncol = K_v)*rho_0*(log(n1)+log(n2))/(2*sqrt(n1*n2))
  P_0_list[[1]]<-PPmat

  AA<-generate_net_dynamic_by_mat(P_0_list[[1]], n1,n2, b, A_last=matrix(0, nrow = n1, ncol = n2)
                                  , membership_u_list[[1]],membership_v_list[[1]])
  A <-as.matrix(AA$net)
  A_list[[1]]<-A
  
  PPmat_last<-matrix(c(0.05, 0.14, 0.12,
                       0.02, 0.09, 0.05,
                       0.14, 0.15, 0.05), nrow = K_u, ncol = K_v)*rho_0*(log(n1)+log(n2))/(2*sqrt(n1*n2))
  
  for (t in 2:(T_)) {
    if (t%%2 == 1) {
      P_0_list[[t]]<-PPmat
    }
    if (t%%2 == 0) {
      P_0_list[[t]]<-PPmat_last
    }
    
    AA<-generate_net_dynamic_by_mat(P_0_list[[t]], n1,n2, b, A_list[[t-1]]
                                    , membership_u_list[[t]],membership_v_list[[t]])
    A <-as.matrix(AA$net)
    A_list[[t]]<-A
    
  }
  result<-list()
  result<-list('A_list'=A_list,'membership_u_list'=membership_u_list,'membership_v_list'=membership_v_list)
  
  return(result)
}

################################################################################
################################################################################
################################################################################

#plot(as.network(A))

stable_embedding<-function(n1,n2,T_,A_list,K_u,K_v,r,epsi,lambda){
  {
    u0_list<-list()
    v0_list<-list()

    {
    for (t in 1:T_) {
      v0_list[[t]]<-svd(A_list[[t]])$v[,1:r]
      u0_list[[t]]<-svd(A_list[[t]])$u[,1:r]
    }

    term1=0
    for (t in 1:T_) {
      term1<-term1+(A_list[[t]] )%*%t(A_list[[t]] )
    }
    u0<-svd(term1)$u[,1:r]

    term1=0
    for (t in 1:T_) {
      term1<-term1+t(A_list[[t]] )%*%(A_list[[t]] )
    }
    v0<-svd(term1)$v[,1:r]
    }
      
    b0=0
    
    #######################################################interitive1
    b1<-0
    u1<-matrix(0,nrow=n1,ncol = r) 
    v1<-matrix(0,nrow=n2,ncol = r) 
    
    v1_list<-list()
    u1_list<-list()
    
    for (t in 1:T_) {
      v1_list[[t]]<-matrix(0,nrow=n2,ncol = r)
      u1_list[[t]]<-matrix(0,nrow=n1,ncol = r)
      
    }
    
    A_list_trans<-list()
    for (t in 1:T_) {
      A_list_trans[[t]]<-t(A_list[[t]])
    }
    
    #u1_list=u0_list
    p0_list_u<-list()
    p0_list_u[[1]]<-Probmat_new(u0, v0_list[[1]], b0, matrix(0, nrow = n1, ncol = n2) )
    for (t in 2:T_) {
      p0_list_u[[t]]<-Probmat_new(u0, v0_list[[t]], b0, A_list[[t-1]] )
    }
    p1_list_u<-list()
    for (t in 1:T_) {
      p1_list_u[[t]]<-matrix(0,nrow = n1, ncol = n2)
    }
    
    p0_list_v<-list()
    p0_list_v[[1]]<-Probmat_new(v0,u0_list[[1]], b0, matrix(0, nrow = n2, ncol = n1) )
    for (t in 2:T_) {
      p0_list_v[[t]]<-Probmat_new(v0,u0_list[[t]], b0, A_list_trans[[t-1]] )
    }
    p1_list_v<-list()
    for (t in 1:T_) {
      p1_list_v[[t]]<-matrix(0,nrow = n2, ncol = n1)
    }
    
    
    
    h=1
    st=50*sqrt(n1*n2) #learning rate
  }
  
  
  {
    
    u0=(u0-mean(u0))
    v0=(v0-mean(v0))
    
    for (t in 1:T_) {
      u0_list[[t]]=(u0_list[[t]]-mean(u0_list[[t]]))
    }
    for (t in 1:T_) {
      v0_list[[t]]=(v0_list[[t]]-mean(v0_list[[t]]))
    }
    
    
    
    l0_u<-negloglikelihood(A_list,p0_list_u, u0, v0_list, b0 ,lambda)
    
    l0_v<-negloglikelihood(A_list_trans,p0_list_v, v0, u0_list, b0 ,lambda)
    
    #(l0_u+l0_v)/2 
  }
  
  ################################################################################
  ################################################################################
  ################################################################################
  
  eta=1
  cumst<-0
  while (TRUE) { #
    
    for (s in 1:n1) { #update u
      delta_u=partiallik_u(s, u0, v0_list, A_list, p0_list_u, lambda)
      u1[s,]<-u0[s,]-st*delta_u/2
    }
    for (t in 1:T_) {
      for (s in 1:n2) { #update v
        delta_vlist=partiallik_v_list(s,u0, v0_list[[t]], A_list[[t]], p0_list_u[[t]], lambda)
        v1_list[[t]][s,]<-v0_list[[t]][s,]-st*delta_vlist/2
      }
    }
    
    
    for (s in 1:n2) { #update u
      delta_v=partiallik_u(s, v0, u0_list, A_list_trans, p0_list_v, lambda)
      v1[s,]<-v0[s,]-st*delta_v/2 
    }
    for (t in 1:T_) {
      for (s in 1:n1) { #update v
        delta_ulist=partiallik_v_list(s,v0,u0_list[[t]], A_list_trans[[t]], p0_list_v[[t]], lambda)
        u1_list[[t]][s,]<-u0_list[[t]][s,]-st*delta_ulist/2
      }
    }
    
    delta_b=(partiallik_b( u0, v0_list, A_list, p0_list_u)+partiallik_b( v0, u0_list, A_list_trans, p0_list_v))/2
    cumst<-max(0.99*cumst+(1-0.99)*((delta_b)^2), cumst)
    b1<-b0-eta*delta_b/(sqrt(cumst ))

    p1_list_u[[1]]<-Probmat_new(u1, v1_list[[1]], b1, matrix(0, nrow = n1, ncol = n2) )
    for (t in 2:T_) {
      p1_list_u[[t]]<-Probmat_new(u1,v1_list[[t]], b1, A_list[[t-1]])
    }
    
    p1_list_v[[1]]<-Probmat_new(v1, u1_list[[1]], b1, matrix(0, nrow = n2, ncol = n1) )
    for (t in 2:T_) {
      p1_list_v[[t]]<-Probmat_new(v1,u1_list[[t]], b1, A_list_trans[[t-1]])
    }
    
    
    
    l1_u<-negloglikelihood(A_list, p1_list_u, u1, v1_list, b1, lambda)
    l1_v<-negloglikelihood(A_list_trans, p1_list_v, v1, u1_list, b1, lambda) 
    error=abs((l0_u+l0_v)/2 - (l1_u+l1_v)/2 ) / ((l0_u+l0_v)/2)
    
    if (h %% 10 == 0){
      print(c(l1_u,l1_v,error,h,st))
    }
    if (error<epsi & h>300) {
      break
    }
    if (((l1_u+l1_v)/2>(l0_u+l0_v)/2)) {
      st=st/2
      h=h+1
    }
    else{
      #W0=W1
      #st=st*1.2
      
      u0=u1
      v0=v1
      
      v0_list=v1_list
      u0_list=u1_list
      
      b0=b1
      
      p0_list_u=p1_list_u
      p0_list_v=p1_list_v
      
      
      l0_u=l1_u
      l0_v=l1_v
      
      h=h+1
    }
  }
  result<-list('b0'=b0,'u0'=u0,'v0'=v0,'v0_list'=v0_list,'u0_list'=u0_list)
  
  
  return(result)
}

########################################################################################
####################################################################################################
####################################################################################################
# K_u=3;K_v=3
unbias_SoS<-function(A_list,membership_u_list,membership_v_list,K_u,K_v,r=3){
  term1=0
  term2=0
  T_=length(A_list)
  
  for (t in 1:T_) {
    term1<-term1+(A_list[[t]] )%*%t(A_list[[t]] )
  }
  for (t in 1:T_) {
    term2<-term2+t(A_list[[t]] )%*%(A_list[[t]] )
  }
  
  clZ<-spectral_cluster_svd(term1-diag(diag(term1)),r,K_u,K_v)
  clY<-spectral_cluster_svd(term2-diag(diag(term2)),r,K_u,K_v)
  
  ave_ari_Z=0
  ave_ari_Y=0
  
    for (t in 1:T_) {
      ave_ari_Z=ave_ari_Z+ adjustedRandIndex(clZ$clustering_u ,c(membership_u_list[[t]]))
      ave_ari_Y=ave_ari_Y+ adjustedRandIndex(clY$clustering_v ,c(membership_v_list[[t]]))
    }
  result<-list()
  result$ari_Z<-ave_ari_Z/T_
  result$ari_Y<-ave_ari_Y/T_
  return(result)
}

embedding_SoS<-function(u0, v0_list,u0_list,v0,membership_u_list,membership_v_list,K_u,K_v,r=3){
  term1=0
  term2=0
  
  T_=length(v0_list)
  
  for (t in 1:T_) {
    term1<-term1+t(u0%*%t(v0_list[[t]]))%*%(u0%*%t(v0_list[[t]])) #t(u1_list[[t]]%*%t(v1))%*%(u1_list[[t]]%*%t(v1))
    term2<-term2+(u0_list[[t]]%*%t(v0))%*%t(u0_list[[t]]%*%t(v0)) #t(u1_list[[t]]%*%t(v1))%*%(u1_list[[t]]%*%t(v1))
  }
  clY<-spectral_cluster_svd(term1,r,K_u,K_v)
  clZ<-spectral_cluster_svd(term2,r,K_u,K_v)
  
  ave_ari_Y=0
  ave_ari_Z=0
  for (t in 1:T_) {
    ave_ari_Y=ave_ari_Y+adjustedRandIndex(clY$clustering_v,c(membership_v_list[[t]]))
    ave_ari_Z=ave_ari_Z+adjustedRandIndex(clZ$clustering_u,c(membership_u_list[[t]]))
    
  }
  result<-list()
  result$ari_Z<-ave_ari_Z/T_
  result$ari_Y<-ave_ari_Y/T_
  return(result)
}

oneside_embedding_SoS<-function(u0, v0_list,u0_list,v0,membership_u_list,membership_v_list,K_u,K_v,r=3){
  T_=length(v0_list)
  ave_ari_Y=0
  ave_ari_Z=0
  for (t in 1:T_) {
    cl_u<-pam(u0_list[[t]],K_u)
    cl_v <- pam(v0_list[[t]],K_v)
    
    ave_ari_Y=ave_ari_Y+adjustedRandIndex(cl_v$clustering,c(membership_v_list[[t]]))
    ave_ari_Z=ave_ari_Z+adjustedRandIndex(cl_u$clustering,c(membership_u_list[[t]]))
    
  }
  result<-list()
  result$ari_Z<-ave_ari_Z/T_
  result$ari_Y<-ave_ari_Y/T_
  return(result)
}


UASE<-function(A_list,membership_u_list,membership_v_list){
  #A_list<-list()
  A_unfold_u<-c()
  A_unfold_v<-c()
  T_=length(A_list)
  
  for (t in 1:T_) {
    A_unfold_v<-rbind(A_unfold_v,A_list[[t]])
    A_unfold_u<-cbind(A_unfold_u,A_list[[t]])
  }
  V<-svd(A_unfold_v)$v[,1:r]
  U<-svd(A_unfold_u)$u[,1:r]
  
  cl_u<-pam(U,K_u)
  cl_v <- pam(V,K_v)
  
  ave_ari_Z=0
  for (t in 1:T_) {
    ave_ari_Z=ave_ari_Z+adjustedRandIndex(cl_u$clustering ,c(membership_u_list[[t]]))
  }
  ave_ari_Y=0
  for (t in 1:T_) {
    ave_ari_Y=ave_ari_Y+adjustedRandIndex(cl_v$clustering ,c(membership_v_list[[t]]))
  }
  result<-list()
  result$ari_Z<-ave_ari_Z/T_
  result$ari_Y<-ave_ari_Y/T_
  return(result)
}


sum_spectral<-function(A_list,membership_u_list,membership_v_list){
  #A_list<-list()
  T_=length(A_list)
  A_sum<-matrix(0,dim(A_list[[1]])[1],dim(A_list[[1]])[2])
  for (t in 1:T_) {
    A_sum<-A_sum+A_list[[t]]
  }
  V<-svd(A_sum)$v[,1:r]
  U<-svd(A_sum)$u[,1:r]
  
  cl_u<-pam(U,K_u)
  cl_v <- pam(V,K_v)
  
  ave_ari_Z=0
  for (t in 1:T_) {
    ave_ari_Z=ave_ari_Z+adjustedRandIndex(cl_u$clustering ,c(membership_u_list[[t]]))
  }
  ave_ari_Y=0
  for (t in 1:T_) {
    ave_ari_Y=ave_ari_Y+adjustedRandIndex(cl_v$clustering ,c(membership_v_list[[t]]))
  }
  result<-list()
  result$ari_Z<-ave_ari_Z/T_
  result$ari_Y<-ave_ari_Y/T_
  return(result)
}


spectral_cluster_svd<-function(A,r,K_u,K_v){ 
  #r=min(K_u,K_v)
  
  rn=dim(A)[1]
  cn=dim(A)[2]

  L_n<-A
  
  UDV<-svd(L_n)
  #U<-UDU$vectors[,1:K]
  U<-UDV$u[,1:r]
  V<-UDV$v[,1:r]
  
  
  cl_u<-pam(U,K_u)
  cl_v <- pam(V,K_v)
  
  result<-list()
  result$clustering_u<-cl_u$clustering
  result$clustering_v<-cl_v$clustering
  
  
  return(result)
}





# longitudinal envelope with dennis
# Jan14th 2021
# Lan Liu
rm(list=ls())
library(Matrix)
library(xtable)

# Function for calculating the product of matrix then take summation in calculation of estimate for alpha
prod_fun=function(Z,X,Y,Sigma_inv){
    for(i in 1:n){
        W_i=kronecker(t(Z[,i]),X)
        if(i==1){
            num=t(W_i)%*%Sigma_inv%*%Y[,i]
            denom=t(W_i)%*%Sigma_inv%*%W_i
        }else{
            num=num+t(W_i)%*%Sigma_inv%*%Y[,i]
            denom=denom+t(W_i)%*%Sigma_inv%*%W_i
        }
        if(i==n){
            return(solve(denom)%*%num)
        }
        
    }    
}


prod2_fun=function(i,Z,X,Y,Sigma_inv){
 #   for(i in 1:n){
        W_i=kronecker(t(Z[,i]),X)
 #       if(i==1){
            num=t(W_i)%*%Sigma_inv%*%Y[,i]
            denom=t(W_i)%*%Sigma_inv%*%W_i
 #       }else{
  #          num=num+t(W_i)%*%Sigma_inv%*%Y[,i]
  #          denom=denom+t(W_i)%*%Sigma_inv%*%W_i
  #      }
  #      if(i==n){
  #          return(solve(denom)%*%num)
 #       }
        
 #   }    
 return(c(as.vector(num),as.vector(denom)))
}


#######################################
# Simulation to compare the variance  #
#######################################
if_env_win=1
# if_env_win=0

library('Renvlp')
library('MASS')
library('Matrix')
set.seed(1)
r=20
p=8
u=6 
 q=15*if_env_win+6*(1-if_env_win) # q=15 env win q=6 tilde beta win
# q=u # optimal choice
 # q1 is the dimension of span(X) intersect with span(Gamma), #q1 should be <= min(colrank(Gamma), colrank(X))=min(u,q)
q1=4 
q2= q-q1 # q2 is the dimension of span(X) intersect with span(Gamma0)
n=200
set.seed(10)
#lambda=50
#lambda0=20
# lambda=0.5
# lambda0=50
#Omega=diag(u)*lambda
#Omega_0=diag(r-u)*lambda0
#Omega=lambda*diag(c(rep(1, u-q1),rep(1,q1)))
#Omega_0=lambda0*diag(c(rep(1,q2),rep(1,r-u-q2)))
# env win
if(if_env_win){
#Omega=diag(rnorm(u,0,5)^2+40)
 #Omega=diag(c(rep(.5, u-q1),rep(0.3,q1)))
  Omega=diag(c(rep(.5, u-q1),rep(1.5,q1)))
# Omega=diag(c(0.5*runif(u-q1,1,2),1.5*runif(q1,1,2)))
 Omega_0=diag(c(rep(50,q2),rep(50,r-u-q2)))
# Omega_0=diag(c(rep(50,q2),rep(0.5,r-u-q2)))
} else{ 
# tilde beta win
 Omega=diag(c(rep(50, u-q1),rep(1.5,q1)))
 Omega_0=diag(c(rep(0.5,q2),rep(0.5,r-u-q2)))
}

temp.m=matrix(rnorm(r*r),nrow=r)
#temp.m_normalized=temp.m%*%solve(t(temp.m)%*%temp.m)%*%t(temp.m)
temp.m_normalized=svd(temp.m)$u
Gamma=temp.m_normalized[,1:u]
Gamma_0=temp.m_normalized[,-(1:u)]
Sigma=Gamma%*%Omega%*%t(Gamma)+Gamma_0%*%Omega_0%*%t(Gamma_0)
Sigma_inv=solve(Sigma)

eigen(Gamma[,(u-q1+1):u]%*%t(Gamma[,(u-q1+1):u]))$value
eigen(Gamma_0[,1:q2]%*%t(Gamma_0[,1:q2]))$value


#phi=rbind(matrix(0,nrow=u-q1,ncol=q),diag(q),matrix(0,nrow=r-u-q2,ncol=q)) 
# eta=matrix(rnorm(u*p),nrow=u)
temp_eta1.m=matrix(rnorm(q1*u),nrow=u)
temp_eta2.m=matrix(rnorm(q1*p),nrow=q1)
#eta=rbind(matrix(0,nrow=u-q1,ncol=p),temp_eta2.m) # gaurantee rank(eta) is no greater than q1, thus the rank(beta) is no greater than q1
eta=temp_eta1.m%*%temp_eta2.m # gaurantee rank(eta) is no greater than q1, thus the rank(beta) is no greater than q1
beta=Gamma%*%eta
temp_sigma_Z.m=matrix(rnorm(p^2),nrow=p)
Sigma_Z=temp_sigma_Z.m%*%t(temp_sigma_Z.m)
Sigma_Z=diag(log(c(1:p))+1)%*%(diag(p)*0.15+0.85)%*%diag(log(c(1:p))+1)
cov2cor(Sigma_Z)
#Sigma_Z=diag(p)
Sigma_Z_inv=solve(Sigma_Z)
Z=t(chol(Sigma_Z))%*%matrix(rnorm(n*p),nrow=p) # Z should be mean 0 as when we estimate Sigma_Z we use Z%*%t(Z)/n without subtracting the mean
phi=as.matrix(bdiag(temp_eta1.m,rbind(diag(q2),matrix(0,nrow=r-u-q2,ncol=q2)))) 
X=cbind(Gamma,Gamma_0)%*%phi
P_X=X%*%solve(t(X)%*%X)%*%t(X)
eigen_info_P_X=eigen(P_X)
X0=eigen_info_P_X$vector[,which(eigen_info_P_X$values<0.001)]
alpha=solve(t(X)%*%X)%*%t(X)%*%beta
summary(as.vector(beta-X%*%alpha))
gamma=as.vector(alpha)

q.mis=2
phi.mis=rbind(diag(q.mis),matrix(0,ncol=q.mis,nrow=r-q.mis))
#temp_phi.mis=phi[,-((u-q1+1):(u-q1+2))]
#q.mis=ncol(temp_phi.mis)
#phi.mis=temp_phi.mis#%*%matrix(rnorm(q.mis*3),nrow=q.mis) # multiply the reduced basis by a radom matrix to construct a X.mis
X.mis=cbind(Gamma,Gamma_0)%*%phi.mis
# Function to generate Y with different random error and calculate estimators and variance estimators
sim_fun=function(i){
    set.seed(i)
    print(i)
    Y=beta%*%Z+t(chol(Sigma))%*%matrix(rnorm(r*n),nrow=r)
    Y_D=solve(t(X)%*%X)%*%t(X)%*%Y
    Y_S=t(X0)%*%Y
    
    # env
    hat_u=u.env(t(Z),t(Y))$u.bic
    hat_u.v[i]<<-hat_u
    env_result=env(t(Z),t(Y),hat_u)
    hat_beta_env=env_result[["beta"]]
    hat_avar_hat_beta_env=(env_result[["asySE"]])^2
    
    # OLS
    hat_Sigma_Z=Z%*%t(Z)/n
    hat_Sigma_Z_inv=solve(hat_Sigma_Z)
    hat_beta_ols=Y%*%t(Z)%*%hat_Sigma_Z_inv/n
    hat_Sigma=cov(t(Y-hat_beta_ols%*%Z))
    hat_Sigma_inv=solve(hat_Sigma)
    hat_avar_hat_beta_ols=diag(kronecker(hat_Sigma_Z_inv,hat_Sigma))
   
    # tilde beta
    hat_gamma=prod_fun(Z,X,Y,hat_Sigma_inv)
    hat_alpha=matrix(hat_gamma,nrow=q)
#    tilde_beta=X%*%hat_alpha
    tilde_beta=X%*%solve(t(X)%*%hat_Sigma_inv%*%X)%*%(t(X)%*%hat_Sigma_inv%*%hat_beta_ols)
#    tilde_beta3=X%*%solve(t(X)%*%Sigma_inv%*%X)%*%(t(X)%*%Sigma_inv%*%hat_beta_ols)
#    tilde_beta4=X%*%solve(t(X)%*%Sigma_inv%*%X)%*%(t(X)%*%Sigma_inv%*%beta)
    # avar_tilde_beta=diag(kronecker(hat_Sigma_Z_inv,X%*%solve(t(X)%*%Sigma_inv%*%X)%*%t(X)))
    hat_avar_tilde_beta=diag(kronecker(hat_Sigma_Z_inv,X%*%solve(t(X)%*%hat_Sigma_inv%*%X)%*%t(X)))

    # tilde beta mis
    hat_gamma.mis=prod_fun(Z,X.mis,Y,hat_Sigma_inv)
    hat_alpha.mis=matrix(hat_gamma.mis,nrow=ncol(X.mis))
    tilde_beta.mis=X.mis%*%hat_alpha.mis
    # avar_tilde_beta.mis=diag(kronecker(Sigma_Z_inv,X.mis%*%solve(t(X.mis)%*%Sigma_inv%*%X.mis)%*%t(X.mis)))
    hat_avar_tilde_beta.mis=diag(kronecker(hat_Sigma_Z_inv,X.mis%*%solve(t(X.mis)%*%hat_Sigma_inv%*%X.mis)%*%t(X.mis)))

    # avar_hat_beta_ols=diag(kronecker(hat_Sigma_Z_inv,Sigma))

    # penv    
    hat_u_partial_env=u.penv(X1=t(Z),X2=t(Y_S),Y=t(Y_D))$u.bic
    hat_u_penv.v[i]<<-hat_u_partial_env
    penv_result=penv(X1=t(Z),X2=t(Y_S),Y=t(Y_D),u=hat_u_partial_env)
    hat_alpha_penv=penv_result[["beta1"]]
    hat_beta_penv=X%*%hat_alpha_penv
    hat_avar_hat_beta_new_env=diag(kronecker(diag(p),X)%*%penv_result$covMatrix[1:length(hat_alpha_penv),1:length(hat_alpha_penv)]%*%t(kronecker(diag(p),X)))
    
    
    return(c(hat_beta_ols,hat_beta_env,tilde_beta,tilde_beta.mis, hat_beta_penv,hat_avar_hat_beta_ols,hat_avar_hat_beta_env,hat_avar_tilde_beta,hat_avar_tilde_beta.mis,hat_avar_hat_beta_new_env))
}


sim_num=1000
hat_u.v=rep(NA,sim_num)
hat_u_penv.v=rep(NA,sim_num)
sim_result=apply(as.matrix(1:sim_num),1,sim_fun)
print(summary(hat_u.v))
print(summary(hat_u_penv.v))
sim_result_mean=matrix(apply(sim_result,1,mean),ncol=10)
sim_result_var=matrix(apply(sim_result,1,var),ncol=10)
MC_hat_avar_hat_beta_ols=sim_result_var[,1]
MC_hat_avar_hat_beta_env=sim_result_var[,2]
MC_hat_avar_tilde_beta=sim_result_var[,3]
MC_hat_avar_tilde_beta.mis=sim_result_var[,4]
MC_hat_avar_hat_beta_penv=sim_result_var[,5]
average_hat_avar_hat_beta_ols=sim_result_mean[,6]
average_hat_avar_hat_beta_env=sim_result_mean[,7]
average_hat_avar_tilde_beta=sim_result_mean[,8]
average_hat_avar_tilde_beta.mis=sim_result_mean[,9]
average_hat_avar_hat_beta_penv=sim_result_mean[,10]
mean(average_hat_avar_hat_beta_ols)
mean(average_hat_avar_hat_beta_env)
mean(average_hat_avar_tilde_beta)
mean(average_hat_avar_hat_beta_penv)
avar_hat_beta_ols=diag(kronecker(Sigma_Z_inv,Sigma))
avar_tilde_beta=diag(kronecker(Sigma_Z_inv,X%*%solve(t(X)%*%Sigma_inv%*%X)%*%t(X)))
avar_tilde_beta.mis=diag(kronecker(Sigma_Z_inv,X.mis%*%solve(t(X.mis)%*%Sigma_inv%*%X.mis)%*%t(X.mis)))
tmp_var=kronecker(t(eta),Gamma_0)
U=kronecker(eta%*%Sigma_Z%*%t(eta),solve(Omega_0))+(kronecker(sqrt(Omega),solve(sqrt(Omega_0)))-kronecker(sqrt(Omega_0),solve(sqrt(Omega))))^2
avar_hat_beta_env=diag(kronecker(Sigma_Z_inv,(Gamma%*%Omega%*%t(Gamma)))+tmp_var%*%ginv(U)%*%t(tmp_var))
avar_hat_beta_penv=diag(kronecker(Sigma_Z_inv,solve(t(X)%*%Sigma_inv%*%X)))
cbind(avar_hat_beta_ols,average_hat_avar_hat_beta_ols,MC_hat_avar_hat_beta_ols*n)[1:10,]
cbind(avar_hat_beta_env,average_hat_avar_hat_beta_env,MC_hat_avar_hat_beta_env*n)[1:10,]
cbind(avar_tilde_beta,average_hat_avar_tilde_beta,MC_hat_avar_tilde_beta*n)[1:10,]
#cbind(avar_hat_beta_penv,average_hat_avar_hat_beta_penv,MC_hat_avar_hat_beta_penv*n)[1:10,] #avar_hat_beta_penv was wrong



xtable(rbind(summary(avar_hat_beta_ols),summary(avar_hat_beta_env),summary(avar_tilde_beta)))
xtable(rbind(summary(average_hat_avar_hat_beta_ols),summary(average_hat_avar_hat_beta_env),summary(average_hat_avar_tilde_beta)))
xtable(rbind(summary(MC_hat_avar_hat_beta_ols*n),summary(MC_hat_avar_hat_beta_env*n),summary(MC_hat_avar_tilde_beta*n)))

setwd('/home/lan/Desktop/research/my_paper/constraint_env/simulation')


# cbind(sim_result_var[,c(1,3)]*n,sim_result_mean[,c(7:8)])[1:10,]
cbind(sim_result_var[,1:4]*n,sim_result_mean[,5:8])[1:10,]
# tilde_beta/env var ratio
summary(MC_hat_avar_tilde_beta/MC_hat_avar_hat_beta_env)

# ols/envelope var ratio
# summary(env_result[["ratio"]]^2) # almost the same as below
summary(MC_hat_avar_hat_beta_ols/MC_hat_avar_hat_beta_env)

# ols/tilde_beta var ratio
summary(MC_hat_avar_hat_beta_ols/MC_hat_avar_tilde_beta)

# tilde_beta.mis/tilde_beta var ratio

summary(MC_hat_avar_tilde_beta.mis/MC_hat_avar_tilde_beta)

hat_beta_ols.m=sim_result[1:(p*r),]
hat_beta_env.m=sim_result[(p*r)+(1:(p*r)),]
hat_tilde_beta.m=sim_result[2*(p*r)+(1:(p*r)),]
hat_tilde_beta.mis.m=sim_result[3*(p*r)+(1:(p*r)),]
hat_beta_penv.m=sim_result[4*(p*r)+(1:(p*r)),]


rk_hat_tilde_beta.v=rep(NA,sim_num)
for(i in 1:sim_num){
    tmp=matrix(hat_tilde_beta.m[,i],nrow=nrow(beta))
    rk_hat_tilde_beta.v[i]=qr(tmp)$rank
}
table(rk_hat_tilde_beta.v)
#boxplot(as.vector(hat_beta_env.m-as.vector(beta)),as.vector(hat_tilde_beta.m-as.vector(beta)),outline=F,xaxt='n')


setwd('/home/lan/Desktop/research/my_paper/constraint_env/simulation')
 #pdf(sprintf('Nov24th2018_boxplot_tilde_beta_win_u%d_q1_%d.pdf',u,q_1),width=10, height=5)
xtable(rbind(summary(MC_hat_avar_hat_beta_ols),summary(MC_hat_avar_hat_beta_env),summary(MC_hat_avar_tilde_beta)),digits=4)

if(if_env_win){ 
    pdf(sprintf('Jan14th2021_boxplot_env_win_n%d_u%d_q1_%d_sim%d.pdf',n,u,q1,sim_num),width=10, height=5)
    par(cex.axis=3) 
    boxplot(as.vector(hat_beta_ols.m-as.vector(beta)),as.vector(hat_beta_env.m-as.vector(beta)),as.vector(hat_tilde_beta.m-as.vector(beta)),as.vector(hat_beta_penv.m-as.vector(beta)),outline=F,xaxt='n')
    axis(1,at=c(1:4),adj=1, tick=F, line=3.2,labels=c(expression(paste(hat(beta)['um'])),
                                                      expression(paste(hat(beta)['em'])),
                                                      expression(paste(hat(beta)['cm'])),
                                                      expression(paste(hat(beta)['ecm'])))) # do not plot mis tilde beta since it is way biased as was thought to be all outliers.
    abline(h=0,col=2)
    dev.off()
    save(list = ls(all.names = TRUE),file=sprintf('Jan14th2021_MSE_env_win_n%d_u%d_q1_%d_sim%d.RData',n,u,q1,sim_num))
}else { 
    pdf(sprintf('Jan14th2021_boxplot_tilde_beta_win_n%d_u%d_q1_%d_sim%d.pdf',n,u,q1,sim_num),width=10, height=5)
    par(cex.axis=3) 
    boxplot(as.vector(hat_beta_ols.m-as.vector(beta)),as.vector(hat_beta_env.m-as.vector(beta)),as.vector(hat_tilde_beta.m-as.vector(beta)),as.vector(hat_beta_penv.m-as.vector(beta)),outline=F,xaxt='n')
    axis(1,at=c(1:4),adj=1, tick=F, line=3.2,labels=c(expression(paste(hat(beta)['um'])),
                                                      expression(paste(hat(beta)['em'])),
                                                      expression(paste(hat(beta)['cm'])),
                                                      expression(paste(hat(beta)['ecm'])))) # do not plot mis tilde beta since it is way biased as was thought to be all outliers.
    abline(h=0,col=2)
    dev.off()
    save(list = ls(all.names = TRUE),file=sprintf('Jan14th2021_MSE_tilde_beta_win_n%d_u%d_q1_%d_sim%d.RData',n,u,q1,sim_num))
}

   boxplot(as.vector(hat_beta_ols.m-as.vector(beta)),as.vector(hat_beta_env.m-as.vector(beta)),as.vector(hat_tilde_beta.m-as.vector(beta)),outline=F,xaxt='n')
    axis(1,at=c(1:3),adj=1, tick=F, line=3,labels=c(expression(paste(hat(beta)^'sat')),expression(paste(hat(beta)^'env')),expression(tilde(beta)))) # do not plot mis tilde beta since it is way 

# Average MSE of 4 estimators
#load('Aug11th2020_MSE_tilde_beta_win_u6_q1_4_sim100.RData')
#load('Aug11th2020_MSE_env_win_u6_q1_4_sim100.RData')
MSE_hat_beta_ols=sum(apply((hat_beta_ols.m-as.vector(beta))^2,1,mean))/(p*r)
MSE_hat_beta_env=sum(apply((hat_beta_env.m-as.vector(beta))^2,1,mean))/(p*r)
MSE_hat_tilde_beta=sum(apply((hat_tilde_beta.m-as.vector(beta))^2,1,mean))/(p*r)
MSE_hat_beta_penv=sum(apply((hat_beta_penv.m-as.vector(beta))^2,1,mean))/(p*r)
sum(apply((hat_tilde_beta.m-as.vector(beta))^2,1,mean))/(p*r)
sum(apply((hat_tilde_beta.mis.m-as.vector(beta))^2,1,mean))/(p*r)



###################################
# Simulation to compare the bias  #
###################################
rm(list=ls())
load("/home/lan/Desktop/research/my_paper/constraint_env/simulation/Jan14th2021_MSE_env_win_n200_u6_q1_4_sim1000.RData")
MSE_hat_beta_ols=sum(apply((hat_beta_ols.m-as.vector(beta))^2,1,mean))/(p*r)
MSE_hat_beta_env=sum(apply((hat_beta_env.m-as.vector(beta))^2,1,mean))/(p*r)

cal_tilde_beta_penv_only_1seed_func=function(i,q.mis){
    set.seed(i)
#    print(i)
    Y=beta%*%Z+t(chol(Sigma))%*%matrix(rnorm(r*n),nrow=r)
    hat_Sigma_Z=Z%*%t(Z)/n
    hat_Sigma_Z_inv=solve(hat_Sigma_Z)
    hat_beta_ols=Y%*%t(Z)%*%hat_Sigma_Z_inv/n
    hat_Sigma=cov(t(Y-hat_beta_ols%*%Z))
    hat_Sigma_inv=solve(hat_Sigma)
    phi.mis=rbind(diag(q.mis),matrix(0,ncol=q.mis,nrow=r-q.mis))
    X.mis=cbind(Gamma,Gamma_0)%*%phi.mis
    hat_gamma.mis=prod_fun(Z,X.mis,Y,hat_Sigma_inv)
    hat_alpha.mis=matrix(hat_gamma.mis,nrow=ncol(X.mis))
    tilde_beta.mis=X.mis%*%hat_alpha.mis
    hat_avar_tilde_beta.mis=diag(kronecker(hat_Sigma_Z_inv,X.mis%*%solve(t(X.mis)%*%hat_Sigma_inv%*%X.mis)%*%t(X.mis)))
    P_X.mis=X.mis%*%solve(t(X.mis)%*%X.mis)%*%t(X.mis)
    eigen_info_P_X.mis=eigen(P_X.mis)
    X0.mis=eigen_info_P_X.mis$vector[,which(eigen_info_P_X.mis$values<0.001)]
    Y_D.mis=solve(t(X.mis)%*%X.mis)%*%t(X.mis)%*%Y
    Y_S.mis=t(X0.mis)%*%Y
    
    # penv    
    hat_u_partial_env.mis=u.penv(X1=t(Z),X2=t(Y_S.mis),Y=t(Y_D.mis))$u.bic
     #hat_u_penv.v[i]<<-hat_u_partial_env
    penv_result.mis=penv(X1=t(Z),X2=t(Y_S.mis),Y=t(Y_D.mis),u=hat_u_partial_env.mis)
    hat_alpha_penv.mis=penv_result.mis[["beta1"]]
    hat_beta_penv.mis=X.mis%*%hat_alpha_penv.mis
    hat_avar_hat_beta_new_env.mis=diag(kronecker(diag(p),X.mis)%*%penv_result.mis$covMatrix[1:length(hat_alpha_penv.mis),1:length(hat_alpha_penv.mis)]%*%t(kronecker(diag(p),X.mis)))
 #   if(q.mis>=6){print(hat_u_partial_env.mis)
  #    print(mean(hat_avar_hat_beta_new_env.mis))}
    
    return(c(tilde_beta.mis,hat_beta_penv.mis,hat_avar_tilde_beta.mis,hat_avar_hat_beta_new_env.mis))
    
    
}

cal_tilde_beta_penv_only_func=function(q.mis){
    cat(q.mis,'\n')
    temp_tilde_beta_penv.m=apply(as.matrix(1:sim_num),1,cal_tilde_beta_penv_only_1seed_func,q.mis)
    mean_MSE_tilde_beta=sum(apply((temp_tilde_beta_penv.m[1:(p*r),]-as.vector(beta))^2,1,mean))/(p*r)
    mean_MSE_penv=sum(apply((temp_tilde_beta_penv.m[(p*r)+1:(p*r),]-as.vector(beta))^2,1,mean))/(p*r)
    return(c(mean_MSE_tilde_beta,mean_MSE_penv))
}

tilde_beta_main_fun=function(){
    MSE_result=apply(as.matrix(1:(r-1)),1,cal_tilde_beta_penv_only_func)
}



setwd('/home/lan/Desktop/research/my_paper/constraint_env/simulation')
tilde_beta_penv_mse_result=tilde_beta_main_fun()
if(if_env_win){
    save(list = ls(all.names = TRUE),file=sprintf('Nov1st2021_MSE_env_win_u%d_q1_%d_sim%d_add_vary_X.RData',u,q1,sim_num))
    pdf(sprintf('Nov1st2021_MSE_env_win_u%d_q1_%d_sim%d.pdf',u,q1,sim_num))
}else{
    save(list = ls(all.names = TRUE),file=sprintf('Nov1st2021_MSE_tilde_beta_win_u%d_q1_%d_sim%d_add_vary_X.RData',u,q1,sim_num))
    pdf(sprintf('Nov1st2021_MSE_tilde_beta_win_u%d_q1_%d_sim%d.pdf',u,q1,sim_num))
}
par(cex.axis=2,mar=c(5.1,5.1,4.1,2.1)) 
# save(list = ls(all.names = TRUE),file='tilde_beta_good.RData')
# pdf("tilde_beta_good.pdf")
#plot(5:r,tilde_beta_mse_result[-(1:4)],type='l',ylab='MSE',xlab='q',ylim=c(min(tilde_beta_mse_result[-(1:4)]),max(tilde_beta_mse_result[-(1:4)])+0.02),cex.lab=1.3,cex.axis=1.5,lwd=2)
plot(1:(r-1),tilde_beta_penv_mse_result[1,],type='l',ylab='MSE',xlab=expression(k),cex.lab=2,lwd=2)
#pdf(sprintf('MSE_Omega_%s_Omega0_%d_zoomed.pdf',lambda,lambda0))
#plot(4:r,tilde_beta_mse_result[-(1:3)],type='l',ylab='MSE',xlab='q',ylim=c(min(tilde_beta_mse_result[-(1:4)]),max(tilde_beta_mse_result[-(1:4)])+0.0001),cex.lab=1.3,cex.axis=1.5,lwd=2)
lines(1:(r-1),tilde_beta_penv_mse_result[2,],lty=4,lwd=6)
abline(h=MSE_hat_beta_env,lty=2,lwd=2)
abline(h=MSE_hat_beta_ols,lty=3,lwd=2)
legend('topright',legend=
       c(expression(paste(hat(beta)['cm'])),
         expression(paste(hat(beta)['em'])),
         expression(paste(hat(beta)['um'])),
         expression(paste(hat(beta)['ecm']))),
lty=c(1:3,4),cex=2,lwd=c(2,2,2,6))
dev.off()

 xtable(matrix(tilde_beta_mse_result,byrow=T,ncol=10),digits=4)
tilde_beta_mse_result
MSE_hat_beta_env
MSE_hat_beta_ols



###############################
# Simulation with a general X #
###############################
setwd('/home/lan/Desktop/research/my_paper/constraint_env/simulation')
library('Renvlp')
library('xtable')
# simulation to compare tilde beta and the envelope on alpha
rm(list=ls())
r=20
p=8
n=200
q=15
u_penv=3
temp_penv.m=matrix(rnorm(q*q),nrow=q)
temp_penv.m_normalized=svd(temp_penv.m)$u
Gamma_penv=temp_penv.m_normalized[,1:u_penv]
Gamma_0_penv=temp_penv.m_normalized[,-(1:u_penv)]
Omega_penv=diag(rep(0.5,u_penv))
Omega_0_penv=diag(rep(50,q-u_penv))
Sigma_D_cond_S_penv=Gamma_penv%*%Omega_penv%*%t(Gamma_penv)+Gamma_0_penv%*%Omega_0_penv%*%t(Gamma_0_penv)
X=matrix(rnorm(r*q),nrow=r) # qr(X)$rank
P_X=X%*%solve(t(X)%*%X)%*%t(X)
eigen_info_P_X=eigen(P_X)
X0=eigen_info_P_X$vector[,which(eigen_info_P_X$values<0.001)]
#Y_D=solve(t(X)%*%X)%*%t(X)%*%Y
#Y_S=t(X0)%*%Y
temp_sigma_Z.m=matrix(rnorm(p^2),nrow=p)
Sigma_Z=temp_sigma_Z.m%*%t(temp_sigma_Z.m)
#Sigma_Z=diag(p)
Sigma_Z_inv=solve(Sigma_Z)
Z=t(chol(Sigma_Z))%*%matrix(rnorm(n*p),nrow=p) # Z should be mean 0 as when we estimate Sigma_Z we use Z%*%t(Z)/n without subtracting the mean
eta_penv=matrix(rnorm(u_penv*p),nrow=u_penv)
#eta=rbind(matrix(0,nrow=u-q1,ncol=p),temp_eta2.m) # gaurantee rank(eta) is no greater than q1, thus the rank(beta) is no greater than q1
alpha=Gamma_penv%*%eta_penv
beta=X%*%alpha
phi=matrix(rnorm(q*(r-q)),nrow=q)
A=cbind(X%*%solve(t(X)%*%X),X0)
A_inv=solve(A)
cal_scale=0
if(cal_scale==T){
    Lambda=diag(c(1, runif(n=q-1,2,10)))
} else {
    Lambda=diag(rep(1,q))
}


sim_penv_fun=function(i){
    set.seed(i)
    print(i)
    Y_S = matrix(rnorm((r-q)*n),nrow=r-q)
    Y_D = solve(Lambda)%*%alpha %*% Z+phi%*%Y_S+matrix(rnorm(q*n),nrow=q)
    Y=t(A_inv)%*%rbind(Y_D,Y_S)
    hat_Sigma_Z=Z%*%t(Z)/n
    hat_Sigma_Z_inv=solve(hat_Sigma_Z)
    hat_beta_ols=Y%*%t(Z)%*%hat_Sigma_Z_inv/n
    hat_Sigma=cov(t(Y-hat_beta_ols%*%Z))
    hat_Sigma_inv=solve(hat_Sigma)

    # tilde beta
    tilde_beta=X%*%solve(t(X)%*%hat_Sigma_inv%*%X)%*%(t(X)%*%hat_Sigma_inv%*%hat_beta_ols)
    # avar_tilde_beta=diag(kronecker(hat_Sigma_Z_inv,X%*%solve(t(X)%*%Sigma_inv%*%X)%*%t(X)))
    hat_avar_tilde_beta=diag(kronecker(hat_Sigma_Z_inv,X%*%solve(t(X)%*%hat_Sigma_inv%*%X)%*%t(X)))

    # new env
    hat_u_partial_env=u.penv(X1=t(Z),X2=t(Y_S),Y=t(Y_D))$u.bic
    hat_u_penv.v[i]<<-hat_u_partial_env
    penv_result=penv(X1=t(Z),X2=t(Y_S),Y=t(Y_D),u=hat_u_partial_env)
    hat_alpha_penv=penv_result[["beta1"]]
    hat_beta_penv=X%*%hat_alpha_penv
    hat_avar_hat_beta_new_env=diag(kronecker(diag(p),X)%*%penv_result$covMatrix[1:length(hat_alpha_penv),1:length(hat_alpha_penv)]%*%t(kronecker(diag(p),X)))
    
    diag(kronecker(hat_Sigma_Z_inv,solve(t(X)%*%hat_Sigma_inv%*%X)))

    return(c(tilde_beta,hat_beta_penv, hat_avar_tilde_beta,hat_avar_hat_beta_new_env))
 
    
}


sim_num_penv=1000
hat_u_penv.v=rep(NA,sim_num_penv)

sim_result_penv=apply(as.matrix(1:sim_num_penv),1,sim_penv_fun)
print(summary(hat_u_penv.v))
sim_result_penv_mean=matrix(apply(sim_result_penv,1,mean),ncol=4)
sim_result_penv_var=matrix(apply(sim_result_penv,1,var),ncol=4)
MC_hat_avar_tilde_beta=sim_result_penv_var[,1]
MC_hat_avar_hat_beta_penv=sim_result_penv_var[,2]
mean(MC_hat_avar_tilde_beta)
mean(MC_hat_avar_hat_beta_penv)
average_hat_avar_tilde_beta_penv=sim_result_penv_mean[,3]
average_hat_avar_hat_beta_penv=sim_result_penv_mean[,4]
xtable(rbind(summary(average_hat_avar_tilde_beta_penv),summary(average_hat_avar_hat_beta_penv)))
xtable(rbind(summary(MC_hat_avar_tilde_beta_penv*n),summary(MC_hat_avar_hat_beta_penv*n)))
tilde_beta_penv.m=sim_result_penv[1:(p*r),]
hat_beta_penv.m=sim_result_penv[(p*r)+(1:(p*r)),]
#avar_tilde_beta=diag(kronecker(Sigma_Z_inv,X%*%solve(t(X)%*%Sigma_inv%*%X)%*%t(X)))

pdf(sprintf('Aug11th2020_boxplot_tilde_beta_newenv_sim%d.pdf',sim_num_penv),width=10, height=5)
par(cex.axis=3) 
boxplot(as.vector(tilde_beta_penv.m-as.vector(beta)),as.vector(hat_beta_penv.m-as.vector(beta)),outline=F,xaxt='n')
axis(1,at=c(1:2),adj=1, tick=F, line=3.2,labels=c(expression(tilde(beta)),expression(paste(hat(beta)^'penv')))) # do not plot mis tilde beta since it is way biased as was thought to be all outliers.
abline(h=0,col=2)
dev.off()


    ylim_upper_tildebetawin=20
    pdf('Aug11th2020_tilde_beta_penv_hatavar_MCavar_hist.pdf')
    par(mfrow=c(2,2),cex.axis=2.5)
 #   hist(avar_tilde_beta_penv,breaks=100,cex.lab=1.5,xlab='',main=expression(paste('avar(',tilde(bold(beta)),')')),ylim=c(0,ylim_upper_tildebetawin));
    hist(average_hat_avar_tilde_beta_penv,breaks=100,cex.lab=1.5,xlab='',main=expression(paste(widehat('avar'),'(',tilde(bold(beta)),')')),ylim=c(0,ylim_upper_tildebetawin));
    hist(MC_hat_avar_tilde_beta_penv*n,breaks=100,cex.lab=1.5,xlab='',main=expression(paste('MC avar(',tilde(bold(beta)),')')),ylim=c(0,ylim_upper_tildebetawin))
#    hist(avar_hat_beta_env,breaks=100,cex.lab=1.5,xlab='',main=expression(paste('avar(',paste(hat(bold(beta))^'env'),')')),ylim=c(0,ylim_upper_tildebetawin));
    hist(average_hat_avar_hat_beta_penv,breaks=100,cex.lab=1.5,xlab='',main=expression(paste(widehat('avar'),'(',paste(hat(bold(beta))^'penv'),')')),ylim=c(0,ylim_upper_tildebetawin));
    hist(MC_hat_avar_hat_beta_penv*n,breaks=100,cex.lab=1.5,xlab='',main=expression(paste('MC avar(',paste(hat(bold(beta))^'penv'),')')),ylim=c(0,ylim_upper_tildebetawin))
    dev.off()
    
MSE_tilde_beta_penv=sum(apply((tilde_beta_penv.m-as.vector(beta))^2,1,mean))/(p*r) # 0.004349856
MSE_hat_beta_penv=sum(apply((hat_beta_penv.m-as.vector(beta))^2,1,mean))/(p*r) # 0.001048832

    
   save(list = ls(all.names = TRUE),file='Aug11th2020_add_new_env.RData')

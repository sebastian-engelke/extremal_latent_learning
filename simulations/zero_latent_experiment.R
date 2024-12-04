library(graphicalExtremes)
library(igraph)
library(tidyverse)
library(latex2exp)
library(glmnet)  
library(egg)
library(cowplot)
library(tictoc)
library(here)
library(clusterGeneration)
library(pracma)
library(matrixcalc)
library(CVXR)



source(here("functions/functions_paper.R"))


set.seed(123)
d<-20 # number of observed variables 

ratio <- 0.7 # k = n^ratio where k is the effective sample size
n <- as.integer(c(exp(log(2000)/ratio)))
rholist <- append(seq(0.0001, 0.001, length.out = 10),seq(0.001, 0.05, length.out = 20)); # this represents sequence of \lambda_n for eglatent and rho for eglearn
alpha = 4 # how much rho gets scaled for eglearn
method="maxstable"
m = 2
gen_model = "BA_fixed" # BA model of degree m=2 
val_set = TRUE
num_iter<-50
lambda_2_vec <- c(1,4,8,20) #gamma for eglatent
BA_model <- generate_BA_model(d = d, m = m)
p=1-n^0.7/n
h=0 # no latent variables

prob_recov_cv <- list()
prob_recov_oracle <- list()
F1score_cv <- list()
F1score_oracle <- list()
rk_oracle <- list()
rk_cv <- list()
F1score_eglearn_cv <- list()
F1score_eglearn_oracle <- list()





F1score_cv_vec <- list()
rk_cv_vec <-list()
F1score_eglearn_cv_vec <- list()
likelihood_vec <- list()
likelihood_eglearn_vec<-list()
sparsity_vec <- list()
sparsity_eglearn_vec <- list()


for (lambda_2_iter in 1:length(lambda_2_vec)){
  lambda_2 <- lambda_2_vec[lambda_2_iter]
  likelihood_vec_temp<-  c()
  F1score_cv_vec_temp <-c()
  rank_cv_vec_temp <-c()
  F1score_eglearn_cv_temp <-c()
  likelihood_eglearn_vec_temp <- c()
  sparsity_vec_temp <- c()
  sparsity_eglearn_vec_temp <- c()
  
  for (iter_ind in 1:num_iter){
   

      # runs both the estimator without latent variables and the estimator with latent variables for range of regularization parameters
      output <- eglatent_path(d = d, n=n, p=p, m = m, h=h, lambda_2 = lambda_2, rholist = rholist, method=method, gen_model = gen_model, val_set = val_set,alpha = 4,BA_model = BA_model)
      
      ##### extracting and analyzing results for the latent model
      likelihood<- output$loglik_latent_refit # likelihood scores of the models learned on validation data of same size
      F1_latent<- output$F1_latent # F1 of graph structure
      rk_output<- output$rk[which.max(likelihood)] # number of latent variables
      sparsity<- output$sparse_latent
      
     
    # score across different lambda_2
    likelihood_vec_temp<-  append(likelihood_vec_temp,max(likelihood))
    F1score_cv_vec_temp <-append(F1score_cv_vec_temp,F1_latent[which.max(likelihood)])
    rank_cv_vec_temp <-append(rank_cv_vec_temp,rk_output)
    sparsity_vec_temp<-append(sparsity_vec_temp,sparsity[which.max(likelihood)])
    
    
    likelihood_eglearn<- output$loglik_eglearn
    F1_eglearn<- output$F1_eglearn
    temp<-which(is.na(likelihood_eglearn)==FALSE)
    likelihood_eglearn_vec_temp<-append(likelihood_eglearn_vec_temp,max(likelihood_eglearn[temp]))
    F1score_eglearn_cv_temp <- append(F1score_eglearn_cv_temp,F1_eglearn[temp[which.max(likelihood_eglearn[temp])]])
    sparsity_eglearn_vec_temp<-append(sparsity_eglearn_vec_temp,output$sparse_eglearn[temp[which.max(likelihood_eglearn[temp])]])
    
  }
  
  F1score_cv_vec <- append(F1score_cv_vec,list(F1score_cv_vec_temp))
  rk_cv_vec <-append(rk_cv_vec,list(rank_cv_vec_temp))
  likelihood_vec <- append(likelihood_vec,list(likelihood_vec_temp))
  sparsity_vec <- append(sparsity_vec,list(sparsity_vec_temp))
  
}
F1score_eglearn_cv_vec <- append(F1score_eglearn_cv_vec,list(F1score_eglearn_cv_temp))
likelihood_eglearn_vec<-append(likelihood_eglearn_vec,list(likelihood_eglearn_vec_temp))
sparsity_eglearn_vec<-append(sparsity_eglearn_vec,list(sparsity_eglearn_vec_temp))

data<- tibble(F1score_cv_vec = F1score_cv_vec, rk_cv_vec = rk_cv_vec, likelihood_vec = likelihood_vec, F1score_eglearn_cv_vec = F1score_eglearn_cv_vec, likelihood_eglearn_vec = likelihood_eglearn_vec)
save(data,file = here("simulations/data/data_zero_latent.Rdata"))              




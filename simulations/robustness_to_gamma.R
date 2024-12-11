## stability to the choice of \gamma


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


set.seed(1242)
rholist <- seq(0.0001, .05, length.out = 20) # this represents sequence of \lambda_n for eglatent and rho for eglearn
alpha =10 # how much rho gets scaled for eglearn
method="maxstable"
gen_model ="latent_cycle" # latent model with a cycle graph
d = 30 # number of observed variables
ratio <- 0.65 # k = n^ratio where k is the effective sample size
nvec <- rev(as.integer(c(exp(log(1000)/ratio))))
hvec <- c(1,2,3) # number of latent variables
m = 1
plot_result = FALSE
lambda_2_vec<- c(2,4,6)# gamma for eglatent
num_iter <- 50 # number of iterations

F1score_cv <- list()
F1score_oracle <- list()
rk_oracle <- list()
rk_cv <- list()
F1score_eglearn_cv <- list()
F1score_eglearn_oracle <- list()
likelihood_cv <- list()
likelihood_eglearn_cv <- list()

for (lambda_2_iter in 1:length(lambda_2_vec)){
  lambda_2 <- lambda_2_vec[lambda_2_iter]
for (n_iter in 1:length(nvec)){
  n <- nvec[n_iter]
  for (h_iter in 1:length(hvec)){
    h <- hvec[h_iter]
    
    prob_recov_cv_vec<- c()
    prob_recov_oracle_vec<-c()
    F1score_cv_vec <- c()
    F1score_oracle_vec <- c()
    rk_cv_vec <-c()
    rk_oracle_vec <- c()
    F1score_eglearn_cv_vec <- c()
    F1score_eglearn_oracle_vec <- c()
    likelihood_vec <- c()
    likelihood_eglearn_vec <- c()
    
    
    for (iter in 1:num_iter){
      print(iter)
      # runs both the estimator without latent variables and the estimator with latent variables for range of regularization parameters
      output <- eglatent_path(d = d, n = n, p=1-n^ratio/n,h = h, m = m, lambda_2 = lambda_2, rholist = rholist, method=method, gen_model =gen_model, val_set = TRUE,plot_result = plot_result,alpha = alpha)
      
      
      ##### extracting and analyzing results for the latent model
      likelihood<-output$loglik_latent_refit # likelihood scores of the models learned on validation data of same size
      F1_latent<-output$F1_latent # F1 of graph structure
      rk<-output$rk # number of latent variables
      
      # F1 score across the 50 trials
      F1score_cv_vec<-append(F1score_cv_vec,F1_latent[which.max(likelihood)])
      F1score_oracle_vec <- append(F1score_oracle_vec,max(F1_latent))
      
      # rank across the 50 trials
      rk_cv_vec <- append(rk_cv_vec,rk[which.max(likelihood)])
      rk_oracle_vec <- append(rk_oracle_vec,rk[which.max(F1_latent)])
      
      likelihood_vec <- append(likelihood_vec,max(likelihood))
      
      
      
      # extracting results for the model without latent variables
      likelihood_eglearn<-output$loglik_eglearn
      F1_eglearn<-output$F1_eglearn
      
      # store results across ten trials 
      temp<-which(is.na(likelihood_eglearn)==FALSE)
      likelihood_eglearn_vec <- append(likelihood_eglearn_vec,max(likelihood_eglearn[temp]))
      F1score_eglearn_cv_vec <- append(F1score_eglearn_cv_vec,F1_eglearn[temp[which.max(likelihood_eglearn[temp])]])
      F1score_eglearn_oracle_vec <- append(F1score_eglearn_oracle_vec,max(F1_eglearn))
      
    }
    F1score_cv <- append(F1score_cv,list(F1score_cv_vec))
    F1score_oracle <- append(F1score_oracle,list(F1score_oracle_vec))
    rk_oracle <- append(rk_oracle,list(rk_oracle_vec))
    rk_cv <- append(rk_cv,list(rk_cv_vec))
    likelihood_cv<-append(likelihood_cv,list(likelihood_vec))
    
    F1score_eglearn_cv <- append(F1score_eglearn_cv,list(F1score_eglearn_cv_vec))
    F1score_eglearn_oracle <- append(F1score_eglearn_oracle,list(F1score_eglearn_oracle_vec))
    likelihood_eglearn_cv<-append(likelihood_eglearn_cv,list(likelihood_eglearn_vec))
    
  }
}

}
data<- tibble(F1score_cv = F1score_cv,  F1score_oracle =  F1score_oracle, rk_oracle = rk_oracle, rk_cv = rk_cv, likelihood_cv = likelihood_cv, F1score_eglearn_cv = F1score_eglearn_cv, F1score_eglearn_oracle = F1score_eglearn_oracle, likelihood_eglearn_cv = likelihood_eglearn_cv)
save(data,file = here("simulations/data/robustness_to_gamma.Rdata"))              




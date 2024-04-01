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
library(ggplot2)
library(ggpubr)

source(here("functions/functions_paper.R"))
figure_dest_folder <- here("figures/plots/")

####################

#################################
##### Regularization Paths ######
#################################



set.seed(333242)

d <- 30 # number of observed nodes
k <- 5000 # number of effective extreme samples
n <- floor(k^{1/.7})
p <- 1 - n^0.7 / n
m <- 2
h <- 2 # number of latent variables
rholist <- seq(0.0001, .05, length.out = 15) # this represents sequence of \lambda_n for eglatent and rho for eglearn
alpha <- 4 # this is how much rholist gets scaled for eglearn
lambda_2 <- 4 # gamma for eglatent
method <- "maxstable"
gen_model <- "latent_cycle" # latent model with a cycle graph
val_set <- TRUE
plot_result <- TRUE 

tbl1 <- eglatent_path(d = d, n = n, p = p, m = m, h = h, lambda_2 = lambda_2, rholist = rholist, method = method, gen_model = gen_model, val_set = val_set,plot_result = plot_result,alpha = alpha,figure_dest_folder=figure_dest_folder)






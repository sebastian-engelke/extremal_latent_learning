# install.packages("devtools")
# devtools::install_github("sebastian-engelke/graphicalExtremes")

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
figure_dest_folder <- "application/plots/"

# Get IATAs from Texas cluster
IATAs <- getFlightDelayData("IATAs", "tcCluster")

MIN_N_FLIGHTS <- 20
flight_graph <- getFlightGraph(IATAs, minNFlights = MIN_N_FLIGHTS)

# Plot airports + connections (indicating number of flights by thickness)

gg0 <- plotFlights(
  IATAs,
  graph = flight_graph,
  useAirportNFlights = TRUE,
  useConnectionNFlights = FALSE,
  returnGGPlot = TRUE,
  clipMap = 1.3,
  xyRatio = 1
)

gg0

save_myplot(
  plt = gg0,
  plt_nm = paste0(figure_dest_folder, "flight_map", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)


# Check whether all dates from the train-test-split are available
# (due to size restrictions, the CRAN version of this package does not contain all dates)

allDatesAvailable <- tryCatch(
  {
    getFlightDelayData("dates", dateFilter = c("tcAll"))
    TRUE
  },
  error = function(...) FALSE
)
cat("All dates available:", allDatesAvailable, "\n")

# Full data set 3603 x 29 (if allDatesAvailable == TRUE)
mat <- rbind(drop(getFlightDelayData("delays", "tcCluster", "tcTrain")), drop(getFlightDelayData("delays", "tcCluster", "tcTest")))

eglatent_test <- 0
eglearn_test <- 0
empirical_test <- 0
flight_test <- 0  

# Probability threshold to define extremes
n <- nrow(mat)
p <- 0.95 #1-n^0.65/n #0.85 
k <- (1-p) * n
d <- ncol(mat)

# Regularization parameter to enforce sparsity
lambda1_range <- seq(0, .1, by = 0.01) # seq(0, 0.1, by = 0.002)
ll <- length(lambda1_range)
# Regularization parameter to enforce low-rank
lambda2_range <- c(4) # 2 # seq(0, 0.3, by = 0.08) #seq(0, 0.3, by = 0.02)
# for eglearn
rholist <- lambda1_range



## cross-validation for model assessment
set.seed(4134999)

splitInd <- floor(nrow(mat) * 1 / 5)

for (cvsplit in 1:5) {

  print(cvsplit)

  # Validation set for this iteration
  ValInd <- c((splitInd * (cvsplit - 1) + 1):(splitInd * cvsplit))
  matVal <- mat[ValInd, ]
  # Training set for this iteration
  matEst <- mat[setdiff(1:nrow(mat), ValInd), ]

  # Normalize data to multivariate Pareto scale
  train_data <- data2mpareto(data = matEst, p = p)
  test_data <- data2mpareto(data = matVal, p = p)
  
  # empirical variogram on training and test sets
  Gamma_train <- emp_vario(data = train_data)
  Gamma_test <- emp_vario(data = test_data)

  # fit latent model
  fit_eglatent <- eglatent(
    Gamma = Gamma_train,
    lam1_list = lambda1_range,
    lam2_list = lambda2_range,
    refit = TRUE
  )

  eglatent_test <- eglatent_test + sapply(1:(ll*length(lambda2_range)), FUN = function(i) {
    loglik_HR(data = test_data, Gamma = fit_eglatent$G_refit[[i]], cens = FALSE)[1]
  })

  flight_test <- flight_test + loglik_HR(data = test_data, Gamma = complete_Gamma(Gamma = Gamma_train, graph = flight_graph), cens = FALSE)[1]

  empirical_test <- empirical_test + loglik_HR(data = test_data, Gamma = Gamma_train, cens = FALSE)[1]

  # fit eglearn model

  fit_eglearn <- eglearn(
    data = train_data,
    rholist = rholist,
    reg_method = "ns",
    complete_Gamma = TRUE
  )

  connected_eglearn <- sapply(1:ll, FUN = function(i) ifelse(is_connected(fit_eglearn$graph[[i]]), 1, 3))

  eglearn_test <- eglearn_test + sapply(1:ll, FUN = 
  function(i) ifelse(connected_eglearn[i] == 1, loglik_HR(data = test_data, Gamma = complete_Gamma(Gamma = Gamma_train, graph = fit_eglearn$graph[[i]]), cens = FALSE)[1], NA))
}



dat <- tibble(
  rholist = lambda1_range,
  eglearn_test = eglearn_test,
  flight_test = rep(flight_test, times = ll), 
  empirical_test = rep(empirical_test, times = ll), 
  eglatent_test = eglatent_test
)


k_latent_max <- which.max(eglatent_test)
k_eglearn_max <- which.max(eglearn_test)



gg2 <- ggplot(dat) +
  geom_line(aes(x = rholist, y = eglatent_test / 5), linetype = "solid") +
  geom_point(aes(x = rholist, y = eglatent_test / 5), shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(yintercept = eglearn_test[1] /5, linetype = "longdash") +
  geom_line(aes(x = rholist, y = eglearn_test / 5), linetype = "dashed") +
  geom_point(aes(x = rholist, y = eglearn_test / 5), shape = 21, size = 3, stroke = 1, fill = "white") +
   xlab(TeX("Regularization parameter $\\lambda_n$")) +
    ylab("Log-likelihood") +
    scale_x_continuous(
      breaks = unique(round(rholist, 2)),
      labels = unique(round(rholist, 2))
    )
gg2






## Fiting the models on the whole data set

lambda2 <- 2

train_data <- data2mpareto(data = mat, p = p)
Gamma_train <- emp_vario(data = train_data)
fit_eglatent <- eglatent(
  Gamma = Gamma_train,
  lam1_list = lambda1_range,
  lam2_list = lambda2,
  refit = TRUE
)
eglatent_sparse <- sapply(1:ll, FUN = function(i) length(E(fit_eglatent$graph[[i]])))
eglatent_rk <- fit_eglatent$rk


fit_eglearn <- eglearn(
  data = train_data,
  rholist = rholist,
  reg_method = "ns",
  complete_Gamma = TRUE
)
connected_eglearn <- sapply(1:ll, FUN = function(i) ifelse(is_connected(fit_eglearn$graph[[i]]), 1, 3))
eglearn_sparse <- sapply(1:ll, FUN = function(i) length(E(fit_eglearn$graph[[i]])))



dat_sparse <- tibble(rholist = rholist, eglatent_edges = eglatent_sparse, eglatent_rk = eglatent_rk, eglearn_edges = eglearn_sparse)


gg1 <- ggplot(dat_sparse) +
  geom_line(aes(x = rholist, y = eglatent_edges), linetype = "solid") +
  geom_point(aes(x = rholist, y = eglatent_edges), shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_line(aes(x = rholist, y = eglearn_edges), linetype = "dashed") +
  geom_point(aes(x = rholist, y = eglearn_edges), shape = 21, size = 3, stroke = 1, fill = "white") +
    xlab(TeX("Regularization parameter $\\lambda_n$")) +
    ylab("Number of edges") +
    scale_x_continuous(
      breaks = unique(round(rholist, 2)),
      labels = unique(round(rholist, 2)),
      sec.axis = sec_axis(
        trans = ~., breaks = rholist,
        labels = eglatent_rk,
        name = ""
      )
    )

gg1

gg_both <- ggarrange(gg1, gg2, nrow=1, align = "hv")

save_myplot(
  plt = gg_both,
  plt_nm = paste0(figure_dest_folder, "sparsity_likelihoods_p95", ".pdf"),
  width = 10,
  height = 5,
  cairo = FALSE
)


gg3 <- plotFlights(
  IATAs,
  graph = fit_eglatent$graph[[k_latent_max]],
  xyRatio = 1,
  clipMap = 1.3,
  returnGGPlot = TRUE,
  useAirportNFlights = TRUE,
)
gg3

nb_neighbor <- sapply(1:d, FUN = function(i) length(neighbors(fit_eglatent$graph[[k_latent_max]], v = i))) 

IATAs[which(nb_neighbor == max(nb_neighbor))]

save_myplot(
  plt = gg3,
  plt_nm = paste0(figure_dest_folder, "flights_latent_p95", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)



gg4 <- plotFlights(
  IATAs,
  graph = fit_eglearn$graph[[k_eglearn_max]],
  xyRatio = 1,
  clipMap = 1.3,
  returnGGPlot = TRUE,
  useAirportNFlights = TRUE
)
gg4

save_myplot(
  plt = gg4,
  plt_nm = paste0(figure_dest_folder, "flights_eglearn_p95", ".pdf"),
  width = 5,
  height = 5,
  cairo = FALSE
)



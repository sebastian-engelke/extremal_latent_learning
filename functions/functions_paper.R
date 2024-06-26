
eglatent_path <- function(d = 30,
                          n = 1000,
                          p = NULL,
                          m = 2,
                          h = 2,
                          lambda_2,
                          rholist,
                          method = c("maxstable", "mpareto"),
                          gen_model = c("BA", "latent_cycle", "latent_random","flights"),
                          val_set = FALSE, plot_result = FALSE,alpha = 10, BA_model = NULL,figure_dest_folder=NULL, ...) {
  if (gen_model == "BA") {
    BA_model <- generate_BA_model(d = d, m = m)
    G <- BA_model$G
    g <- BA_model$graph
  }
  if (gen_model == "BA_fixed") {
    G <- BA_model$G
    g <- BA_model$graph
  }
  if (gen_model == "latent_cycle") {
    latent_model <- generate_latent_model_cycle(p = d, h)
    G <- 10 * latent_model$Gamma[1:d, 1:d]
    g <- subgraph(latent_model$graph, 1:d)
  }
  if (gen_model == "latent_random") {
    latent_model <- generate_latent_model_random(p = d, h)
    G <- 10 * latent_model$Gamma[1:d, 1:d]
    g <- subgraph(latent_model$graph, 1:d)
  }
  
  
  # perform simulation
  if (method == "maxstable") X <- rmstable(n = n, d = d, model = "HR", par = G)
  if (method == "mpareto") X <- rmpareto(n = n, d = d, model = "HR", par = G)
  
  G_emp <- emp_vario(data = X, p = p)
  
  
  # fit latent model
  fit_latent <- eglatent(Gamma = G_emp, lam1_list = rholist, lam2_list = lambda_2, refit = TRUE)
  
  F1_latent <- sapply(1:length(rholist), FUN = function(i) F1_score(g = g, gest = fit_latent$graph[[i]]))
  sparse_latent <- sapply(1:length(rholist), FUN = function(i) length(E(fit_latent$graph[[i]])) / (d * (d - 1) / 2))
  rk <- fit_latent$rk
  
  # fit eglearn model
  fit_eglearn <- eglearn(data = X, p = p, rholist = alpha* rholist, reg_method = "ns")
  F1_eglearn <- sapply(1:length(rholist), FUN = function(i) F1_score(g = g, gest = fit_eglearn$graph[[i]]))
  connected_eglearn <- sapply(1:length(rholist), FUN = function(i) ifelse(is_connected(fit_eglearn$graph[[i]]), 21, 4))
  sparse_eglearn <- sapply(1:length(rholist), FUN = function(i) length(E(fit_eglearn$graph[[i]])) / (d * (d - 1) / 2))
  
  
  dat_sparse <- tibble(rholist = rholist, F1_latent = F1_latent, F1_eglearn = F1_eglearn, sparse_latent = sparse_latent, rk = rk, sparse_eglearn = sparse_eglearn)
  
  if (plot_result){
    gg1 <- ggplot(dat_sparse) +
      geom_line(aes(x = rholist, y = F1_latent), linetype = "solid") +
      geom_point(aes(x = rholist, y = F1_latent), shape = 21, size = 3, stroke = 1, fill = "white") +
      geom_line(aes(x = rholist, y = F1_eglearn), linetype = "dashed") +
      geom_point(aes(x = rholist, y = F1_eglearn), shape = connected_eglearn, size = 3, stroke = 1, fill = "white") +
      xlab(TeX("Regularization parameter $\\lambda_n$")) +
      ylab("F-score") +
      ylim(0, 1) +
      scale_x_continuous(
        breaks = unique(round(rholist, 2)),
        labels = unique(round(rholist, 2)),
        sec.axis = sec_axis(
          trans = ~., breaks = rholist,
          labels = rk,
          name = ""
        )
      )
    gg1
  }
  
  
  if (val_set == TRUE) {
    if (method == "maxstable") {
      X_val_tmp <- rmstable(n = n, d = d, model = "HR", par = G)
      X_val <- data2mpareto(X_val_tmp, p = p)
    }
    if (method == "mpareto") X_val <- rmpareto(n = n, d = d, model = "HR", par = G)
    
    loglik_latent <- sapply(1:length(rholist), FUN = function(i) {
      loglik_HR(data = X_val, Gamma = fit_latent$G_est[[i]], cens = FALSE)[1]
    })
    
    loglik_latent_refit <- sapply(1:length(rholist), FUN = function(i) {
      loglik_HR(data = X_val, Gamma = fit_latent$G_refit[[i]], cens = FALSE)[1]
    })
    
    loglik_eglearn <- sapply(1:length(rholist), FUN = function(i) ifelse(connected_eglearn[i] == 21, loglik_HR(data = X_val, Gamma = complete_Gamma(Gamma = G_emp, graph = fit_eglearn$graph[[i]]), cens = FALSE)[1], NA))
    
    
    dat <- tibble(rholist = rholist, loglik_latent_refit = loglik_latent_refit, loglik_eglearn = loglik_eglearn)
    
    if (plot_result){
      gg2 <- ggplot(dat) +
        geom_line(aes(x = rholist, y = loglik_latent_refit), linetype = "solid") +
        geom_point(aes(x = rholist, y = loglik_latent_refit), shape = 21, size = 3, stroke = 1, fill = "white") +
        geom_line(aes(x = rholist, y = loglik_eglearn), linetype = "dashed") +
        geom_point(aes(x = rholist, y = loglik_eglearn), shape = connected_eglearn, size = 3, stroke = 1, fill = "white") +
        xlab(TeX("Regularization parameter $\\lambda_n$")) +
        ylab("Log-likelihood") +
        coord_cartesian(ylim=c(min(loglik_latent_refit), max(loglik_eglearn, loglik_latent_refit))) +
        # ylim(min(loglik_latent_refit), max(loglik_eglearn, loglik_latent_refit)) + # - .5*(max(loglik_latent_refit) - min(loglik_latent_refit))
        scale_x_continuous(
          breaks = unique(round(rholist, 2)),
          labels = unique(round(rholist, 2)),
          sec.axis = sec_axis(
            trans = ~., breaks = rholist,
            labels = sapply(
              fit_latent$graph,
              igraph::gsize
            )
          )
        )
      
      gg2
      gg_both <- ggarrange(gg1, gg2, nrow=1, align = "hv")
      
      figure_dest_folder <- 
      
      save_myplot(
        plt = gg_both,
        plt_nm = paste0(figure_dest_folder, "F1_likelihood", ".pdf"),
        width = 10,
        height = 5,
        cairo = FALSE
      )
    }
  }
  return(list(F1_latent=F1_latent,rk=rk,loglik_latent_refit=loglik_latent_refit,F1_eglearn=F1_eglearn,loglik_eglearn=loglik_eglearn))
}


generate_latent_model_random <- function(p, h) {
  
  max_degree <- 5
  min_degree <- 1
  non_pos_def<- TRUE
  not_connected <- TRUE
  while (max_degree > 3 || min_degree <1 || non_pos_def){
    W <- matrix(1, p + h, p + h)
    
    S<-as.matrix(as.matrix(as.matrix(erdos.renyi.game(p,0.08))))
    not_connected <- !(is_connected(graph_from_adjacency_matrix(S)))
    max_degree <- max(rowSums(S))
    min_degree <- min(rowSums(S))
    
    
    W[1:p, 1:p] <- S
    diag(W) <- 0
    L <- matrix(runif((p + h)^2, 2, 2), nrow = p + h)
    L[(p + 1):(p + h), (p + 1):(p + h)] <- diag(h)
    
    
    z <- 30/ (sqrt(as.integer(p / h))) * matrix(runif(as.integer(p / h), 1, 2), nrow = as.integer(p / h))
    
    for (i in 1:h) {
      L[1:p, (p + i):(p + i)] <- 0
      L[(p + i):(p + i), 1:p] <- 0
      F<-seq(i, p, h)
      L[F[1:length(z)], (p + i):(p + i)] <- z
      L[(p + i):(p + i), F[1:length(z)]] <- z
    }
    
    W <- W * L
    W[lower.tri(W)] <- t(W)[lower.tri(W)]
    Theta <- diag(rowSums(W)) - W
    non_pos_def <- (min(eig(Theta))< -10^(-6))
  }
  Lst <- (Theta[1:p, ((p + 1):(p + h))]) %*% solve(Theta[(p + 1):(p + h), (p + 1):(p + h)]) %*% t((Theta[1:p, ((p + 1):(p + h))]))
  inc <- max(diag(svd(Lst)$u[, 1:h] %*% t(svd(Lst)$u[, 1:h])))
  G <- Theta2Gamma(Theta)  
  Lst <- (Theta[1:p, ((p + 1):(p + h))]) %*% solve(Theta[(p + 1):(p + h), (p + 1):(p + h)]) %*% t((Theta[1:p, ((p + 1):(p + h))]))
  return(list(Gamma = G, graph = Gamma2graph(G), Lst = Lst))
}


generate_latent_model_cycle <- function(p, h) {
  W <- matrix(1, p + h, p + h)
  S <- matrix(0, p, p)
  for (i in 1:(p - 1)) {
    S[i, i + 1] <- 1
  }
  S[1, p] <- 1
  S <- S + t(S)
  diag(S) <- 0
  
  W[1:p, 1:p] <- S
  diag(W) <- 0
  L <- matrix(runif((p + h)^2, 2, 2), nrow = p + h)
  if(h!=0){
  L[(p + 1):(p + h), (p + 1):(p + h)] <- diag(h)
  
  
  z <- 30 / (sqrt(as.integer(p / h))) * matrix(runif(as.integer(p / h), 1, 2), nrow = as.integer(p / h))
  
  for (i in 1:h) {
    L[1:p, (p + i):(p + i)] <- 0
    L[(p + i):(p + i), 1:p] <- 0
    L[seq(i, p, h), (p + i):(p + i)] <- z
    L[(p + i):(p + i), seq(i, p, h)] <- z
  }
  }
  W <- W * L
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  Theta <- diag(rowSums(W)) - W
  G <- Theta2Gamma(Theta)
  return(list(Gamma = G, graph = Gamma2graph(G), Lst = NULL))
}




mychol <- function(M){
  d <- nrow(M)
  n <- rankMatrix(M)
  if(n==d) R <- chol(M)
  else{
    R1 <- chol(M[1:n, 1:n])
    R <- cbind(R1, solve(t(R1)) %*% M[1:n, (n+1):d])
  }
  return(R)
}






sim_study_latent <- function(d = 5, 
                          n = 100,
                          p = NULL,
                          method = c("maxstable", "mpareto"), 
                          m = 2,
                          h = 2, 
                          lambda_2,
                          gen_model = c("BA", "latent"),
                          reg_method = c("eglearn", "MTP2", "eglatent"),
                          rhostring = "seq(0.01,0.15,length=20)",
                          rng = NULL){
  ## perform a simulation study to measure performance of the EMTP2 block descent algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - n: number of samples.
  ##    - p: threshold probability.
  ##    - gen_method: data generation method.
  ##    - m: the number of edges added in each setp of the Barabasi-Albert model (m=1 is a tree)
  ##    - reg_method: regression method used to estimate the extremal graphical structure.
  ##    - rholist: the list of penality parameters; must be given a string.
  ## Returns:
  ##  a tibble with F scores
  
  # check arguments
  
  rholist <-  eval(parse(text = rhostring))
  
  F1 <- numeric(length(rholist))
  F1_val <- NA
  loglik_max <- NA
  rk_max <- NA
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  
  if(gen_model=="BA"){
    BA_model <- generate_BA_model(d = d, m = m)
    G <- BA_model$G
    g <- BA_model$graph
  }
  if (gen_model == "latent") {
    latent_model <- generate_latent_model(p = d, h)
    G <- 10 * latent_model$Gamma[1:d, 1:d]
    g <- subgraph(latent_model$graph, 1:d)
  }
      
  # perform simulation
  if(method=="maxstable")  X <- rmstable(n=2*n, d=d, model="HR", par=G)
  if(method=="mpareto")  X <- rmpareto(n=2*n, d=d, model="HR", par=G)
  
  X_train <- X[1:n,]
  X_val <- X[(n+1):(2*n),]
  G_emp <- emp_vario(data = X_train, p = p)


  if(reg_method=="MTP2"){
    ptm <- proc.time()[1]
    G_emtp2 <- emtp2(G_emp, tol=1e-6,verbose = FALSE)$G_emtp2 
    time <- proc.time()[1] - ptm
    adj_emtp2 <- (abs(Gamma2Theta(G_emtp2)) >= 1e-4) 
    graph_emtp2 <- igraph::graph_from_adjacency_matrix(adj_emtp2, mode = "undirected", diag = FALSE)
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=graph_emtp2))
    loglik_val <- rep(loglik_HR(data = X_val, p=p, Gamma = G_emtp2, cens = FALSE)[1], times = length(rholist))
  }
  else if(reg_method=="eglearn"){
    ptm <- proc.time()[1]
    fit <- eglearn(data = X_train, p=p, rholist = 7*rholist, reg_method = "ns", complete_Gamma = FALSE)
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph[[i]]))  
    connected_eglearn <- sapply(1:length(rholist), FUN = function(i) ifelse(is_connected(fit$graph[[i]]), 1,0))
    loglik_val <- sapply(1:length(rholist), FUN = function(i) ifelse(connected_eglearn[i]==1, loglik_HR(data = X_val, p=p, Gamma = complete_Gamma(Gamma = G_emp, graph = fit$graph[[i]]), cens = FALSE)[1], NA))
  }
  else if(reg_method=="eglatent"){
    ptm <- proc.time()[1]
    fit_latent <- eglatent(Gamma = G_emp, lam1_list = rholist, lam2_list = lambda_2, refit = TRUE)    
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g = g, gest = fit_latent$graph[[i]]))
    rk <- fit_latent$rk
    loglik_val <- sapply(1:length(rholist), FUN = function(i) {loglik_HR(data = X_val, p=p, Gamma = fit_latent$G_refit[[i]], cens = FALSE)[1]})
  }



  F1_max <- max(F1)
  idx_val_max <- which.max(loglik_val)
  F1_val <- F1[idx_val_max]
  loglik_max <- max(loglik_val, na.rm = TRUE)
  if(reg_method=="eglatent") rk_max <- rk[idx_val_max] 

  
  tbl <- tibble(type = paste0("time"), 
                value = time) %>% 
    bind_rows(tibble(type = paste0("F1_score", 1:length(rholist)),
                     value =  F1)) %>% 
     bind_rows(tibble(type = c("F1_max"),
                      value =  F1_max)) %>% 
    bind_rows(tibble(type = c("F1_val"),
                     value =  F1_val)) %>% 
    bind_rows(tibble(type = c("loglik_max"),
                     value =  loglik_max)) %>% 
    bind_rows(tibble(type = c("rk_max"),
                     value =  rk_max))
  
  return(tbl)
}



F1_score <- function(g, gest) {
  (2*ecount(intersection(gest, g)))/( 2*ecount(intersection(gest, g)) + 
                                        ecount(intersection(complementer(g), gest)) +
                                        ecount(intersection(g, complementer(gest))))
}


wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Ahttps://uniart1.wixsite.com/uni-artrgs:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_fn, args = sim_fn_args[i, ]) %>% #ML: fixed the name of the last argument. It used to be fun_args
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}


#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



generate_BA_model <- function(d,m){
  g <- sample_pa(n=d, m=m, zero.appeal=1,directed=FALSE)
  W <- as_adj(g, sparse=F) * matrix(runif(d^2,2, 5), nrow=d) #matrix(2 + rexp(d^2, rate = 1), nrow=d) #matrix(runif(d^2,2, 5), nrow=d) #  # 
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  O <- diag(rowSums(W)) - W
  G <- Theta2Gamma(O)
  return(list(G = G, graph = Gamma2graph(G)))
}



generate_block_model <- function(ncliques, clique_size, alphad = 1){
  kk <- clique_size
  GG <- matrix(NA, ncliques*(kk-1) + 1, ncliques*(kk-1) + 1)
  for(i in 1:ncliques){
    bigS <- rcorrmatrix(kk, alphad = alphad)
    G1 <- Sigma2Gamma(bigS, full = TRUE)
    if(i==1) GG[1:kk, 1:kk] <- G1
    else GG[(kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1), (kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1)] <- G1
  }
  G <- complete_Gamma(GG)
  round(Gamma2Theta(G),2)
  sum(sum(round(Gamma2Theta(G),2) > 0)) - nrow(G)
  sum(sum(round(Gamma2Theta(G),2) < 0))
  
  return(list(G=G, graph = Gamma2graph(G,to_plot = FALSE)))
}


list2tibble <- function(param_grid, robj, rname){ 
  ## tibble list character -> tibble
  ## convert list to tibble
  
  purrr::map_dfr(seq_len(nrow(param_grid)), function(i){
    tibble(
      k = param_grid[i, ]$k,
      h = param_grid[i, ]$h,
      val = robj[[i]],
      reg_method = rname
    )
  })
}



save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}



lst_methods <- list("eglatent_cv" = "eglatent_cv",
                    "eglatent_oracle" = "eglatent_oracle",
                    "eglearn_cv" = "eglearn_cv",
                    "eglearn_oracle" = "eglearn_oracle",
                    "rk_cv" = "rk_cv",
                    "rk_oracle" = "rk_oracle"
)

my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)



my_palette_methods <- list(
  c("reg_method" = lst_methods$eglatent_cv, "color" = my_palette$red, "fill" = "white"),
  c("reg_method" = lst_methods$eglatent_oracle, "color" = my_palette$blue, "fill" = "white"),
  c("reg_method" = lst_methods$eglearn_cv, "color" = my_palette$green, "fill" = "white"),
  c("reg_method" = lst_methods$eglearn_oracle, "color" = my_palette$yellow, "fill" = "white"),
  c("reg_method" = lst_methods$rk_cv, "color" = my_palette$red, "fill" = "white"),
  c("reg_method" = lst_methods$rk_oracle, "color" = my_palette$blue, "fill" = "white")
) %>% 
  purrr::transpose() %>% 
  as_tibble() %>% 
  unnest(cols = c(reg_method, color, fill)) 

my_col <-  my_palette_methods %>% 
  dplyr::select(reg_method, color) %>% 
  deframe()

my_fill <-  my_palette_methods %>% 
  dplyr::select(reg_method, fill) %>% 
  deframe()


refactor_methods <- function(methods, lst_methods){
  ## character_vector list with mapping -> factor
  ## refactor column with methods
  
  lst_methods <- lst_methods
  
  
  unique_methods <- unique(methods)
  
  new_levels <- names(lst_methods)
  new_labels <- lst_methods %>% unlist() %>% unname()
  
  factor(methods,
         levels = new_levels,
         labels = new_labels)
}


theme_fct <- function(font_size1=11,  font_size2=7.5){
  theme_set(theme_bw() +
              theme(
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.background = element_blank(),
                strip.background = element_rect(fill = "white"),
                plot.caption=element_text(size=font_size2, hjust=0, 
                                          margin=margin(t=15)),
                text = element_text(size = font_size1),
                axis.ticks = element_blank(),
                axis.text = element_text(size = font_size1),
                panel.grid.major = element_line(linewidth = 0.25)
              ) 
            #+
            # theme_cowplot(font_size = 11)
  )
}

theme_fct()


create_palette_levels <- function(reg_method_n, palette_tbl){
  
  str_spl <- strsplit(reg_method_n, "__")
  
  my_tbl <- tibble(
    reg_method = purrr::map_chr(str_spl, function(el){el[1]}),
    level =  purrr::map_chr(str_spl, function(el){el[2]})
  ) %>%
    left_join(palette_tbl, by = "reg_method") %>%
    mutate(fill = if_else(level == "2", color, fill)) %>%
    mutate(reg_method_lev = paste(reg_method, level, sep = "__"))
  
  my_col <-  my_tbl %>%
    dplyr::select(reg_method_lev, color) %>%
    deframe()
  
  my_fill <-  my_tbl %>%
    dplyr::select(reg_method_lev, fill) %>%
    deframe()
  
  
  list(cols = my_col, fills = my_fill)
  
}

texify_column <- function(column, letter) {
  ## vector character -> factor
  ## paste latex formula of the form "$letter = column$"
  factor(column,
         levels = unique(column),
         labels = TeX(paste0(
           "$", letter, " = ",
           unique(column), "$"
         ))
  )
}
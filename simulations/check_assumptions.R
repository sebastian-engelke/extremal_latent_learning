latent_strength<-2
edge_strength<-0.2
p <- 30
h <- 1
set.seed(1)

set.seed(1);return1 <- check_assumptions(p, h,latent_strength,edge_strength,0.005,0.1);
set.seed(1);return2 <- check_assumptions(p, h,latent_strength,edge_strength,0.001,0.1)




library(tidyverse)
library(latex2exp)
library(egg)
library(cowplot)
library(here)
library(ggpubr)

figure_dest_folder <- here("figures/plots/")

source(here("functions/functions_paper.R"))

load(file=here("simulations/data/data_zero_latent.Rdata"))


dat_lik <- bind_rows(
  tibble(method = "gamma=1", val=data$likelihood_vec[[1]], nb_edges = mean(data$sparsity_vec[[1]])),
  tibble(method = "gamma=2", val=data$likelihood_vec[[2]], nb_edges = mean(data$sparsity_vec[[2]])),
  tibble(method = "gamma=3", val=data$likelihood_vec[[3]], nb_edges = mean(data$sparsity_vec[[3]])),
  tibble(method = "gamma=4", val=data$likelihood_vec[[4]], nb_edges = mean(data$sparsity_vec[[4]])),  
  tibble(method = "eglearn", val=data$likelihood_eglearn_vec[[1]], nb_edges = mean(data$sparsity_eglearn_vec[[1]]))
)


dat_lik <- dat_lik %>%
  mutate(method = factor(method, levels = c(
    "gamma=1",
    "gamma=2",
    "gamma=3",
    "gamma=4",
    "eglearn"
  )))

dat_lik_means <- dat_lik %>%
  group_by(method) %>% 
  summarise(mean_nb_edges = mean(nb_edges))


gg_lik <- ggplot(dat_lik) +
  geom_boxplot(aes(x = method, y = val)) +
  xlab("") +
  ylab("Validation log-likelihood") +
  scale_x_discrete(
    breaks = unique(dat_lik$method),
    labels = c(
      TeX("$\\gamma=1$"),
      TeX("$\\gamma=4$"),
      TeX("$\\gamma=8$"),
      TeX("$\\gamma=20$"),
      "eglearn"
    )
  ) +
  geom_text(
    data = dat_lik_means,
    aes(
      x = method, y = -46500, # !!! specify the height for y ticks on top
      label = mean_nb_edges
    ),
    size = 3, color = "black"
  ) +
  coord_cartesian(ylim = c(-52000, -46900), clip = "off") + # !!! specify where plot frame stops
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines"))

gg_lik



dat_F1 <- bind_rows(
  tibble(method = "gamma=1", val=data$F1score_cv_vec[[1]], rank = mean(data$rk_cv_vec[[1]])),
  tibble(method = "gamma=2", val=data$F1score_cv_vec[[2]], rank = mean(data$rk_cv_vec[[2]])),
  tibble(method = "gamma=3", val=data$F1score_cv_vec[[3]], rank = mean(data$rk_cv_vec[[3]])),
  tibble(method = "gamma=4", val=data$F1score_cv_vec[[4]], rank = mean(data$rk_cv_vec[[4]])),  
  tibble(method = "eglearn", val=data$F1score_eglearn_cv_vec[[1]], rank = NaN)
)


dat_F1 <- dat_F1 %>%
  mutate(method = factor(method, levels = c(
    "gamma=1",
    "gamma=2",
    "gamma=3",
    "gamma=4",
    "eglearn"
  )))

dat_F1_means <- dat_F1 %>%
  group_by(method) %>% 
  summarise(mean_rank = mean(rank))


gg_F1 <- ggplot(dat_F1) +
  geom_boxplot(aes(x = method, y = val)) +
  xlab("") +
  ylab(TeX("$F$-score$")) +
  scale_x_discrete(
    breaks = unique(dat_F1$method),
    labels = c(
      TeX("$\\gamma=1$"),
      TeX("$\\gamma=4$"),
      TeX("$\\gamma=8$"),
      TeX("$\\gamma=20$"),
      "eglearn"
    )
  ) +
  geom_text(
    data = dat_F1_means,
    aes(
      x = method, y = 1.05, # !!! specify the height for y ticks on top
      label = mean_rank
    ),
    size = 3, color = "black"
  ) +
  coord_cartesian(ylim = c(.4, 1), clip = "off") + # !!! specify where plot frame stops
  theme(plot.margin = unit(c(2, 1, 1, 1), "lines"))

gg_F1


gg_both <- ggarrange(gg_F1, NULL, gg_lik, nrow = 1, align = "v", widths = c(1, 0.1, 1))

save_myplot(gg_both, plt_nm = here("figures/plots/zero_latent_both.pdf"), width = 11, height = 5)















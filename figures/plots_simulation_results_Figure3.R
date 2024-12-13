library(tidyverse)
library(latex2exp)
library(egg)
library(cowplot)
library(here)
library(ggpubr)


source(here("functions/functions_paper.R"))
load(here("simulations/data/random_graph_data.Rdata"))



F1score_cv<-data$F1score_cv
F1score_oracle<-data$F1score_oracle
rk_oracle <- data$rk_oracle
rk_cv <-data$rk_cv
likelihood_cv <- data$likelihood_cv
F1score_eglearn_cv<-data$F1score_eglearn_cv
F1score_eglearn_oracle <- data$F1score_eglearn_oracle
likelihood_eglearn_cv <- data$likelihood_eglearn_cv


likelihood_cv_diff <- mapply(function(x, y) x-y, likelihood_cv, likelihood_eglearn_cv, SIMPLIFY = FALSE)



param_grid <- expand_grid(k = c(5e3,1e3,2e2), h = c(1, 2, 3))


## F score
dat_fscore <- bind_rows(
  list2tibble(param_grid, F1score_cv, "eglatent_cv"),
  list2tibble(param_grid, F1score_oracle, "eglatent_oracle"),
  list2tibble(param_grid, F1score_eglearn_cv, "eglearn_cv"),
  list2tibble(param_grid, F1score_eglearn_oracle, "eglearn_oracle")
) %>% 
  mutate(h = texify_column(h, "h"))



gg1 <- ggplot(dat_fscore,  aes(x=factor(k), y=val, fill= reg_method, col=reg_method)) + 
  facet_grid(~ h, scales = "free_x", labeller = label_parsed) +
  geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
  xlab("k") +
  ylab("F-score") +
  ylim(0,1) +
  theme(legend.position = c(0.5, -0.23), legend.direction = "horizontal",
        # strip.text.x = element_blank()
        ) +
  labs(color="Method:", fill = "Method:") +
  scale_fill_manual(values = my_fill) +
  scale_colour_manual(values=my_col)
 
  
# save_myplot(gg1, plt_nm = here("figures/plots/simu_fscore.pdf"), width = 3, height = 3)

  
## Rank
dat_rank <- bind_rows(
  list2tibble(param_grid, rk_cv, "rk_cv"),
  list2tibble(param_grid, rk_oracle, "rk_oracle")
) %>%  
  mutate(h = texify_column(h, "h"))


gg2 <- ggplot(dat_rank,  aes(x=factor(k), y=val, fill= reg_method, 
                             col=reg_method)) + 
  facet_grid(~ h, scales = "free_x", labeller = label_parsed) +
  geom_boxplot(alpha=.8, size=.35, outlier.size = .7) +
  xlab("k") +
  ylab("Rank") +
  ylim(0,7) +
  theme(legend.position = c(0.5, -0.23), legend.direction = "horizontal",
        # strip.text.x = element_blank()
  ) +
  labs(color="Method:", fill = "Method:") +
  scale_fill_manual(values = my_fill) +
  scale_colour_manual(values=my_col)





## test likelihood
dat_likelihood <- bind_rows(
  list2tibble(param_grid, likelihood_cv_diff, "eglatent_cv"),
) %>% 
  mutate(h = texify_column(h, "h"))


gg3 <- ggplot(dat_likelihood,  aes(x=factor(k), y=val, fill= reg_method, col=reg_method)) + 
  facet_grid(~ h, scales = "free_x", labeller = label_parsed) +
  geom_boxplot(alpha=.8, size=.35, outlier.size = .7, show.legend = FALSE) +
  xlab("k") +
  ylab("Log-likelihood difference") +
  geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  #ylim(0,1) +
  #theme(legend.position = c(0.5, -0.23), legend.direction = "horizontal",
        # strip.text.x = element_blank()
      #  ) +
  #labs(color="Method:", fill = "Method:") +
  scale_fill_manual(values = "white") +
  scale_colour_manual(values="black")

gg3

  
gg_all <- ggarrange(gg1,NULL,gg2, NULL, gg3, nrow=5, align = "hv",heights = c(1, 0.1, 1, 0.1, 1))

save_myplot(gg_all, plt_nm = here("figures/plots/simu_combined.pdf"), width = 9, height = 10)























library(tibble)

# Multinomial parameters
zeta1 <- c(4.9, -0.5, -0.1)
zeta2 <- c(-5, -0.5, 0.08)
zeta3 <- c(4.6, 0.5, -0.1)
zeta4 <- c(-3.9, 0.5, 0.06)

param_tbl <- tibble(
  setting = c("1", "2", "3", "X"),
  
  theta = list(
    rbind(
      c(-6, 1.5, .07),
      c(-7.18, 1.5, .082),
      c(-9, 1.5, .082),
      c(-8.5, 1.1, .083),
      c(-9.58, 1.15, .085),
      c(-8.86, 1.5, .09)
    ),
    rbind(
      c(-6, 1.5, .07),
      c(-7, 1.5, 0.080),
      c(-7.42, 1.5, .081),
      c(-8.7, 1.1, .088),
      c(-9.1, 1.15, .1),
      c(-8.52, 1.5, .095)
    ),
    rbind(
      c(-7.2, 1.5, .06),
      c(-8.5, 1.5, 0.068),
      c(-8.2, 1.2, .068),
      c(-8.7, 1.1, .072),
      c(-10.6, 1.15, 0.085),
      c(-8.5, 1.5, .072)
    ),
    rbind(
      c(-7.8, 1.5, .06),
      c(-7.4, 1.5, .06),
      c(-7.2, 1.5, .06),
      c(-7.9, 1.55, .06),
      c(-7.4, 1.45, .06),
      c(-9, 1.5, 0.068),
      c(-8.7, 1.2, .068),
      c(-10.2, 1, .072),
      c(-10.6, 1.15, 0.085),
      c(-10.2, 1.5, .068)
    )
  ),
  
  gamma = list(
    rbind(
      c(-5.5, 1.5, .1),
      c(-1.135, 1, .08),
      c(-0.908, 1.5, .09),
      c(-1.33, 1, .09),
      c(-3.8, 1.15, .1),
      c(-2.197, 1.5, .1)
    ),
    rbind(
      c(-5.5, 1.5, .1),
      c(-2.971, 1, .08),
      c(-3.775, 1.5, .081),
      c(-0.718, 1, .09),
      c(-2.818, 1.15, .1),
      c(-1.668, 1.5, .1)
    ),
    rbind(
      c(-5.5, 1.5, .1),
      c(-3.5, 1, .08),
      c(-3.775, 1.5, .081),
      c(-0.718, 1, .09),
      c(-2.818, 1.15, .1),
      c(-1.668, 1.5, .1)
    ),
    rbind(
      c(-5.5, 1.5, .1),
      c(-5.45, 1.5, .1),
      c(-5.45, 1.5, .1),
      c(-5.55, 1.5, .1),
      c(-5.55, 1.5, .1),
      c(-2.971, 1, .08),
      c(-3.775, 1.5, .081),
      c(-0.718, 1, .09),
      c(-2.818, 1.15, .1),
      c(-1.668, 1.5, .1)
    )
  ),
  
  zeta = list(rbind(zeta1, zeta2, zeta3, zeta4))
)

# Tibble that contains all scenarios studied.
dgm_param_tbl <- tibble(
  setting = c("1", "2", "3", "X"),
  # Number of trials in each simulation.
  n_trials = 5
) %>%
  expand_grid(tibble(
    # Case-cohort sampling indicator. If TRUE, the data are generated according to a
    # case-cohort sampling design. If FALSE, the data are generated according to a
    # full-cohort design.
    CC_sampling = c(FALSE, TRUE)
  )) %>%
  expand_grid(tibble(
    # Number of subjects in each trial.
    n_t = c(3000, 8000)
  )) %>%
  # joint with parameter values for each setting.
  left_join(param_tbl, by = "setting")

# For the first two settings, we only consider n_t = 3000 since there are many
# clinical events in those settings. For the third setting, we only consider n_t
# = 8000 since there are few clinical events in that setting.
dgm_param_tbl <- dgm_param_tbl %>%
  filter(
    (setting %in% c("1", "2") & n_t == 3000) |
      (setting == "3" & n_t == 8000) |
      (setting == "X" & n_t %in% c(3000))
  )

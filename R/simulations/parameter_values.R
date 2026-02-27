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
      c(-9, 1, .083),
      c(-9.58, 1.15, .085),
      c(-8.86, 1.5, .09)
    ),
    rbind(
      c(-6, 1.5, .07),
      c(-7, 1.5, 0.080),
      c(-7.42, 1.5, .081),
      c(-9.7, 1, .088),
      c(-9.1, 1.15, .1),
      c(-8.52, 1.5, .095)
    ),
    rbind(
      c(-7.8 1.5, .06),
      c(-9, 1.5, 0.068),
      c(-8.7, 1.2, .068),
      c(-10.2, 1, .072),
      c(-10.6, 1.15, 0.085),
      c(-10.2, 1.5, .068)
    ),
    rbind(
      c(-6, 1.5, .1),
      c(-9.01, 1.5, .082),
      c(-8.99, 1.5, .082),
      c(-8.5, 1, .088),
      c(-9, 1.15, .065),
      c(-10, 1.5, .085)
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
      c(-2.971, 1, .08),
      c(-3.775, 1.5, .081),
      c(-0.718, 1, .09),
      c(-2.818, 1.15, .1),
      c(-1.668, 1.5, .1)
    ),
    rbind(
      c(-5.5, 1.5, .1),
      c(-1, 1, .08),
      c(1.1, 1.5, .081),
      c(-2.9, 1, .09),
      c(-2.5, 1.15, .1),
      c(2, 1.5, .1)
    )
  ),
  
  zeta = list(rbind(zeta1, zeta2, zeta3, zeta4))
)

# Tibble that contains all scenarios studied.
dgm_param_tbl <- tibble(
  setting = c("1", "2", "3"),
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
    n_t = c(2000, 5000)
  )) %>%
  # joint with parameter values for each setting.
  left_join(param_tbl, by = "setting")

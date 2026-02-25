zeta <- c(4.9, -0.5, -0.1) #p5
eta <- c(-5, -0.5, 0.08) #p4
lambda <- c(4.6, 0.5, -0.1) #p3
xi <- c(-3.9, 0.5, 0.06) #p2

# zeta <- c(-0.17, 1, 0) #p5
# eta <- c(1.2, 1.5, -0.025) #p4
# lambda <- c(2.35, 2, -0.05) #p3
# xi <- c(3.02, 3, -.075) #p2

# high event rates, ~zero correlation
if(setting == 1){
  theta_0 <- matrix(c(-6, 1.5, .1), ncol = 1)
  theta_1 <- matrix(c(-8.18, 1.5, .082), ncol = 1)
  theta_2 <- matrix(c(-10, 1.5, .082), ncol = 1)
  theta_3 <- matrix(c(-9.5, 1, .088), ncol = 1)
  theta_4 <- matrix(c(-10.08, 1.15, .09), ncol = 1)
  theta_5 <- matrix(c(-9.36, 1.5, .09), ncol = 1)
  theta.set = list(theta_0, theta_1, theta_2, theta_3, theta_4, theta_5)
  
  gamma_0 <- matrix(c(-5.5, 1.5, .1), ncol = 1)
  gamma_1 <- matrix(c(-1.135, 1, .08), ncol = 1)
  gamma_2 <- matrix(c(-0.908, 1.5, .081), ncol = 1)
  gamma_3 <- matrix(c(-1.33, 1, .09), ncol = 1)
  gamma_4 <- matrix(c(-2.356, 1.15, .1), ncol = 1)
  gamma_5 <- matrix(c(-1.197, 1.5, .1), ncol = 1)
  gamma.set = list(gamma_0, gamma_1, gamma_2, gamma_3, gamma_4, gamma_5)
}

# high event rates, positive correlation rho_s=rho_p=.9
if(setting == 2){
  theta_0 <-      matrix(c(-6, 1.5, .1), ncol = 1)
  theta_1 <-      matrix(c(-6.69, 1.5, .1), ncol = 1)
  theta_2 <-      matrix(c(-6.42, 1.5, .081), ncol = 1)
  theta_3 <-      matrix(c(-9.7, 1, .088), ncol = 1)
  theta_4 <-      matrix(c(-8.1, 1.15, .1), ncol = 1)
  theta_5 <-      matrix(c(-9.52, 1.5, .095), ncol = 1)
  theta.set = list(theta_0, theta_1, theta_2, theta_3, theta_4, theta_5)
  
  gamma_0 <-      matrix(c(-5.5, 1.5, .1), ncol = 1)
  gamma_1 <-      matrix(c(-2.971, 1, .08), ncol = 1)
  gamma_2 <-      matrix(c(-3.775, 1.5, .081), ncol = 1)
  gamma_3 <-      matrix(c(-0.718, 1, .09), ncol = 1)
  gamma_4 <-      matrix(c(-2.818, 1.15, .1), ncol = 1)
  gamma_5 <-      matrix(c(-1.668, 1.5, .1), ncol = 1)
  gamma.set = list(gamma_0, gamma_1, gamma_2, gamma_3, gamma_4, gamma_5)
}

if(setting == 3){
  theta_0 <- matrix(c(-10.38, 1.5, .1), ncol = 1)
  theta_1 <- matrix(c(-12.5, 1.5, .081), ncol = 1)
  theta_5 <- matrix(c(-10.72, 1.3, .081), ncol = 1)
  theta_3 <- matrix(c(-10.82, 1, .09), ncol = 1)
  theta_4 <- matrix(c(-10.9, 1.15, .1), ncol = 1)
  theta_2 <- matrix(c(-10.54, 1.1, .099), ncol = 1)
  theta.set = list(theta_0, theta_1, theta_2, theta_3, theta_4, theta_5)
  
  gamma_0 <- matrix(c(-5.849, 1.5, .1), ncol = 1)
  gamma_1 <- matrix(c(-0.61, 1, .08), ncol = 1)
  gamma_5 <- matrix(c(-1.152, 1.5, .081), ncol = 1)
  gamma_3 <- matrix(c(-2.494, 1, .09), ncol = 1)
  gamma_4 <- matrix(c(-4.404, 1.15, .1), ncol = 1)
  gamma_2 <- matrix(c(-0.715, 1.5, .1), ncol = 1)
  gamma.set = list(gamma_0, gamma_1, gamma_2, gamma_3, gamma_4, gamma_5)
}


# 
if(setting == "X"){
  theta_0 <- matrix(c(-6, 1.5, .1), ncol = 1)
  theta_1 <-
    matrix(c(-9.01, 1.5, .082), ncol = 1)
  theta_2 <-
    matrix(c(-8.99, 1.5, .082), ncol = 1)
  theta_3 <-
    matrix(c(-8.5, 1, .088), ncol = 1)
  theta_4 <-
    matrix(c(-9, 1.15, .09), ncol = 1)
  theta_5 <-
    matrix(c(-10, 1.5, .09), ncol = 1)
  theta.set = list(theta_0, theta_1, theta_2, theta_3, theta_4, theta_5)
  
  gamma_0 <-
    matrix(c(-5.5, 1.5, .1), ncol = 1)
  gamma_1 <-
    matrix(c(-1, 1, .08), ncol = 1)
  gamma_2 <-
    matrix(c(1.1, 1.5, .081), ncol = 1)
  gamma_3 <-
    matrix(c(-2.9, 1, .09), ncol = 1)
  gamma_4 <-
    matrix(c(-2.5, 1.15, .1), ncol = 1)
  gamma_5 <-
    matrix(c(2, 1.5, .1), ncol = 1)
  gamma.set = list(gamma_0, gamma_1, gamma_2, gamma_3, gamma_4, gamma_5)
}
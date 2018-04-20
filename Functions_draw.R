### Code to plot function drawn from a Gaussian Process with mean zero and squared exponential covariance function

# Required libraries:
library(MASS)
library(ggplot2)
library(reshape)
library(dplyr)

# Squared exponential covariance function:
RBF <- function(xi,xj){
  exp(-.5*(xi-xj)^2)
}

# Function for drawing prior functions from the prior
## priors_No is the number of functions to draw
## mean is a logical parameter to check wheter the mean of the draws is shown

plot_priors <- function(priors_No,mean){
set.seed(116687)
  X_vec <- seq(from=-5,to=5,length.out=100)
  n <- length(X_vec)
  
  
  
  KernMat <- outer(X_vec, X_vec , FUN="RBF")
  y <- mvrnorm(n=priors_No, mu=rep(0,n), Sigma=KernMat, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  if (priors_No > 1) y <- t(y)
  y <- data.frame(y)
  if (mean==T){
    
    y$mean <- apply(y,1,mean)
    functions_prior <- data.frame(X=X_vec,y)  %>% melt(.,id='X', variable_name='Function')  
    ggplot(functions_prior,aes(x=X,y=value,color=Function)) + geom_line() + 
      stat_summary(fun.y = "mean", colour = "red", size = 2, geom = "point") +
      theme(plot.title = element_text(hjust = 0.5),legend.position="none") + ggtitle('Prior Distribution Draws') + ylab('f')
  } else {
    functions_prior <- data.frame(X=X_vec,y)  %>% melt(.,id='X', variable_name='Function')  
    ggplot(functions_prior,aes(x=X,y=value,color=Function)) + geom_line() +
      theme(plot.title = element_text(hjust = 0.5)) + ggtitle('Prior Distribution Draws') + ylab('f')
  }
}



# Function for plotting functions drawn from the posterior
## post_No is the number of functions to draw
## X_obs is the values of x for which there are observations
# s2 is the variance of the error distribution 
        #if s2=0 the observations are noise-less, s2>0 the functions will be shown as points instead of lines
## Note: the true function is sin(-X_obs)+cos(-X_obs) 
plot_posterior <- function(post_No,X_obs,s2){
  set.seed(116687)
  X_vec <- seq(from=-5,to=5,by=.1)
  n <- length(X_vec)
  noise <- mvrnorm(n=1, mu=rep(0,n), Sigma=s2*diag(n), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  y_true <- sin(-X_obs)+cos(-X_obs)  
  y_obs <- y_true + noise[which(X_vec %in% round(X_obs,digits=1))]
    
    
  KernMat <- outer(X_vec, X_vec , FUN="RBF") + s2*diag(n)
  KernX_obs <- outer(X_obs, X_obs , FUN="RBF")
  KernX_vecX_obs <- outer(X_vec, X_obs , FUN="RBF")
  KernX_obsX_vec <- outer(X_obs, X_vec , FUN="RBF")
  KernX_obs_inv <- solve(KernX_obs)
  
  if (s2 == 0 ){
    mu_post <- KernX_vecX_obs%*%KernX_obs_inv%*%y_true
  } else {
    mu_post <- KernX_vecX_obs%*%KernX_obs_inv%*%y_obs
  }
  Sigma_post <- KernMat-KernX_vecX_obs%*%KernX_obs_inv%*%KernX_obsX_vec
  
  post_f <- mvrnorm(n=post_No, mu=mu_post, Sigma=Sigma_post, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  if (post_No > 1) post_f <- t(post_f)
  functions_posterior <- data.frame(X=X_vec,f=post_f) %>% melt(.,id='X', variable_name='Function')  
  
  if (s2 == 0){
    ggplot(functions_posterior,aes(x=X,y=value,color=Function)) + geom_line() +
      stat_summary(fun.y = "mean", colour = "red", size = 1.5, geom = "point") +
      theme(plot.title = element_text(hjust = 0.5),legend.position="none") + ggtitle('Posterior Distribution Draws') +
      geom_point(data=data.frame(X=X_obs,y=y_true), aes(x=X, y=y), colour="black", size=3) + ylab('f')
  } else {
    ggplot(functions_posterior,aes(x=X,y=value,color=Function)) + geom_point() +
      stat_summary(fun.y = "mean", colour = "red", size = 3, geom = "point") +
      theme(plot.title = element_text(hjust = 0.5),legend.position="none") + ggtitle('Posterior Distribution Draws') +
      geom_line(data=data.frame(X=X_vec,y=sin(-X_vec)+cos(-X_vec)), aes(x=X, y=y), colour="black", size=1) + ylab('f') +
      geom_point(data=data.frame(X=X_obs,y=y_obs), aes(x=X, y=y), colour="black", size=3)
  }
}

# Draws functions from prior
plot_priors(priors_No=3,mean=F)
plot_priors(priors_No=100,mean=T)
# Draws the covariance between f(0) and f(x)
xseq <- seq(from=-5,to=5,length.out=100)
data.frame(X=xseq,Cov=RBF(0,xseq)) %>% ggplot(.,aes(x=X,y=Cov)) + geom_line() + ylab(expression(Cov(f(0),f(x)))) +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none") + ggtitle(expression(paste('Covariance Decay Around ',x=0)))


# Draws functions from posterior:
# Noise-free observations
s2=0
plot_posterior(post_No=30,X_obs = c(-.5),s2)
plot_posterior(post_No=30,X_obs = c(-4,-.5),s2)
plot_posterior(post_No=30,X_obs = c(-4,-3,-.5),s2)
plot_posterior(post_No=30,X_obs = c(-4,-3,-.5,3.5),s2)
plot_posterior(post_No=30,X_obs = c(-4,-3,-2.5,-1,-.5,.5,3.5,4),s2)

# Noisy observations
s2=0.01
plot_posterior(post_No=3,X_obs = c(-.5),s2)
plot_posterior(post_No=3,X_obs = c(-4,-.5),s2)
plot_posterior(post_No=3,X_obs = c(-4,-3,-.5),s2)
plot_posterior(post_No=3,X_obs = c(-4,-3,-.5,3.5),s2)
plot_posterior(post_No=3,X_obs = c(-4,-3,-2.5,-1,-.5,.5,3.5,4),s2)
plot_posterior(post_No=30,X_obs = c(-4,-3,-2.5,-1,-.5,.5,3.5,4),s2)

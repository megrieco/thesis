#################################################
### Example R Code for Data Simulation        ###
### Megan Grieco                              ###
### 2023/06/28                    					  ###
#################################################

###################################
### Load Libraries / Functions  ###
###################################
library(MASS)
library(tidyverse)
library(broom)
library(janitor)
library(corrplot)
library(gWQS)
library(qgcomp)
library(bkmr)
library(bkmrhat)
library(coda)

#function to convert correlation matrix to covariance matrix using standard deviations
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}

##################################
### Generate Covariance Matrix ###
##################################
# Note: valid matrix always given by constant or exponential decay 

set.seed(1111)  #Set seed for reproduceability
nexp=10         #Number of exposures to simulate

#Make correlation matrix of weakly/positively correlated (0.1-0.35)
R <- matrix(runif(nexp^2,min=0.5,max=0.75), ncol=nexp) 
R[lower.tri(R)] = t(R)[lower.tri(R)]
#Set diagonals of matrix = 1
diag(R)<-1
#vector of standard deviations = all 1
S <- rep(1,10)

#Calculate covariance matrix based on correlations and standard deviations
sample_covariance_matrix <- cor2cov(R,S)

#Mean value for each exposure pulled from normal distribution
sample_means <- rnorm(nexp,5,1)

##########################
### Generate Datasets  ###
##########################

datasets <- list()
N=100          #Define total number of simulated datasets
n=1000           #Define number of observations per dataset

#Set beta coefficients (weights)
betas <- data.frame(b1=.1, 
                    b2=.25,
                    b3=.65,
                    b4=0,
                    b5=0,
                    b6=0,
                    b7=0,
                    b8=0,
                    b9=0,
                    b10=0)
#Note: for WQS - need to make sure betas sum to 1


#loop through to create N simulated datasets
for (i in 1:N) {
  # Generate exposures (Multivariate Normal)
  set.seed(i)
  simulation <- mvrnorm(n = n,mu = sample_means, Sigma = sample_covariance_matrix)
  #Sample exposures from multivariate normal distribution
  
  # exponentiate to get Multivariate log-normal distribution
  x <- simulation %>% as.data.frame(.) %>% exp(.)
  #change names to be x1, ..., x10
  names(x) <- gsub(x=names(x),pattern="V",replacement="x")
  
  # Generate outcome variable
  df <- x %>% mutate(
    mu=betas$b1*x1 + betas$b2*x2 + betas$b3*x3 + betas$b4*x4 + betas$b5*x5 
    + betas$b6*x6 + betas$b7*x7 + betas$b8*x8 + betas$b9*x9 + betas$b10*x10)
  
  df$y = rnorm(n,df$mu,10)  # Y_i ~ N(mu_i, sd=10)
  
  ### save whole dataset
  datasets[[i]] <- df
}

##########################
### Run MLR            ###
##########################

#output for beta estimates
beta_est <- as.data.frame(matrix(nrow=0,ncol=2))
names(beta_est) <- c("exp","estimate")

#loop through N datasets
for (i in 1:N) {
  df <- datasets[[i]]
  
  model <- df %>% dplyr::select(-mu) %>% 
    #Run MLR (Y ~ x1 + x2 + ...)
    
    lm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10 , data = .) %>%
    summary()
  
  #retain estimates for each beta coefficient
  beta1s <- model$coefficients[-1,1] %>% as.data.frame() %>% rownames_to_column()
  names(beta1s) <- c("exp","estimate")
  
  #append to list of beta estimates for all datasets
  beta_est <- rbind(beta_est,beta1s)
  
}


beta_est %>% 
  #convert to numeric columns
  mutate(estimate=as.numeric(estimate)) %>% 
  #Take average of betas for each exposure 
  group_by(exp) %>% summarize(mean=mean(estimate),
                              sd=sd(estimate)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3)
#     exp mean    sd lower upper
# 1   x1  0.1 0.008 0.100 0.101
# 2  x10  1.0 0.013 0.999 1.001
# 3   x2  0.2 0.008 0.199 0.200
# 4   x3  0.3 0.006 0.299 0.300
# 5   x4  0.4 0.017 0.399 0.401
# 6   x5  0.5 0.003 0.500 0.500
# 7   x6  0.6 0.004 0.600 0.600
# 8   x7  0.7 0.008 0.700 0.701
# 9   x8  0.8 0.005 0.799 0.800
# 10  x9  0.9 0.001 0.900 0.900

##############################
### Generate WQS Datasets  ###
##############################

datasets <- list()
N=100           #Define total number of simulated datasets
n=1000          #Define number of observations per dataset
q=10            #Define number of quantiles for WQS
betas=c(0.1,0.25,0.65,rep(0,7))
#Set beta coefficients (weights)
betas <- data.frame(b1=0.1, 
                    b2=0.25,
                    b3=0.65,
                    b4=0,
                    b5=0,
                    b6=0,
                    b7=0,
                    b8=0,
                    b9=0,
                    b10=0)
#Note: for WQS - need to make sure betas sum to 1


#loop through to create N simulated datasets
for (i in 1:N) {
  # Generate exposures (Multivariate Normal)
  set.seed(i)
  simulation <- mvrnorm(n = n,mu = sample_means, Sigma = sample_covariance_matrix)
  #Sample exposures from multivariate normal distribution
  
  # exponentiate to get Multivariate log-normal distribution
  x <- simulation %>% as.data.frame(.) %>% exp(.)
  #change names to be x1, ..., x10
  names(x) <- gsub(x=names(x),pattern="V",replacement="x")
  
  ### transform x's into deciles before calculating y
  trans_decile <- function (x) ntile(x, q)
  x_trans <- x
  x_trans <- plyr::colwise(trans_decile)(x_trans)
  
  # Generate outcome variable
  ## can also add scaling factor to mu to change overall effect
  x_trans <- x_trans %>% mutate(
    mu=betas$b1*x1 + betas$b2*x2 + betas$b3*x3 + betas$b4*x4 + betas$b5*x5 
    + betas$b6*x6 + betas$b7*x7 + betas$b8*x8 + betas$b9*x9 + betas$b10*x10)
  
  df <- cbind(x, mu=x_trans$mu)
  df$y = rnorm(n,df$mu,1)  # Y_i ~ N(mu_i, sd=10)
  
  ### save whole dataset
  datasets[[i]] <- df
}

##########################
### Run WQS            ###
##########################

Xs <- names(x)

#output for WQS weights
weights_all <- as.data.frame(matrix(nrow=0,ncol=2))
names(weights_all) <- c("exp","estimate")

#loop through N datasets
for (i in 1:N) {
  df <- datasets[[i]]
  
  # we run the model and save the results in the variable "results"
  results <- gwqs(y ~ wqs, mix_name = Xs, data = df, 
                  q = q, validation = 0.6, b = 100, b1_pos = T, 
                  b1_constr = T, family = "gaussian", seed = 2016)
                  #Error if I set b1_constr=F: There are no positive b1 in the bootstrapped models
  
  df_weights <- results$final_weights %>% rename(exp=mix_name, estimate=mean_weight) %>% remove_rownames()
  
  #add overall estimate
  df_weights$exp <- as.character(df_weights$exp)
  df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
  
  #combine with beta estimates of other datasets
  weights_all <- rbind(weights_all,df_weights)
  
}


weights_all %>% 
  #convert to numeric columns
  mutate(weight=as.numeric(estimate)) %>% 
  #Take average of betas for each x 
  group_by(exp) %>% summarize(mean=mean(weight),
                              sd=sd(weight)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  arrange(-mean)

#With B1=0.3, B2=0.7, everything else=0
#        exp   mean    sd  lower  upper
# 1  Overall 1.031 0.027 1.026 1.037 N
# 2       x2 0.667 0.022 0.663 0.671 N
# 3       x1 0.278 0.017 0.274 0.281 N
# 4      x10 0.007 0.008 0.006 0.009
# 5       x3 0.007 0.008 0.005 0.008
# 6       x4 0.007 0.009 0.005 0.009
# 7       x5 0.007 0.008 0.005 0.009
# 8       x6 0.007 0.008 0.006 0.009
# 9       x7 0.007 0.008 0.006 0.009
# 10      x8 0.007 0.009 0.005 0.009
# 11      x9 0.006 0.006 0.005 0.007

#With B1=0.1, B2=0.25, B3=0.65, everything else=0
#        exp   mean    sd  lower  upper
# 1  Overall 1.022 0.029 1.016 1.028 N
# 2       x3 0.625 0.023 0.620 0.629 N
# 3       x2 0.234 0.018 0.231 0.238 N
# 4       x1 0.085 0.019 0.081 0.089 N
# 5       x7 0.009 0.009 0.008 0.011
# 6       x4 0.008 0.009 0.006 0.010
# 7       x5 0.008 0.009 0.007 0.010
# 8       x6 0.008 0.009 0.007 0.010
# 9       x8 0.008 0.009 0.006 0.010
# 10     x10 0.007 0.008 0.006 0.009
# 11      x9 0.007 0.007 0.005 0.008

#With B1=0.1, B2=0.25, B3=0.65, everything else=0 AND increasing correlation (0.5-0.75)
#        exp   mean    sd  lower  upper
# 1  Overall 1.002 0.020 0.999 1.006
# 2       x3 0.620 0.029 0.614 0.626
# 3       x2 0.221 0.028 0.215 0.227
# 4       x1 0.063 0.026 0.058 0.068
# 5       x7 0.020 0.017 0.017 0.024
# 6       x5 0.018 0.019 0.015 0.022
# 7       x4 0.014 0.016 0.011 0.018
# 8       x8 0.014 0.015 0.011 0.017
# 9       x6 0.013 0.014 0.011 0.016
# 10      x9 0.010 0.010 0.007 0.012
# 11     x10 0.006 0.010 0.005 0.008

#With betas = seq(1:10)/sum(seq(1:10))
#        exp   mean    sd  lower  upper
# 1  Overall 0.991 0.028 0.985 0.996
# 2      x10 0.183 0.020 0.179 0.187
# 3       x9 0.161 0.018 0.157 0.165
# 4       x8 0.145 0.021 0.140 0.149
# 5       x7 0.127 0.021 0.123 0.132
# 6       x6 0.109 0.021 0.105 0.113
# 7       x5 0.088 0.020 0.084 0.092
# 8       x4 0.069 0.021 0.064 0.073
# 9       x3 0.057 0.021 0.053 0.061
# 10      x2 0.039 0.018 0.036 0.043
# 11      x1 0.023 0.016 0.020 0.026

#With betas = seq(1:10)/sum(seq(1:10)) AND increasing correlation (0.2-0.45)
#        exp   mean    sd  lower  upper
# 1  Overall 0.993 0.024 0.988 0.998
# 2      x10 0.183 0.021 0.179 0.188
# 3       x9 0.160 0.021 0.156 0.164
# 4       x8 0.144 0.022 0.140 0.149
# 5       x7 0.126 0.020 0.122 0.130
# 6       x6 0.107 0.021 0.103 0.112
# 7       x5 0.089 0.024 0.084 0.093
# 8       x4 0.070 0.022 0.065 0.074
# 9       x3 0.058 0.021 0.054 0.062
# 10      x2 0.039 0.021 0.035 0.043
# 11      x1 0.023 0.017 0.020 0.027


#With betas = seq(1:10)/sum(seq(1:10)) AND increasing correlation (0.3-0.55)
#        exp   mean    sd  lower  upper
# 1  Overall 0.994 0.022 0.990 0.998
# 2      x10 0.183 0.023 0.179 0.188
# 3       x9 0.159 0.023 0.155 0.164
# 4       x8 0.144 0.025 0.139 0.149
# 5       x7 0.126 0.023 0.121 0.131
# 6       x6 0.107 0.023 0.102 0.112
# 7       x5 0.088 0.025 0.083 0.093
# 8       x4 0.069 0.024 0.064 0.074
# 9       x3 0.059 0.022 0.054 0.063
# 10      x2 0.041 0.023 0.036 0.045
# 11      x1 0.025 0.020 0.021 0.028

#With betas = seq(1:10)/sum(seq(1:10)) AND increasing correlation (0.4-0.65)
#        exp   mean    sd  lower  upper
# 1  Overall 0.996 0.020 0.992 1.000 Y
# 2      x10 0.185 0.023 0.180 0.190 Y
# 3       x9 0.157 0.029 0.152 0.163 N
# 4       x8 0.144 0.028 0.138 0.150 Y
# 5       x7 0.123 0.027 0.118 0.128 Y
# 6       x6 0.107 0.025 0.102 0.112 Y
# 7       x5 0.088 0.031 0.082 0.094 Y
# 8       x4 0.067 0.026 0.062 0.072 N
# 9       x3 0.061 0.024 0.056 0.065 N
# 10      x2 0.042 0.024 0.037 0.046 N
# 11      x1 0.027 0.022 0.022 0.031 N

#With betas = seq(1:10)/sum(seq(1:10)) AND increasing correlation (0.5-0.75)
#        exp   mean    sd  lower  upper
# 1  Overall 0.996 0.019 0.993 1.000 Y
# 2      x10 0.186 0.031 0.180 0.192 Y
# 3       x9 0.157 0.031 0.151 0.163 N
# 4       x8 0.141 0.033 0.135 0.148 Y
# 5       x7 0.119 0.030 0.113 0.125 N
# 6       x6 0.105 0.030 0.099 0.111 Y
# 7       x5 0.083 0.036 0.076 0.090 Y
# 8       x4 0.067 0.030 0.061 0.073 Y
# 9       x3 0.065 0.029 0.059 0.071 N
# 10      x2 0.045 0.028 0.040 0.051 N
# 11      x1 0.033 0.027 0.027 0.038 N

##########################
### Run qgcomp         ###
##########################

#output for qgcomp weights
weights_all <- as.data.frame(matrix(nrow=0,ncol=2))
names(weights_all) <- c("exp","estimate")

#loop through N datasets
for (i in 1:N) {
  df <- datasets[[i]]
  
  # we run the model and save the results in the variable "results"
  results <- qgcomp.boot(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                         expnms=Xs,
                         df, family=gaussian(), q=q, B=100, seed=2016)
  #Error if I set b1_constr=F: There are no positive b1 in the bootstrapped models
  
  df_weights <- results[["fit"]][["coefficients"]][-1] %>% as.data.frame() %>% rownames_to_column(var = "exp") 
  names(df_weights) <- c("exp","estimate")
  
  #add overall estimate
  df_weights$exp <- as.character(df_weights$exp)
  df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
  
  #combine with beta estimates of other datasets
  weights_all <- rbind(weights_all,df_weights)
  
}


weights_all %>% 
  #convert to numeric columns
  mutate(weight=as.numeric(estimate)) %>% 
  #Take average of betas for each x 
  group_by(exp) %>% summarize(mean=mean(weight),
                              sd=sd(weight)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3) %>% arrange(-mean)


#With betas = seq(1:10)/sum(seq(1:10)), standard correlation (0.1-0.35)
#        exp   mean    sd  lower  upper
# 1  Overall 0.998 0.024 0.994 1.003 Y
# 2      x10 0.183 0.014 0.180 0.185 Y
# 3       x9 0.162 0.013 0.160 0.165 Y
# 4       x8 0.146 0.013 0.143 0.148 Y
# 5       x7 0.126 0.014 0.123 0.129 Y
# 6       x6 0.110 0.011 0.107 0.112 Y
# 7       x5 0.090 0.013 0.088 0.093 Y
# 8       x4 0.072 0.013 0.069 0.074 Y
# 9       x3 0.056 0.013 0.053 0.058 Y
# 10      x2 0.038 0.012 0.035 0.040 Y
# 11      x1 0.017 0.014 0.014 0.020 Y

##############################
### Generate BKMR Datasets ###
##############################

nexp=5         #Number of exposures to simulate

#Make correlation matrix of weakly/positively correlated (0.1-0.35)
R <- matrix(runif(nexp^2,min=0.1,max=0.35), ncol=nexp) 
R[lower.tri(R)] = t(R)[lower.tri(R)]
#Set diagonals of matrix = 1
diag(R)<-1
#vector of standard deviations = all 1
S <- rep(1,nexp)

#Calculate covariance matrix based on correlations and standard deviations
sample_covariance_matrix <- cor2cov(R,S)

#Mean value for each exposure pulled from normal distribution
sample_means <- rnorm(nexp,5,1)

datasets <- list()
N=100           #Define total number of simulated datasets
n=1000          #Define number of observations per dataset

#Set beta coefficients (weights)
## try with 2 weights to start, plot two variables to check association
## later on can make a b6 that is an interaction term
betas <- data.frame(b1=.1, 
                    b2=.25,
                    b3=.65,
                    b4=0,
                    b5=0)
#Note: for WQS - need to make sure betas sum to 1


#loop through to create N simulated datasets
for (i in 1:N) {
  # Generate exposures (Multivariate Normal)
  set.seed(i)
  simulation <- mvrnorm(n = n,mu = sample_means, Sigma = sample_covariance_matrix)
  #Sample exposures from multivariate normal distribution
  
  # exponentiate to get Multivariate log-normal distribution
  x <- simulation %>% as.data.frame(.) %>% exp(.)
  #change names to be x1, ..., x10
  names(x) <- gsub(x=names(x),pattern="V",replacement="x")

  # Generate outcome variable
  df <- x %>% mutate(mu=as.matrix(x) %*% t(betas))
  df$y = rnorm(n,df$mu,10)  # Y_i ~ N(mu_i, sd=10)
  
  ### save whole dataset
  datasets[[i]] <- df
}

##########################
### Run BKMR           ###
##########################

#output for BKMR PIPs
pips_all <- as.data.frame(matrix(nrow=0,ncol=2))
names(pips_all) <- c("variable","PIP")

#loop through N datasets
for (i in 1:N) {
  print(paste0("Running BKMR on dataset ",i))
  df <- datasets[[i]]

  y <- df$y
  Z <-df[,1:nexp]
  
  set.seed(111)
  
  #varsel = T: perform variable selection
  fitkm <- kmbayes(y = y, Z = Z, iter = 100, verbose = FALSE, varsel = T)
  
  #Estimated posterior inclusion probabilities
  pips <- ExtractPIPs(fitkm)
  
  #combine with PIP estimates of other datasets
  pips_all <- rbind(pips_all,pips)
  
}

pips_all %>% 
  #convert to numeric columns
  mutate(PIP=as.numeric(PIP)) %>% 
  #Take average of betas for each x 
  group_by(variable) %>% summarize(mean=mean(PIP),
                              sd=sd(PIP)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3) %>% arrange(-mean)


#with B1=3, B2=7, B3=B4=B5=0
#   variable  mean    sd lower upper
# 1       x2 0.611 0.166 0.579 0.644
# 2       x5 0.551 0.117 0.528 0.574
# 3       x4 0.451 0.101 0.431 0.471
# 4       x1 0.379 0.121 0.355 0.403
# 5       x3 0.340 0.117 0.317 0.364

#with B1=.1, B2=.25, B3=0.65, B4=B5=0 (100 iterations)
#   variable  mean    sd lower upper
# variable  mean    sd lower upper
# 1       x3 0.736 0.422 0.652 0.820
# 2       x4 0.441 0.240 0.393 0.489
# 3       x2 0.249 0.375 0.175 0.324
# 4       x5 0.219 0.270 0.166 0.273
# 5       x1 0.097 0.065 0.084 0.110

#with B1=.3, B2=.7, B3=B4=B5=0 (only 10 iterations)
#   variable  mean    sd lower upper
# 1       x2 0.618 0.155 0.587 0.648
# 2       x5 0.547 0.108 0.525 0.568
# 3       x4 0.453 0.100 0.433 0.472
# 4       x1 0.374 0.110 0.352 0.396
# 5       x3 0.340 0.117 0.317 0.364

#with B1=.1, B2=.25, B3=.65, B4=B5=0 (only 50 iterations)
#   variable  mean    sd lower upper
# 1       x2 0.986 0.076 0.971 1.001
# 2       x4 0.905 0.248 0.856 0.954
# 3       x5 0.785 0.399 0.706 0.864
# 4       x3 0.710 0.454 0.620 0.800
# 5       x1 0.698 0.460 0.607 0.789

pips_all %>% mutate(d=rep(1:N, times=1, each=nexp)) %>%
  pivot_wider(id_cols="d",names_from = "variable",values_from = "PIP") %>%
  ggplot(aes(x=x1,y=x2)) +
  geom_point(alpha=0.5) +
  geom_smooth(se=F)
  labs(x= "X1", y="X2")

############################################################  
#### Without any correlation                          ######
############################################################
#Note: no significant difference from 0.1-0.35 range of weakly
#     correlated exposures

#Make Sigma matrix of 0's
R <- matrix(rep(0,nexp^2), ncol=nexp) 
R[lower.tri(R)] = t(R)[lower.tri(R)]
#Set diagonals of matrix = 1
diag(R)<-1

#Calculate covariance matrix based on correlations and standard deviations
sample_covariance_matrix <- cor2cov(R,S)

##############################
### Generate WQS Datasets  ###
##############################

datasets <- list()
N=100          #Define total number of simulated datasets
n=1000          #Define number of observations per dataset

#Set beta coefficients (weights)
betas <- data.frame(b1=0.01818182, 
                    b2=0.0363636,
                    b3=0.05454545,
                    b4=0.07272727,
                    b5=0.09090909,
                    b6=0.10909091,
                    b7=0.12727273,
                    b8=0.14545455,
                    b9=0.16363636,
                    b10=0.18181818)
#Note: for WQS - need to make sure betas sum to 1


#loop through to create N simulated datasets
for (i in 1:N) {
  # Generate exposures (Multivariate Normal)
  set.seed(i)
  simulation <- mvrnorm(n = n,mu = sample_means, Sigma = sample_covariance_matrix)
  #Sample exposures from multivariate normal distribution
  
  # exponentiate to get Multivariate log-normal distribution
  x <- simulation %>% as.data.frame(.) %>% exp(.)
  #change names to be x1, ..., x10
  names(x) <- gsub(x=names(x),pattern="V",replacement="x")
  
  # Generate outcome variable
  df <- x %>% mutate(
    mu=betas$b1*x1 + betas$b2*x2 + betas$b3*x3 + betas$b4*x4 + betas$b5*x5 
    + betas$b6*x6 + betas$b7*x7 + betas$b8*x8 + betas$b9*x9 + betas$b10*x10)
  
  df$y = rnorm(n,df$mu,10)  # Y_i ~ N(mu_i, sd=10)
  
  ### save whole dataset
  datasets[[i]] <- df
}

##########################
### Run WQS            ###
##########################

Xs <- names(x)

#output for WQS weights
weights_all <- as.data.frame(matrix(nrow=0,ncol=2))
names(weights_all) <- c("exp","estimate")

#loop through N datasets
for (i in 1:N) {
  df <- datasets[[i]]
  
  # we run the model and save the results in the variable "results"
  results <- gwqs(y ~ wqs, mix_name = Xs, data = df, 
                  q = 10, validation = 0.6, b = 100, b1_pos = T, 
                  b1_constr = T, family = "gaussian", seed = 2016)
  #Error if I set b1_constr=F: There are no positive b1 in the bootstrapped models
  
  df_weights <- results$final_weights %>% rename(exp=mix_name, estimate=mean_weight) %>% remove_rownames()
  
  #add overall estimate
  df_weights$exp <- as.character(df_weights$exp)
  df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
  
  #combine with beta estimates of other datasets
  weights_all <- rbind(weights_all,df_weights)
  
}


weights_all %>% 
  #convert to numeric columns
  mutate(weight=as.numeric(estimate)) %>% 
  #Take average of betas for each x 
  group_by(exp) %>% summarize(mean=mean(weight),
                              sd=sd(weight)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  arrange(mean)

#With B1=0.3, B2=0.7, everything else=0
#        exp   mean    sd  lower  upper
# 1  Overall 46.161 3.840 45.922 46.399
# 2       x1  0.239 0.034  0.237  0.241
# 3      x10  0.017 0.017  0.016  0.018
# 4       x2  0.621 0.041  0.618  0.623
# 5       x3  0.017 0.017  0.016  0.018
# 6       x4  0.018 0.017  0.017  0.019
# 7       x5  0.017 0.017  0.016  0.018
# 8       x6  0.018 0.017  0.017  0.019
# 9       x7  0.017 0.017  0.016  0.018
# 10      x8  0.017 0.017  0.016  0.018
# 11      x9  0.018 0.017  0.017  0.019

#With betas = seq(1:10)/sum(seq(1:10))
#        exp   mean    sd  lower  upper
# 1       x4   0.024 0.018   0.020   0.027
# 2       x1   0.028 0.019   0.024   0.032
# 3       x2   0.031 0.019   0.027   0.034
# 4       x7   0.048 0.023   0.043   0.052
# 5       x3   0.049 0.022   0.045   0.054
# 6      x10   0.060 0.023   0.056   0.065
# 7       x8   0.086 0.024   0.082   0.091
# 8       x6   0.091 0.024   0.086   0.096
# 9       x5   0.124 0.023   0.119   0.128
# 10      x9   0.460 0.029   0.454   0.466
# 11 Overall 109.204 6.407 107.933 110.475

##########################
### Run qgcomp         ###
##########################

#output for qgcomp weights
weights_all <- as.data.frame(matrix(nrow=0,ncol=2))
names(weights_all) <- c("exp","estimate")

#loop through N datasets
for (i in 1:N) {
  df <- datasets[[i]]
  
  # we run the model and save the results in the variable "results"
  results <- qgcomp.boot(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                         expnms=Xs,
                         df, family=gaussian(), q=10, B=100, seed=2016)
  
  df_weights <- results[["fit"]][["coefficients"]][-1] %>% as.data.frame() %>% rownames_to_column(var = "exp") 
  names(df_weights) <- c("exp","estimate")
  
  #add overall estimate
  df_weights$exp <- as.character(df_weights$exp)
  df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
  
  #combine with beta estimates of other datasets
  weights_all <- rbind(weights_all,df_weights)
  
}


weights_all %>% 
  #convert to numeric columns
  mutate(weight=as.numeric(estimate)) %>% 
  #Take average of betas for each x 
  group_by(exp) %>% summarize(mean=mean(weight),
                              sd=sd(weight)) %>% 
  mutate(lower=mean - (qt(0.975,df=N-1)*sd/sqrt(N)),
         upper=mean + (qt(0.975,df=N-1)*sd/sqrt(N))) %>%
  as.data.frame() %>% 
  mutate_if(is.numeric, round, digits=3) 
#        exp   mean    sd  lower  upper
# 1  Overall 44.152 3.055 43.962 44.341
# 2       x1 10.094 1.308 10.013 10.175
# 3      x10  0.161 1.177  0.088  0.234
# 4       x2 27.863 1.637 27.762 27.965
# 5       x3 -0.046 1.185 -0.119  0.028
# 6       x4  1.426 1.120  1.356  1.495
# 7       x5  1.158 1.210  1.083  1.234
# 8       x6  1.075 1.116  1.006  1.144
# 9       x7  0.802 1.156  0.731  0.874
# 10      x8  0.909 1.186  0.836  0.983
# 11      x9  0.708 1.151  0.637  0.780


#################################################
### BKMR Simulation for running in background ###
### Megan Grieco                              ###
### 2023/07/21                    					  ###
#################################################

###################################
### Load Libraries / Functions  ###
###################################
library(MASS)
library(tidyverse)
library(broom)
library(janitor)
library(corrplot)

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
betas <- data.frame(b1=.3, 
                    b2=.7,
                    b3=0,
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
  fitkm <- kmbayes(y = y, Z = Z, iter = 1000, verbose = FALSE, varsel = T)
  
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


#with B1=3, B2=7, B3=B4=B5=0 (100 iterations)
#   variable  mean    sd lower upper
# 1       x2 0.611 0.166 0.579 0.644
# 2       x5 0.551 0.117 0.528 0.574
# 3       x4 0.451 0.101 0.431 0.471
# 4       x1 0.379 0.121 0.355 0.403
# 5       x3 0.340 0.117 0.317 0.364

#with B1=3, B2=7, B3=B4=B5=0 (10000 iterations)
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

# pips_all %>% mutate(d=rep(1:N, times=1, each=nexp)) %>%
#   pivot_wider(id_cols="d",names_from = "variable",values_from = "PIP") %>%
#   ggplot(aes(x=x1,y=x2)) +
#   geom_point(alpha=0.5) +
#   geom_smooth(se=F)
# labs(x= "X1", y="X2")
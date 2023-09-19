source("~/OneDrive - Emory University/Courses/Thesis/simulate_power_typeI_source.R") #functions

### Set simulation parameters
nexp <- 10 #number of exposures
sample_means <- rnorm(nexp,5,1) #Mean value for each exposure pulled from normal distribution
S <- rep(1,nexp) #vector of standard deviations = all 1


effect=1      #Define overall effect size
N=100           #Define total number of simulated datasets
n=1000          #Define number of observations per dataset
q=4            #Define number of quantiles for WQS
betas=c(0.1,0.25,0.65,rep(0,7)) #define weights of exposures
covar_betas <- rep(-1,5) # define covariate effects
  #parity -> more births decreases gestational birth age
  #bmi ->  bmi decreases gestational birth age relative to "normal" (3 dummy vars)
  #alcohol use -> yes decreases gestational birth age

#generate covariance matrix based on correlations between exposures
sample_covariance_matrix <- create_covmatrix(min=0.5,max=0.75,n=nexp,S=S)

#simulate N datasets without confounding covariates 
datasets <- generate_datasets(N=N, n=n, q=q,betas=betas, effect=effect, sample_means=sample_means, sample_covariance_matrix=sample_covariance_matrix, sd=1)

#Run WQS and qgcomp on datasets
print("Running WQS and qgcomp without covariates")
simulate_analysis(datasets=datasets,nexp=nexp,q=q,b=b,betas=betas,effect=1,covariates=F)

#simulate N datasets wit confounding covariates 
datasets_covars <- generate_datasets(N=N, n=n, q=q,betas=betas, effect=effect, sample_means=sample_means, sample_covariance_matrix=sample_covariance_matrix, sd=1, covariates=T,covar_betas=covar_betas)

#Run WQS and qgcomp on datasets with covariates included
print("Running WQS and qgcomp with covariates")
simulate_analysis(datasets=datasets_covars,nexp=nexp,q=q,b=b,betas=betas,effect=1,covariates=T)

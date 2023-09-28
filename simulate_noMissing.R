project_dir <- "~/OneDrive - Emory University/Courses/Thesis/"

source(paste0(project_dir,"simulate_power_typeI_source.R"))
set.seed(12345)

# Inputted simulation parameters
correlation=0.5    #Define correlation between exposures
n=1000          #Define number of observations per dataset
q=4            #Define number of quantiles

# Other simulation parameters
nexp <- 10 #number of exposures
sample_means <- rnorm(nexp,5,0.4) #Mean value for each exposure pulled from normal distribution
S <- rep(1,nexp) #vector of standard deviations = all 1
effect=1      #Define overall effect size
N=100           #Define total number of simulated datasets
b=100           #Define number of bootstrap samples
betas=c(0.1,0.25,0.65,rep(0,7)) #define weights of exposures

#generate covariance matrix based on correlations between exposures
sample_covariance_matrix <- create_covmatrix(min=correlation,max=correlation,n=nexp,S=S)

#simulate N datasets
datasets <- generate_datasets(N=N, n=n, q=q,betas=betas, effect=effect, sample_means=sample_means, sample_covariance_matrix=sample_covariance_matrix, sd=1)

sim_results <- simulate_analysis(datasets=datasets,nexp=nexp,q=q,b=b,betas=betas,effect=effect,covariates=F) %>%
  bind_rows(.id="method") %>% mutate(correlation=correlation,
                                     number_missing=number_missing,
                                     sample_size=n,
                                     q=q)

out_name <- paste0(project_dir,"Results/No_Missing/results_",
                   correlation,"corr_",
                   n,"n_",
                   q,"quantiles.Rds")

saveRDS(object=sim_results,file=out_name)
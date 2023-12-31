project_dir <- "/home/mgrieco/Thesis/Test-231003/"

source(paste0(project_dir,"simulate_power_typeI_source.R"))
set.seed(123)

args = commandArgs(trailingOnly=TRUE)

# Inputted simulation parameters
correlation=as.numeric(args[1])   #Define correlation between exposures
n=as.numeric(args[2])           #Define number of observations per dataset
q=as.numeric(args[3])           #Define number of quantiles

print("Simulation with no missing data")
print(paste0("Correlation: ", correlation))
print(paste0("Sample size: ", n))
print(paste0("Quantiles: ", q))

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
datasets <- generate_datasets(N=N, n=n, q=q,betas=betas, effect=effect, sample_means=sample_means, sample_covariance_matrix=sample_covariance_matrix, SD=1)

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
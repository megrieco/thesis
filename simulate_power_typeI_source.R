library(MASS)
library(tidyverse)
library(broom)
library(janitor)
library(corrplot)
library(gWQS)
library(qgcomp)

#imputation packages
library(mice)
library(missForest)
library(missRanger)
library(VIM)

#function to convert correlation matrix to covariance matrix using standard deviations
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}

#creates covariance matrix based on specified min and max correlation between n variables
### min: minimum possible correlation to be generated between variables
### max: maximum possible correlation to be generated between variables
### n: number of variables
### S: vector of standard deviations for n variables
create_covmatrix <- function(min,max,n,S, seed=1111){
  set.seed(seed)
  
  #Generate correlation matrix
  R <- matrix(runif(n^2,min=min,max=max), ncol=n) 
  R[lower.tri(R)] = t(R)[lower.tri(R)]
  #Set diagonals of matrix = 1
  diag(R)<-1
  
  #Calculate covariance matrix based on correlations and standard deviations
  sample_covariance_matrix <- cor2cov(R,S)
  
  return(sample_covariance_matrix)
}

#creates and returns N datasets based on specified exposure weights and overall effect
### N: number of simulated datasets to produce
### n: number of observations in each dataset
### q: number of quartiles to transform exposures into
### betas: vector of weights for exposure variables
### effect: overall effect size of exposure variables on outcome
### sample_means: vector of means for each exposure used to generate multivariate normal data
### sample_covariance_matrix: covariance matrix for exposures used to generate multivariate normal data
### sd: standard deviation of outcome variable
### covariates: whether or not covariates will be included in the model
### covar_betas: if covariates = T, then the weights for each of the covariates (parity, bmi under, bmi_over, bmi_obese, alcohol use)
generate_datasets <- function(N, n, q=10,betas, effect, sample_means, sample_covariance_matrix, sd, covariates=F,covar_betas=NULL){
  datasets <- list()
  
  #loop through to create N simulated datasets
  for (i in 1:N) {
    # Generate exposures (Multivariate Normal)
    set.seed(i)
    simulation <- mvrnorm(n = n,mu = sample_means, Sigma = sample_covariance_matrix)
    
    # exponentiate to get Multivariate log-normal distribution
    x <- simulation %>% as.data.frame(.) %>% exp(.)
    #change names to be x1, ..., x10
    names(x) <- gsub(x=names(x),pattern="V",replacement="x")
    
    ### transform x's into deciles before calculating y
    trans_decile <- function (x) ntile(x, q)
    x_trans <- x
    x_trans <- plyr::colwise(trans_decile)(x_trans)
    
    if(covariates){
      # generate covariates
      covars <- data.frame(cbind)
      parity=rnorm(n=n, mean=2,sd=1) #continuous
      bmi=rchisq(n = n, df = 2) 
      bmi=pmin(bmi * 2 + 18, 58)
      bmi_under=ifelse(bmi < 18.5,1,0) #ordinal - reference is normal
      bmi_over=ifelse(bmi >= 25 & bmi < 30,1,0)
      bmi_ob=ifelse(bmi >= 30,1,0)
      
      alc_use=rbinom(n=n,p=0.5, size=1) #binary - update probability as necessary
      #https://stats.oarc.ucla.edu/r/codefragments/mesimulation/
      
      covars <- data.frame(cbind(parity,bmi_under,bmi_over,bmi_ob,alc_use))
      
      # Generate outcome variable with covariates
      df <- x %>% mutate(mu=effect*(as.matrix(x_trans) %*% betas)+(as.matrix(covars) %*% covar_betas))
      df <- cbind(df,covars)
      
    } else{
      # Generate outcome variable with no covariates
      df <- x %>% mutate(mu=effect*(as.matrix(x_trans) %*% betas))
      
    }
    df$y = rnorm(n,df$mu,sd)  # Y_i ~ N(mu_i, sd)
    
    ### save whole dataset
    datasets[[i]] <- df
  }
  return(datasets)
}

#Run WQS on dataframe and return estimate, std error, and CI for each exposure and overall effect
#nexp: number of exposures
#q: number of quantiles that exposures will be transformed to
#b: number of bootstrap iterations
#covariates: whether to include covariates in the analysis
#Xs: vector with names of exposure variables
#i: ith dataset being run in simulation (for seed setting)
simulate_wqs <- function(df, nexp=10, q=10 , b=100,covariates=F,Xs,i=1 ){
  
  if(covariates){
      results <- gwqs(y ~ wqs+parity+bmi_under+bmi_over+bmi_ob+alc_use, mix_name = Xs, data = df, 
                      q = q, validation = 0.6, b = b, b1_pos = T, 
                      b1_constr = T, family = "gaussian", seed = i)
    } else{
      # we run the model and save the results in the variable "results"
      results <- gwqs(y ~ wqs, mix_name = Xs, data = df, 
                    q = q, validation = 0.6, b = b, b1_pos = T, 
                    b1_constr = T, family = "gaussian", seed = i)
    #Error if I set b1_constr=F: There are no positive b1 in the bootstrapped models
    
    }
    #final weights for exposure variables
    df_weights <- results$final_weights %>% dplyr::rename(exp=mix_name, estimate=mean_weight) %>%
      remove_rownames()
    
    #add overall estimate 
    df_weights$exp <- as.character(df_weights$exp)
    df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1][1]))
    
    # 2.5th/97.5th quartiles for exposures
    conf_band <- apply(results$bres[,1:nexp],2,quantile, probs = c(0.025,0.975)) %>% t() %>% as.data.frame() %>% rownames_to_column(var="exp")
    
    #add standard error for exposures
    conf_band$stde <- sapply(results$bres[1:nexp], function(x){sd(x)/sqrt(length((x)))})
    
    #add 95% CI and standard error (z dist'n) for overall
    overall_est <-summary(results)$coefficients[-1,1][1]
    overall_stde <- summary(results)$coefficients[-1,2][1]
    conf_band <- rbind(conf_band, c("Overall",
                                    overall_est - (1.96*overall_stde),
                                    overall_est + (1.96*overall_stde),
                                    overall_stde))
    
    #combine estimates and confidance band
    df_weights <- df_weights %>% inner_join(conf_band,by=c("exp")) %>%
      dplyr::rename(Lower=`2.5%`,Upper=`97.5%`)
  
  return(df_weights)
}

#Run qgcomp on dataframe and return estimate and CI for each exposure and overall effect
#q: number of quantiles that exposures will be transformed to
#b: number of bootstrap iterations
#Xs: vector with names of exposure variables
#i: ith dataset being run in simulation (for seed setting)
simulate_qgcomp <- function(df,q=10,b=100,Xs,i=1){

   # we run the model and save the results in the variable "results"
    results <- qgcomp.boot(y ~ .,
                           expnms=Xs,
                           df, family=gaussian(), q=q, B=b, seed=i)
    df_weights <- results[["fit"]][["coefficients"]][-1] %>% as.data.frame() 
    #select only exposure variables
    df_weights <- df_weights %>%
      filter(row.names(df_weights) %in% Xs) %>%
      #rename columns
      rownames_to_column(var = "exp") 
    names(df_weights) <- c("exp","estimate")
    
    #add overall estimate and standard deviations
    df_weights$exp <- as.character(df_weights$exp)
    df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
    
    stde_weights <- summary(results[["fit"]])$coefficients[Xs,2]
    
    df_weights$stde <- c(stde_weights ,Overall=summary(results)$coefficients[-1,2])
    
    ###use normal distribution to calculate confidence interval for all exposures/overall effect
    df_weights <- df_weights %>% 
      mutate(Lower=as.numeric(estimate) - (1.96*stde),
             Upper=as.numeric(estimate) + (1.96*stde))
    
    return(df_weights)
}

#function that calculates whether the confidence interval for each exposure estimate and overall estimate contains 0 and true rate
### df_weights: data frame with 5 columns containing the exposure name, stde, lower CI bound, upper CI bound
### betas_df: dataframe where first column is exposure names and second column are the true betas
### threshold: threshold for lower bound of CI to be considered approximately 0
### weights_all: data frame with 3 additional columns (true beta, contains0, containsT) with results from all simulations thus far
calculate_CI_contains <- function(df_weights,betas_df,threshold,weights_all){
  #Calculate whether interval contains 0 and contains true beta
  df_weights <- df_weights %>% left_join(betas_df,c("exp"="Xs")) %>% as.data.frame %>%
    mutate(containsT=case_when(
      betas >= as.numeric(Lower) & betas <= as.numeric(Upper) ~1, TRUE ~ 0),
      contains0=case_when(
        as.numeric(Lower) < threshold ~ 1, #threshold to consider as null
        TRUE ~ 0))
  
  weights_all <- rbind(weights_all, df_weights)
  
  return(weights_all)
}

#Calculate Power, Type I error, and whether true value is contained for each exposure and overall effect
#weights_all: dataframe containing results from all simulations with the following columns:
  #estimate: estimate from the simulation
  #betas: true value used to generate datasets
  #contains0: whether the CI for the estimate contains 0
  #containsT: whether the CI for the estimate contains the true value (betas)
  #exp: exposure (or "Overall") corresponding to each row
#Xs: vector with name of exposure variables
calculate_power_typeI <- function(weights_all,Xs){
  output <- weights_all %>% 
    #convert to numeric columns
    mutate(weight=as.numeric(estimate),
           stde=as.numeric(stde)) %>% 
    #Take average of betas for each x 
    reframe(true=mean(betas),
            mean=mean(weight),
            average_stde = mean(stde),
            Power=case_when(betas !=0 ~ sum(contains0==0)/n()),
            TypeI=case_when(betas ==0 ~ sum(contains0==0)/n()),
            containsTrate=sum(containsT==1)/n(), .by = exp) %>% unique() %>%
    as.data.frame() %>% 
    mutate_if(is.numeric, round, digits=3) %>%
    arrange(-mean)
  
  #convert exposure column to rownmaes
  output <- output %>% column_to_rownames(var="exp")
  #sort by X1, X2, ...
  output <- output[c("Overall",Xs),]
  
  #revert rownames back to column
  output <- output %>% rownames_to_column(var="exp")
  
  return(output)
}

#Wrapper function that runs WQS and qgcomp on each simulated dataset and outputs Type I error, Power, and contains True rate for each exposure and overall effect
#datasets: list of simulated dataframes
#nexp: number of exposure variables
#q: number of quantiles that exposures will be transformed to
#b: number of bootstrap iterations
#betas: true weights of exposure variables used to generate datasets
#effect: overall effect size used to generate datasets
#covariates: whether to include covariates in the model
simulate_analysis <- function(datasets, nexp=10, q=10 , b=100, betas, effect=1,covariates=F){
  #names of exposure columns
  Xs <- datasets[[1]] %>% dplyr::select(contains("x")) %>% names(.)
  
  #output for WQS weights
  weights_all_wqs <- as.data.frame(matrix(nrow=0,ncol=8))
  names(weights_all_wqs) <- c("exp","estimate","lower","upper","stder","containsT","contains0")
  
  #output for qgcomp weights
  weights_all_qgcomp <- weights_all_wqs
  
  #loop through N datasets
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
  
    #remove mu before running
    df <- df %>% select(-mu)
    
    #Run WQS on dataframe
    df_weights_wqs <- simulate_wqs(df, nexp=nexp, q=q , b=b, covariates=covariates,Xs=Xs,i=i)
    #Run qgcomp on dataframe
    df_weights_qgcomp <- simulate_qgcomp(df, q=q , b=b,Xs=Xs,i=i)
    
    #Combine true betas and overall effect into one dataframe
    betas_df <- as.data.frame(cbind(Xs,betas))
    betas_df <- rbind(betas_df,c("Overall",effect))
    betas_df$betas <- as.numeric(betas_df$betas)
    
    #Calculate whether interval contains 0 and contains true beta for this dataset 
    weights_df_wqs <- calculate_CI_contains(df_weights=df_weights_wqs,
                                             weights_all=weights_all_wqs,
                                             betas_df=betas_df,threshold=.001)
    weights_df_qgcomp <- calculate_CI_contains(df_weights=df_weights_qgcomp,
                                                weights_all=weights_all_qgcomp,
                                                betas_df=betas_df,threshold=.001)
    weights_all_wqs <- rbind(weights_all_wqs,weights_df_wqs)
    weights_all_qgcomp <- rbind(weights_all_qgcomp, weights_df_qgcomp)
    }
  
  out_wqs <- calculate_power_typeI(weights_all_wqs,Xs)
  out_qgcomp <- calculate_power_typeI(weights_all_qgcomp,Xs)
  
  return(list(WQS=out_wqs,qgcomp=out_qgcomp))
}

#randomly drop % of observations for each exposure
###datasets = list of dataframes from which to drop observations
### names: names of exposures (or covariates) to drop
###percent = percent of observations to make missing for each exposure
random_missing <- function(datasets, names, percent){
  set.seed(1111)
  
  datasets_incomplete <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    df_retain <- df %>% dplyr::select(!one_of(names))
    x_incomplete <- df %>% dplyr::select(one_of(names)) %>%
      lapply(., function(x){replace(x, sample(length(x), percent*length(x)/100), NA)}) %>% as.data.frame()
    
    df_incomplete <- cbind(x_incomplete, df_retain)
    names(df_incomplete) <- c(names(x_incomplete),names(df_retain))
    
    datasets_incomplete[[i]] <- df_incomplete
  }
  
  return(datasets_incomplete)
}

#Drop observations that are below the limit of detection in each dataset
### datasets: list of dataframes
### names: names of exposures (or covariates) to drop
### lods: vector with limit of detection for each exposure
missing_below_lod <- function(datasets, names, lods=NULL){
  set.seed(1111)
  
  datasets_incomplete <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    x_incomplete <- df %>% dplyr::select(one_of(names))
    df_retain <- df %>% dplyr::select(!one_of(names))

    #drop based on being below lod
    x_missing <- lapply(seq_along(x_incomplete), function(j){replace(x_incomplete[,j], which(x_incomplete[,j]<lods[j]), NA)}) %>% as.data.frame()

    df_incomplete <- cbind(x_missing, df_retain)
    names(df_incomplete) <- c(names(x_incomplete),names(df_retain))
    
    datasets_incomplete[[i]] <- df_incomplete
  }
  
  return(datasets_incomplete)
}

#impute missing values with LOD / sqrt(2) for each exposure
### datasets: list of dataframes with missing exposure values
### names: names of exposures to impute
### lods: vector with limit of detection for each exposure (same order as columns in datasets)
# requires that exposure names contain "x"
single_impute <- function(datasets_missing, names, lods){

  datasets_imputed <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets_missing[[i]]
    #separate exposure variables from non-exposure variables
    x_missing <- df %>% dplyr::select(one_of(names))
    df_retain <- df %>% dplyr::select(!one_of(names))
    
    #impute based on LOD/sqrt(2) for each exposure
    x_imputed <- lapply(seq_along(x_missing), function(j){replace(x_missing[,j], is.na(x_missing[,j]), lods[j]/sqrt(2))}) %>% as.data.frame()
    
    #recombine with non-missing columns
    df_imputed <- cbind(x_imputed, df_retain)
    names(df_imputed) <- c(names(x_missing),names(df_retain))
    
    #add to list of dataframes to output
    datasets_imputed[[i]] <- df_imputed
  }
  
  return(datasets_imputed)
}

#Drop observations that are below the limit of detection in each dataset
### datasets: list of dataframes
### lods: vector with limit of detection for each exposure
### pct: percent to drop lowest values of each exposure (for simulation purposes instead of lods)
drop_one_exposure <- function(datasets,exp,pct){
  set.seed(1111)
  
  datasets_incomplete <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    df[,c(exp)] <- replace(df[,c(exp)], which(df[,c(exp)]< quantile(df[,c(exp)], pct/100)), NA)
    
    datasets_incomplete[[i]] <- df
  }
  
  return(datasets_incomplete)
}

#drop whole class of chemicals for certain % of samples
###datasets = list of dataframes from which to drop observations
###names = vector of the names of exposures in the chemical class to drop
###percent = percent of observations to make missing for whole chemical class
drop_one_class <- function(datasets, names, pct){
  set.seed(1111)
  
  datasets_incomplete <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    x_incomplete <- df %>% dplyr::select(one_of(names)) 
    df_retain <- df %>% dplyr::select(!one_of(names)) 
    
    to_drop <- sample(nrow(x_incomplete), pct*nrow(x_incomplete)/100)
    x_missing <- lapply(x_incomplete, function(x){replace(x, to_drop, NA)}) %>% as.data.frame()
    
    df_incomplete <- cbind(x_missing, df_retain)
    names(df_incomplete) <- c(names(x_incomplete),names(df_retain))
    
    datasets_incomplete[[i]] <- df_incomplete
  }
  
  return(datasets_incomplete)
}

#takes in a list of dataframes with missing observations and returns a list of the dataframes with complete cases only
missing_completecases <- function(datasets_missing){
  complete_cases_df <- list()
  
  for (i in 1:length(datasets)) {
    df_incomplete <- datasets_missing[[i]]
    df_incomplete_drop <- df_incomplete[complete.cases(df_incomplete), ]
    
    complete_cases_df[[i]] <- df_incomplete_drop
  }
  
  return(complete_cases_df)
}

### Function that imputes missing values in a list of dataframes
# datasets_missing: list of dataframes that have missing values
# vars_to_use: vector of variables to use to impute (includes missing vars)
# missing_vars: vector of variables that have missing values
# method: imputation method to use (one of "mice","missForest","missRanger",
#         "hotdeck","irmi", or "knn")
# m: the number of imputations to perform (for mice imputation)
impute <- function(datasets_missing, vars_to_use, missing_vars, method,m=1){
  set.seed(1111)
  
  datasets_imputed <- list()
  
  for (i in 1:length(datasets_missing)) {
    
    df <- datasets_missing[[i]]
    df_to_use <- df %>% dplyr::select(one_of(vars_to_use)) 
    df_to_not_use <- df %>% dplyr::select(!one_of(vars_to_use))
    
    #log transform columns with missing values
    df_to_use[,missing_vars] <- log(df_to_use[,missing_vars])
    
    if(method=="mice"){
      #imputed_mice <- mice(df_to_use, m=1, maxit = 50, method = 'pmm', printFlag=F)
      imputed_mice <- mice(df_to_use,m=m)
    
      #get first imputation
      missing_imputed <- complete(imputed_mice, 1)
    } else if(method=="missForest"){
      missing_imputed <- missForest(df_to_use)$ximp
      
    } else if(method=="missRanger"){
      missing_imputed <- missRanger(df_to_use, num.trees = 100, verbose = 0)
      
    } else if(method == "hotdeck"){
      missing_imputed <- hotdeck(df_to_use,variable=missing_vars, imp_var = F)
      
    } else if(method == "irmi"){
      #remove any rows that have all NA
      df_to_use <- df_to_use[rowSums(is.na(df_to_use)) != ncol(df_to_use), ]
      
      #eps=5, maxit=50
      missing_imputed <- irmi(df_to_use, imp_var=F)
      
    } else if(method=="knn"){
      missing_imputed <- kNN(df_to_use,variable=missing_vars, imp_var = F)
      
    } else{
      stop("Error: Imputation method not valid")
    }
    #exponentiate to reverse transformation (ensures positive values)
    missing_imputed[,missing_vars] <- exp(missing_imputed[,missing_vars])
    
    df_imputed <- cbind(missing_imputed,df_to_not_use)
    names(df_imputed) <- c(names(missing_imputed),names(df_to_not_use))
    
    datasets_imputed[[i]] <- df_imputed
  }
  
  return(datasets_imputed)
}





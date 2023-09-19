library(MASS)
library(tidyverse)
library(broom)
library(janitor)
library(corrplot)
library(gWQS)
library(qgcomp)

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
### betas: vector of weights for exposure variables
### effect: overall effect size of exposure variables on outcome
### sample_means: vector of means for each exposure used to generate multivariate normal data
### sample_covariance_matrix: covariance matrix for exposures used to generate multivariate normal data
### sd: standard deviation of outcome variable
### covariates: whether or not covariates will be included in the model
### covar_betas: if covariates = T, then the weights for each of the covariates (parity, bmi under, bmi_over, bmi_obese, alcohol use)
generate_datasets <- function(N, n, betas, effect, sample_means, sample_covariance_matrix, sd, covariates=F,covar_betas=NULL){
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
      parity=rnorm(n=n, mu=2,sd=1) #continuous
      bmi=rchisq(n = 10, df = 5.5) 
      bmi=pmin(bmi * 2 + 18, 58)
      bmi_under=ifelse(bmi < 18.5,1,0) #ordinal - reference is normal
      bmi_over=ifelse(bmi >= 25 & bmi < 30,1,0)
      bmi_ob=ifelse(bmi >= 30,1,0)
      
      alc_use=rbinom(n=n,p=0.5) #binary - update probability as necessary
      #https://stats.oarc.ucla.edu/r/codefragments/mesimulation/
      
      covars <- data.frame(cbind(parity,bmi_under,bmi_over,bmi_ob,alc_use))
      
      # Generate outcome variable with covariates
      df <- x_trans %>% mutate(mu=effect*(as.matrix(x_trans) %*% betas)+(as.matrix(covars) %*% covar_betas))
      
      
    } else{
      # Generate outcome variable with no covariates
      df <- x_trans %>% mutate(mu=effect*(as.matrix(x_trans) %*% betas))
      
    }
    df$y = rnorm(n,df$mu,sd)  # Y_i ~ N(mu_i, sd)
    
    ### save whole dataset
    datasets[[i]] <- df
  }
  return(datasets)
}

simulate_wqs <- function(datasets, nexp=10, q=10 , b=100, betas, effect=1 ){
  #names of exposure columns
  Xs <- datasets[[1]] %>% dplyr::select(-y,-mu) %>% names(.)
  
  #output for WQS weights
  weights_all <- as.data.frame(matrix(nrow=0,ncol=8))
  names(weights_all) <- c("exp","estimate","stde","lower","upper","containsT","contains0")
  
  #loop through N datasets
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    # we run the model and save the results in the variable "results"
    results <- gwqs(y ~ wqs, mix_name = Xs, data = df, 
                    q = q, validation = 0.6, b = b, b1_pos = T, 
                    b1_constr = T, family = "gaussian", seed = 2016)
    #Error if I set b1_constr=F: There are no positive b1 in the bootstrapped models
    
    #final weights for exposure variables
    df_weights <- results$final_weights %>% dplyr::rename(exp=mix_name, estimate=mean_weight) %>%
      remove_rownames()
    
    #add overall estimate 
    df_weights$exp <- as.character(df_weights$exp)
    df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
    
    # 2.5th/97.5th quartiles for exposures
    conf_band <- apply(results$bres[,1:nexp],2,quantile, probs = c(0.025,0.975)) %>% t() %>% as.data.frame() %>% rownames_to_column(var="exp")
    
    #add 95% CI (z dist'n) for overall
    overall_est <-summary(results)$coefficients[-1,1]
    overall_stde <- summary(results)$coefficients[-1,2]
    conf_band <- rbind(conf_band, c("Overall",
                                    overall_est - (1.96*overall_stde),
                                    overall_est + (1.96*overall_stde)))
    
    #combine estimates and confidance band
    df_weights <- df_weights %>% inner_join(conf_band,by=c("exp")) %>%
      dplyr::rename(Lower=`2.5%`,Upper=`97.5%`)
    
    #add true value, and if true value falls in CI
    betas_df <- as.data.frame(cbind(Xs,betas))
    betas_df <- rbind(betas_df,c("Overall",effect))
    betas_df$betas <- as.numeric(betas_df$betas)
    
    df_weights <- df_weights %>% left_join(betas_df,c("exp"="Xs")) %>% as.data.frame %>%
      mutate(containsT=case_when(
        betas >= as.numeric(Lower) & betas <= as.numeric(Upper) ~1, TRUE ~ 0),
        contains0=case_when(
          as.numeric(Lower) < 0.001 ~ 1, #threshold to consider as null
          TRUE ~ 0))
    
    #combine with beta estimates of other datasets
    weights_all <- rbind(weights_all,df_weights)
    
    
  }
  
  output <- weights_all %>% 
    #convert to numeric columns
    mutate(weight=as.numeric(estimate)) %>% 
    #Take average of betas for each x 
    reframe(true=mean(betas),
            mean=mean(weight),
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
  
  return(output)
}

simulate_qgcomp <- function(datasets,betas,q=10,B=100, effect=1){
  #names of exposure columns
  Xs <- datasets[[1]] %>% dplyr::select(-y,-mu) %>% names(.)
  
  #output for qgcomp weights
  weights_all <- as.data.frame(matrix(nrow=0,ncol=8))
  names(weights_all) <- c("exp","estimate","stde","lower","upper","containsT","contains0")
  
  #loop through N datasets
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    # we run the model and save the results in the variable "results"
    results <- qgcomp.boot(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
                           expnms=Xs,
                           df, family=gaussian(), q=q, B=B, seed=2016)
    
    df_weights <- results[["fit"]][["coefficients"]][-1] %>% as.data.frame() %>%
      rownames_to_column(var = "exp") 
    names(df_weights) <- c("exp","estimate")
    
    #add overall estimate and standard deviations
    df_weights$exp <- as.character(df_weights$exp)
    df_weights <- rbind(df_weights,c("Overall",summary(results)$coefficients[-1,1]))
    
    stde_weights <- sqrt(diag(results$cov.yhat))
    
    df_weights$stde <- c(stde_weights ,summary(results)$coefficients[-1,2])
    
    #add confidence interval 
    
    ###use z-score rather than t
    df_weights <- df_weights %>% 
      mutate(Lower=as.numeric(estimate) - (1.96*stde),
             Upper=as.numeric(estimate) + (1.96*stde))
    
    #add true value, and if true value falls in CI
    betas_df <- as.data.frame(cbind(Xs,betas))
    betas_df <- rbind(betas_df,c("Overall",effect))
    betas_df$betas <- as.numeric(betas_df$betas)
    
    df_weights <- df_weights %>% left_join(betas_df,c("exp"="Xs")) %>% as.data.frame %>%
      mutate(containsT=case_when(
        betas >= as.numeric(Lower) & betas <= as.numeric(Upper) ~1, TRUE ~ 0),
        contains0=case_when(
          as.numeric(Lower) <= 0 & as.numeric(Upper) >= 0 ~ 1,
          #as.numeric(Lower) < 0.001 ~ 1, #threshold to consider as null
          TRUE ~ 0))
    
    #combine with beta estimates of other datasets
    weights_all <- rbind(weights_all,df_weights)
    
    
  }
  
  output <- weights_all %>% 
    #convert to numeric columns
    mutate(weight=as.numeric(estimate)) %>% 
    #Take average of betas for each x 
    reframe(true=mean(betas),
            mean=mean(weight),
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
  
  return(output)
}


#randomly drop % of observations for each exposure
###datasets = list of dataframes from which to drop observations
###percent = percent of observations to make missing for each exposure
random_missing <- function(datasets, percent){
  set.seed(1111)
  
  datasets_incomplete <- list()
  
  for (i in 1:length(datasets)) {
    df <- datasets[[i]]
    x_incomplete <- df %>% dplyr::select(-y,-mu) %>%
      lapply(., function(x){replace(x, sample(length(x), percent*length(x)/100), NA)}) %>% as.data.frame()
    
    df_incomplete <- cbind(x_incomplete, df$mu, df$y)
    names(df_incomplete) <- c(names(x_incomplete),"mu","y")
    df_incomplete_drop <- df_incomplete[complete.cases(df_incomplete), ]
    
    datasets_incomplete[[i]] <- df_incomplete
  }
  
  return(datasets_incomplete)
}

missing_completecases <- function(datasets_missing){
  complete_cases_df <- list()
  
  for (i in 1:length(datasets)) {
    df_incomplete <- datasets_missing[[i]]
    df_incomplete_drop <- df_incomplete[complete.cases(df_incomplete), ]
    
    complete_cases_df[[i]] <- df_incomplete_drop
  }
  
  return(complete_cases_df)
}

mice_impute <- function(datasets_missing){
  set.seed(1111)
  
  datasets_imputed <- list()
  
  for (i in 1:length(datasets_missing)) {
    df <- datasets_missing[[i]]
    x_missing <- df %>% dplyr::select(-y,-mu) 
    y_mu <- df %>% dplyr::select(y,mu)
    
    #imputed_mice <- mice(x_missing, m=5, maxit = 50, method = 'pmm', printFlag=F)
    imputed_mice <- mice(x_missing)
    
    #add y and mu back in
    imputed_mice <- cbind(imputed_mice, y_mu)
    datasets_imputed[[i]] <- complete(imputed_mice)
  }
  
  return(datasets_imputed)
  }
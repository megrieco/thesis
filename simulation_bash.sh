#!/usr/bin/env bash
#SBATCH --job-name=simulate_missing                            # Job name
#SBATCH --mail-type=END,FAIL                          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=megan.grieco@emory.edu         # Where to send mail
#SBATCH --nodes=1                                    # Run all processes on a single node
#SBATCH --ntasks-per-node=10                        # Number of CPU cores per task
#SBATCH --mem=128000                               # Job memory request
#SBATCH --array=1-6                        # job array with index values 1, 2, ... 8 for each correlation
#SBATCH --error=job.%A_%a.err
#SBATCH --output=job.%A_%a.out
#SBATCH --partition=week-long-cpu

#input arguments to script: sample size and number of quantiles
n=$1 #sample size
q=$2 #number of quantiles


correlations=(50 55 60 65 70 75) #correlation between exposures
imputation_methods=("none" "single" "mice" "missForest" "missRanger" "hotdeck" "irmi" "knn") #imputation method for missing values
number_missings=(2 3 4 5) #number of exposures missing
missing_percents=(5 10 20 30 50) #percent missingness for each exposure

# correlation based on job array
correlation =${correlations[$SLURM_ARRAY_TASK_ID]}

module load R/4.2.2

#Run no missing case
Rscript simulate_noMissing.R $correlation $n $q

#loop through each imputation method, number missing, and missing percent for LOD case and random missingness case
for imputation_method in ${imputation_methods[@]}; do
    for number_missing in ${number_missings[@]}; do
        for missing_percent in ${missing_percents[@]}; do
            #Run random missing case
            Rscript missing_random_sim.R $correlation $number_missing $n $q $imputation_method $missing_percent

            #Run LOD missing case
            Rscript missing_lod_sim.R $correlation $number_missing $n $q $imputation_method $missing_percent
        done
    done
done

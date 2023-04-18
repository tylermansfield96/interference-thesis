# Header -----------------------------------------------------------------------
# Proj: Randomization Test for the Spillover Effect in Two-Stage Observational Data
# Author: Tyler Mansfield (tyler.mansfield96@gmail.com)
# Desc: Code to run repeated simulations. The main function, run_sims, allows 
#       various arguments to give flexibility in the applied method (such as
#       which outcome function to use, whether or not to match on cluster-level
#       covariates, etc). The actual simulations used in the paper are called 
#       in the final section of the script, which takes around 48 hours to run
#       on my laptop, but which could easily be parallelized.


### Load the functions needed to generate data and apply the method -----------

source("01_simulation_setup.R")

### Function to run entire simulation -----------------------------------------

# This function runs repeated simulations that each output a given p-value for
# the procedure dictated in the thesis. First, a dataset is generated using the
# outcome_fxn (This process is repeated n_data times). For each dataset, the 
# randomization test is ran n_runs times. The resulting p-values are stored in
# a n_data x n_runs matrix and is returned
# Other arguments:
# - lambda_vals: If there are multiple values in this vector, the process above
#         will repeat for each value of lambda
# - match_on_cluster_cov: If this is TRUE, then within a fixed propensity score,
#         we will also match on cluster level covariates and size
# - data_gen_name: A string used when saving the results to identify the simulation 
#         name
# - save: Whether or not to save the results to a csv file
# - more_int_more_outcome: a Boolean value that specifies whether having more 
#         interference is expected to result in a higher outcome value (since we 
#         are doing a one-sided p-value). If interference causes control units
#         in treated clusters to have higher outcome values, this argument should 
#         be TRUE
# - use_all_psuedocontrols: If FALSE, we use the method and test statistic
#         dictated in Proposition 1. Otherwise, use the methodology dictated in
#         Proposition 2.
run_sims <- function(outcome_fxn, 
                     lambda_vals, 
                     n_data,
                     n_runs, 
                     match_on_cluster_cov = FALSE,
                     data_gen_name = NULL, 
                     more_int_more_outcome = NULL,
                     save = T, 
                     use_all_psuedocontrols = NULL)
{
  assert_that(is.null(data_gen_name) == F | save == F)
  assert_that(is.logical(match_on_cluster_cov),
              is.logical(more_int_more_outcome),
              is.logical(use_all_psuedocontrols),
              is.number(n_data),
              is.number(n_runs),
              is.numeric(lambda_vals),
              is.function(outcome_fxn))
  
  for (lambda in lambda_vals)
  {
    message("Lambda ", lambda)
    set.seed(1)
    
    # Create a matrix to store p-values
    pval_matrix <- matrix(nrow = n_data, ncol = n_runs)
    
    for (i in 1:n_data)
    {
      if (i %% 100 == 0 | n_runs > 100)
      {
        message("Dataset ", i)
      }
      
      # Generate a dataset
      dataset <- generate_data(lambda = lambda, outcome_fxn = outcome_fxn)
      
      p_val_dist <- rep(NA, n_runs)
      
      for (j in 1:n_runs)
      {
        if (j %% 100 == 0)
        {
          message("Run ", j)
        }
        
        # Draw virtual treatment, match clusters
        results <- randomization_test_part1(dataset,
                                            match_on_cluster_cov = match_on_cluster_cov)
        dataset_aug <- results$dataset_aug
        matches <- results$matches
        
        # Select focal units (Proposition 1), compute p-value 
        p_val_dist[j] <- randomization_test_part2(dataset_aug = dataset_aug, 
                                                  our_matches = matches,
                                                  more_int_more_outcome = more_int_more_outcome,
                                                  use_all_psuedocontrols = use_all_psuedocontrols)
      }
      
      # Save results
      pval_matrix[i,] <- p_val_dist
    }
    if (save == T)
    {
      # Create file name
      file_name <- paste0("../data/simulation_pvals_sim",
                          data_gen_name,
                          "_",
                          n_data,
                          "datasets_",
                          n_runs,
                          "runs_",
                          ifelse(use_all_psuedocontrols == FALSE,
                                 "prop1_",
                                 "prop2_"),
                          ifelse(match_on_cluster_cov == FALSE,
                                 "",
                                 "clustercovmatch_"),
                          lambda,
                          "lambda.csv")
      
      write.csv(pval_matrix, file = file_name)
    }
  }
  return(pval_matrix)
}

### Run simulations ------------------------------------------------------------


# For each value of lambda, outcome scenario (A-C), and method (Proposition 1 
# and 2), we simulate 1000 datasets and run the randomization test exactly once 
# for each dataset

for (use_all_psuedocontrols in c(F, T)) # Iterating between Proposition 1 and 2
{
  run_sims(outcome_fxn = outcome_simulationA,
           lambda_vals = c(0,1,2.5,5,10),
           n_data = 1000,
           n_runs = 1,
           match_on_cluster_cov = FALSE,
           data_gen_name = "A",
           more_int_more_outcome = TRUE,
           save = TRUE,
           use_all_psuedocontrols = use_all_psuedocontrols)

  run_sims(outcome_fxn = outcome_simulationB,
           lambda_vals = c(0,0.1,0.25,0.5,1),
           n_data = 1000,
           n_runs = 1,
           match_on_cluster_cov = FALSE,
           data_gen_name = "B",
           more_int_more_outcome = FALSE,
           save = TRUE,
           use_all_psuedocontrols = use_all_psuedocontrols)

  run_sims(outcome_fxn = outcome_simulationC,
           lambda_vals = c(0,0.25,0.5,1,2.5),
           n_data = 1000,
           n_runs = 1,
           match_on_cluster_cov = FALSE,
           data_gen_name = "C",
           more_int_more_outcome = TRUE,
           save = TRUE,
           use_all_psuedocontrols = use_all_psuedocontrols)
}

# For each value of lambda, we simulate 8 datasets using outcome setting B and 
# run the randomization test (Using proposition 2) 1000 times for each dataset

run_sims(outcome_fxn = outcome_simulationB, 
         lambda_vals = c(0,0.1,0.25,0.5,1), 
         n_data = 8,
         n_runs = 1000,
         match_on_cluster_cov = FALSE,
         data_gen_name = "B", 
         more_int_more_outcome = FALSE,
         save = TRUE, 
         use_all_psuedocontrols = TRUE)

# For lambda = 0 and outcome scenario C, we simulate 1000 datasets and run the 
# randomization test using Proposition 2 exactly once for each dataset while 
# matching on cluster covariates within fixed propensity score strata

run_sims(outcome_fxn = outcome_simulationC,
         lambda_vals = c(0),
         n_data = 1000,
         n_runs = 1,
         match_on_cluster_cov = TRUE,
         data_gen_name = "C",
         more_int_more_outcome = TRUE,
         save = TRUE,
         use_all_psuedocontrols = TRUE)

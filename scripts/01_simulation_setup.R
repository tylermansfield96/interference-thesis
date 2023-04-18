# Header -----------------------------------------------------------------------
# Proj: Randomization Test for the Spillover Effect in Two-Stage Observational Data
# Author: Tyler Mansfield (tyler.mansfield96@gmail.com)
# Desc: Functions to generate data from various simulation distributions and run the
#       randomization test for spillover effect outlined in Propositions 1 and 2. For
#       an example of how to use the functions listed here, see end of script

# Load relevant libraries --------------------------------------------------

library(data.table)
library(magrittr)
library(ggplot2)
library(assertthat)
library(MatchIt)


# Data generating functions --------------------------------------------------


# outcome_fxn should be of the form
# function(unit_treatment, unit_cov1, unit_cov2, lambda, mean_others_treated)
# and should be vectorized (e.g. return multiple outcomes if the arguments
# are vectors). See below for examples
generate_data <- function(lambda = 1, 
                          outcome_fxn)
{
  ## Constants ----------------------------------------------------------
  n_clusters <- 125
  # P(n_i = cluster_sizes_support[i]) = cluster_sizes_probs[i]
  cluster_sizes_support <- c(8,22,40) 
  cluster_sizes_probs <- c(0.4, 0.35, 0.25)
  
  # P(Z_{ij} = 1) = plogis(treatment_beta_vals[1] * 1 +
  #                        treatment_beta_vals[2] * unit_cov1 +
  #                        treatment_beta_vals[3] * unit_cov2 +
  #                        treatment_beta_vals[4] * cluster_cov1 +
  #                        treatment_beta_vals[5] * cluster_cov2)
  treatment_beta_vals <- c(0.75, -0.015, -0.025, 0, 1)
  
  # Determine cluster size, covariates, and propensity score -----------
  multinomial_draw <- rmultinom(1, n_clusters, prob = cluster_sizes_probs) |> 
    as.numeric()
  cluster_sizes <- mapply(rep, cluster_sizes_support, multinomial_draw) |> 
    unlist() |>
    sample()
  
  cluster_cov1 <- rnorm(n = n_clusters, mean = 6, sd = 1)
  cluster_cov2 <- rnorm(n = n_clusters, mean = 0, sd = sqrt(0.75))
  
  # P(theta = 0.3) = 0.25 + 0.25*I(cluster_cov1 <= 5.5)
  # P(theta = 0.5) = 0.25 + 0.25*I(5.5 < cluster_cov1 < 6.5)
  # P(theta = 0.7) = 0.25 + 0.25*I(cluster_cov1 >= 6.5)
  cluster_thetas <- rep(NA, n_clusters)
  theta_0 <- runif(n = n_clusters) # Use theta_0 to determine propensity score (theta)
  cluster_thetas[theta_0 < 0.25] <- 0.3
  cluster_thetas[theta_0 >= 0.25 & theta_0 < 0.50] <- 0.5
  cluster_thetas[theta_0 >= 0.50 & theta_0 < 0.75] <- 0.7
  cluster_thetas[theta_0 >= 0.75 & cluster_cov1 < 5.5] <- 0.3
  cluster_thetas[theta_0 >= 0.75 & cluster_cov1 >= 5.5 & cluster_cov1 < 6.5] <- 0.5
  cluster_thetas[theta_0 >= 0.75 & cluster_cov1 >= 6.5] <- 0.7
  
  dataset <- data.table()
  
  for (i in 1:n_clusters)
  {
    n_i <- cluster_sizes[i]
    
    # Generate cluster level treatment and covariates -------------------------
    cluster_treatment <- rbinom(1, 1, prob = cluster_thetas[i])
    
    unit_cov1 <- rnorm(n_i, 40, sqrt(5))
    unit_cov2 <- rnorm(n_i, cluster_cov1[i], sqrt(0.2))
    
    # Generate individual-level treatment status ------------------------------
    # Log odds of treatment is beta_1 + beta_2 * unit_cov1 + beta_3 * unit_cov2 +
    #                          beta_4 * cluster_cov1 + beta_5 * cluster_cov2 
    # Compute via matrix multiplication
    unit_covariate_matrx <- matrix(data = c(rep(1, n_i),
                                            unit_cov1,
                                            unit_cov2,
                                            rep(cluster_cov1[i], n_i),
                                            rep(cluster_cov2[i], n_i)),
                                   nrow = n_i)
    
    unit_log_odds_treat <- unit_covariate_matrx %*% treatment_beta_vals |>
      as.numeric()
    
    # Convert log-odds to probabilities
    unit_pi <- plogis(unit_log_odds_treat) # Pi is the unit level propensity score
    
    # Unit level treatment, if the cluster is treated
    unit_treatment <- rbinom(n_i, 1, prob = unit_pi) * cluster_treatment
    
    # Number of others treated 
    mean_others_treated <- (sum(unit_treatment) - unit_treatment) / (n_i - 1)
    
    
    # Generate individual-level outcome -----------------------------------------
    # Outcome function should be of the form
    # function(unit_treatment, unit_cov1, unit_cov2, lambda, mean_others_treated)
    unit_outcome <- outcome_fxn(unit_treatment, 
                                unit_cov1, 
                                unit_cov2, 
                                lambda, 
                                mean_others_treated)
    
    # Combine individual and cluster level info into one data.table -------------
    cluster_data <- data.table(cluster_id = i, 
                               cluster_treatment = cluster_treatment,
                               cluster_theta = cluster_thetas[i],
                               cluster_size = n_i,
                               cluster_cov1 = cluster_cov1[i],
                               cluster_cov2 = cluster_cov2[i],
                               unit_id = 1:n_i,
                               unit_cov1 = unit_cov1,
                               unit_cov2 = unit_cov2,
                               unit_treatment = unit_treatment,
                               unit_pi = unit_pi,
                               unit_outcome = unit_outcome)
    
    dataset <- rbind(dataset, cluster_data)
  }
  
  dataset
}

outcome_simulationA <- function(unit_treatment, 
                                unit_cov1, 
                                unit_cov2, 
                                lambda, 
                                mean_others_treated)
{
  # Make sure arguments have same length
  vecs_same_len <- list(unit_treatment, unit_cov1, unit_cov2, mean_others_treated)
  vec_lens <- lapply(vecs_same_len, length) |> unlist()
  assert_that(length(unique(vec_lens)) == 1)
  
  n_i <- length(unit_treatment)
  
  # Outcome log odds
  outcome_log_odds <- 0.1 - 0.05 * unit_cov1 + 0.5 * unit_cov2 - 0.5 * unit_treatment +
    lambda * mean_others_treated * (0.2 - 0.25 * unit_treatment)
  
  unit_prob_outcome <- plogis(outcome_log_odds)
  
  unit_outcome <- rbinom(n_i, 1, prob = unit_prob_outcome)
  
  unit_outcome
}

outcome_simulationB <- function(unit_treatment, 
                                unit_cov1, 
                                unit_cov2, 
                                lambda, 
                                mean_others_treated)
{
  # Make sure arguments have same length
  vecs_same_len <- list(unit_treatment, unit_cov1, unit_cov2, mean_others_treated)
  vec_lens <- lapply(vecs_same_len, length) |> unlist()
  assert_that(length(unique(vec_lens)) == 1)
  
  n_i <- length(unit_treatment)
  
  # Outcome log odds
  outcome_log_odds_no_interfere <- 0.1 - 0.05 * unit_cov1 + 0.5 * unit_cov2 - 0.5 * unit_treatment
  
  unit_prob_outcome <- plogis(outcome_log_odds_no_interfere) * (1 - (mean_others_treated * lambda))
  
  unit_outcome <- rbinom(n_i, 1, prob = unit_prob_outcome)
  
  unit_outcome
}

outcome_simulationC <- function(unit_treatment, 
                                unit_cov1, 
                                unit_cov2, 
                                lambda, 
                                mean_others_treated)
{
  # Make sure arguments have same length
  vecs_same_len <- list(unit_treatment, unit_cov1, unit_cov2, mean_others_treated)
  vec_lens <- lapply(vecs_same_len, length) |> unlist()
  assert_that(length(unique(vec_lens)) == 1)
  
  n_i <- length(unit_treatment)
  
  # Random noise
  epsilon <- rnorm(n_i, mean = 0, sd = 0.3)
  
  unit_outcome <- 0.1 - 0.05 * unit_cov1 + 0.5 * unit_cov2 - 0.5 * unit_treatment +
    lambda * mean_others_treated * (0.2 - 0.25 * unit_treatment) + epsilon
  
  unit_outcome
}



# Randomization test functions -------------------------------------------------

# Part 1 of the randomization test (which is the same for Propositions 1 and 2)
# creates the virtual treatment assignments, discards clusters without any control
# units, and performs cluster matching exactly on propensity scores and (potentially)
# on cluster-level covariates within exact propensity score stratum

randomization_test_part1 <- function(dataset, match_on_cluster_cov = NULL)
{
  assert_that(is.logical(match_on_cluster_cov))
  
  n <- nrow(dataset)
  
  # Create virtual treatment
  dataset$virtual_treatment <- rbinom(n, 1, prob = dataset$unit_pi)
  dataset[cluster_treatment == 1, virtual_treatment := unit_treatment]

  # Discard any clusters with no controls
  dataset %<>% .[, n_controls := sum(virtual_treatment == 0), by = cluster_id] %>%
    .[n_controls > 0] %>%
    .[, !"n_controls"]
  
  # Match clusters on propensity scores
  cluster_cols <- grep("cluster", names(dataset), value=T)
  cluster_data <- dataset %>% 
    .[, .SD[1], by = cluster_id] %>%
    .[, cluster_cols, with = FALSE]
  
  matches <- data.table()
  
  for (i in unique(cluster_data$cluster_theta))
  {
    cluster_data_subset <- cluster_data[cluster_theta == i]
    
    if(match_on_cluster_cov == T)
    {
      formula <- as.formula("cluster_treatment ~ cluster_cov1 + cluster_cov2 + cluster_size")
    } else {
      formula <- as.formula("cluster_treatment ~ cluster_size")
    }
    
    suppressWarnings({
      m.out <- matchit(formula=formula,
                       data=cluster_data_subset,
                       distance='mahalanobis')
      
      cluster_matches <- get_matches(m.out,data=cluster_data_subset) %>%
        as.data.table() %>%
        .[order(subclass, cluster_treatment), .(match_id = as.numeric(subclass) + nrow(matches)/2, 
                                                cluster_id, 
                                                cluster_treatment)] 
      
      
      
      })
    
    matches %<>% rbind(cluster_matches)
  }
  
  return(list(matches = matches, dataset_aug = dataset))
}

# Part 2 of the randomization test is different for Propositions 1 and 2
# For Proposition 1, (use_all_psuedocontrols = F) we uniformly chooses a single 
# pseudo-control unit from each cluster to be the focal unit. For Proposition 2,
# (use_all_psuedocontrols = T), we use all pseudo-controls in each cluster.
# We then calculate the test statistic, and then computes the p-value. If the 
# outcome is binary, the p-value can be computed directly, otherwise it will be 
# estimated using Monte Carlo simulation.
# The first two arguments are the outputs from the part 1 function. The
# argument more_int_more_outcome is a boolean value that specifies whether 
# having more interference is expected to result in a higher outcome value
# (since we are doing a one-sided p-value). If interference causes control units
# in treated clusters to have higher outcome values, this argument should be TRUE

randomization_test_part2 <- function(dataset_aug, 
                                     our_matches, 
                                     more_int_more_outcome = NULL,
                                     use_all_psuedocontrols = NULL)
{
  assert_that(is.logical(more_int_more_outcome), is.logical(use_all_psuedocontrols))
  
  if (use_all_psuedocontrols == F)
  {
    # Choose focal units 
    focal_units_potential <- dataset_aug[virtual_treatment == 0]
    focal_units <- focal_units_potential[sample(1:nrow(focal_units_potential))] %>%
      .[, .SD[1], by = cluster_id] %>%
      .[, .(cluster_id, focal_unit_id = unit_id, focal_unit_outcome = unit_outcome)]
    
    our_matches %<>% merge(focal_units, by = "cluster_id", all.x = T)
    
    # Pivot wider
    matches_wide <- dcast(our_matches, match_id ~ cluster_treatment, 
                          value.var = "focal_unit_outcome") %>%
      setnames(c("0", "1"), c("Control_Outcome", "Treated_Outcome")) %>%
      .[, difference := Treated_Outcome - Control_Outcome]
  } else {
    
    # Choose focal units 
    focal_units <- dataset_aug[virtual_treatment == 0]
    
    # Take mean of focal units
    focal_units_mean <- focal_units[, .(mean_focal_unit_outcome = mean(unit_outcome)), 
                                    by = cluster_id]
    
    our_matches %<>% merge(focal_units_mean, by = "cluster_id", all.x = T)
    
    # Pivot wider
    matches_wide <- dcast(our_matches, match_id ~ cluster_treatment, 
                          value.var = "mean_focal_unit_outcome") %>%
      setnames(c("0", "1"), c("Control_Outcome", "Treated_Outcome")) %>%
      .[, difference := Treated_Outcome - Control_Outcome]
    
  }
  
  # P-value for binary outcomes is easy
  if (all(matches_wide$difference %in% c(-1,0,1)))
  {
    # The only important difference to permute are the non-zero values
    matches_wide <- matches_wide[difference != 0]
    n_ones <- nrow(matches_wide[difference == 1])
    n_negativeones <- nrow(matches_wide[difference == -1])
    
    if (more_int_more_outcome == T)
    {
      p_val <- pbinom(n_negativeones, n_ones + n_negativeones, prob = 0.5)
    } else {
      p_val <- pbinom(n_ones, n_ones + n_negativeones, prob = 0.5)
    }
  } else # Estimate p-value by repeatedly permuting which cluster in each pair was treated
  {
    B <- 100000
    d <- matches_wide$difference
    m0 <- mean(d) # Our observed test statistic
    
    # This is our randomization distribution under the null (H0: mu_d = 0)
    rndmdist <- replicate(B,mean((rbinom(length(d),1,.5)*2-1)*d))
    
    if (more_int_more_outcome == T) #H1 is mu_d > 0
    {
      p_val <- sum(rndmdist >= m0)/length(rndmdist)
    } else { #H1 is mu_d < 0
      p_val <- sum(rndmdist <= m0)/length(rndmdist)
    }
  }
  
  p_val
}

# Example of how to use the code above -----------------------------------------

# dataset <- generate_data(lambda = 0, outcome_fxn = outcome_simulationA)
# results <- randomization_test_part1(dataset, match_on_cluster_cov = F)
# dataset_aug <- results$dataset_aug
# our_matches <- results$matches
# pval <- randomization_test_part2(dataset_aug, 
#                                  our_matches, 
#                                  more_int_more_outcome = T, 
#                                  use_all_psuedocontrols = T)
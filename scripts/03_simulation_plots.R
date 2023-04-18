# Header -----------------------------------------------------------------------
# Proj: Randomization Test for the Spillover Effect in Two-Stage Observational Data
# Author: Tyler Mansfield (tyler.mansfield96@gmail.com)
# Desc: Code to generate the plots used in the thesis write-up. These functions
#       make use of the data saved from running `02_run_simulation.R`


### Load relevant libraries --------------------------------------------------

library(tidyverse)
library(latex2exp)
library(data.table)
library(stringr)


### Function to extract simulation details from data filename ----------------

extract_metadata <- function(file_names)
{
  results <- data.frame(simulation = c(),
                        lambda = c(),
                        test_stat = c(),
                        match = c()) 
  
  for (i in 1:length(file_names))
  {
    file <- file_names[i]
    
    sim_num <- gsub("_.+", "", gsub(".+_sim", "", file))
    lambda  <- gsub(".+_","",gsub("lambda.+", "", file))
    test_stat <- ifelse(grepl("prop1", file),
                        "Proposition 2", 
                        "Proposition 1")
    match <- ifelse(grepl("clustercovmatch", file),
                    "Covariate Match", 
                    "No Covariate Match")
    n_datasets <- gsub(".+_","",gsub("datasets.+", "", file))
    n_runs <- gsub(".+_","",gsub("runs.+", "", file))
    
    
    results <- rbind(results,
                     data.frame(simulation = sim_num,
                                n_datasets = n_datasets,
                                n_runs = n_runs,
                                lambda = lambda,
                                test_stat = test_stat,
                                match = match,
                                file_name = file))
    
  }
  
  results
}


### Figure 1: Validity of Method -----------------------------------------------

target_files <- list.files(path = "../data", pattern = "sim[A-Z]") %>%
  extract_metadata() %>%
  filter(n_runs == "1", match == "No Covariate Match") %>%
  pull(file_name)

# Input all p-values in target_files into a single dataframe
p_val_df <- data.frame()

for (file in target_files)
{
  # Upload a single one
  pvals <- fread(paste0("../data/",file)) %>%
    .[,-1]
  
  data_matrix <- as.matrix(pvals)
  p_val_df_i <- cbind(extract_metadata(file), t(data_matrix))
  p_val_df <- rbind(p_val_df, p_val_df_i)
}

# Now pivot the dataframe to be in long format
value_columns <- names(p_val_df)[grepl("\\d", names(p_val_df))]

p_val_df_long <- pivot_longer(p_val_df, cols = value_columns) 

# Pad lambda values for ease of reading, rename columns, make columns factors
p_val_df_long$lambda <- sprintf(as.numeric(p_val_df_long$lambda), fmt = '%#.2f') %>%
  str_pad(width = 5, side = "left", pad = "0")
cols_to_factor <- names(p_val_df_long)[!grepl("value", names(p_val_df_long))]
p_val_df_long[cols_to_factor] <- lapply(p_val_df_long[cols_to_factor] , factor)
p_val_df_long <- rename(p_val_df_long, dataset_id = name)
summary(p_val_df_long)

# ECDFs of p-value distribution across various propositions and simulations 
ggplot(p_val_df_long, aes(x = value, group = lambda)) +
  stat_ecdf(aes(col = lambda), linewidth = 1) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), lty = "dotted", linewidth = 1) +
  labs(x = "x",
       y = TeX("$P(pval \\leq x)$"),
       col = TeX("$\\lambda$ value"),
       lty = TeX("$\\lambda$ value")) +
  facet_grid(test_stat~simulation) +
  #scale_color_brewer(palette = "Dark2") +
  theme_minimal() +
  theme(panel.spacing = unit(1.5, "lines"),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8))

ggsave(paste0("../images/ecdf_contrast.png"), 
       height = 3.5,
       width = 6,
       dpi = 300)

# Clean up environment. Not necessary, but names are reused below
rm(p_val_df, p_val_df_i, p_val_df_long, pvals, data_matrix, cols_to_factor,
   file, target_files, value_columns)


### Figure 2a: Plot multiple p-value histograms together ----------------------

target_files <- list.files(path = "../data", pattern = "sim[A-Z]") %>%
  extract_metadata() %>%
  mutate(n_runs = as.numeric(n_runs)) %>%
  filter(n_runs > 1) %>%
  pull(file_name)

p_val_df <- data.frame()

for (file in target_files)
{
  # Upload a file
  pvals <- fread(paste0("../data/",file)) %>%
    .[,-1]
  
  data_matrix <- as.matrix(pvals)
  
  # Add the p-values for each unique dataset one-by-one, adding a unique dataset
  # identifier
  for (i in 1:nrow(data_matrix))
  {
    relevent_metadata <- extract_metadata(file) %>%
      select(simulation, lambda) %>%
      mutate(dataset = i)
    
    p_val_df_i <- cbind(relevent_metadata, t(data_matrix[i, ]))
    
    # Pivot long
    value_columns <- names(p_val_df_i)[grepl("\\d", names(p_val_df_i))]
    p_val_df_i_long <- pivot_longer(p_val_df_i, 
                                  cols = value_columns, 
                                  names_to = "rep", 
                                  values_to = "p_val") 
    
    
    p_val_df <- rbind(p_val_df, p_val_df_i_long)
  }
  
  rm(p_val_df_i, p_val_df_i_long, data_matrix, pvals, relevent_metadata, file, i,
     value_columns)
}

p_val_df$lambda <- sprintf(as.numeric(p_val_df$lambda), fmt = '%#.2f') %>%
  str_pad(width = 5, side = "left", pad = "0")
p_val_df$dataset_lambda <- factor(paste0(p_val_df$lambda,
                                         p_val_df$dataset))

ggplot(p_val_df, aes(y = p_val, x = dataset_lambda)) +
  geom_boxplot(aes(group = dataset_lambda, col = lambda), linewidth = 0.5) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2") +
  labs(x = "",
       y = "Distribution of P-Value",
       title = "Various Datasets from Simulation B",
       col = TeX("$\\lambda$ Value")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))

ggsave(paste0("../images/single_dataset_contrast_dist.png"), dpi = 300)


### Figure 2b: Proportion of times reject at 0.05 ------------------------------

p_val_df$significant <- as.numeric(p_val_df$p_val < 0.05)
p_val_df_agg <- p_val_df %>%
  group_by(dataset_lambda, lambda) %>%
  summarise(prop_reject = mean(significant))

ggplot(p_val_df_agg, aes(y = prop_reject, x = dataset_lambda)) +
  geom_bar(aes(fill = lambda), linewidth = 0.5, stat = "identity") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "",
       y = TeX("Proportion of $\\textit{p}-values < 0.05$"),
       title = "Various Datasets from Simulation B",
       fill = TeX("$\\lambda$ Value")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12))


ggsave(paste0("../images/single_dataset_contrast_reject.png"), dpi = 300)


### Figure 3: Check for z-dependence in Simulation 3 with covariate matching --

target_files <- list.files(path = "../data", pattern = "sim[A-Z]") %>%
  extract_metadata() %>%
  filter(simulation == "C", match == "Covariate Match") %>%
  pull(file_name)

# Input all p-values in target_files into a single dataframe
p_val_df <- data.frame()

for (file in target_files)
{
  # Upload a single one
  pvals <- fread(paste0("../data/",file)) %>%
    .[,-1]
  
  data_matrix <- as.matrix(pvals)
  p_val_df_i <- cbind(extract_metadata(file), t(data_matrix))
  p_val_df <- rbind(p_val_df, p_val_df_i)
}

# Now pivot the dataframe to be in long format
value_columns <- names(p_val_df)[grepl("\\d", names(p_val_df))]

p_val_df_long <- pivot_longer(p_val_df, cols = value_columns) 

ggplot(p_val_df_long) +
  geom_histogram(aes(x = value), breaks = seq(0,1,by = 0.05)) +
  labs(x = TeX("Distribution of p-value")) +
  theme_minimal()

ggsave(paste0("../images/ecdf_covariate_match.png"), 
       height = 3,
       width = 6,
       dpi = 300)

ggplot(p_val_df_long) +
  geom_density(aes(x = value)) +
  theme_minimal()

hist(p_val_df_long$value)
ks.test(p_val_df_long$value,"punif")
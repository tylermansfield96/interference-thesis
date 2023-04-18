# interference-thesis

## Overview

This project was originally proposed by Sam Pimentel. This repository represents the work done by Tyler Mansfield (tyler.mansfield96@gmail.com) in the Spring 2023 as part of the completion of the Masters' thesis project. The manuscript of the thesis is contained in main directory, while the rest of the folders contain the relevant code, simulated data, results, images, and some of the literature used in the write-up.

## Running the code

There are three main R scripts used in this project. The first, `01_simulation_setup.R`, generates simulated data from a data-generating process defined by Chapter 3 in the manuscript. The script also contains code that implements the methodology proposed in Chapter 2. The second script, `02_run_simulations.R`, uses the code in the first script to repeatedly generate data and perform the randomization test for identifying interference (e.g. spillover), with the main output being p-value(s). Results from the script are written in the `data` folder. The final script `03_simulation_plots.R` uses the output of the second script stored in the `data` folder to generate the figures used in the manuscript. 

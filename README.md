# Wright-Malecot Model in R for Spatial Population Genetic Inference
Welcome! This repository includes the code and example data for running the Wright-Malecot model described in Hancock et al. (2023). This model will estimate a range of parameters, including Wright's neighborhood size (Nbhd), mutation rate (m), "s" (which is 1 - pi in the collecting phase), etc. 

You can find the preprint for this model here: https://doi.org/10.1101/2023.03.10.532094.

# What's in the repo

At the highest level you'll find the scripts for running the model locally: you'll need exe_WM.R, wm_lib.R, and wmModel_plots.R. These will call the model from the folder "models". The folder "hpcc_processing_slim" includes slurm scripts for running the model and processing output on a computing cluster (will need to be modified to fit your needs if you want to run them!). The "empirical_data" folder includes the dataset of Amphiprion bicinctus and how to process that dataset specifically. The "slim_recipes" folder includes the demographic models we ran for the manuscript. 

# How to run the model

Example datasets have been provided so that you can see how the model is performed and what the output looks like. The example datasets are example_locs.txt and example_pwp.csv, the first is geographic locations of samples in lat/long and the second is a matrix of pairwise pi. 

The file exe_WM_local.R will walk you through the steps - you will first source to the libraries you will need (all of which are in this repo), you will then convert the locations file to a pairwise distance matrix, and then you will build the STAN datablock. 

The STAN datablock has a few tricky things. First, you'll need to have an estimate of L, which is the number of independent loci in the dataset. Second, you'll need to provide a value of k - this is kappa in the manuscript, and represents the minimum distance between two samples upon which random mating can be assumed. The model is sensitive to arbitrarily high values, but not low - so, generally start with small numbers and increase if needed.

The model will build an Robj file that includes the posterior distribution of all of your estimates, and we also include a data visualization script that will automatically be called and produce a pdf that plots the empirical IBD curve and the model fit across each chain. 

# Questions?

If you have any questions about running this model, please contact Zach at hancockz (at) umich.edu. 
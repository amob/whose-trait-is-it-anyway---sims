# whose-trait-is-it-anyway---sims

This repository includes scripts associated with a manuscript. Detailed explanation of goals and output is available there. 
Script and final figure files (pdfs) are included. Some intermediate files (.txt, .csv) are included, but scripts will have to be re-run fully if changing any parameters.
In general, R scripts do the work, shell scripts are examples of how to run on a cluster, as some take a long time to run. 
  When R scripts have a paired shell script, they have matching names, differing only in ending.

Explanation of scripts

Scripts in series 01 provide simulation functions.
  Only host-micr-fitconfl_01_simfunction.R is used for analyses
  host-micr-fitconfl_01_simmultipopwithmigrate.R may be of interest to some readers, but is beyond the scope of the manuscript and remains a work in progress

Scripts in series 02 and 05 run simulations over different conditions. 
  host-micr-fitconfl_02_evoquestions_shortexfig.R runs the example simulations and makes most figures associated with these simulations
  host-micr-fitconfl_02_evoquestions_parametereffs.R explores effects of all parameters singly, (and host & microbial locus number in combination)
  host-micr-fitconfl_05_specialparametereffs.R explores effects of parameters altering evolutionary scenarios in combination (factorially)

Scripts in series 03 and 06 process simulation results
  e.g. calculate summary matrices. 
  Names link them to either exploration of all parameters ("parametereffs"), or evolutionary scenario parameters ("specialparametereffs").

Scripts 04 and 07 make figures for the exploration of results across parameter space.  
  Names link them to either exploration of all parameters ("testparamE"), or evolutionary scenario parameters ("specialparamE").

Scripts in series 08 may be of interest to some readers, but are beyond the scope of the manuscript
  These make use of host-micr-fitconfl_01_simmultipopwithmigrate.R and likewise remain a work in progress
 
Scripts in series 10 write the functions, and run analyses associated with GWAS results in manuscript
  host-micr-fitconfl_10_gwasfunctionsonly.R sets out the functions used in host-micr-fitconfl_10_gwasdemos.R
  host-micr-fitconfl_10_gwasdemos.R processeses simulation output from "02_shortexfig" and prepares for GWAS run in host-micr-fitconfl_10_gwasdemos.sh

Script in secies 11 processes GWAS output and makes associated figures in manuscript



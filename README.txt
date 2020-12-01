# whose-trait-is-it-anyway---sims

This repository includes scripts associated with the following article:

O’Brien, AM, C Jack, ML Friesen, and ME Frederickson.  Whose trait is it anyways?: Coevolution of joint phenotypes and genetic architecture in mutualisms. Proceedings of the Royal Society B: Biological Sciences. DOI: 10.1098/rspb.2020.2483 

Detailed explanation of goals and output is available therein. 
Script and final figure files (pdfs) are included. Some intermediate files (.txt, .csv) are included, but scripts will have to be re-run fully if changing any parameters.
In general, R scripts do the work, shell scripts are examples of how to run on a cluster, as some take a long time to run. 
  When R scripts have a paired shell script, they have matching names, differing only in ending.

EXPLANATION OF SCRIPTS

FitnessAlignmentConceptFig.R
  Produces panels for figure 1.
  Uses data extracted from Haney et al 2015 (please see manuscript for full citation), and available in fit_align.csv


Numbered scripts pertain to the simulations. 
**IMPORTANT EXPLANATION OF PARAMETER NOTATION. Some notation differs from manuscript to code
  - "theta" is written in the code as "Zopt" or "Zo" or similar
  - Code subscript for "host" is often "p" or "plant"; these should be read as equivalent to "host" or subscript "H"
  - "L" is "nL" for number of loci
  - greek letters are often written out by name (e.g. "Lambda") or substituted with a visually similar english alphabet letter
  - further abbreviations exist and are defined/described in code


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
 
Scripts in series 10 write the functions, and run analyses associated with GWAS results in manuscript
  host-micr-fitconfl_10_gwasfunctionsonly.R sets out the functions used in host-micr-fitconfl_10_gwasdemos.R
  host-micr-fitconfl_10_gwasdemos.R processeses simulation output from "02_shortexfig" and prepares for GWAS run in host-micr-fitconfl_10_gwasdemos.sh
  The output files from gemma are included in this directory ("assoc.txt" and .csv files named with pattern "micrloci" or "plantloci"), so that these very last steps can be run by interested parties, without having to re-run simulations (as these take awhile on a personal computer)

Script in series 11 processes GWAS output and makes associated figures in manuscript


ADDITIONAL FILES
The output files from gemma (10 series .R scripts) are included in this directory ("assoc.txt" and .csv files named with pattern "micrloci" or "plantloci")


###########
#Libraries#
###########
#Problems to fix
#1) Need to split/rename lipid objectives into brown and white
#2) Need to re-run results on white/brown after step 1)
#3) Change names in biomass barplot from reaction abbreviations to common metabolites
#4) Break up subscripts to be self-contained?

#Each of these sections is designed to be modular and standalone.

#Note, it's assumped the OS is linux
#Set the working directory.  Every output file/directory is always relative to this directory
wd <- c('/usr3/graduate/akram/metabolic_modeling/metabolic_modeling')
od <- c('./modeling_output_files/')
dir.create(od)
setwd(wd)

library(lattice)  #lattice must be explicitly called now to prevent a downstream error. Not sure why.
library(cplexAPI)
library(sybil)
library(sybilSBML)
library(sybilcycleFreeFlux)
library(parallel)
library(sqldf)
library(ggplot2)
library(reshape2)

###########
#Functions#
###########
source('~/metabolic_modeling/metabolic_modeling/metabolic_modeling_functions_linux.R')
source('~/metabolic_modeling/metabolic_modeling/sysBiolAlg_easyConstraintfvClass.R')
source('~/metabolic_modeling/metabolic_modeling/ACHR_efficient.R')
source('~/metabolic_modeling/metabolic_modeling/ACHR_warmup.R')
source('~/metabolic_modeling/metabolic_modeling/model_validation_functions/model_ATP_yield.R')

##########################
#Data and global settings#
##########################
options(stringsAsFactors=F) #Force of habit

#Change some global settings.
SYBIL_SETTINGS(c("SOLVER"), c("cplexAPI"))
cores <- 16

#Bring in the model
#Skip to loading the pread model (line 173 if not interested in making the original model)
#Note: This step can take up to 20 min for very large models (~70000 reactions)
Recon21x_sybil <- readSBMLmod(paste0(wd,"/Recon2_1x.xml"))
Recon2_sybil <- Recon21x_sybil

#Bring in some annotation and constraint data
annotation_reactions <- read.delim(file="sybil_Recon21x_model_react.tsv") 
DMEM_media_constraints <- read.csv(file="Recon2_constraints_media_adjusted_nonzero.csv")
DMEM_names <- read.csv(file=paste0("./output_files/DMEM_names.csv"))
blocked_reactions <- read.csv(file="Recon2_blocked_reactions.csv")
adipocyte_data <- read.csv(file="BA_WA_input_data.csv")
carbon_inputs <- read.csv(file="atp_theoretical_yields.csv")
escher_reactions <- read.csv(file="escher_reactions.csv")
nova_data <- read.csv(file="nova_inputs_model.csv")
lipid_objective <- read.csv(file="lipid_objective_coefficients.csv")

########################################################################
#Model Validation Of Recon 2.1X and Creating + Validating of Recon 2.1A#
########################################################################
#Running time: ~20 hrs with 16 cores.  Don't you just love FVA?
source('./Publication_model_validation.R')

#Load the SBML to skip the model validation and refinement
#This loaded model already has the media constraints on by default
pread_model <- readSBMLmod(paste0(wd,"/Recon_21A.xml"))

###################################################################
#FVA to reveal which fluxes are altered with seahorse measurements#
###################################################################
#Running time: ~4 hrs
source('./Publication_model_seahorse_base')

#################################################################################################
#Search for differences via sampling + MTF in brown and white adipocytes using the Seahorse data#
#################################################################################################
#Running time: ~40 min each with 16 cores
#Brown Adipocytes
source('./Publication_model_seahorse_brown_predictions.R')

#White Adipocytes
source('./Publication_model_seahorse_white_predictions.R')

#############################################################################################################
#Search for differences via sampling + MTF in brown and white adipocytes using the seahorse data + nova data#
#############################################################################################################
#Running time: ~40 min each with 16 cores
#Brown Adipocytes
source('./Publication_model_seahorse_nova_brown_predictions.R')

#White Adipocytes
source('./Publication_model_seahorse_nova_white_predictions.R')

########################################################################################
#Search for differences via FVA in brown and white Adipocytes with seahorse + nova data#
########################################################################################
#Running time: ~4 hrs each with 16 cores
#Brown Adipocytes
source('./Publication_model_seahorse_nova_brown_FVA_predictions.R')

#White Adipocytes
source('./Publication_model_seahorse_nova_white_FVA_predictions.R')

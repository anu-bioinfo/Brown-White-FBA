###########
#Libraries#
###########
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
pread_model <- readSBMLmod(paste0(wd,"/Recon_21A.xml"))

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

#############################
#Apply the media constraints#
#############################
#Note: The table was created for Recon2, not Recon2.1x, hence the indirect way of setting up the bounds
#Originally, the 24 hr approximation was used.  However, Brown adipocytes experimentally has higher flux
#For glucose uptake, thus the 3 hr approximation is used.
Recon2_constraints <- DMEM_media_constraints
DMEM_ids <- paste0(Recon2_constraints[,"Reaction_name"],"in")
DMEM_ids <- c(DMEM_ids, "carbon_uptake")
DMEM_lb <- rep(0, length(DMEM_ids))
DMEM_ub <- c(DMEM_media_constraints[,"Theoretical_max_flux..mmol.gDw.hr."],1000) #3hr approx
pread_model <- changeBounds(pread_model, DMEM_ids, lb=DMEM_lb, ub=DMEM_ub)

#########################
#Define the measurements#
#########################
exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "biomass_reaction", 
               "EX_lac_L(e)in", "EX_lac_L(e)ex")

adipocyte_data_original <- adipocyte_data
nsamples <- 150
pread_model <- rmReact(pread_model, react=blocked_reactions[,"abbreviation"]) #Reduces running time significantly
secretory_reactions <- react_id(pread_model)
ba_output <- matrix(nrow=length(secretory_reactions), ncol=nsamples)

#Sample
library(foreach)
library(doMC)
registerDoMC(cores=cores)

ba_output <- foreach(i=1:nsamples, .combine=cbind) %dopar% {
  tissue <- "hBA"
  atps4m_ba <- -1
  atp_leak_ba <- -1
  lpsolution <- 0
  
  #Before MTF, we need to make sure the sampled inputs and FBA make sense
  while(atps4m_ba <= 0 | atp_leak_ba <= 0 | lpsolution != 1) {
    adipocyte_data <- sample_seahorse_tseng(adipocyte_data_original)
    ocr_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"]
    
    #atps4m_ba must be postive
    atps4m_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_flux"] 
    
    atp_leak_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_leak"]
    ocr_mito_min_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    ocr_mito_max_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_fccp"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    h_sec <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"PPR_basal"]
    biomass_ba <- log(2)/76.8
    
    ba_lb <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_min_ba, biomass_ba,0, h_sec)
    ba_ub <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_max_ba, biomass_ba,0, h_sec)
    ba_model <- changeBounds(pread_model, react= exp_coefs, lb=ba_lb, ub=ba_ub)
    ba_model <- changeObjFunc(ba_model, "biomass_reaction", obj_coef=1)
    
    ba_test <- changeObjFunc(ba_model, "biomass_reaction", obj_coef=1)
    ba_test <- changeBounds(ba_test, c("PALFATPtc", "RTOTAL2FATPc_pmt", "RTOTALFATPc_pmt", "RTOTAL3FATPc_pmt"), lb=rep(0,4), ub=rep(0,4))
    ba_test_flux <- optimizeProb(ba_test, algorithm="fba", lpdir="max")
    
    lpsolution <- ba_test_flux@lp_stat
  }
  
  ba_fluxes <- optimizeProb(ba_model, alg="mtf", mtfobj=mod_obj(ba_test_flux))
  ba_fluxes_df <- getFluxDist(ba_fluxes)
  ba_output <- ba_fluxes_df
}
rownames(ba_output) <- secretory_reactions
write.csv(ba_output, file=paste0(od, "/model_seahorse_brown_predictions_v2.csv"))
save(ba_output, file=paste0(od, "/model_seahorse_brown_predictions.RData"))

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
library(parallel)
library(sqldf)

###########
#Functions#
###########
source('~/metabolic_modeling/metabolic_modeling/metabolic_modeling_functions_linux.R')

##########################
#Data and global settings#
##########################
options(stringsAsFactors=F) #Force of habit

#Change some global settings.
SYBIL_SETTINGS(c("SOLVER"), c("cplexAPI"))
SYBIL_SETTINGS(c("TOLERANCE"), 1e-9)
cores <- 16

#Bring in the model
#Note: This step can take up to 20 min for very large models (~70000 reactions)
pread_model <- readSBMLmod(paste0(wd,"/Recon_21A.xml"))

#Bring in some annotation and constraint data
annotation_reactions <- read.delim(file="sybil_Recon21x_model_react.tsv") 
DMEM_media_constraints <- read.csv(file="Recon2_media_constraints_revised.csv")
blocked_reactions <- read.csv(file="Recon2_blocked_reactions.csv")
adipocyte_data <- read.csv(file="BA_WA_seahorse_data_v2.csv")

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
               "EX_h(e)in", "EX_h(e)ex")

adipocyte_data_original <- adipocyte_data
nsamples <- 150
pread_model <- rmReact(pread_model, react=blocked_reactions[,"abbreviation"]) #Reduces running time significantly
secretory_reactions <- react_id(pread_model)

#Sample
library(foreach)
library(doMC)
registerDoMC(cores=cores)

ba_output <- foreach(i=1:nsamples, .combine=cbind) %dopar% {
  tissue <- "hWA"
  atps4m_ba <- -1
  atp_leak_ba <- -1
  lpsolution <- 0
  
  #Before MTF, we need to make sure the sampled inputs and FBA make sense (ie are the constraints feasible.  If not, we sample again)
  while(atps4m_ba <= 0 | atp_leak_ba <= 0 | lpsolution != 1) {
    adipocyte_data <- sample_seahorse_tseng(adipocyte_data_original)
    ocr_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"]
    
    #atps4m_ba must be postive
    atps4m_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_flux"] 
    
    atp_leak_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_leak"]
    ocr_mito_min_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    ocr_mito_max_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_fccp"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    h_sec <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"PPR_basal"]
    biomass_ba <- 0 #Used to be nonzero, but now zero
    
    ba_lb <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_min_ba, biomass_ba,0, h_sec)
    ba_ub <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_max_ba, biomass_ba,0, h_sec)
    ba_model <- changeBounds(pread_model, react= exp_coefs, lb=ba_lb, ub=ba_ub)
    ba_model <- changeObjFunc(ba_model, "ATPS4m", obj_coef=1)
    
    ba_test <- changeObjFunc(ba_model, "ATPS4m", obj_coef=1)
    ba_test <- changeBounds(ba_test, c("PALFATPtc", "RTOTAL2FATPc_pmt", "RTOTALFATPc_pmt", "RTOTAL3FATPc_pmt"), lb=rep(0,4), ub=rep(0,4))
    ba_test_flux <- optimizeProb(ba_test, algorithm="fba", lpdir="max")
    
    lpsolution <- ba_test_flux@lp_stat
  }
  
  ba_fluxes <- optimizeProb(ba_model, alg="mtf", mtfobj=mod_obj(ba_test_flux))
  ba_fluxes_df <- getFluxDist(ba_fluxes)
  ba_output <- ba_fluxes_df
}
rownames(ba_output) <- secretory_reactions
write.csv(ba_output, file=paste0(od, "/model_seahorse_white_predictions_v6.csv"))

###################
#Predict palmitate#
###################
#To make a prediction for palmitate, it must be added to the media.
pread_model <- changeBounds(pread_model, c("EX_hdca(e)in"), lb=0, ub=1000)

ba_output2 <- foreach(i=1:nsamples, .combine=cbind) %dopar% {
  tissue <- "hWA"
  atps4m_ba <- -1
  atp_leak_ba <- -1
  lpsolution <- 0
  
  #Before MTF, we need to make sure the sampled inputs and FBA make sense (ie are the constraints feasible.  If not, we sample again)
  while(atps4m_ba <= 0 | atp_leak_ba <= 0 | lpsolution != 1) {
    adipocyte_data <- sample_seahorse_tseng(adipocyte_data_original)
    ocr_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"]
    
    #atps4m_ba must be postive
    atps4m_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_flux"] 
    
    atp_leak_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_leak"]
    ocr_mito_min_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    ocr_mito_max_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_fccp"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    h_sec <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"PPR_basal"]
    biomass_ba <- 0 #Used to be nonzero, but now zero
    
    ba_lb <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_min_ba, biomass_ba,0, h_sec)
    ba_ub <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_max_ba, biomass_ba,0, h_sec)
    ba_model <- changeBounds(pread_model, react= exp_coefs, lb=ba_lb, ub=ba_ub)
    ba_model <- changeObjFunc(ba_model, "ATPS4m", obj_coef=1)
    
    ba_test <- changeObjFunc(ba_model, "ATPS4m", obj_coef=1)
    ba_test <- changeBounds(ba_test, c("PALFATPtc", "RTOTAL2FATPc_pmt", "RTOTALFATPc_pmt", "RTOTAL3FATPc_pmt"), lb=rep(0,4), ub=rep(0,4))
    ba_test_flux <- optimizeProb(ba_test, algorithm="fba", lpdir="max")
    
    lpsolution <- ba_test_flux@lp_stat
  }
  
  ba_fluxes <- optimizeProb(ba_model, alg="mtf", mtfobj=mod_obj(ba_test_flux))
  ba_fluxes_df <- getFluxDist(ba_fluxes)
  ba_output <- ba_fluxes_df
}
rownames(ba_output2) <- secretory_reactions
write.csv(ba_output2, file=paste0(od, "/model_seahorse_white_predictions_v6pmt.csv"))

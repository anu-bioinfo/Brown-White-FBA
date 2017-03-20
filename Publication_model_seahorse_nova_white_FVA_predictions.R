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
exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","EX_lac_L(e)in","EX_lac_L(e)ex", "ATPS4m", "DM_atp_m_", "O2tm", "biomass_reaction", 
               "EX_glc(e)in", "EX_glc(e)ex", "EX_hdca(e)in", "EX_hdca(e)ex", "EX_gln_L(e)in", "EX_gln_L(e)ex", 
               "EX_nh4(e)in", "EX_nh4(e)ex", "EX_glu_L(e)in", "EX_glu_L(e)ex","EX_h(e)in", "EX_h(e)ex", "EX_Rtotal2(e)ex_pmt",
               "EX_Rtotal(e)ex_pmt")

adipocyte_data_original <- adipocyte_data
nova_data_original <- nova_data

#the nsamples = 1 is to take advantage of the existing loop for FVA.
nsamples <- 1
pread_model <- rmReact(pread_model, react=blocked_reactions[,"abbreviation"])

secretory_reactions <- react_id(pread_model)

ba_output <- matrix(nrow=length(secretory_reactions), ncol=nsamples)
ba_upper <- matrix(nrow=length(secretory_reactions), ncol=nsamples)
ba_lower <- matrix(nrow=length(secretory_reactions), ncol=nsamples)
ba_input <- matrix(nrow=length(exp_coefs), ncol=nsamples)
tag_obj <- annotation_reactions[grep("^EX_tag_hs\\(e\\)ex", annotation_reactions[,"abbreviation"]),"abbreviation"]
tag_coef <- rep(1, length(tag_obj))

#White Adipocytes
i <- 1
while(i <= nsamples){
  tissue <- "hWA"
  atps4m_ba <- -1
  atp_leak_ba <- -1
  lpsolution <- 0
  
  while(atps4m_ba <= 0 | atp_leak_ba <= 0 | lpsolution != 1) {
    adipocyte_data <- adipocyte_data_original
    nova_data_s <- nova_data_original
    ocr_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"]
    
    #Lactate is tricky since BA may be postive or negative
    ppr_ba <- nova_data_s[nova_data_s[,"Tissue"] == tissue,"Lac_mean"]
    if(ppr_ba >= 0){
      ppr_ba_ex <- ppr_ba
      ppr_ba_in <- 0
    } else{
      ppr_ba_ex <- 0
      ppr_ba_in <- abs(ppr_ba)
    }
    
    #atps4m_ba must be postive
    atps4m_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_flux"] #Again, note the removal -ppr
    
    atp_leak_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_leak"]
    ocr_mito_min_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    ocr_mito_max_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_fccp"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
    glc_uptake <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"glucose_uptake"]
    pmt_uptake <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"Palmitate_uptake"]
    nh4_sec <- abs(nova_data_s[nova_data_s[,"Tissue"] == tissue,"NH4._mean"])
    gln_uptake <- abs(nova_data_s[nova_data_s[,"Tissue"] == tissue,"Gln_mean"])
    glu_sec <- abs(nova_data_s[nova_data_s[,"Tissue"] == tissue,"Glu_mean"])
    na1_uptake <- abs(nova_data_s[nova_data_s[,"Tissue"] == tissue,"Na._mean"])
    k_uptake <- abs(nova_data_s[nova_data_s[,"Tissue"] == tissue,"K._mean"])
    h_sec <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"PPR_basal"]
    
    biomass_ba <- log(2)/76.8
    ba_lb <- c(ocr_ba,0,ppr_ba_in,ppr_ba_ex,atps4m_ba, atp_leak_ba, ocr_mito_min_ba, biomass_ba, glc_uptake,0,pmt_uptake, 0,gln_uptake,0,0,nh4_sec,
               0, glu_sec, 0, h_sec,0,0)
    ba_ub <- c(ocr_ba,0,ppr_ba_in,ppr_ba_ex,atps4m_ba, atp_leak_ba, ocr_mito_max_ba, 1000, glc_uptake,0,pmt_uptake, 0, gln_uptake,0,0,nh4_sec,
               0, glu_sec, 0, h_sec,0,0)
    ba_model <- changeBounds(pread_model, react= exp_coefs, lb=ba_lb, ub=ba_ub)
    ba_model <- changeObjFunc(ba_model, "lipid_reaction_white", obj_coef=1)
    
    ba_test <- changeObjFunc(ba_model, "lipid_reaction_white", obj_coef=1)
    ba_test_flux <- optimizeProb(ba_test, algorithm="fba")
    
    lpsolution <- ba_test_flux@lp_stat
  }
  
  ba_fluxes <- multiDel(ba_model, nProc=cores,del1=secretory_reactions, todo="fluxVar")
  ba_fluxes_df <- lapply(ba_fluxes, fluxvar_rewrite)
  ba_fluxes_df <- do.call("rbind", ba_fluxes_df)
  ba_mean <- apply(ba_fluxes_df,1,mean)
  
  ba_output[,i] <- ba_mean
  ba_input[,i] <- ba_lb
  ba_upper[,i] <- ba_fluxes_df[,2] 
  ba_lower[,i] <- ba_fluxes_df[,1]
  i <- i + 1
}
rownames(ba_output) <- secretory_reactions
rownames(ba_upper) <- secretory_reactions
rownames(ba_lower) <- secretory_reactions
rownames(ba_input) <- exp_coefs
write.csv(ba_output, file=paste0(od, "/model_seahorse_nova_white_FVA_predictions_mean_matrix.csv"))
write.csv(ba_upper, file=paste0(od, "/model_seahorse_nova_white_FVA_predictions_ub_matrix.csv"))
write.csv(ba_lower, file=paste0(od, "/model_seahorse_nova_white_FVA_predictions_lb_matrix.csv"))
write.csv(ba_input, file=paste0(od, "/model_seahorse_nova_white_FVA_predictions_input_matrix.csv"))

save(ba_output, ba_input, ba_upper, ba_lower,
     file=paste0(od, "/model_seahorse_nova_white_FVA_predictions.RData"))

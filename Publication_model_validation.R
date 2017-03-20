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


##############################
#Shut off all exchange fluxes#
##############################
#First we restrict all the exchange input reactions to 0.
#Note:In this case, findExchReact will not find all exchange reactions due to the unique exchange definition in Recon 2.1x
#since there are boundary metabolites and carbon uptake.  
#A simple way to get the exchange reactions is to use the annotation file
#The annotation file was generated from the Recon2.1x SBML via the function model2tsv
closed_exchange_ids <- annotation_reactions[grepl("Exchange", annotation_reactions[,"subsystem"]),]
closed_exchange_ids <- closed_exchange_ids[closed_exchange_ids[,"compartment"] == "e",]
closed_exchange_ids <- closed_exchange_ids[grep("in|carbon_uptake", closed_exchange_ids[,"abbreviation"]),]
closed_exchange_ids_names <- closed_exchange_ids[,"abbreviation"]
closed_exchange_lb <- rep(0, length(closed_exchange_ids_names))
closed_exchange_ub <- rep(0, length(closed_exchange_ids_names))
Recon2_sybil_closed <- changeBounds(Recon2_sybil, closed_exchange_ids_names,lb=closed_exchange_lb, ub=closed_exchange_ub)

##################################
#Quantitative tests on Recon 2.1X#
##################################
#The goal here is to iterate through each carbon source and checking that we get the proper amount of cytosolic ATP
#First, open up the necessary minimal media constraints
bm <- c("ca2", "cl", "fe2", "fe3", "h2o", "k", "nh4","na1", "so4", "pi", "o2", "h")
bm <- c(paste0("EX_", bm, "(e)in"))
bm_ub <- rep(1000, length=length(bm))
bm_lb <- rep(0, length=length(bm))
Recon_21atp <- changeBounds(Recon2_sybil_closed, bm, lb=bm_lb, ub=bm_ub)

atp_yields <- model_ATP_yield(Recon_21atp, carbon_sources = carbon_inputs[,"carbon_source"])
yield_ratio <- round(atp_yields[,"Aerobic_ATP_yield"]/carbon_inputs[,"food_yield"], 2)
atp_results <- data.frame(carbon_inputs, round(atp_yields,2), yield_ratio)
write.csv(atp_results, file=paste0(od,"ATP_yield_recon21x.csv"), row.names=F)

#Calculate a correlation and make a plot
cor(atp_results[,"food_yield"], atp_results[,"Aerobic_ATP_yield"], method="pearson", use="pairwise.complete.obs")
atp_fit <- lm(atp_results[,"food_yield"] ~ atp_results[,"Aerobic_ATP_yield"])

pdf(file=paste0(od,"ATP_recon21x_correlation_plot.pdf"), height=7, width=7)
plot(atp_results[,"Aerobic_ATP_yield"], atp_results[,"food_yield"], xlim=c(0,1000), ylim=c(0,1000), pch=16,
     xlab="Model ATP Yield (mmol/gDw/hr)", ylab="Theoretical ATP Yield (mmol/gDw/hr)")
abline(a=0, b=1, col="red", lty=2)
dev.off()

#Test the every element in DMEM to determine which ones are essential.
Recon2_constraints <- DMEM_media_constraints
DMEM_ids <- paste0(Recon2_constraints[,"Reaction_name"],"in")
DMEM_ids <- c(DMEM_ids, "carbon_uptake")
DMEM_lb <- rep(0, length(DMEM_ids))
DMEM_ub <- rep(1000, length(DMEM_ids))
Recon2_sybil_mc <- changeBounds(Recon2_sybil_closed, DMEM_ids, lb=DMEM_lb, ub=DMEM_ub)

essential_reactions <- DMEM_ids
essential_results <- essential_test(Recon2_sybil_mc, essential_reactions)
essential_results <- sort(essential_results)

pdf(file=paste0(od,"Essential_DMEM_constituents_recon21x.pdf"), width=10, height=10)
par(mar=c(10,4,4,4))
barplot(essential_results, col="black", las=2, ylab="biomass (mmol/gDw/hr)")
dev.off()

#Note: The next two parallelized FVA steps have significant running time (~8 hrs each)
#With zero input, check the range of all fluxes (FVA)
reacts_t0 <- react_id(Recon2_sybil_closed)
fluxvar_no_input <- multiDel(Recon2_sybil_closed, del1=reacts_t0, nProc=cores, todo="fluxVar")
fluxvar_no_input_df <- lapply(fluxvar_no_input, fluxvar_rewrite)
fluxvar_no_input_df <- do.call("rbind", fluxvar_no_input_df)
write.csv(fluxvar_no_input_df, file=paste0(od,"FVA_no_media_recon21x.csv"))

#Conversely, let's test for blocked reactions by opening every exchange reaction (input and output)
open_exchange_ids <- annotation_reactions[grepl("Exchange|carbon_uptake", annotation_reactions[,"subsystem"]),]
open_exchange_ids <- open_exchange_ids[open_exchange_ids[,"compartment"] == "e",]
open_exchange_ids_names <- open_exchange_ids[,"abbreviation"]
open_exchange_lb <- rep(0, length(open_exchange_ids_names))
open_exchange_ub <- rep(1000, length(open_exchange_ids_names))
Recon2_sybil_open <- changeBounds(Recon2_sybil_closed, open_exchange_ids_names, lb=open_exchange_lb, ub=open_exchange_ub)

blocked_reactions <- multiDel(Recon2_sybil_open, del1=react_id(Recon2_sybil_open), nProc=cores, todo="fluxVar")
blocked_reactions_df <- lapply(blocked_reactions, fluxvar_rewrite)
blocked_reactions_df <- do.call("rbind", blocked_reactions_df)
write.csv(blocked_reactions_df, file=paste0(od,"FVA_full_media_recon21x.csv"))

#####################
#Creating Recon 2.1A#
#####################
#Reactions to shutoff to give reasonable ATP yields per substrate through the textbook metabolic reactions
or <- c("r0081","r0122", "r0153", "r0165", "r0173", "r0202", "r0280", "r0354", "r0355", "r0413", 
        "r0509",  "r2520", "r0885", "ICDHyrm", "r0016", "r1453", "P450SCC1m", "BAAT3x", "SUCOASm", "FUMtm",
        "PALFATPtc", "RTOTAL2FATPc_pmt", "RTOTALFATPc_pmt", "RTOTAL3FATPc_pmt")
or <- c(or, annotation_reactions[grep("DNDPt", annotation_reactions[,"abbreviation"], ignore.case=T),"abbreviation"])
or_lb <- rep(0,length(or))
or_ub <- rep(0,length(or))
Recon2_sybil_closed <- changeBounds(Recon2_sybil_closed, react=or, lb=or_lb, ub=or_ub)

ir <- c("carbon_uptake", "GLUt2m", "ACACt2m", "ACETONEt2m", "PIt2m")
ir_ub <- rep(1000, length=length(ir))
ir_lb <- rep(0, length=length(ir))
Recon2_sybil_closed <- changeBounds(Recon2_sybil_closed, ir, lb=ir_lb, ub=ir_ub)

#Now let's fix the biomass-associated mistakes (just phenylalanine):
Recon2_sybil_closed <- changeBounds(Recon2_sybil_closed, react="r0399", lb=0, ub=1000)

##################################################
#Adding the missing reactions and lipid objective#
##################################################
#Add a mitochondrial ATP demand to simulate proton leak and add the bounds
pread_model <- addReact(Recon2_sybil_closed, id="DM_atp_m_", met=c("atp[m]","h2o[m]","adp[m]","h[m]","pi[m]"),
                        Scoef=c(-1,-1,1,1,1), ub=0, lb=0, subSystem="Exchange/demand reaction", metComp=rep("m", 5))

#Let's add the necessary reactions for lipid mass generation
#Triglycerides to generic
tags <- annotation_reactions[grep("^EX_tag_hs\\(e\\)ex", annotation_reactions[,"abbreviation"]),"abbreviation"]
for(i in 1:length(tags)){
  tags[i] <- names(get_metabolites(pread_model, tags[i]))
}
tags_id <- gsub("\\[e\\]", "", tags)
tags <- gsub("\\[e\\]", "\\[c\\]", tags)

for(i in 1:length(tags)){
  pread_model <- addReact(pread_model, id = paste0("Generic_", tags_id[i]), met = c(tags[i], "tag_hs[c]"), lb=0, ub=1000, obj=0,
                          Scoef = c(-1,1), subSystem="Generic TAG", reactName="Generic TAG reaction", metComp=rep("c", 2)) 
}

#Diglycerides to generic
dags <- annotation_reactions[grep("^EX_dag_hs\\(e\\)ex", annotation_reactions[,"abbreviation"]),"abbreviation"]
for(i in 1:length(dags)){
  dags[i] <- names(get_metabolites(pread_model, dags[i]))
}
dags_id <- gsub("\\[e\\]", "", dags)
dags <- gsub("\\[e\\]", "\\[c\\]", dags)

for(i in 1:length(dags)){
  pread_model <- addReact(pread_model, id = paste0("Generic_", dags_id[i]), met = c(dags[i], "dag_hs[c]"), lb=0, ub=1000, obj=0,
                          Scoef = c(-1,1), subSystem="Generic DAG", reactName="Generic DAG reaction", metComp=rep("c", 2)) 
}

#Ceramides to generic
crms <- annotation_reactions[grep("^EX_crm_hs\\(e\\)ex", annotation_reactions[,"abbreviation"]),"abbreviation"]
for(i in 1:length(crms)){
  crms[i] <- names(get_metabolites(pread_model, crms[i]))
}

crms_id <- gsub("\\[e\\]", "", crms)
crms <- gsub("\\[e\\]", "\\[c\\]", crms)

for(i in 1:length(crms)){
  pread_model <- addReact(pread_model, id = paste0("Generic_", crms_id[i]), met = c(crms[i], "crm_hs[c]"), lb=0, ub=1000, obj=0,
                          Scoef = c(-1,1), subSystem="Generic CRM", reactName="Generic Ceramide reaction", metComp=rep("c", 2)) 
}

#FFAs to generic
ffas <- c("td[c]", "pmt[c]", "hd[c]", "st[c]", "ode[c]", "lnlc[c]", "lnlnca[c]", "lnlncg[c]")
ffas_id <- gsub("\\[c\\]", "", ffas)

for(i in 1:length(crms)){
  pread_model <- addReact(pread_model, id = paste0("Generic_", ffas_id[i]), met = c(ffas[i], "ffa_hs[c]"), lb=0, ub=1000, obj=0,
                          Scoef = c(-1,1), subSystem="Generic CRM", reactName="Generic Ceramide reaction", metComp=rep("c", 2)) 
}

#And finally, the lipid mass reaction
pread_model <- addReact(pread_model, id="lipid_reaction_white", met=lipid_objective[,"metabolites"], lb=0, ub=1000, obj=0, 
                        Scoef=-lipid_objective[,"sat_comp"], subSystem="Lipid_mass", reactName="Lipid mass reaction white", 
                        metComp=rep("c", length(lipid_objective[,"metabolites"])))

pread_model <- addReact(pread_model, id="lipid_reaction_brown", met=lipid_objective[,"metabolites"], lb=0, ub=1000, obj=0, 
                        Scoef=-lipid_objective[,"bat_comp"], subSystem="Lipid_mass", reactName="Lipid mass reaction brown", 
                        metComp=rep("c", length(lipid_objective[,"metabolites"])))

#Change some annotation
pread_model@mod_name <- c("Recon 2.1A")
save(pread_model, file="Recon21A_robject.RData")

#Now that all's said and done,save the SBML
#Note 1: the function fails for level 3 annotation.
#Note 2: the function also fails for printNotes=printAnnos=T
#Note 3: I figured it it out: for level=3, printNotes=printAnnos=T, you must set fbcLevel=2.  The manual lies about
#fbcLevel=2.  
#Note 4: Nope, it failed again.  The actual problem is in printAnnos=T when NAs exit in those dataframes, which
#eventually produces a segfault.  The "fix" is to replace those NAs with the relevant annotation, or just "" out of laziness.
#pread_model is recon2.1A
pread_model@react_attr[is.na(pread_model@react_attr)] <- ""
pread_model@met_attr[is.na(pread_model@met_attr)] <- ""
pread_model@met_attr[,"charge"] <- "0"

writeSBML(morg=pread_model, filename="Recon_21A.xml", level=3, version=1, printNotes=T, printAnnos=T, validation=F, fbcLevel=2)

##################################
#Quantitative tests on Recon 2.1A#
##################################
#The goal here is to iterate through each carbon source and checking that we get the proper amount of cytosolic ATP
#First, open up the necessary minimal media constraints
Recon21A_sybil_closed <- changeBounds(pread_model, closed_exchange_ids_names,lb=closed_exchange_lb, ub=closed_exchange_ub)
bm <- c("ca2", "cl", "fe2", "fe3", "h2o", "k", "nh4","na1", "so4", "pi", "o2", "h")
bm <- c(paste0("EX_", bm, "(e)in"))
bm <- c(bm, "carbon_uptake")
bm_ub <- rep(1000, length=length(bm))
bm_lb <- rep(0, length=length(bm))
Recon_21atp <- changeBounds(Recon21A_sybil_closed, bm, lb=bm_lb, ub=bm_ub)

atp_yields <- model_ATP_yield(Recon_21atp, carbon_sources = carbon_inputs[,"carbon_source"])
yield_ratio <- round(atp_yields[,"Aerobic_ATP_yield"]/carbon_inputs[,"food_yield"], 2)
atp_results <- data.frame(carbon_inputs, round(atp_yields,2), yield_ratio)
write.csv(atp_results, file=paste0(od,"ATP_yield_recon21A.csv"), row.names=F)

#Calculate a correlation and make a plot
cor(atp_results[,"food_yield"], atp_results[,"Aerobic_ATP_yield"], method="pearson", use="pairwise.complete.obs")
atp_fit <- lm(atp_results[,"food_yield"] ~ atp_results[,"Aerobic_ATP_yield"])

pdf(file=paste0(od,"ATP_recon21A_correlation_plot.pdf"), height=7, width=7)
plot(atp_results[,"Aerobic_ATP_yield"], atp_results[,"food_yield"], xlim=c(0,150), ylim=c(0,150), pch=16,
     xlab="Model ATP Yield (mmol/gDw/hr)", ylab="Theoretical ATP Yield (mmol/gDw/hr)")
abline(a=0, b=1, col="red", lty=2)
abline(a= atp_fit[[1]][[1]], b=atp_fit[[1]][[2]])
dev.off()

#Test the every element in DMEM to determine which ones are essential.
Recon2_constraints <- DMEM_media_constraints
DMEM_ids <- paste0(Recon2_constraints[,"Reaction_name"],"in")
DMEM_ids <- c(DMEM_ids, "carbon_uptake")
DMEM_lb <- rep(0, length(DMEM_ids))
DMEM_ub <- rep(1000, length(DMEM_ids))
Recon2_sybil_mc <- changeBounds(Recon21A_sybil_closed, DMEM_ids, lb=DMEM_lb, ub=DMEM_ub)

essential_reactions <- DMEM_ids
essential_results <- essential_test(Recon2_sybil_mc, essential_reactions)
essential_results <- sort(essential_results)

pdf(file=paste0(od,"Essential_DMEM_constituents_recon21A.pdf"), width=10, height=10)
par(mar=c(10,4,4,4))
barplot(essential_results, col="black", las=2, ylab="biomass (mmol/gDw/hr)")
dev.off()

#Note: The next two parallelized FVA steps have significant running time (~8 hrs each)
#With zero input, check the range of all fluxes (FVA)
reacts_t0 <- react_id(Recon21A_sybil_closed)
fluxvar_no_input <- multiDel(Recon21A_sybil_closed, del1=reacts_t0, nProc=cores, todo="fluxVar")
fluxvar_no_input_df <- lapply(fluxvar_no_input, fluxvar_rewrite)
fluxvar_no_input_df <- do.call("rbind", fluxvar_no_input_df)
write.csv(fluxvar_no_input_df, file=paste0(od,"FVA_no_media_recon21A.csv"))

#Conversely, let's test for blocked reactions by opening every exchange reaction (input and output)
open_exchange_ids <- annotation_reactions[grepl("Exchange|carbon_uptake", annotation_reactions[,"subsystem"]),]
open_exchange_ids <- open_exchange_ids[open_exchange_ids[,"compartment"] == "e",]
open_exchange_ids_names <- open_exchange_ids[,"abbreviation"]
open_exchange_lb <- rep(0, length(open_exchange_ids_names))
open_exchange_ub <- rep(1000, length(open_exchange_ids_names))
Recon2_sybil_open <- changeBounds(Recon21A_sybil_closed, open_exchange_ids_names, lb=open_exchange_lb, ub=open_exchange_ub)

blocked_reactions <- multiDel(Recon2_sybil_open, del1=react_id(Recon2_sybil_open), nProc=cores, todo="fluxVar")
blocked_reactions_df <- lapply(blocked_reactions, fluxvar_rewrite)
blocked_reactions_df <- do.call("rbind", blocked_reactions_df)
write.csv(blocked_reactions_df, file=paste0(od,"FVA_full_media_recon21A.csv"))

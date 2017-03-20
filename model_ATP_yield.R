#This function iterates through each provided carbon source and returns a dataframe of the results
model_ATP_yield <- function(model, carbon_sources){
  #Check for carbon sources not found in the model
  no_source <- carbon_sources[!(carbon_sources %in% react_id(model))]
    if(length(no_source) > 0){
      print(paste("The following metabolites were not found in the model:", no_source))
    }
    
  atp_yield <- vector("numeric", length=length(carbon_sources))
  atp_yield_an <- vector("numeric", length=length(carbon_sources))
  
  names(atp_yield) <- carbon_sources
  model_atp <- changeObjFunc(model,react="DM_atp_c_", obj_coef=1)
  
  #Aerobic conditions
  for(i in 1:length(carbon_sources)){
    if(carbon_sources[i] %in% no_source){
      atp_yield[i] <- NA
      next
    }
    model_temp <- changeBounds(model_atp, carbon_sources[i], lb=0, ub=1)
    atp_flux <- optimizeProb(model_temp, algorithm="fba")
    atp_yield[i] <- mod_obj(atp_flux)
    print(paste0("Finished ",carbon_sources[i]))
  }
  
  #Anaerobic conditions
  model_atp <- changeBounds(model_atp, "EX_o2(e)in", lb=0, ub=0)
  for(i in 1:length(carbon_sources)){
    if(carbon_sources[i] %in% no_source){
      atp_yield_an[i] <- NA
      next
    }
    model_temp <- changeBounds(model_atp, carbon_sources[i], lb=0, ub=1)
    atp_flux <- optimizeProb(model_temp, algorithm="fba")
    atp_yield_an[i] <- mod_obj(atp_flux)
    print(paste0("Finished ",carbon_sources[i]))
  }
  final_atp <- cbind(atp_yield, atp_yield_an)
  colnames(final_atp) <- c("Aerobic_ATP_yield", "Anaerobic_ATP_yield")
  final_atp
}
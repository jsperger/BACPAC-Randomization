
#' Calculates study-wide maximum treatment imbalance, 
#' study-wide maximum treatment imbalance by factor level,
#' site-level treatment imbalance, and
#' site by covariate level imbalance
#'
#' @param indf data frame with study data. Site must be in a column named Site 
#' and treatment must be in a column named A1
#' @param

AnalyzeImbalances <- function(indf, vars.to.check){
  imbalance_across_study <- vector("numeric", length = length(vars.to.check))
  
  ImbalanceHelper <- function(indf, cur.var){
    allocation_table <- indf %>% group_by(!!sym(cur.var)) %>% group_map( .f = ~table(.$A1))
    imbalance_by_level <- map_dbl(allocation_table, ~diff(range(.)))
    
    return(imbalance_by_level)
  }
  
  
  study_level_imbalances <- map(vars.to.check, ~ImbalanceHelper(indf, .)) %>% unlist
  
  cluster_level_imbalances <- indf %>%  group_by(Site) %>% 
    group_map(., .f = function(x, ...) map(vars.to.check, ~ImbalanceHelper(indf = x, .) %>% unlist)) %>% 
    unlist
  
  cluster_treatment_imbalances <- indf %>% group_by(Site) %>% 
    group_map(., ~max(diff(range(table(.$A1))))) %>% 
    unlist
  
  imbalance_summary <- c("MaxTreatmentImbalance" = max(diff(range(table(indf$A1)))),
                         "MaxStudywideCovImbalance" = max(study_level_imbalances),
                         "MaxSiteTreatmentImbalance" = max(cluster_treatment_imbalances),
                         "MaxSiteCovImbalance" = max(cluster_level_imbalances))  
  return(imbalance_summary)
}

#' Using data.table because map is slow
#' 
#' @param data.w.assignments
#' @param vars.to.check
AnalyzeImbalancesDT <- function(data.w.assignments, vars.to.check){
  overall_imbalances <- data.w.assignments[,.(OverallImbalance = diff(range(table(A1)))), by = .(Setting, Replication)]
  site_imbalances <- data.w.assignments[,.(SiteImbalance = diff(range(table(A1)))), by = .(Setting, Replication, Site)][,. (OverallSiteImbalance = max(SiteImbalance)), by = .(Setting, Replication)]
  
  
  covar_imbalances <- .DTCalcImbalanceByCovar(data.w.assignments, vars.to.check[1])
  
  for(cur.var in vars.to.check[-1]){
    temp_imbalance_df <- .DTCalcImbalanceByCovar(data.w.assignments, cur.var)
    
    covar_imbalances <- merge.data.table(covar_imbalances, temp_imbalance_df, all = TRUE, key = "Setting, Replication")
  }
  
  largest_covar_imbalances <- covar_imbalances[, .(LargestCovarImbalance = max(.SD)), by = .(Setting, Replication)]
  imbalance_summaries <- merge.data.table(overall_imbalances, site_imbalances, all = TRUE, key = "Setting, Replication")
  imbalance_summaries <- merge.data.table(overall_imbalances, covar_imbalances, all = TRUE, key = "Setting, Replication")
  
  imbalance_summaries <- merge.data.table(overall_imbalances, largest_covar_imbalances, all = TRUE, key = "Setting, Replication")
  
  return(imbalance_summaries)
}

.DTCalcOverallImbalance <- function(sim.dt){
  overall_imbalances <- sim.dt[,.(OverallImbalance = diff(range(table(A1)))), by = .(Setting, Replication)]
  
  return(overall_imbalances)
}

.DTCalcSiteImbalance <- function(sim.dt){
  site_imbalances <- sim.dt[,.(SiteImbalance = diff(range(table(A1)))), 
                            by = .(Setting, Replication, Site)][,
                                                                .(OverallSiteImbalance = max(SiteImbalance)), by = .(Setting, Replication)]
  
  
  return(site_imbalances)
}

.DTCalcImbalanceByCovar <- function(sim.dt, var.name){
  imbalance_var_name <- paste0("LargestFactorImb", var.name)
  
  overall_imbalance_by_covar <- sim.dt[,
                                       .(FactorImbalance = diff(range(table(A1)))), 
                                       by = .(Setting, Replication, get(var.name))][
                                         ,.(LargestFactorImb = max(FactorImbalance)), 
                                         by = .(Setting, Replication)]
  
  setnames(overall_imbalance_by_covar, "LargestFactorImb", imbalance_var_name)
  
  return(overall_imbalance_by_covar)
}


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
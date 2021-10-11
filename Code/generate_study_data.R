#' Generate covariate data for a study. Does not include treatment assignments
#' @param i iteration number. A workaround parameter to facilitate working with map functions
#' @param n.subj number of subjects
#' @param k.arms number of arms
#' @param n.sites number of sites
#' @param bin.props vector of probabilities of successes for each binary covariate. 
#' The length of the vector determines the number of binary covariates
#' @param bin.iccs vector of intra-cluster correlation coefficients the same length as bin.props, or a single constant that will be used for all covariates
#' 
#' @return a tibble with the study data
GenerateStudyCovariateData <- function(i = 1, 
                                       n.subj = 740, 
                                       k.arms = 4, 
                                       n.sites = 10, 
                                       bin.props = runif(5, min = .1, max = .9), 
                                       bin.iccs = .1,
                              ...) {
  
  study_data <- tibble(ID = 1:n.subj,
                       Site = sample(1:n.sites, size = n.subj, replace = TRUE)) %>%
    MakeCorrelatedBinaryCovars(indf = ., bin.props = bin.props, bin.iccs = bin.iccs)
  
  covar_matrix <- study_data %>% 
    select(Site, starts_with("X_")) %>% 
    mutate_all(factor) %>%  mutate_all(as.numeric)  %>% #balanceRandomize function is expecting factor levels starting at 1 not 0 or -1
    as.matrix
  
  study_data <- AddContraVar(study.data = study_data, ...)
  
  return(study_data)
}

#' Add assignments to study data created by \code{GenerateStudyCovariateData}
#' 
#' @param study.data data frame containing study data
#' @param n.sites integer specifying the number of study sites
#' @param k.arms  integer specifying the number of arms in the study
#' @param stratify.by.site logical indicating whether to stratify by site before applying minimization algorithm
GenerateAssignments <- function(study.data, n.sites, k.arms, stratify.by.site = FALSE, ...){
  covar_matrix <- study.data %>% 
    select(Site, starts_with("X_")) %>% 
    mutate_all(factor) %>%  mutate_all(as.numeric)  %>% #balanceRandomize function is expecting factor levels starting at 1 not 0 or -1
    as.matrix
  
  if(stratify.by.site == TRUE){
    study.data$A1 <- NA
    for(site in 1:n.sites){
      site_participant_lgl <- study.data$Site == site
      site_assignments <- MinimRandomize(covariate.mat = covar_matrix[site_participant_lgl,], 
                                         contraindications = study.data$Contra[site_participant_lgl], 
                                         k.arms = k.arms,  
                                         ...)$Assignments
      study.data$A1[site_participant_lgl] <- site_assignments
    }
  } else{
    study.data$A1 <- MinimRandomize(covariate.mat = covar_matrix, 
                                    contraindications = study.data$Contra, 
                                    k.arms = k.arms,  
                                    ...)$Assignments
  }
  
  return(study.data)
}


#' Create cluster-correlated binary covariates using the \code{fabricatr} package. 
#' @param indf the current data frame
#' @param bin.props vector of probs; probabilities (proportions) with binary
#' covariate = 1
#' @param bin.iccs vector of intra-cluster correlation coefficients the same length as bin.props, or a single constant that will be used for all covariates
#' 
#' @return the input data frame with the covariates appended
MakeCorrelatedBinaryCovars <- function(indf, bin.props, bin.iccs){
  stopifnot("Site" %in% colnames(indf))
  num_covars <- length(bin.props)
  
  if(num_covars == 0) return(indf)
  
  # Recycle bin.iccs so its length matches num_covars
  if(length(bin.iccs) < num_covars){
    if(length(bin.iccs) > 1) warning("Length of bin.iccs and bin.props do not match and bin.iccs is not a constant. Recycling bin.iccs, verify this mismatch is intended")
    bin.iccs <- rep(bin.iccs, length.out = num_covars)
  }
  
  # Create covariates
  bin_covar_mat <- matrix(nrow = nrow(indf), ncol = num_covars)
  
  for(i in 1:num_covars){
    bin_covar_mat[,i] <- fabricatr::draw_binary_icc(prob = bin.props[i],
                                                    clusters = indf$Site,
                                                    ICC = bin.iccs[i])
  }
  
  # Naming
  covar_names <- paste("X", 1:num_covars, sep = "_")
  colnames(bin_covar_mat) <- covar_names
  
  df_with_bin_covars <- bind_cols(indf, as_tibble(bin_covar_mat))
  
  return(df_with_bin_covars)
}

#' Generate contraindications and append the variable to an existing data frame
#' @param study.data data frame with the current study data
#' @param contra.probs vector of contraindication probabilities. Should equal the number of arms in the study.
#' Or if there are no contraindications a single 0 is acceptable
#' 
#' @return original data frame with a column containing the contraindications appended
AddContraVar <- function(study.data, contra.probs, ...){
  if(length(contra.probs) == 1) {
    study.data$Contra <- rbinom(n = nrow(study.data), size = 1, p = contra.probs)
    return(study.data)
  }
  contraindications_draw <- rmultinom(n = nrow(study.data), size = 1, prob = contra.probs)
  # Extract the actual result and make it start at zero (indicating no contraindications)
  contraindications <- apply(contraindications_draw, 2, which.max) - 1
  
  study.data$Contra <- contraindications
  
  return(study.data)
}
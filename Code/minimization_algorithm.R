#' Randomize patients using Pocock's Minimization method with a biased coin
#' 
#' 
#' @param covariate.mat
#' @param k.arms
#' @param imbalance.func logical indicating whether to use a formula which includes a term for the overall treatment imbalance or not
#' @param use.overall.imbalance
#' @param covar.weights
#' @param coin.bias
#' @param other.arm.p.func
#' 
#' 
MinimRandomize <- function(covariate.mat, contraindications, k.arms, imbalance.func =  function(x) diff(range(x)), 
                           use.overall.imbalance = FALSE, 
                           covar.weights = 1,
                                    coin.bias = 2/3, other.arm.p.func = .MinimEqualProb, ...){
  
  .MinimRandomizeCheckInputs(covariate.mat = covariate.mat,
                             k.arms = k.arms,
                                 coin.bias = coin.bias)
  
  # Preprocess inputs into desired format
  if(is.data.frame(covariate.mat) == TRUE) covariate.mat <- covariate.mat %>% mutate_all(as.numeric) %>% as.matrix
  if(sum(covar.weights != 1)) covar.weights <- covar.weights / sum(covar.weights)
  
  # Define useful constants
  n_subj <- nrow(covariate.mat)
  i_factors <- ncol(covariate.mat)
  max_fac_levels <- apply(covariate.mat, 2, max)
  arm.vec <- 1:k.arms
    
  # Pre-allocate resources
  allocation_vec <- vector("integer", length = n_subj)
  assignment_probs_before_rand <- matrix(nrow = n_subj, ncol = k.arms)
  assignments_by_covars <- array(data = NA, dim = c(i_factors, max(max_fac_levels), k.arms )) # gives number of subjects with level j of factor i in group k
  
  
  for (i in 1:i_factors) assignments_by_covars[i, 1:max_fac_levels[i],] <- 0 #0 for feasible levels of factor, NA otherwise
  
  hypothetical_imbalances <- vector("integer", length = i_factors)
  weighted_imbalance_scores <- vector("numeric", length = k.arms)

  for(cur_subj in 1:n_subj){
    
    subj_covars <- covariate.mat[cur_subj,]
    eligible_arms <- setdiff(arm.vec, contraindications[cur_subj])
    assignment_probs <- rep(0, times = k.arms)
    
    for(arm in 1:k.arms){
      hypo_factor_counts <- assignments_by_covars
      
      for(covariate in 1:i_factors){
      #  print(hypo_factor_counts[covariate, subj_covars[covariate], arm])
        
        hypo_factor_counts[covariate, subj_covars[covariate], arm] <- hypo_factor_counts[covariate, subj_covars[covariate], arm] + 1
        
    #    print(hypo_factor_counts[covariate, subj_covars[covariate], arm])
     #   print(imbalance.func(hypo_factor_counts[covariate, subj_covars[covariate],]))
        hypothetical_imbalances[covariate] <- imbalance.func(hypo_factor_counts[covariate, subj_covars[covariate],])
      }
      
      if(use.overall.imbalance == TRUE){
        # TODO: Better workaround to ensure that 1:k.arms isn't empty. 
        #Adding pseudovalues doesn't change range-based imbalance functions but could change variance-based ones
        
        current_allocations <- table(c(allocation_vec[1:(cur_subj-1)], 1:k.arms))
        
        if(cur_subj %% 100 == 0) print(current_allocations)
        
        hypothetical_allocations <- current_allocations
        hypothetical_allocations[arm] <- hypothetical_allocations[arm] + 1
        
        if(cur_subj %% 100 == 0) print(hypothetical_allocations)
        
        
        overall_imbalance <- imbalance.func(hypothetical_allocations)
        
        if(cur_subj %% 100 == 0) print(overall_imbalance)
        
      }
      
      weighted_imbalance_scores[arm] <- ifelse(use.overall.imbalance == TRUE,
                                               sum(covar.weights * c(overall_imbalance, hypothetical_imbalances)),
                                               sum(covar.weights * hypothetical_imbalances))
    }
    print(weighted_imbalance_scores)
    # Remove arms that are contraindicated
    feasible_imbalance_scores <- weighted_imbalance_scores[eligible_arms]
    # Which arm(s) minimizes the weighted imbalance score
    minimizing_arms <- which(feasible_imbalance_scores == min(feasible_imbalance_scores))
    
    #Break ties randomly
    preferred_arm <- ifelse(length(minimizing_arms == 1), minimizing_arms,
                            sample(x = minimizing_arms, size = 1)) #if x has length 1 sample samples from 1:x which is not desired
      
    
    other_arm_probs <- other.arm.p.func(feasible_imbalance_scores, coin.bias)
    
    # Map back to arm IDs
    preferred_arm_id <- eligible_arms[preferred_arm]
    other_arm_ids <- eligible_arms[-preferred_arm]
    
    assignment_probs[preferred_arm_id] <- coin.bias
    assignment_probs[other_arm_ids] <- other_arm_probs
    
    assignment_probs_before_rand[cur_subj,] <- assignment_probs
    
    assigned_arm <- sample(1:k.arms, size = 1, prob = assignment_probs)
    allocation_vec[cur_subj] <- assigned_arm
    
    for(covariate in 1:i_factors){
      assignments_by_covars[covariate, subj_covars[covariate], assigned_arm] <- assignments_by_covars[covariate, subj_covars[covariate], assigned_arm] + 1
    }
  }
  
  allocation_list <- list(Assignments = allocation_vec,
                          CovarInfo = assignments_by_covars,
                          Probs = assignment_probs_before_rand)
  
  return(allocation_list)
}

################################################################################
## Helper Functions
##
################################################################################


#' Helper function to check the inputs to \code{MinimRandomize}
.MinimRandomizeCheckInputs <- function(covariate.mat, k.arms,
                                           coin.bias){
  
  stopifnot(coin.bias <= 1 & coin.bias >= 0)
  stopifnot(is.numeric(k.arms) & k.arms %% 1 == 0)
  
  
  return(TRUE)
}


StratifyThenMinimizeRand <- function(covariate.mat, contraindications, duloxetine.grp.num = 2, k.arms, ...){
  covariate_mat_w_contras <- cbind(contraindications, covariate.mat)
  
  duloxetine_eligible_subjects <- covariate_mat_w_contras[contraindications != duloxetine.grp.num,]
  duloxetine_ineligible_subjects <- covariate_mat_w_contras[contraindications == duloxetine.grp.num,]
  
  dulox_eligible_assignments <- MinimRandomize(covariate.mat = duloxetine_eligible_subjects, contraindications = contraindications[contraindications != 1],
                                               k.arms = k.arms, 
                                               ...)
  
  dulox_ineligible_assignments <- MinimRandomize(covariate.mat = duloxetine_ineligible_subjects, 
                                                 contraindications = contraindications[contraindications == 1],
                                                 k.arms = k.arms,
                                               ...)
  
  # Make the ineligible contraindication covariate numbering match 
  dulox_inel_counts <- c(table(dulox_ineligible_assignments$Assignments))
  dulox_ineligible_assignments$CovarInfo[1,1:k.arms,] <- 0
  dulox_ineligible_assignments$CovarInfo[1,duloxetine.grp.num,] <- dulox_inel_counts
  
  
  combined_assignments <- list(Assignments = c(dulox_eligible_assignments$Assignments, dulox_ineligible_assignments$Assignments),
                               CovarInfo = dulox_eligible_assignments$CovarInfo + dulox_ineligible_assignments$CovarInfo,
                               EligibleCovarInfo = dulox_eligible_assignments$CovarInfo,
                               IneligibleCovarInfo = dulox_ineligible_assignments$CovarInfo)
  
  return(combined_assignments)
}



################################################################################
## Assignment Probability Functions
##
################################################################################

.MinimEqualProb <- function(elig.imbalance.scores, coin.bias){
  assignment_probs <- (1 - coin.bias) / (length(elig.imbalance.scores) - 1)
  
  return(assignment_probs)
}

covar_matrix <- sim_data_list_method_a[[1]][[1]] %>% 
  select(Site, starts_with("X_")) %>% 
  mutate_all(factor) %>%  mutate_all(as.numeric)  %>% #balanceRandomize function is expecting factor levels starting at 1 not 0 or -1
  as.matrix


contraindications_draw <- rmultinom(n = nrow(covar_matrix), size = 1, prob = c(.55, 0, .4, .025, .025))

# Extract the actual result and make it start at zero (indicating no contraindications)
contraindications <- apply(contraindications_draw, 2, which.max) - 1

covar_matrix_w_contras <- cbind(covar_matrix, contraindications)

set.seed(42)
test_alloc <- MinimRandomize(covariate.mat = covar_matrix, 
                             contraindications = rep(0, times = nrow(covar_matrix)), k.arms = 4, 
                             imbalance.func = function(x) diff(range(x)), 
                             use.overall.imbalance = FALSE,
                             covar.weights = 1,
                             coin.bias = 2/3, 
                             other.arm.p.func = .MinimEqualProb)

set.seed(42)
test_alloc_overall <- MinimRandomize(covariate.mat = covar_matrix, 
                                     contraindications = rep(0, times = nrow(covar_matrix)), k.arms = 4, 
                                     imbalance.func = function(x) diff(range(x)), 
                                     use.overall.imbalance = TRUE,
                                     covar.weights = 1,
                                     coin.bias = 2/3, 
                                     other.arm.p.func = .MinimEqualProb)

set.seed(42)
test_alloc_stratify <- StratifyThenMinimizeRand(covariate.mat = covar_matrix, 
                                     contraindications = contraindications, 
                                     k.arms = 4, 
                                     imbalance.func = function(x) diff(range(x)), 
                                     use.overall.imbalance = FALSE,
                                     covar.weights = 1,
                                     coin.bias = 2/3, 
                                     other.arm.p.func = .MinimEqualProb)

set.seed(42)
test_alloc_contra_min <- MinimRandomize(covariate.mat = cbind(contraindications, covar_matrix), 
                                     contraindications = contraindications, 
                                     k.arms = 4, 
                                     imbalance.func = function(x) diff(range(x)), 
                                     use.overall.imbalance = FALSE,
                                     covar.weights = 1,
                                     coin.bias = 2/3, 
                                     other.arm.p.func = .MinimEqualProb)

test_alloc_ignore_contra <- MinimRandomize(covariate.mat = covar_matrix, 
                                        contraindications = contraindications, 
                                        k.arms = 4, 
                                        imbalance.func = function(x) diff(range(x)), 
                                        use.overall.imbalance = FALSE,
                                        covar.weights = 1,
                                        coin.bias = 2/3, 
                                        other.arm.p.func = .MinimEqualProb)

set.seed(42)
josh_alloc <- balanceRandomize(Z = covar_matrix, N = 4, bias = 2/3)
josh_alloc1 <- balanceRandomize(Z = covar_matrix, N = 4, bias = 1)
josh_alloc2 <- balanceRandomize(Z = covar_matrix, N = 4, bias = .3)


#' 
#' @param n.subj number of subjects
#' @param k.arms number of arms
#' @param n.sites number of sites
#' @param bin.props vector of probabilities of successes for each binary covariate. 
#' The length of the vector determines the number of binary covariates
#' @param bin.iccs vector of intra-cluster correlation coefficients the same length as bin.props, or a single constant that will be used for all covariates
#' 
#' 
GenerateStudyData <- function(i = 1, n.subj = 600, k.arms = 4, n.sites = 10, bin.props = runif(5, min = .1, max = .9), bin.iccs = .1,
                              ...) {
  study_data <- tibble(ID = 1:n.subj,
                       Site = sample(1:n.sites, size = n.subj, replace = TRUE)) %>%
    MakeCorrelatedBinaryCovars(indf = ., bin.props = bin.props, bin.iccs = bin.iccs)
  
  covar_matrix <- study_data %>% 
    select(Site, starts_with("X_")) %>% 
    mutate_all(factor) %>%  mutate_all(as.numeric)  %>% #balanceRandomize function is expecting factor levels starting at 1 not 0 or -1
    as.matrix
  
  study_data$A1Josh <- balanceRandomize(Z = covar_matrix, 
                                        N = k.arms, 
                                        ...)$assignment
  
  study_data$A1John <- MinimRandomize(covariate.mat = covar_matrix, 
                                      contraindications = rep(0, times = nrow(covar_matrix)), k.arms = 4, 
                                      imbalance.func = function(x) diff(range(x)), 
                                      covar.weights = 1,
                                      coin.bias = 2/3, 
                                      other.arm.p.func = .MinimEqualProb)$assignment
  
  return(study_data)
}



####################################################################
############
####################################################################
dulox_eligible_assignments <- MinimRandomize(covariate.mat = duloxetine_eligible_subjects, 
                                             contraindications = contraindications[contraindications != 1],
                                             k.arms = 4,
                                             imbalance.func = function(x) diff(range(x)), 
                                             use.overall.imbalance = FALSE,
                                             covar.weights = 1,
                                             coin.bias = 2/3, 
                                             other.arm.p.func = .MinimEqualProb)

dulox_ineligible_assignments <- MinimRandomize(covariate.mat = duloxetine_ineligible_subjects, 
                                               contraindications = contraindications[contraindications == 1],
                                               k.arms = 4,
                                               imbalance.func = function(x) diff(range(x)), 
                                               use.overall.imbalance = FALSE,
                                               covar.weights = 1,
                                               coin.bias = 2/3, 
                                               other.arm.p.func = .MinimEqualProb)

if (FALSE){
  covariate.mat = duloxetine_ineligible_subjects
  covariate.mat = duloxetine_ineligible_subjects
  
  contraindications = contraindications[contraindications == 1]
  k.arms = 4
  imbalance.func = function(x) diff(range(x))
  use.overall.imbalance = FALSE
  covar.weights = 1
  coin.bias = 2/3
  other.arm.p.func = .MinimEqualProb
}


contra_asgn <- MinimRandomize(covariate.mat = covar_matrix_w_contras, 
                              contraindications = covar_matrix_w_contras[,1],
                              k.arms = 4,
                              imbalance.func = function(x) diff(range(x)), 
                              use.overall.imbalance = TRUE,
                              covar.weights = 1,
                              coin.bias = 2/3, 
                              other.arm.p.func = .MinimEqualProb)

if (FALSE) {
  covar_matrix_w_contras[,1] <- covar_matrix_w_contras[,1] + 1
  covariate.mat = covar_matrix_w_contras
  contraindications = covar_matrix_w_contras[,1]
  k.arms = 4
  imbalance.func = function(x) diff(range(x))
  use.overall.imbalance = FALSE
  covar.weights = 1
  coin.bias = 2/3
  other.arm.p.func = .MinimEqualProb
}
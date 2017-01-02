source('modules/Classes/Estimators.R')

runEstimators <- function(){
  my.computer <- 1
  expand.grid(
    dirichlet.alpha = c(.01, 1),
    num.theta = 2,
    num.actions = c(2, 6),
    num.states = 
      seq(100, 110, length.out = 1000) %>%
      as.integer
  ) %>% 
    mt(random.seed = 10^4 * my.computer + sample(n())) %>%
    d_('random.seed', function(s){
      set.seed(s$random.seed)
      
      Estimators$new(
        beta = .999,
        dirichlet.alpha = s$dirichlet.alpha,
        num.theta = s$num.theta,
        num.actions = s$num.actions,
        num.states = s$num.states, 
        num.obs = 300
      ) %>% {
        .$initializeEstimators() 
        # if(s$num.states < 3000) .$runNFXP_timed()
        .$runNPL_timed()
        # .$runSNFXP_timed()
        .$runSNPL_timed()
        list(
          experiment = s,
          times = .$times,
          theta = .$theta
        )
      } %>%
        se(s$random.seed, loc = experimentSave)
    }, .parallel = TRUE)
}

counterfactualAnalysis <- function(){
  my.trial <- 9
  
  expand.grid(
    dirichlet.alpha = c(.01, 1),
    num.theta = 2,
    num.actions = c(2, 6),
    num.states = 48000
  ) %>% 
    ag(num.states) %>%
    mt(random.seed = row_number()) %>% 
    d_('random.seed', function(s){
      set.seed(s$random.seed)
      
      Estimators$new(
        beta = .999,
        dirichlet.alpha = s$dirichlet.alpha,
        num.theta = s$num.theta,
        num.actions = s$num.actions,
        num.states = s$num.states, 
        num.obs = 300
      ) %>% {
        .$initializeEstimators() 
        # .$solve_DP_full_timed()
        .$solve_DP_shape_timed()
        list(
          experiment = s,
          times = .$times
        )
      } %>%
        se(paste(my.trial, s$random.seed, sep ='_'), loc = counterfactualSave)
    }, .parallel = FALSE)
}

memoryTests <- function(){
  df(
    dirichlet.alpha = .01,
    num.theta = 2,
    num.actions = 2,
    num.states = 30000
  ) %>% {
    Estimators$new(
      beta = .999,
      dirichlet.alpha = .$dirichlet.alpha,
      num.theta = .$num.theta,
      num.actions = .$num.actions,
      num.states = .$num.states, 
      num.obs = 300
    )
  } %>% {
    .$initializeEstimators() 
    .$runSNPL(num.VI.iters = 1)
    .$runNPL()
  } 
}

discountFactorMispecificationTests <- function(){
  expand.grid(
    dirichlet.alpha = .01,
    num.theta = 2,
    num.actions = 2,
    num.states = 1000,
    beta.real = c(0.99, 1.01),
    beta.used = c(0.99, 1.01),
    num.obs = seq(500, 50000, 250)
  ) %>%
    mt(
      random.seed = num.obs,
      trial = row_number()
    ) %>% 
    d_('trial', function(s){
      set.seed(s$random.seed)
      
      Estimators$new(
        beta = s$beta.real,
        dirichlet.alpha = s$dirichlet.alpha,
        num.theta = s$num.theta,
        num.actions = s$num.actions,
        num.states = s$num.states, 
        num.obs = s$num.obs
      ) %>% {
        .$initializeEstimators() 
        .$resetBeta(s$beta.used)
        .$runSNPL(num.VI.iters = 1)
        
        list(
          experiment = s,
          theta = .$theta
        )
      } %>%
        se(s$trial, loc = mispecificationSave)
    }, .parallel = TRUE)
}
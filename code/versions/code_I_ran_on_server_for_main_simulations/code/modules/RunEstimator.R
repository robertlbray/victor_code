source('modules/Classes/Estimators.R')

runEstimators <- function(){
  my.computer <- 4
  expand.grid(
    dirichlet.alpha = c(.01, 1),
    num.theta = 2,
    num.actions = c(2, 6),
    num.states = 
      seq(100, 10000, length.out = 200) %>%
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
        if(s$num.states < 2000) .$runNFXP_timed()
        .$runNPL_timed()
        .$runNNFXP_timed()
        .$runNNPL_timed()
        
        list(
          experiment = s,
          times = .$times,
          theta = .$theta
        )
      } %>%
        se(s$random.seed, loc = experimentSave)
    }, .parallel = TRUE)
}



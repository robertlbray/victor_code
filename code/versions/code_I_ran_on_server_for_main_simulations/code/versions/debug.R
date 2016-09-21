  
119%>%
  ll(function(l){
    print(l)
    set.seed(l)
    df(
      dirichlet.alpha = .01,
      num.theta = 2,
      num.actions = 2,
      num.states = 2839
    )  ->s
    Estimators$new(
      beta = .9,
      dirichlet.alpha = s$dirichlet.alpha,
      num.theta = s$num.theta,
      num.actions = s$num.actions,
      num.states = s$num.states, 
      num.obs = 300
    ) %T>% {
      .$initializeEstimators() 
      # .$runNPL_timed()
      .$runNNFXP_timed()
    } %>%
    se(l, loc = experimentSave)
  }, .parallel = TRUE) 
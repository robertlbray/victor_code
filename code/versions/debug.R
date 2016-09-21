dp <<- 
  DPEngine$new(
    beta = .99,
    num.states = 20000, 
    num.actions = 2, 
    num.theta = 2,
    dirichlet.alpha = 1
  ) %T>% {
    .$initializeDP()
  }

theta <- rnorm(2)



system.time(dp$calcDP_policy(theta = theta, num.value.updates = 1))

dp$resetSystem()

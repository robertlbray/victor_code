Transition <- 
  setRefClass(
    'Transition',
    fields = list(
      num.states = 'numeric',
      num.actions = 'numeric',
      dirichlet.alpha = 'numeric',
      tran.mat = 'list'
    ),
    methods = list(
      calc_transition = function(){
        #Reads: num.states, num.actions, dirichlet.alpha
        #Writes: tran.mat
        library('gtools')
        tran.mat <<- 
          seq(num.actions) %>%
          map(~ rdirichlet(n = num.states, alpha = rep(dirichlet.alpha, num.states)))
      }
    )
  )
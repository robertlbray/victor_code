PayoffFunction <- 
  setRefClass(
    'PayoffFunction',
    fields = list(
      payoff.stat = 'list', #list of matricies
      utility = 'list' #payoff.stat * diag(theta)
    ),
    methods = list(
      createPayoffStats = function(num.states, num.actions, num.theta){
        #Reads: num.states, num.actions, num.theta (type: integer)
        #Writes: payoff.stat
        payoff.stat <<-
          seq(num.actions) %>%
          map(
            ~ rnorm(num.states * num.theta) %>%
                matrix(num.states)
          )
      },
      
      updatePayoffs = function(theta){
        #Payoff vector \pi_q * \theta
        #Reads: c.stat, and theta (type: data.frame)
        #Writes: theta, cash.flows
        utility <<-
          payoff.stat %>%
          map(~ . %*% theta)
      }
    )
  )
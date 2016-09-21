source('modules/Classes/Transition.R')
source('modules/Classes/PayoffFunction.R')

DPEngine <- 
  setRefClass(
    'DPEngine',
    fields = list(
      iter.stop.threshold = 'numeric',
      beta = 'numeric',
      num.states = 'numeric',
      num.actions = 'numeric',
      num.theta = 'numeric',
      dirichlet.alpha = 'numeric',
      tran = 'Transition',
      payoff = 'PayoffFunction',
      CCP = 'list',
      CCP.saved = 'list', #Last "saved" CCPs
      value.fn = 'numeric',
      policy.iteration.mat = 'matrix'
    ),
    methods = list(
      initializeDP = function(){
        iter.stop.threshold <<- 10^-10
        
        tran <<-
          Transition$new(
            num.states = num.states,
            num.actions = num.actions,
            dirichlet.alpha = dirichlet.alpha) %T>% {
              .$calc_transition()
            }
        
        payoff <<-
          PayoffFunction$new() %T>% {
            .$createPayoffStats(num.states, num.actions, num.theta)
          }
        
        resetSystem()
      },
      
      updateValues_VI = function(){
        #Maps CCPs and values to values
        value.fn <<-
          list(
            CCP,
            payoff$utility,
            tran$tran.mat
          ) %>%
          pmap(function(ccp, payoff, mat){
            ccp * (payoff - log(ccp) + beta * mat %*% value.fn)
          }) %>%
          Reduce('+', .) %>% 
          c
      },
      
      updateValues_PI = function(){
        #Maps CCPs and payoffs to values
        value.fn <<-
          map2(
            CCP,
            payoff$utility,
            ~ .x * (.y - log(.x))
          ) %>%
          Reduce('+', .) %>% {
            policy.iteration.mat %*% .
          } %>%
          c
      },
      
      updateCCPs = function(){
        #Maps values to CCPs
        CCP <<- 
          map2(
            payoff$utility,
            tran$tran.mat,
            ~ .x + beta * .y %*% value.fn
          ) %>% {
            max.CSV <- Reduce(pmax, .) 
            map(., ~ . - max.CSV)
          } %>%
          map(~ exp(.)) %>% {
            denominator <-
              Reduce('+', .)
            map(., ~ c(./denominator))
          } %>%
          map(~ pmax(10^-100, .))
      },
      
      operator_VI = function(){
        updateCCPs()
        updateValues_VI()
      },
      
      operator_PI = function(){
        updateValues_PI()
        updateCCPs()
      },
      
      calcPolicyIterationMat = function(){
        policy.iteration.mat <<-
          map2(
            CCP,
            tran$tran.mat,
            ~ .x * .y
          ) %>%
          Reduce('+', .) %>% {
            solve(diag(ncol(.)) - beta * .)
          } 
      },
      
      estimateInitialCCPs = function(panel) {
        CCP <<-
          seq(payoff$payoff.stat) %>%
          ld(function(l){
            payoff$payoff.stat[[l]] %>%
              ad %>%
              mt(
                action = l,
                state = row_number()
              ) %>%
              ml(c('action', 'state'), variable.name = 'stat')
          }) %>% 
          dc(state ~ stat + action) %>% {
            new.data <- .
            my.formula <- 
              sl(., -state) %>% 
              names %>%
              paste(collapse = '+') %>% {
                paste('action', ., sep = '~')
              } %>% 
              as.formula
            
            ij(., panel) %>%
              nnet::multinom(
                my.formula,
                data = .,
                weights = n
              ) %>%
              predict(., type = 'probs', newdata = new.data) 
          } %>% {
            if(class(.) == 'matrix') {
              alply(., 2, function(l) l)
            } else {
              list(1 - ., .)
            }
          }
          
      },
      
      calcDP_policy = function(theta){
        payoff$updatePayoffs(theta)
        delta <- 1
        print(h(value.fn))
        while(delta > iter.stop.threshold){
          print('value iteration')
          CCP.prior <- CCP
          operator_VI()
          delta <-
            map2(
              CCP.prior,
              CCP,
              ~ max(abs(.x - .y))
            ) %>%
            rd(max)
        }
      },
      
      calcDP_value = function(theta){
        #Rust's method for finding value function: value interations then Newton steps
        payoff$updatePayoffs(theta)
        delta.prior <- 0
        delta <- 1
        newton <- FALSE
        newton.iterations <- 0 
        while(delta > iter.stop.threshold){
          value.fn.prior <- value.fn
          if(newton) {
            print('Newton policy iteration')
            newton.iterations <- newton.iterations + 1
            calcPolicyIterationMat()
            operator_PI()
          } else {
            print('Newton value iteration')
            operator_VI()
          }
          print(delta)
          delta <- max(abs(value.fn.prior - value.fn))
          if((delta < 10^-3) & (abs(delta/delta.prior - beta) < 10^-4)) newton <- TRUE
          if(newton.iterations == 3) delta <- 0 #Avoid infinite loops due to lack of precision
          delta.prior <- delta
        }
      },
      
      resetSystem = function(){
        value.fn <<- rep(0, num.states)
        CCP <<- #set all CCPs  1/num.actions
          seq(num.actions) %>%
          map(~ rep(1/num.actions, num.states))
        saveCCPs()
      },
      
      saveCCPs = function(){
        CCP.saved <<- CCP
      },
      
      loadCCPs = function(){
        CCP <<- CCP.saved
      }
    )
  )

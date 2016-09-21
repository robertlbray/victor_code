source('modules/Classes/Transition.R')
source('modules/Classes/PayoffFunction.R')

PolicyIteration <- 
  setRefClass(
    'PolicyIteration',
    fields = list(
      beta = 'numeric',
      iter.stop.threshold = 'numeric',
      num.states = 'integer',
      num.actions = 'integer',
      num.theta = 'integer',
      dirichlet.alpha = 'numeric',
      tran = 'Transition',
      payoff = 'PayoffFunction',
      CCP = 'list', #Working CCP value
      CCP.1 = 'list', #Tentative next-iteration CCP
      CCP.0 = 'list', #initial CCP estimate
      CSV = 'list',
      value.fn = 'numeric',
      policy.iteration.mat.NPL = 'matrix',
      policy.iteration.mat.SANPL = 'matrix',
      policy.iteration.mat.ANNPL = 'matrix'
    ),
    methods = list(
      initializePolIter = function(){
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
        
        resetCCPs()
      },
      
      calcPolicyIterationMat_NPL = function(){
        policy.iteration.mat.NPL <<-
          map2(
            CCP,
            tran$tran.mat,
            ~ .x * .y
          ) %>%
          Reduce('+', .) %>% {
            solve(diag(ncol(.)) - beta * .)
          } 
      },
      
      calcPolicyIterationMat_SANPL = function(){
        policy.iteration.mat.SANPL <<-
          tran$tran.mat[[1]] %>% {
            solve(diag(ncol(.)) - beta * .)
          } 
      },
      
      calcPolicyIterationMat_ANNPL = function(){
        policy.iteration.mat.ANNPL <<-
          map2(
            CCP.0,
            tran$tran.mat,
            ~ .x * .y
          ) %>%
          Reduce('+', .) %>% {
            solve(diag(ncol(.)) - beta * .)
          } 
      },
      
      calcValueFn_NPL = function(){
        value.fn <<-
          map2(
            CCP,
            payoff$utility,
            ~ .x * (.y - log(.x))
          ) %>%
          Reduce('+', .) %>% {
            policy.iteration.mat.NPL %*% .
          } %>%
          c
      },
      
      calcValueFn_SANPL = function(){
        value.fn <<- { 
          policy.iteration.mat.SANPL %*% (payoff$utility[[1]] - log(CCP[[1]])) 
        } %>%
          c
      },
      
      calcValueFn_ANNPL = function(){
        value.fn <<- { 
          policy.iteration.mat.ANNPL %*% (payoff$utility[[1]] - log(CCP[[1]])) 
        } %>%
          c
      },
      
      calcCSV = function(){
        CSV <<-
          map2(
            payoff$utility,
            tran$tran.mat,
            ~ .x + beta * .y %*% value.fn
          )
      },
      
      calcCCP = function(){
        CCP.1 <<- 
          CSV %>%
          map(~ . - Reduce(pmax, CSV)) %>%
          map(~ exp(.)) %>% {
            denominator <-
              Reduce('+', .)
            
            map(., ~ c(./denominator))
          } %>%
          map(~ pmax(10^-100, .))
      },
      
      policyIterationOperator_NPL = function(){
        calcValueFn_NPL()
        calcCSV()
        calcCCP()
      },
      
      policyIterationOperator_SANPL = function(){
        calcValueFn_SANPL()
        calcCSV()
        calcCCP()
      },
      
      policyIterationOperator_ANNPL = function(){
        calcValueFn_ANNPL()
        calcCSV()
        calcCCP()
      },
      
      iterateToConvergence_NPL = function(theta){
        pol.iter$payoff$updatePayoffs(theta)
        delta <- 1
        while(delta > iter.stop.threshold){
          print(delta)
          calcPolicyIterationMat_NPL()
          policyIterationOperator_NPL()
          delta <-
            map2(
              CCP,
              CCP.1,
              ~ max(.x - .y)
            ) %>%
            rd(max)
          ccpUpdate()
        }
      },
      
      iterateToConvergence_SANPL = function(theta){
        # browser()
        pol.iter$payoff$updatePayoffs(theta)
        delta <- 1
        calcPolicyIterationMat_SANPL()
        while(delta > iter.stop.threshold){
          print(delta)
          policyIterationOperator_SANPL()
          delta <-
            map2(
              CCP,
              CCP.1,
              ~ max(.x - .y)
            ) %>%
            rd(max)
          ccpUpdate()
        }
      },
      
      estimateInitialCCPs = function(panel) {
        CCP.0 <<-
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
      
      resetCCPs = function(){
        if(length(CCP.0) == 0) {
          CCP <<- #set all CCPs  1/num.actions
            seq(num.actions) %>%
            map(~ rep(1/num.actions, num.states))
        } else {
          CCP <<- CCP.0
        }
      },
      
      ccpUpdate = function(){
        CCP <<- CCP.1
      }
    )
  )

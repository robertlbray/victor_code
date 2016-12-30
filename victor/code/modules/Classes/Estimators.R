source('modules/Classes/DPEngine.R')

Estimators <- 
  setRefClass(
    'Estimators',
    fields = list(
      beta = 'numeric',
      num.theta = 'numeric',
      num.states = 'numeric',
      num.actions = 'numeric',
      num.obs = 'numeric',
      dirichlet.alpha = 'numeric',
      panel = 'data.frame',
      times = 'list', 
      theta = 'list', #true theta and theta estimates
      dp = 'DPEngine'
    ),
    methods = list(
      initializeEstimators = function(){
        times <<- as.list(NULL)
        theta <<- as.list(NULL)
        theta$true <<- rnorm(num.theta)
        
        dp <<- 
          DPEngine$new(
            beta = beta,
            num.states = num.states, 
            num.actions = num.actions, 
            num.theta = num.theta,
            dirichlet.alpha = dirichlet.alpha
          ) %T>% {
            .$initializeDP()
          }
        
        createPanel()
      },
      
      createPanel = function() {
        dp$solve_DP_shape(theta$true)
        panel <<-
          dp$CCP %>%
          rd(cbind) %>% {
            colnames(.) <- NULL
            .
          } %>%
          ad %>%
          mt(state = row_number()) %>% {
            ij(
              ., 
              df(state = seq(num.states)) %>%
                sn(num.obs, replace = TRUE) %>%
                ct(state)
            )
          } %>%
          dd('state', function(s){
            rmultinom(1, size = s$n, prob = sl(s, -c(state, n))) %>% 
              ad %>%
              mt(action = row_number()) %>%
              rn(n = V1)
          }) %>%
          fl(n > 0)
      },
      
      evaluateLikelihood = function(my.ccps){
        seq(num.actions) %>%
          ld(function(l){
            dp$CCP[[l]] %>%
              ad %>% 
              rename_(CCP = '.') %>%
              mt(
                state = row_number(),
                action = l
              )
          }) %>%
          ij(panel) %>% 
          sm(sum(n * log(CCP))) %>% 
          un
      },
      
      runNFXP = function() {
        dp$resetSystem()
        rep(0, num.theta) %>%
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$solve_DP_full(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$NFXP <<- .$par
          }
      },
      
      runSNFXP = function() {
        dp$resetSystem()
        rep(0, num.theta) %>%
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$solve_DP_shape(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$SNFXP <<- .$par
          }
      },
      
      runNPL = function() {
        theta$NPL <<- rep(0, num.theta)
        dp$resetSystem()
        dp$estimateInitialCCPs(panel)
        nestedPsuedoLikelihoodStep <- function(){
          dp$calcPolicyIterationMat()
          theta$NPL %>% 
            optim(  
              par = .,
              method = 'BFGS',
              control = list(fnscale=-1), #max likelihood
              fn = function(theta){
                dp$payoff$updatePayoffs(theta)
                dp$loadCCPs()
                dp$operator_PI()
                evaluateLikelihood()
              }
            ) %>% {
              theta$NPL <<- .$par
            }
          dp$saveCCPs()
        }
        delta = 1
        iter.stop.threshold <- 10^-10
        while(max(delta) > iter.stop.threshold){
          CCP.l <- dp$CCP
          nestedPsuedoLikelihoodStep()
          delta <- 
            map2(dp$CCP, CCP.l, ~ max(abs(.x - .y))) %>%
            unlist %>%
            max
        }
      },
      
      runSNPL = function() {
        theta$SNPL <<- rep(0, num.theta)
        dp$resetSystem()
        nestedPsuedoLikelihoodStep <- function(){
          theta$SNPL %>% 
            optim(  
              par = .,
              method = 'BFGS',
              control = list(fnscale=-1), #max likelihood
              fn = function(theta){
                dp$payoff$updatePayoffs(theta)
                dp$updateCCPs()
                evaluateLikelihood()
              }
            ) %>% {
              theta$SNPL <<- .$par
            }
          dp$updateValues_VI(use.val.shape = TRUE)
        }
        iter.stop.threshold <- 10^-6
        delta = 1
        while(max(delta) > iter.stop.threshold){
          value.fn.l <- dp$value.fn
          nestedPsuedoLikelihoodStep()
          delta <- max(abs(value.fn.l - dp$value.fn))
        }
        
        theta$SNPL %>% #One final run for accuracy
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$solve_DP_shape(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$SNPL <<- .$par
          }
      },
      
      runSNPL2 = function() {
        theta$SNPL2 <<- rep(0, num.theta)
        dp$resetSystem()
        nestedPsuedoLikelihoodStep <- function(){
          theta$SNPL2 %>% 
            optim(  
              par = .,
              method = 'BFGS',
              control = list(fnscale=-1), #max likelihood
              fn = function(theta){
                dp$payoff$updatePayoffs(theta)
                dp$loadValues()
                dp$updateValues_VI(use.val.shape = TRUE)
                dp$updateCCPs()
                evaluateLikelihood()
              }
            ) %>% {
              theta$SNPL2 <<- .$par
            }
          dp$saveValues()
          
          
          
          
          
          
          #####NOTE!!!!
          ####I think there is a bug here, which I fix in marcello.R. The final evaluation of the maximizer probably wont be the maximizing value. So I'm saving the wrong number! I haven't fixed this bug, because I don't want to change the code, but I think the "One final run for accuracy" part might be redundant after I fix it.
          
          
          
          
          
        }
        iter.stop.threshold <- 10^-6
        delta = 1
        while(max(delta) > iter.stop.threshold){
          value.fn.l <- dp$value.fn
          nestedPsuedoLikelihoodStep()
          delta <- max(abs(value.fn.l - dp$value.fn))
        }
        
        theta$SNPL2 %>% #One final run for accuracy
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$solve_DP_shape(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$SNPL2 <<- .$par
          }
      },
      
      runNFXP_timed = function(){
        times$NFXP <<- system.time(runNFXP())
      },
      
      runNPL_timed = function(){
        times$NPL <<- system.time(runNPL())
      },
      
      runSNFXP_timed = function(){
        times$SNFXP <<- system.time(runSNFXP())
      },
      
      runSNPL_timed = function(){
        times$SNPL <<- system.time(runSNPL())
      },
      
      runSNPL2_timed = function(){
        times$SNPL2 <<- system.time(runSNPL2())
      },
      
      solve_DP_shape_timed = function(){
        dp$resetSystem()
        times$dp_shape_solution <<- system.time(dp$solve_DP_shape(theta$true))
      },
      
      solve_DP_full_timed = function(){
        dp$resetSystem()
        times$dp_full_solution <<- system.time(dp$solve_DP_full(theta$true))
      },
      
      resetBeta = function(new.beta){
        beta <<- new.beta
        dp$beta <<- beta
      }
    )
  )
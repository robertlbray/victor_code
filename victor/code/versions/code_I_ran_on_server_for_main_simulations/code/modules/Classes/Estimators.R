source('modules/Classes/DPEngine.R')

Estimators <- 
  setRefClass(
    'Estimators',
    fields = list(
      iter.stop.threshold = 'numeric',
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
        iter.stop.threshold <<- 10^-9
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
        dp$calcDP_policy(theta$true)
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
              dp$calcDP_value(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$NFXP <<- .$par
          }
      },
      
      runNNFXP = function() {
        dp$resetSystem()
        rep(0, num.theta) %>%
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$calcDP_policy(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$NNFXP <<- .$par
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
        while(max(delta) > iter.stop.threshold){
          theta.l <- theta$NPL
          nestedPsuedoLikelihoodStep()
          delta <- max(abs(theta.l - theta$NPL))
        }
      },
      
      runNNPL = function() {
        theta$NNPL <<- rep(0, num.theta)
        dp$resetSystem()
        nestedPsuedoLikelihoodStep <- function(){
          theta$NNPL %>% 
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
              theta$NNPL <<- .$par
            }
          dp$updateValues_VI()
        }
        delta = 1
        while(max(delta) > iter.stop.threshold){
          theta.l <- theta$NNPL
          nestedPsuedoLikelihoodStep()
          delta <- max(abs(theta.l - theta$NNPL))
        }
        
        theta$NNPL %>% #One final run for accuracy
          optim(  
            par = .,
            method = 'BFGS',
            control = list(fnscale=-1), #max likelihood
            fn = function(theta){
              dp$calcDP_policy(theta)
              evaluateLikelihood()
            }
          ) %>% {
            theta$NNPL <<- .$par
          }
      },
      
      runNFXP_timed = function(){
        times$NFXP <<- system.time(runNFXP())
      },
      
      runNPL_timed = function(){
        times$NPL <<- system.time(runNPL())
      },
      
      runNNFXP_timed = function(){
        times$NNFXP <<- system.time(runNNFXP())
      },
      
      runNNPL_timed = function(){
        times$NNPL <<- system.time(runNNPL())
      }
    )
  )
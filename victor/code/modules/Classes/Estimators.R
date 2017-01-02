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
        times <<- as.list(NULL)
        theta <<- as.list(NULL)
        theta$true <<- rnorm(num.theta)
        
        iter.stop.threshold <<- 10^-9
        
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
          ij(panel, by = c('state', 'action')) %>% 
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
        
        ####NOTE: The monte carlo simulation ran from a previous version. That version is in versions/Econmetrica_submission. It didn't include the final NPL_step(theta$NPL) call before the dp$saveCCPs() line. 
        ####Also, I canged the definition of SNPL since then. But SNPL doesn't ever converge to something nice, so I'm dropping it.
        
        NPL_step <- function(theta) {
          dp$payoff$updatePayoffs(theta)
          dp$loadCCPs()
          dp$operator_PI()
        }
        
        theta$NPL <<- rep(0, num.theta)
        dp$resetSystem()
        dp$estimateInitialCCPs(panel)
        delta = 1
        while(delta > iter.stop.threshold){
          value.fn.l <- dp$value.fn
          
          dp$calcPolicyIterationMat()
          theta$NPL %>% 
            optim(  
              par = .,
              method = 'BFGS',
              control = list(fnscale=-1), #max likelihood
              fn = function(theta){
                NPL_step(theta)
                evaluateLikelihood()
              }
            ) %>% {
              theta$NPL <<- .$par
            }
          
          NPL_step(theta$NPL)
          dp$saveCCPs()
          delta <- max(abs(value.fn.l - dp$value.fn))
        }
      },
      
      runSNPL = function() {
        SNPL_step <- function(theta) {
          dp$payoff$updatePayoffs(theta)
          dp$loadValues()
          dp$loadCCPs()
          updateValues_VI(use.val.shape)
          updateCCPs()
        }
        
        theta$SNPL <<- rep(0, num.theta)
        dp$resetSystem()
        delta = 1
        while(delta > iter.stop.threshold){
          value.fn.l <- dp$value.fn
          theta$SNPL %>% 
            optim(  
              par = .,
              method = 'BFGS',
              control = list(fnscale=-1), #max likelihood
              fn = function(theta){
                SNPL_step(theta)
                evaluateLikelihood()
              }
            ) %>% {
              theta$SNPL <<- .$par
            }
          SNPL_step(theta$SNPL)
          dp$saveValues()
          dp$saveCCPs()
          delta <- max(abs(value.fn.l - dp$value.fn))
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
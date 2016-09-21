consolidateResults <- function() {
  dir(experimentSave, recursive = TRUE) %>%
    map(~ le(., loc = experimentSave)) %>%
    ld(function(l){
      rbind(
        l$theta %>%
          ld(function(m) max(abs(m - l$theta$true))) %>%
          rn(
            estimator = .id,
            error = V1
          ) %>%
          ml('estimator') %>%
          fl(estimator != 'true'),
        l$times %>% 
          ld(function(m) m[1]) %>%
          rn(
            estimator = .id,
            time = user.self
          ) %>%
          ml('estimator')
      ) %>%
        cbind(l$experiment)
    }) %>%
    se('experiments.rds')
}

tabErrors <- function() {
  le('experiments.rds') %>% 
    fl(variable == 'error') %>% 
    dc(
      random.seed ~ estimator, 
      value.var = 'value'
    ) %>%
    fl(!is.na(NNFXP + NNPL + NPL)) %>% 
    mt(large.prob = is.na(NFXP)) %>% 
    ml(c('random.seed', 'large.prob')) %>%
    na.omit %>%
    dd(
      c('large.prob', 'variable'), 
      function(s){
        quantile(s$value, probs = c(.1, .25, .5, .75, .9)) %>%
          ad %>% 
          add_rownames
      }
    ) %>% 
    mt(
      rowname = str_replace(rowname, '%', '\\\\%'),
      large.prob = setLevels(
        large.prob, 
        c('FALSE', 'TRUE'),
        c('Fewer Than 3,000 States', ' More Than 3,000 States')
      ),
      variable = setLevels(
        variable,
        c('NFXP', 'NNFXP', 'NPL', 'NNPL'),
        c('NFXP', 'SNFXP', 'NPL', 'SNPL')
      )
    ) %>% 
    make.table( 
      d = .,
      formula = variable ~ large.prob + rowname,
      value = '.',
      table.type = 'Summary', 
      out.file = paste0(data.out, 'tables/estimateErrors.tex'), 
      hline = TRUE,
      round.digits = 5,
      nsmall = 5,
      align.cols = 'lccccc|ccccc'
    )
}

plotMispecifiedBeta <- function() {
  dir(mispecificationSave) %>%
    map(~ le(., loc = mispecificationSave)) %>%
    map_df(~ cbind(
      .$experiment, 
      error = max(abs(.$theta$true - .$theta$SNPL)))
    ) %>%
    mt(beta.real = setLevels(
      beta.real,
      c(.99, 1.01),
      paste('Actual Discount Factor:', c(.99, 1.01))
    ),
    error = wins(error, c(0, .95))
    ) %>% {
      ggplot(
        data = ., 
        aes(x = num.obs, y = error, colour = as.factor(beta.used))
      ) + 
        geom_point(size = .75) +
        geom_smooth(size = .75, se = FALSE) +
        facet_wrap(~ beta.real) +
        labs(
          x = 'Number of Observations',
          y = 'Estimation Error'
        ) +
        scale_colour_grey(start = 0, end = .8, name = 'Calibrated\nDiscount\nFactor') +
        theme(
          strip.text.x=element_text(size=13.5),
          strip.text.y=element_text(size=13.5),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=15),
          legend.text=element_text(size=13.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm")
        )
    } %>% 
    ggsave(plot = ., filename = p0(data.out, "plots/misspecification.pdf"), width=9, height=6)
}

plotMemoryRequirements <- function() {
  paste0(data.loc, 'intermediate/memoryExperiments/memory.csv') %>%
    read_csv %>%
    mt(
      States = 1000 * States,
      Estimator = setLevels(
        Estimator, 
        c('NPL', 'NNPL'),
        c('NPL', 'SNPL')
      )
    ) %>% {
      ggplot(
        data = ., 
        aes(x = States, y = GiB, colour = Estimator)
      ) + 
        geom_line(size = .75) +
        labs(
          x = 'Number of States',
          y = 'GiB of Memory'
        ) +
        scale_colour_grey(start = .7, end = 0) +
        theme(
          strip.text.x=element_text(size=13.5),
          strip.text.y=element_text(size=13.5),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=13.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %>% 
    ggsave(plot = ., filename = p0(data.out, "plots/memory.pdf"), width=5, height=4)
}

plotCounterfactuals <- function() {
  dir(counterfactualSave, recursive = TRUE) %>%
    map_df(
      ~ le(., loc = counterfactualSave) %>% {
        merge(
          .$experiment,
          df(
            full.time = ifelse(is.null(.$times$dp_full_solution[[1]]), NA, .$times$dp_full_solution[[1]]), 
            shape.time = .$times$dp_shape_solution[[1]]
          )
        )
      }
    ) %>%
    ml(measure.vars = c('full.time', 'shape.time')) %>%
    na.omit %>%
    mt(
      value = value / 60,
      num.actions = setLevels(num.actions, c(2, 6), c('Two Actions', 'Six Actions')),
      dirichlet.alpha = setLevels(dirichlet.alpha, c(.01, 1), c('Sparse Matrix', 'Dense Matrix')),
      variable = setLevels(
        variable,
        c('full.time', 'shape.time'),
        c('Without Strong Convergence', 'With Strong Convergence')
      )
    ) %>% {
      ggplot(
        data = ., 
        aes(x = num.states, y = value, colour = variable)
      ) + 
        geom_point(size = .75) +
        geom_smooth(size = .75, se = FALSE) +
        geom_hline(yintercept = 0, colour = 'black', size = .25) +
        facet_wrap(~ num.actions + dirichlet.alpha, nrow = 1 ) +
        labs(
          x = 'Number of States',
          y = 'Number of Minutes'
        ) +
        scale_colour_grey(start = .7, end = 0) +
        scale_x_continuous(breaks = c(0, 25000, 50000)) +
        theme(
          strip.text.x=element_text(size=13.5),
          strip.text.y=element_text(size=13.5),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=13.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.25, "cm"),
          legend.position = 'top'
        )
    } %T>% {
      ggsave(plot = ., filename = p0(data.out, "plots/counterfactual.pdf"), width=9, height=12)
    } %>% {
      ggsave(plot = ., filename = p0(data.out, "pres/counterfactual.png"), width=10, height=7)
    }
}

plotTimes <- function() {
  le('experiments.rds') %>% 
    fl(variable == 'time') %>% 
    mt(
      newtonless = str_sub(estimator, 1, 2) == 'NN',
      rust = str_detect(estimator, 'X'),
      value = value / 60,
      dirichlet.alpha = setLevels(
        dirichlet.alpha, 
        c(0.01, 1),
        c('Sparse Matrix', 'Dense Matrix')
      ),
      num.actions = setLevels(
        num.actions, 
        c(2, 6),
        c('Two Actions', 'Six Actions')
      ),
      rust = setLevels(
        rust,
        c(TRUE, FALSE),
        c('NFXP', 'NPL')
      ),
      newtonless = setLevels(
        newtonless,
        c(FALSE, TRUE),
        c('Without Strong Convergence', 'With Strong Convergence')
      )
    ) %>% 
    gb(rust, newtonless) %>%
    mt(value = wins(value, c(0, .99))) %>% {
      ggplot(
        data = ., 
        aes(x = num.states, y = value, colour = newtonless)
      ) + 
        geom_point(size = .75) +
        geom_smooth(size = .75, se = FALSE) +
        geom_hline(yintercept = 0, colour = 'black', size = .25) +
        facet_grid(rust ~ num.actions + dirichlet.alpha, scale = 'free_y') +
        labs(
          x = 'Number of States',
          y = 'Number of Minutes'
        ) +
        scale_colour_grey(start = .7, end = 0) +
        theme(
          strip.text.x=element_text(size=13.5),
          strip.text.y=element_text(size=13.5),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=11),
          axis.text.y=element_text(size=11),
          legend.title=element_text(size=0),
          legend.text=element_text(size=13.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %T>% {
      ggsave(plot = ., filename = p0(data.out, "plots/times.pdf"), width=9.5, height=10)
    }  %>%
    ggsave(plot = ., filename = p0(data.out, "pres/times.png"), width=10, height=7)
}

plotTimeRatio <- function(variables) {
  le('experiments.rds') %>%
    fl(variable == 'time') %>% 
    mt(
      newtonless = str_sub(estimator, 1, 2) == 'NN',
      rust = str_detect(estimator, 'X'),
      num.states = pmin(max(num.states - 1), num.states),
      num.states = num.states - num.states %% 200,
      dirichlet.alpha = setLevels(
        dirichlet.alpha, 
        c(0.01, 1),
        c('Sparse Transition Matrix', 'Dense Transition Matrix')
      ),
      num.actions = setLevels(
        num.actions, 
        c(2, 6),
        c('Two Actions', 'Six Actions')
      ),
      rust = setLevels(
        rust,
        c(TRUE, FALSE),
        c('NFXP', 'NPL')
      )
    ) %>%
    gb(rust, newtonless, dirichlet.alpha, num.actions, num.states) %>%
    sm(value = mean(value)) %>% 
    dc(... ~ newtonless, value.var = 'value') %>% 
    mt(ratio = `TRUE`/`FALSE`) %>%  
    na.omit  %>% {
      ggplot(
        data = ., 
        aes(
          x = num.states, 
          y = ratio
        )
      ) + 
        geom_point(size = 1) +
        geom_line(size = .75) +
        labs(
          x = 'Number of States',
          y = 'Ratio of Mean Times'
        ) +
        facet_grid(dirichlet.alpha + num.actions ~ rust, scales = 'free_x') +
        scale_colour_grey(start = 0, end = .7) +
        scale_y_log10() +
        theme(
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=14),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_line(colour = "grey95"), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %>%
    ggsave(plot = ., filename = p0(data.out, "plots/ratio.pdf"),  dpi=300, width=8.5, height=10)
}

plotRustExample <- function() {
  num.actions <- 2
  num.theta <- 1
  dirichlet.alpha <- 1
  num.states <- 12
  lambda <- 1
  iter.stop.threshold <- 10^-10
  theta <- c(1, 10)
  
  dp <- 
    DPEngine$new(
      num.states = num.states, 
      num.actions = num.actions, 
      num.theta = num.theta,
      dirichlet.alpha = dirichlet.alpha
    ) %T>% {
      .$initializeDP()
    }
  
  dp$tran$tran.mat <- 
    list(
      cbind(
        rep(0, num.states - 1),
        diag(num.states - 1)
      ) %>%
        rbind(c(rep(0, num.states - 1), 1)),
      cbind(
        rep(1, num.states), 
        matrix(rep(0, num.states^2 - num.states), num.states)
      )
    )
  
  dp$payoff$payoff.stat <-
    list(
      -matrix(theta[1] * seq(num.states)),
      -matrix(theta[2] * rep(1, num.states))
    )
  
  dp$payoff$updatePayoffs(1)
  
  c(.9, .99, .999, .9999) %>%
    ld(function(my.beta){
      dp$beta <- my.beta
      dp$resetSystem()
      delta.value <- c()
      delta.policy <- c()
      delta.value.shape <- c()
      delta <- 1
      while(delta > iter.stop.threshold){
        print(delta)
        value.fn.prior <- dp$value.fn
        CCP.prior <- dp$CCP[[1]]
        dp$operator_VI()
        delta <- max(abs(value.fn.prior - dp$value.fn))
        delta.value <- c(delta.value, delta)
        delta.value.shape <- c(delta.value.shape, max(abs(value.fn.prior - value.fn.prior[1] - dp$value.fn + dp$value.fn[1])))
        delta.policy <- c(delta.policy, max(abs(CCP.prior - dp$CCP[[1]])))
      } 
      
      (dp$CCP[[1]]  * dp$tran$tran.mat[[1]] + dp$CCP[[2]]  * dp$tran$tran.mat[[2]]) %>%
        eigen %>% {
          .$values
        } %>% 
        Mod %>% {
          .[2]
        } %>% {
          print(c('eigenvalue:', .))
        }
      
      cbind(
        delta.value,
        delta.value.shape,
        delta.policy
      ) %>% 
        ad %>% 
        slice(-1) %>%
        mt(
          Iteration = row_number(),
          my.beta = my.beta
        ) %>% 
        ml(c('Iteration', 'my.beta')) %>% 
        mt(
          variable =
            setLevels(
              variable, 
              c('delta.policy', 'delta.value.shape', 'delta.value'),
              c('Policy Function', 'Differenced Value Function', 'Value Function')
            ),
          my.beta = paste(expression('beta'), '==', my.beta)
        ) %>%
        fl(value > 10^-10)
    }) %>% {
      ggplot(
        data = ., 
        aes(
          x = Iteration, 
          y = value,
          colour = variable
        )
      ) + 
        geom_line(size = 1) +
        facet_wrap(
          ~ my.beta,
          labeller = label_parsed,
          nrow = 2
        ) +
        labs(
          x = 'Value Iteration Algorithm Step Number',
          y = 'Function Change'
        ) +
        scale_colour_grey(start = .7, end = 0) +
        scale_y_log10() +
        scale_x_log10() +
        theme(
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=14),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_line(colour = "grey95"), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %T>%
    ggsave(plot = ., filename = p0(data.out, "plots/rustModel.pdf"),  dpi=300, width=9, height=6) %>% 
    ggsave(plot = ., filename = p0(data.out, "pres/rustModel.png"),  dpi=300, width=9, height=6)
}

plotRustExample2 <- function() {
  library('expm')
  
  num.actions <- 2
  num.theta <- 1
  dirichlet.alpha <- 1
  num.states <- 150
  my.beta <- .999
  lambda <- 1
  iter.stop.threshold <- 10^-12
  theta <- c(1, 10)
  
  dp <- 
    DPEngine$new(
      num.states = num.states, 
      num.actions = num.actions, 
      num.theta = num.theta,
      dirichlet.alpha = dirichlet.alpha,
      beta = my.beta
    ) %T>% {
      .$initializeDP()
    }
  
  dp$tran$tran.mat <- 
    list(
      cbind(
        rep(0, num.states - 1),
        diag(num.states - 1)
      ) %>%
        rbind(c(rep(0, num.states - 1), 1)),
      cbind(
        rep(1, num.states), 
        matrix(rep(0, num.states^2 - num.states), num.states)
      )
    )
  
  dp$payoff$payoff.stat <-
    list(
      -matrix(theta[1] * seq(num.states)),
      -matrix(theta[2] * rep(1, num.states))
    )
  
  dp$payoff$updatePayoffs(1)
  
  #Solve policy (I can't use solve_DP_shape method, because it updates the payments)
  delta <- 1
  while(delta > iter.stop.threshold){
    value.fn.prior <- dp$value.fn
    dp$operator_VI(use.val.shape = TRUE)
    delta <- max(abs(value.fn.prior - dp$value.fn))
  }
  
  tran.mat <-
    map2(
      dp$CCP,
      dp$tran$tran.mat,
      ~ .x * .y
    ) %>%
    Reduce('+', .)
  
  c(0, 2^(0:6)) %>%
    map_df(~ 
             df(
               state = seq(num.states),
               time = .,
               `New Engine` = (tran.mat %^% .)[1,],
               `Used Engine` = (tran.mat %^% .)[3, ]
             )
    ) %>%
    ml(c('state', 'time')) %>%
    fl(state < 8) %>%
    mt(
      state = (state - 1) * 10^4,
      time = setLevels(time, unique(time), c('t', p0('t + ', unique(time)[-1])))
    ) %>% {
      ggplot(
        data = ., 
        aes(
          x = state, 
          y = value
        )
      ) + 
        geom_bar(stat='identity') +
        facet_grid(
          time ~ variable
        ) +
        labs(
          x = 'Engine Mileage',
          y = 'Density'
        ) +
        scale_y_continuous(breaks = c(0, .5, 1)) +
        theme(
          strip.text.x=element_text(size=14),
          strip.text.y=element_text(size=14),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_line(colour = "grey95"), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.y = grid::unit(.4, "cm"),
          legend.position = 'top'
        )
    } %T>% {
      ggsave(plot = ., filename = p0(data.out, "plots/rustModel2.pdf"),  dpi=300, width=7, height=9) 
    } %>% {
      ggsave(plot = ., filename = p0(data.out, "pres/rustModel2.png"),  dpi=300, width=7, height=6.5) 
    }
}

plotRustExample3 <- function() {
  library('gridExtra')
  
  num.actions <- 2
  num.theta <- 1
  dirichlet.alpha <- 1
  num.states <- 12
  lambda <- 1
  iter.stop.threshold <- 10^-10
  theta <- c(1, 10)
  my.beta <- .999
  
  dp <- 
    DPEngine$new(
      num.states = num.states, 
      num.actions = num.actions, 
      num.theta = num.theta,
      dirichlet.alpha = dirichlet.alpha,
      beta = my.beta
    ) %T>% {
      .$initializeDP()
    }
  
  dp$tran$tran.mat <- 
    list(
      cbind(
        rep(0, num.states - 1),
        diag(num.states - 1)
      ) %>%
        rbind(c(rep(0, num.states - 1), 1)),
      cbind(
        rep(1, num.states), 
        matrix(rep(0, num.states^2 - num.states), num.states)
      )
    )
  
  dp$payoff$payoff.stat <- 
    list(
      -matrix(theta[1] * (seq(num.states) - 1)),
      -matrix(theta[2] * rep(1, num.states))
    )
  
  dp$payoff$updatePayoffs(1)
  
  dp$value.fn <- rnorm(num.states, 0, sd = 1)
  
  df(
    iteration = c(0, 2^(seq(0, 6))),
    num.to.do = iteration - lag(iteration, default = 0)
  ) %>%
    dd('iteration', function(s){
      iter.num <- 0
      
      while(s$num.to.do > iter.num){
        dp$operator_VI()
        iter.num <- iter.num + 1
      } 
      
      df(value.fn = -dp$value.fn) %>%
        mt(state = 10000 * (row_number() - 1))
    }) %>% 
    mt(
      iteration = as.factor(iteration),
      iteration = setLevels(
        iteration,
        levels(iteration),
        paste('Iteration', levels(iteration))
      )
    ) %>%
    fl(state <= 100000) %>% {
    ggplot(
      data = ., 
        aes(
          x = state, 
          y = value.fn
        )
      ) + 
        geom_line(size = 1) +
        facet_wrap(
          ~ iteration,
          nrow = 2,
          scales = 'free_y'
        ) +
        labs(
          x = 'Engine Mileage',
          y = 'Expected Discounted Cost'
        ) +
        scale_colour_grey(start = 0, end = .8) +
        scale_x_continuous(breaks = c(25000, 75000)) + 
        theme(
          strip.text.x=element_text(size=14),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_line(colour = "grey95"), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %T>% 
    ggsave(
      plot = ., 
      filename = p0(data.out, "plots/rustModel3.pdf"),
      width=9,
      height=6
    ) %>%
    ggsave(
      plot = .,
      filename = p0(data.out, "pres/rustModel3.png"),
      dpi=300,
      width=8, 
      height=5
    )
    }
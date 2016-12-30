regAdditiveModel <- function() {
  list(
    le('base.R') %>% 
      me(
        funs(wins(., q = c(0, .99))),
        comp, dist, subcomponents, employees, sales, assets, peak, patents, trademarks
      ),
    le('probRegs.rds') %>%
      mt(match.dummy = 1),
    expand.grid(
      dep.var = 'comp ~ ',
      indep.var =
        'dist',
      controls = 
        c(
          '',
          ' + luxury + truck + gen + solvency + subcomponents + sales + peak + patents + trademarks + assets + employees' 
        ),
      fixed.effects = 
        c(
          '',
          '+ as.factor(year) + as.factor(brand) + as.factor(part) + as.factor(country.m) + as.factor(country.s) - 1'
        ),
      heckman.correction = c('', '+ I(-digamma(1)-log(choice.prob))'),
      match.dummy = 1,
      stringsAsFactors = FALSE
    )
  ) %>% 
    rd(ij) %>%
    runRegs %>%
    se('regAdditiveModel.rds')
}

regCompetition <- function() {
  list(
    le('base.R'),
    le('probRegs.rds') %>%
      mt(match.dummy = 1),
    le('supplierDistances.rds') %>% 
      gb(part, plant.m) %>%
      sm(
        first = nth(distance, 1, order_by = distance),
        second = nth(distance, 2, order_by = distance),
        third = nth(distance, 3, order_by = distance),
        fourth = nth(distance, 4, order_by = distance)
      ) %>%
      ug %>%
      me(funs(log(. + 1)), first, second, third, fourth),
    expand.grid(
      dep.var = 'I(log(comp + 1)) ~ ',
      indep.var = 'I(log(dist + 1))',
      controls = 
        c(
          '',
          ' + first + second + third + fourth + luxury + truck + gen + solvency + I(log(subcomponents)) + I(log(sales)) + I(log(peak)) + I(log(1 + patents)) + I(log(1 + trademarks)) + I(log(assets)) + I(log(employees))'
        ),
      fixed.effects = 
        c(
          '',
          '+ as.factor(year) + as.factor(brand) + as.factor(part) + as.factor(country.m) + as.factor(country.s) - 1'
        ),
      heckman.correction = c('', '+ I(-digamma(1)-log(choice.prob))'),
      match.dummy = 1,
      stringsAsFactors = FALSE
    )
  ) %>% 
    rd(ij) %>%
    runRegs %>%
    se('competitionRegs.rds')
}
***********

plotUSAMap <- function() {
  library('ggmap')
  bounds <- c(-105, 23, -65, 50)
  le('base.R') %>%
    fl(
      samp == 0,
      assembler %in% c('FORD', 'GENERAL MOTORS', 'CHRYSLER'),
      part %in% c('ENGINES & ENGINE COMPONENTS', 'TRANSMISSION', 'BRAKING')
    ) %>%
    mt(
      assembler =
        setLevels(
          assembler,
          c('CHRYSLER', 'FORD', 'GENERAL MOTORS'),
          c('Chrysler', 'Ford', 'General Motors')
        ),
      part =
        setLevels(
          part,
          c('ENGINES & ENGINE COMPONENTS', 'TRANSMISSION', 'BRAKING'),
          c('Engines', 'Transmissions', 'Brakes')
        )
    ) %>%
    sl(assembler, part, latMF, lonMF, latSF, lonSF) %>% 
    mt(
      lonMF = pmin(lonMF, bounds[3]),
      lonSF = pmin(lonSF, bounds[3]),
      lonMF = pmax(lonMF, bounds[1]),
      lonSF = pmax(lonSF, bounds[1]),
      latMF = pmin(latMF, bounds[4]),
      latSF = pmin(latSF, bounds[4]),
      latMF = pmax(latMF, bounds[2]),
      latSF = pmax(latSF, bounds[2])
    ) %>% 
    ds %>% {
      qmap(
        location = bounds,
        maptype = 'toner-background'
      ) +
        geom_segment(
          data = .,
          aes(x=lonSF, y=latSF, xend=lonMF, yend=latMF)
        ) +
        facet_grid(part ~ assembler) +
        coord_cartesian(
          xlim = c(bounds[1] + 5, bounds[3] - 5),
          ylim = c(bounds[2] + 3, bounds[4] - 3)
        ) +
        theme(
          strip.text.x=element_text(size=15),      
          strip.text.y=element_text(size=15),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14),
          axis.title.y=element_text(size=18),
          axis.title.x=element_text(size=18),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(fill='white'),
          panel.margin.x = grid::unit(1, "cm")
        ) +
        labs(x = 'Longitude', y = 'Latitude')
    } %>%
    ggsave(
      plot = ., 
      filename = p0(figureSave, 'USAmap.png'),
      dpi=300, 
      width=11,
      height=9
    )
}


regChoiceProbs <- function() {
  library(geosphere)
  library(mlogit)
  
  ij(
    le('base.R') %>% 
      sl(part, plant.m, latMF, lonMF, country.m) %>%
      ds,
    le('supplierLocs.R') %>%
      ds 
  ) %>% {#Plant distances
    mt(
      .,
      distance = 
        distVincentyEllipsoid(
          cbind(lonMF, latMF) %>% as.matrix,
          cbind(lonSF, latSF) %>% as.matrix
        ),
      distance = distance / 1000,
      outsource = country.s != country.m
    ) %>% 
      gb(part, plant.m, supplier) %>% 
      fl(distance == min(distance)) %>% 
      fl(row_number() == 1) %>%
      ug %>%
      sl(part, plant.m, supplier, distance, outsource)
  } %>% {#Contract info
    ij(
      .,
      le('base.R') %>%
        sl(samp, nboot, part, plant.m, supplier) %>%
        rn(supplier.win = supplier)
    ) %>%
      mt(contract.win = supplier == supplier.win) %>%
      sl(-supplier.win)
  } %>% {#Supplier info
    ij(
      .,
      le('base.R') %>%
        sl(supplier, solvency, employees, assets, patents, trademarks) %>%
        ds
    )
  } %>%
    ug %>%
    ag(samp, part, nboot) %>%
    dd(c('samp', 'part'), function(s){ #Regress by part
      try(
        mlogit(
          data = s,
          formula = contract.win ~ I(log(distance + 1)) + outsource + solvency + I(log(1 + patents)) + I(log(1 + trademarks)) + I(log(assets)) + I(log(employees)) | -1, 
          alt.var = 'supplier',
          chid.var = 'nboot',
          shape = 'long'
        ),
        silent = TRUE
      ) %>% {#Try a simpiler regression if first failed
        if(class(.) == "try-error") {
          try(
            mlogit(
              data = s,
              formula = contract.win ~ I(log(distance + 1)) + outsource | -1, 
              alt.var = 'supplier',
              chid.var = 'nboot',
              shape = 'long'
            ),
            silent = TRUE
          )
        } else {
          .
        }
      } %>% {
        if(class(.) == "try-error") {
          return() 
        } else {
          cbind(
            sl(s, nboot) %>%
              ds,
            .$probabilities
          ) %>%
            ad %>% 
            ml('nboot') %>% 
            rn(
              supplier = variable,
              choice.prob = value
            ) %>%
            mt(supplier = as.character(supplier))
        }
      }
    }) %>%
    se('probRegs.rds')
}

regBase <- function() {
  list(
    le('base.R') %>% 
      me(
        funs(wins(., q = c(0, .99))),
        comp, dist, subcomponents, employees, sales, assets, peak, patents, trademarks
      ),
    le('probRegs.rds') %>%
      mt(match.dummy = 1),
    expand.grid(
      coreIndep =
        c(
          'I(log(dist + 1))',
          'outsource'
        ),
      interactors = 
        c(
          '',
          ' * I(log(dist + 1))',
          ' * I(log(subcomponents))',
          ' * luxury',
          ' * truck',
          ' * gen', 
          ' * solvency', 
          ' * I(log(peak))'
        ), 
      controls = 
        ' + luxury + truck + gen + solvency + I(log(subcomponents)) + I(log(sales)) + I(log(peak)) + I(log(1 + patents)) + I(log(1 + trademarks)) + I(log(assets)) + I(log(employees))',
      fixed.effects = 
        c(
          '',
          '+ as.factor(brand) + as.factor(platform) + as.factor(year) + as.factor(part) + as.factor(supplier) + as.factor(country.m) + as.factor(country.s) - 1'
        ),
      heckman.correction = c('', '+ I(-log(choice.prob))'),
      match.dummy = 1,
      stringsAsFactors = FALSE
    )
  ) %>% 
    rd(ij) %>%
    fl(coreIndep != 'I(log(dist + 1))' | interactors != ' * I(log(dist + 1))') %>%
    mt(formula = p0('I(100 * log(comp + 1)) ~ ', coreIndep, interactors, fixed.effects, controls, heckman.correction)) %>% 
    dd(
      c('samp', 'coreIndep', 'interactors', 'controls', 'fixed.effects', 'heckman.correction'),
      function(s){
        lm(first(s$formula), s) %>% 
          Sm %>% {
            .$coefficients
          } %>%
          as.data.frame %>%
          add_rownames
      }
    ) %>%
    se('baseRegs.rds')
}





tabBaseRegressions <- function(){
  le('baseRegs.rds') %T>% {
    fl(., str_detect(rowname, 'factor')) %>%
      sl(rowname) %>% 
      ds %>%
      nrow %>% 
      p0('number of dummies: ', .) %>%
      print()
  } %>% 
    fl(!str_detect(rowname, 'factor')) %>% 
    rn(ste = `Std. Error`) %>% 
    mt(
      controls = ifelse(controls == '', 'No Controls', 'Controls'),
      controls = setLevels(controls, c('No Controls', 'Controls')),
      fixed.effects = ifelse(fixed.effects == '', 'No F.E.', 'F.E.'),
      fixed.effects = setLevels(fixed.effects, c('No F.E.', 'F.E.')),
      model = setLevels(
        p0(coreIndep, interactors),
        c(
          'I(log(dist + 1))',
          'outsource',
          'outsource * I(log(dist + 1))',
          'I(log(dist + 1)) * gen',
          'I(log(dist + 1)) * solvency',
          'I(log(dist + 1)) * luxury',
          'I(log(dist + 1)) * truck',
          'I(log(dist + 1)) * I(log(subcomponents))',
          'I(log(dist + 1)) * I(log(peak))',
          'outsource * gen',
          'outsource * solvency',
          'outsource * luxury',
          'outsource * truck',
          'outsource * I(log(subcomponents))',
          'outsource * I(log(peak))'
        ),
        paste('Model', seq(13))
      ), 
      rowname = setLevels(
        rowname,
        c(
          '(Intercept)',
          'I(log(dist + 1))',
          'outsource',
          'outsource:I(log(dist + 1))',
          'I(log(dist + 1)):gen',
          'outsource:gen',
          'I(log(dist + 1)):solvency',
          'outsource:solvency',
          'I(log(dist + 1)):luxury',
          'outsource:luxury',
          'I(log(dist + 1)):truck',
          'outsource:truck',
          'I(log(dist + 1)):I(log(subcomponents))',
          'outsource:I(log(subcomponents))',
          'I(log(dist + 1)):I(log(peak))',
          'outsource:I(log(peak))',
          'gen',
          'solvency',
          'luxury',
          'truck',
          'I(log(subcomponents))',
          'I(log(peak))',
          'I(log(sales))',
          'I(log(assets))',
          'I(log(employees))',
          'I(log(1 + patents))',
          'I(log(1 + trademarks))'
        ),
        c(
          'Intercept',
          'Log(Distance + 1)',
          'Outsourced',
          'Log(Distance + 1) $ \\cdot $ Outsourced',
          'Log(Distance + 1) $ \\cdot $ Generation',
          'Outsourced $ \\cdot $ Generation',
          'Log(Distance + 1) $ \\cdot $ Solvency',
          'Outsourced $ \\cdot $ Solvency',
          'Log(Distance + 1) $ \\cdot $ Luxury',
          'Outsourced $ \\cdot $ Luxury',
          'Log(Distance + 1) $ \\cdot $ Truck',
          'Outsourced $ \\cdot $ Truck',
          'Log(Distance + 1) $ \\cdot $ Log(Subcomponents)',
          'Outsourced $ \\cdot $ Log(Subcomponents)',
          'Log(Distance + 1) $ \\cdot $ Log(Peak Production)',
          'Outsourced $ \\cdot $ Log(Peak Production)',
          'Generation',
          'Solvency',
          'Luxury',
          'Truck',
          'Log(Subcomponents)',
          'Log(Peak Production)',
          'Log(Sales)',
          'Log(Assets)',
          'Log(Employees)',
          'Log(Patents + 1)',
          'Log(Trademarks + 1)'
        )
      )
    ) %T>% {
      fl(., model %in% paste('Model', seq(3))) %>%
        make.table(
          rowname ~ model + controls + fixed.effects,
          value = 'Estimate',
          table.type = 'Regression', 
          se = 'ste',
          out.file = p0(tableSave, 'univariateDist1.tex'), 
          # align.cols = 'llgggg',
          nsmall = 2,
          hline = TRUE
        )
    } %T>% {
      fl(., model %in% paste('Model', seq(4, 8))) %>%
        make.table(
          rowname ~ model + controls + fixed.effects,
          value = 'Estimate',
          table.type = 'Regression', 
          se = 'ste',
          out.file = p0(tableSave, 'univariateDist2.tex'), 
          # align.cols = 'llgggg',
          nsmall = 2,
          hline = TRUE
        )
    } %>% {
      fl(., model %in% paste('Model', seq(9, 13))) %>%
        make.table(
          rowname ~ model + controls + fixed.effects,
          value = 'Estimate',
          table.type = 'Regression', 
          se = 'ste',
          out.file = p0(tableSave, 'univariateDist3.tex'), 
          # align.cols = 'llgggg',
          nsmall = 2,
          hline = TRUE
        )
    }
  
}

,

runNNPL = function() {
  theta$NNPL <<- rep(0, num.theta)
  dp$resetSystem()
  nestedPsuedoLikelihoodStep <- function(){
    theta$NNPL %>% 
      optim(  
        par = .,
        method = 'BFGS',
        control = list(fnscale=-1, factr = iter.stop.threshold.theta), #max likelihood
        fn = function(theta){
          dp$payoff$updatePayoffs(theta)
          dp$updateCCPs()
          evaluateLikelihood()
        }
      ) %>% {
        theta$NNPL <<- .$par
      }
    dp$updateValues_VI(use.val.shape = TRUE)
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
      control = list(fnscale=-1, factr = iter.stop.threshold.theta), #max likelihood
      fn = function(theta){
        dp$calcDP_value_shape(theta)
        evaluateLikelihood()
      }
    ) %>% {
      theta$NNPL <<- .$par
    }
},

runNNPL2 = function() {
  theta$NNPL2 <<- rep(0, num.theta)
  dp$resetSystem()
  nestedPsuedoLikelihoodStep <- function(){
    theta$NNPL %>% 
      optim(  
        par = .,
        method = 'BFGS',
        control = list(fnscale=-1, factr = iter.stop.threshold.theta), #max likelihood
        fn = function(theta){
          dp$payoff$updatePayoffs(theta)
          dp$updateCCPs()
          evaluateLikelihood()
        }
      ) %>% {
        theta$NNPL2 <<- .$par
      }
    dp$updateValues_VI(use.val.shape = TRUE)
  }
  delta = 1
  while(max(delta) > iter.stop.threshold){
    theta.l <- theta$NNPL
    value.fn.l <- dp$value.fn
    nestedPsuedoLikelihoodStep()
    delta <- max(abs(value.fn.l - dp$value.fn))
  }
}

********
runNNPL = function() {
  theta$NNPL <<- rep(0, num.theta)
  dp$resetSystem()
  nestedPsuedoLikelihoodStep <- function(v){
    theta$NNPL %>% 
      optim(  
        par = .,
        method = 'BFGS',
        control = list(fnscale=-1), #max likelihood
        fn = function(theta){
          dp$payoff$updatePayoffs(theta)
          dp$value.fn <<- v
          dp$updateCCPs()
          dp$updateValues_VI(use.val.shape = TRUE)
          dp$updateCCPs()
          evaluateLikelihood()
        }
      ) %>% {
        theta$NNPL <<- .$par
      }
    dp$updateValues_VI(use.val.shape = TRUE)
  }
  delta = 1
  test = 0
  while(max(delta) > iter.stop.threshold){
    theta.l <- theta$NNPL
    value.fn.l <- dp$value.fn
    nestedPsuedoLikelihoodStep(dp$value.fn)
    print(max(abs(value.fn.l - dp$value.fn)))
    print(max(abs(theta.l - theta$NNPL)))
    print(h(dp$value.fn))
    print(test); test = test+1
    delta <- max(c(abs(theta.l - theta$NNPL), abs(value.fn.l - dp$value.fn)))
  }
}


**************************
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
    dp$updateValues_VI(use.val.shape = TRUE)
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
        dp$calcDP_value_shape(theta)
        evaluateLikelihood()
      }
    ) %>% {
      theta$NNPL <<- .$par
    }
},

********
dp$tran$tran.mat <- 
  list(
    seq(num.states) %>%
      laply(function(l){
        c(rep(0, l - 1), dpois(0:(num.states - l), lambda))
      }) %>% {
        .[, num.states] <- .[, num.states] + 1 - rowSums(.)
        .
      },
    cbind(
      rep(1, num.states), 
      matrix(rep(0, num.states^2 - num.states), num.states)
    )
  )

********
memoryTests <- function(){
  my.computer <- 1
  set.seed(1)
  expand.grid(
    dirichlet.alpha = 1,
    num.theta = 2,
    num.actions = 2,
    # num.states = seq(5000, 20000, 5000)
    num.states = seq(100, 200, 50)
  ) %>% 
    dd('num.states', function(s){
      Estimators$new(
        beta = .999,
        dirichlet.alpha = s$dirichlet.alpha,
        num.theta = s$num.theta,
        num.actions = s$num.actions,
        num.states = s$num.states, 
        num.obs = 300
      ) %>% {
        .$initializeEstimators()
        .$setCasualEstimation()
        
        Rprof(
          filename = p0(varSave, 'memoryTestNPL'),
          memory.profiling=TRUE
        )
        x$runNPL()
        Rprof(NULL)
        
        Rprof(
          filename = p0(varSave, 'memoryTestNNPL'),
          memory.profiling=TRUE
        )
        x$runNNPL()
        Rprof(NULL)
        
        rbind(
          summaryRprof(
            filename = p0(varSave, 'memoryTestNPL'),
            memory = 'stats'
          ) %>% 
            un %>% 
            ad %>% 
            add_rownames %>%
            mt(estimator = 'NPL'),
          summaryRprof(
            filename = p0(varSave, 'memoryTestNNPL'),
            memory = 'stats'
          ) %>% 
            un %>% 
            ad %>% 
            add_rownames %>%
            mt(estimator = 'NNPL')
        ) %>% 
          mutate_(value = '.') %>%
          sl(rowname, estimator, value)
      } 
    }) %>%
    se('memoryTest.rds')


*********
  
  plotEstimationError <- function() {
  le('experiments.rds') %>%
    ml(
      id.vars = c('num.actions', 'num.theta', 'dirichlet.alpha'),
      measure.vars = c('theta.NPL.error', 'theta.SANPL.error', 'theta.ANNPL.error')
    ) %>%
    mt(value = wins(value, c(0, .995))) %>% 
    mt(
      dirichlet.alpha = setLevels(
        dirichlet.alpha, 
        c(0.001, 1),
        c('Sparse Transition Matrix', 'Dense Transition Matrix')
      ),
      num.theta = setLevels(
        num.theta, 
        c(2, 6),
        c('Two Parameters', 'Six Parameters')
      ),
      num.actions = setLevels(
        num.actions, 
        c(2, 6),
        c('Two Actions', 'Six Actions')
      ),
      variable = setLevels(
        variable,
        c('theta.NPL.error', 'theta.SANPL.error', 'theta.ANNPL.error'),
        c('NPL', 'SANPL', 'QNNPL')
      )
    ) %>% {
      ggplot(
        data = ., 
        aes(x = value, colour = variable)
      ) + 
        geom_density(size = .75, adjust = .5) +
        facet_grid(dirichlet.alpha ~ num.theta + num.actions) +
        labs(
          x = 'Estimation Error',
          y = 'Density'
        ) +
        scale_colour_grey(start = 0, end = .8) +
        theme(
          strip.text.x=element_text(size=12),
          strip.text.y=element_text(size=12),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        )
    } %>%
    ggsave(plot = ., filename = p0(data.out, "plots/error.png"),  dpi=300, width=8, height=4)
}


plotIter <- function() {
  le('experiments.rds') %>% 
    ml(
      id.vars = c('num.actions', 'num.theta', 'dirichlet.alpha'),
      measure.vars = c('iter.NPL', 'iter.SANPL', 'iter.ANNPL')
    ) %>%
    mt(value = wins(value, c(0, .98))) %>% 
    mt(
      dirichlet.alpha = setLevels(
        dirichlet.alpha, 
        c(0.001, 1),
        c('Sparse Transition Matrix', 'Dense Transition Matrix')
      ),
      num.theta = setLevels(
        num.theta, 
        c(2, 6),
        c('Two Parameters', 'Six Parameters')
      ),
      num.actions = setLevels(
        num.actions, 
        c(2, 6),
        c('Two Actions', 'Six Actions')
      ),
      variable = setLevels(
        variable,
        c('iter.NPL', 'iter.SANPL', 'iter.ANNPL'),
        c('NPL', 'SANPL', 'QNNPL')
      )
    ) %>% {
      ggplot(
        data = ., 
        aes(x = value, colour = variable)
      ) + 
        geom_density(size = .75) +
        facet_grid(num.theta + num.actions ~ dirichlet.alpha) +
        labs(
          x = 'Number of Iterations',
          y = 'Density'
        ) +
        theme(
          strip.text.x=element_text(size=12),
          strip.text.y=element_text(size=12),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.text.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          legend.title=element_text(size=0),
          legend.text=element_text(size=15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background=element_rect(fill='white', colour='black'),
          strip.background = element_rect(colour="white", fill="white"),
          panel.margin.x = grid::unit(.5, "cm"),
          legend.position = 'top'
        ) +
        scale_colour_grey(start = 0, end = .8)
    } %>%
    ggsave(plot = ., filename = p0(data.out, "plots/iter.png"),  dpi=300, width=7, height=9)
}
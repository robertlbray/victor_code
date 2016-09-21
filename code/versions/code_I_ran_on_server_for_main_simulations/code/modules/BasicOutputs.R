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

basicStats <- function() {
  #Number of iterations
  le('experiments.rds') %>% 
    ml(measure.vars = c('iter.NPL', 'iter.SANPL', 'iter.ANNPL')) %>%
    dd('dirichlet.alpha', function(s){
      lm(value ~ variable - 1, data = s) %>% {
        summary(.)$coefficients
      }  
    })
  
  #Number of iterations
  le('experiments.rds') %>% 
    ml(measure.vars = c('iter.NPL', 'iter.SANPL', 'iter.ANNPL')) %>%
    dd('dirichlet.alpha', function(s){
      lm(value ~ variable - 1, data = s) %>% {
        summary(.)$coefficients %>%
          ad %>%
          add_rownames
      }  
    })
  
  #In sparse matricies, quasi-newton iterates a statistically significantly fewer number of times
  le('experiments.rds') %>% 
    fl(dirichlet.alpha == .001) %>%
    ml(measure.vars = c('iter.SANPL', 'iter.ANNPL')) %>%
    lm(value ~ variable, data = .) %>% {
      summary(.)$coefficients
    }
  
  #Average times
  le('experiments.rds') %>% 
    ml(measure.vars = c('time.NPL', 'time.SANPL', 'time.ANNPL')) %>%
    lm(value ~ variable - 1, data = .) %>% {
      summary(.)$coefficients
    }
  
  #Slopes
  le('experiments.rds') %>% 
    ml(measure.vars = c('time.NPL', 'time.SANPL', 'time.ANNPL')) %>% 
    dd('variable', function(s){
      lm(value ~ num.states, data = s) %>% {
        summary(.)$coefficients
      } %>%
        ad %>%
        add_rownames
    })
  
  #Quasi-newton is 568 seconds faster in the hardest case,  0.001/6/6, which is statistically significant.
  le('experiments.rds') %>% 
    ml(measure.vars = c('time.SANPL', 'time.ANNPL')) %>%
    dd(c('dirichlet.alpha', 'num.actions', 'num.theta'), function(s){
      lm(value ~ variable, data = s) %>% {
        summary(.)$coefficients %>%
          ad %>%
          add_rownames
      }
    }) %>%
    ag(rowname)
  
  #Quasi-newton is significantly worse, but SANPL isn't
  le('experiments.rds') %>% 
    ml(measure.vars = c('theta.NPL.error', 'theta.SANPL.error', 'theta.ANNPL.error')) %>%
    gb(variable) %>%
    mt(value = wins(value, c(0, .95))) %>%
    lm(value ~ variable, data = .) %>% {
      summary(.)$coefficients
    }
}

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

plotTimes <- function() {
  le('experiments.rds') %>% 
    ml(
      id.vars = c('dirichlet.alpha', 'num.theta', 'num.actions', 'num.states'),
      measure.vars = c('time.NPL', 'time.SANPL', 'time.ANNPL')
    ) %>% 
    mt(
      value = value / 60,
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
        c('time.NPL', 'time.SANPL', 'time.ANNPL'),
        c('NPL', 'SANPL', 'QNNPL')
      )
    ) %>%
    gb(dirichlet.alpha, variable) %>%
    mt(value = wins(value, c(0, .97))) %>% {
      ggplot(
        data = ., 
        aes(x = num.states, y = value, colour = variable)
      ) + 
        geom_point(size = .75) +
        geom_smooth(size = .75, se = FALSE) +
        geom_hline(yintercept = 0, colour = 'black', size = .25) +
        facet_grid(dirichlet.alpha ~ num.theta + num.actions, scale = 'free_y') +
        labs(
          x = 'Number of States',
          y = 'Number of Minutes'
        ) +
        scale_colour_grey(start = 0, end = .8) +
        theme(
          strip.text.x=element_text(size=13.5),
          strip.text.y=element_text(size=13.5),
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
    ggsave(plot = ., filename = p0(data.out, "plots/times.png"),  dpi=300, width=10, height=8)
}

plotTimeRatio <- function(variables) {
  le('experiments.rds') %>%
    mt(
      num.states = pmin(8999, num.states),
      num.states = num.states - num.states %% 250
    ) %>% 
    ml(
      id.vars = c('num.states', 'time.NPL', 'time.SANPL', 'time.ANNPL'),
      measure.vars = c('num.theta', 'num.actions', 'dirichlet.alpha')
    ) %>% 
    unite(stat, variable, value) %>%
    gb(num.states, stat) %>%
    sm(
      SANPL = median(time.SANPL)/median(time.NPL),
      QNNPL = median(time.ANNPL)/median(time.NPL)
    ) %>%
    ml(measure.vars = c('SANPL', 'QNNPL')) %>%
    mt(
      stat = setLevels(
        stat,
        c('num.actions_2', 'num.actions_6', 'num.theta_2', 'num.theta_6', 'dirichlet.alpha_0.001', 'dirichlet.alpha_1'),
        c('Two Actions', 'Six Actions', 'Two Parameters', 'Six Parameters', 'Sparse Matrix', 'Dense Matrix')
      )
    ) %>% {
      ggplot(
        data = ., 
        aes(
          x = num.states, 
          y = value, 
          colour = variable
        )
      ) + 
        geom_point(size = 1) +
        geom_line(size = .75) +
        geom_hline(yintercept = 0, colour = 'white') +
        labs(
          x = 'Number of States',
          y = 'Ratio of Median Times'
        ) +
        facet_wrap( ~ stat, nrow = 1) +
        scale_colour_grey(start = 0, end = .7) +
        scale_x_continuous(breaks = c(4000, 6000, 8000)) +
        theme(
          strip.text.x=element_text(size=13.5),
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
    ggsave(plot = ., filename = p0(data.out, "plots/ratio.png"),  dpi=300, width=10, height=5)
}
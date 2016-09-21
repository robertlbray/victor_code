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
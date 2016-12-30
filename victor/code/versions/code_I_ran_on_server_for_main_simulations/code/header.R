options(error=NULL)
# options(error=utils::recover)

laptop <- getwd() == '/Users/robertbray/Dropbox/code/victor/code'
if(laptop){
  data.loc <- '../../../data/victor/'
  funs.loc <- '../../R/funs/'
} else {
  data.loc <- '../data/'
  funs.loc <- '../funs/'
  .libPaths('/home/hailphoebus/R/x86_64-redhat-linux-gnu-library/3.2/')
} 

library('plyr')
l_ply(c('reshape2', 'ggplot2', 'dplyr', 'tidyr', 'stringr', 'doParallel', 'magrittr', 'purrr', 'readr', 'lubridate', 'zoo'), function(l) library(l, character.only=TRUE))

data.in <- paste0(data.loc, 'input/')
varSave <- paste0(data.loc, 'intermediate/variables/')
experimentSave <- paste0(data.loc, 'intermediate/experiments/')
data.out <- paste0(data.loc, 'output/')

l_ply(dir(funs.loc), function(l) source(paste(funs.loc, l, sep="")))
l_ply(dir('modules/', recursive = TRUE), function(l) source(paste('modules/', l, sep="")))
registerDoParallel(cores=32)


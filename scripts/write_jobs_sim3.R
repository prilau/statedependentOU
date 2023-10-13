cat("writing Rev scripts for simulation 3.\n")

# settings
num_tips   = 500
exp_change = c(5, 10, 20)
num_dtraits = 4
num_ctraits = 4


dir.create("jobs", showWarnings=FALSE)

grid = expand.grid(exp_change = exp_change, dtrait=1:num_dtraits,
                   ctrait=1:num_ctraits, stringsAsFactors=FALSE)

# read the stub
lines = readLines("scripts/mcmc_sdOU_simulation_template.Rev")

# make all the job scripts
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {

  this_row         = grid[i,]
  this_exp_change  = this_row[[1]]
  this_num_dtraits = this_row[[2]]
  this_num_ctraits = this_row[[3]]
  
  this_disc_file = paste0('/n500t1r', this_exp_change,
                          'd', this_num_dtraits,
                          '_Discrete.nex")')
  
  this_cont_file = paste0('/n500t1r', this_exp_change,
                          'd', this_num_dtraits,
                          'c', this_num_ctraits)
  
  these_lines = lines

  these_lines[20]  = paste0('tree <- readTrees("data/simulation_3_sameTree_diffRate/n500/t1/tree.tre")[1]')
  these_lines[25]  = paste0('cont <- readContinuousCharacterData("data/simulation_3_sameTree_diffRate/n500/t1/r',
                            this_exp_change, '/d', this_num_dtraits, this_cont_file, '_Continuous.nex")')
  these_lines[30] = paste0('disc <- readDiscreteCharacterData("data/simulation_3_sameTree_diffRate/n500/t1/r',
                           this_exp_change, '/d', this_num_dtraits, this_disc_file)
  these_lines[112] = paste0('monitors.append( mnModel(filename="output/sdOU_simulation3_n500t1r',
                            this_exp_change, 'd', this_num_dtraits,
                            'c', this_num_ctraits, '.log", printgen=10) )')

  this_file = paste0('jobs/sim3_n500t1r', this_exp_change,
                     'd', this_num_dtraits,
                     'c', this_num_ctraits,'.Rev')

  cat(these_lines, file = this_file, sep="\n")

  setTxtProgressBar(bar, i / nrow(grid))

}

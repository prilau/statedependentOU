cat("writing Rev scripts for simulation 1.\n")

# settings
num_tips   = c(100, 250, 500)
num_dtraits = 5
num_ctraits = 5
reps       = 10

dir.create("jobs", showWarnings=FALSE)

grid = expand.grid(num_tips=num_tips, tree=1:reps, dtrait=1:num_dtraits,
                   ctrait=1:num_ctraits, stringsAsFactors=FALSE)

# read the stub
lines = readLines("scripts/mcmc_sdOU_simulation_template.Rev")

# make all the job scripts
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {

  this_row = grid[i,]
  this_num_tips    = this_row[[1]]
  this_tree        = this_row[[2]]
  this_num_dtraits = this_row[[3]]
  this_num_ctraits = this_row[[4]]
  
  this_disc_file = paste0('/n', this_num_tips,
                             't', this_tree, 'd', this_num_dtraits,
                             '_Discrete.nex")')
  
  this_cont_file = paste0('/n', this_num_tips,
                          't', this_tree, 'd', this_num_dtraits,
                          'c', this_num_ctraits)
  
  these_lines = lines

  these_lines[20]  = paste0('tree <- readTrees("data/n', this_num_tips, '/t', this_tree, '/tree.tre")[1]')
  these_lines[25]  = paste0('cont <- readContinuousCharacterData("data/n', this_num_tips,
                            '/t', this_tree, '/d', this_num_dtraits, this_cont_file, '_Continuous.nex")')
  these_lines[30] = paste0('disc <- readDiscreteCharacterData("data/n', this_num_tips,
                           '/t', this_tree, '/d', this_num_dtraits, this_disc_file)
  these_lines[112] = paste0('monitors.append( mnModel(filename="output/sdOU_simulation_n',
                            this_num_tips, 't', this_tree, 'd', this_num_dtraits,
                            'c', this_num_ctraits, '.log", printgen=10) )')

  this_file = paste0('jobs/job_n', this_num_tips,
                     't', this_tree, 'd', this_num_dtraits,
                     'c', this_num_ctraits,'.Rev')

  cat(these_lines, file = this_file, sep="\n")

  setTxtProgressBar(bar, i / nrow(grid))

}

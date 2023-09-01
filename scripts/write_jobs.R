cat("writing Rev scripts for simulation 1.\n")

# settings
num_tips   = c(25, 50, 100)
num_traits = c(1, 2, 4, 8)
true_gamma = c(1, 2, 4, 8)
reps       = 100
runs       = 2
state_dep  = "state_dependent"
noise      = "constant_rate"

dir.create("analyses/simulation_study/simulation_1/jobs", showWarnings=FALSE)

grid = expand.grid(tree=1:reps, num_tips=num_tips, num_traits=num_traits,
                   gamma=true_gamma, state_dep=state_dep, run = 1:runs,
                   noise=noise, stringsAsFactors=FALSE)

# read the stub
lines = readLines("analyses/simulation_study/simulation_1/src/template.Rev")

# make all the job scripts
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {

  this_model = grid[i,]
  these_lines = lines

  these_lines[8]  = paste0("num_tips    = ",   this_model[[2]])
  these_lines[9]  = paste0("num_traits  = ",   this_model[[3]])
  these_lines[10] = paste0("ratio       = ",   this_model[[4]])
  these_lines[11] = paste0("dataset     = ",   this_model[[1]])
  these_lines[12] = paste0("state_model = \"", this_model[[5]],"\"")
  these_lines[13] = paste0("noise_model = \"", this_model[[7]],"\"")
  these_lines[14] = paste0("run_ID      = ",   this_model[[6]])
  these_lines[19] = paste0("seed(",paste0(sample.int(9, replace=TRUE), collapse=""),")")

  this_file = paste0("analyses/simulation_study/simulation_1/jobs/job_",i,".Rev")

  cat(these_lines, file = this_file, sep="\n")

  setTxtProgressBar(bar, i / nrow(grid))

}

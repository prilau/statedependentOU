# statedependentOU
State-dependent OU model for RevBayes. 

## Done
1. Obtained test dataset
    a. Birds (Shiomi 2021, Wright et al. 2022)
    b. Palms (Kissling et al. 2019)
2. Subset the data (in R)
3. First test-run with state-dependent theta and state-independent sigma and alpha
    a. Results: very small alpha for bird data (very large phylogenetic half life)

## To-do
1. Validate that the model is implemented correctly
    a. Compare with state-dependent BM model (set alpha in state-dependent OU model to 0)
    b. Use simulation data
    c. 
2. Testing with empirical data
3. Testing with large dataset
4. Account for intraspecific variations
5. Write simulation function
    a. Define (true) parameter values
    b. Produce n sets of simulated values
    c. Estimate the parameter(s) of the model
    d. Compare the estimated parameter values with true parameter values

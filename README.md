# statedependentOU
State-dependent OU model for RevBayes. 

## Done
1. Obtained test dataset
    a. Birds (Shiomi 2021, Wright et al. 2022)
    b. Palms (Kissling et al. 2019)
2. Subset the data (in R)
3. First test-run with state-dependent theta and state-independent sigma and alpha
    a. Results: very small alpha for bird data (very large phylogenetic half life)

16/08/2023 Updates:
1. Compared with state-dependent BM model - werid
2. Compared with stateless BM model - still weird
3. Compared with stateless OU model - not as bad as above

## To-do
16/08/2023 Updates:
1. Validate that the multiple optima OU part of the model by simluation
    a. Define (true) parameter values
    b. Produce n sets of simulated values
    c. Estimate the parameter(s) of the model
    d. Compare the estimated parameter values with true parameter values
2. Compare the source script for BM and our model at alpha < 1e-20
3. Start drafting the methods part of report (documentating what I have done so far)

Depends on time constraints:
1. Testing with empirical data
2. Testing with large dataset
3. Account for intraspecific variations


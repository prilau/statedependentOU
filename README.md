# statedependentOU
State-dependent OU model for RevBayes. 

## New! Project Overview
1. Validation of model implementation
2. Testing the power of the model<ul>
    <li>same discrete trait-rates, different tree sizes</li>
    <li>different discrete trait-rates, same tree sizes</li>
    <li>discrete trait or continuous trait matters more?</li>
    </ul>
3. Empirical studies





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

31/08/2023 Updates:
1. Completed script to simulate of continuous trait evolution
2. Ran RevBayes on 5 simulations (100-tip * 3, 250-tip * 1, 500-tip * 1)


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

31/08/2023 Updates:
1. Plot current results on R (RevGadgets/ggplot2)
2. Expand simulation datasets
    a. 10 trees per num_tips
    b. 5 discrete trait simulation per tree
    c. 5 continuous trait simulation (same parameters) per discrete trait per tree
    d. 1000-burnin, 10000 iterations, 4 runs per ReyBayes sdOU model fitting
3. Start writing :D


## Side thoughts / Optionals
1. Simulate non-binary discrete trait evolution
2. What if the 250/100-tip trees are subsample from 500-tip trees
    2.1 How does biased sampling affect the inference?
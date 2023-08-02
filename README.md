# statedependentOU
State-dependent OU model for RevBayes. 

## Done
1. Obtained test dataset (birds)
2. Subset the data - 2 binary discrete traits and 1 continuous trait
3. First test-run with state-dependent theta and state-independent sigma and alpha
    a. Results: very small alpha (very large phylogenetic half life)

## To-do
1. Find different datasets to fit the model.
2. Fit the current dataset to different models
    a. Change flexibility of different parameters
    b. Try combined discrete trait (flapper and migrant, non-binary)
    c. Compare MAP (maximum a posteriori) and random character history (aka what we do right now)
3. Size scaling for bird body measurement dataset
4. Account for intraspecific variations
5. Write simulation function
    a. Define (true) parameter values
    b. Produce n sets of simulated values
    c. Estimate the parameter(s) of the model
    d. Compare the estimated parameter values with true parameter values

## Timeline (tentative)
early August 2023   To-do #2-3
late August 2023    To-do #5
September 2023      Let's see...
October 2023

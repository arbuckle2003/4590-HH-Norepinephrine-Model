**What each file does:**
`readFiles.m` Reads the long square current-step sweeps from an NWB file 
`extractFeatures.m`  Measures rheobase, f I slope, adaptation, spike shape, and passive properties
`hhModel.m`  The model: simulates 1 second current steps for one or many parameter sets 
`paramsPyr.m`, `paramsFS.m`  Baseline parameters for each cell 
`runModel.m`  Compares a baseline model with its recording
`runNE.m`  Preliminary result: AHP block of 0–75% on the baseline pyramidal model 
`popFit.m` Intermediate 1: random-walk search for every parameter set that fits each cell 
`plotPopulation.m` Intermediate 1 figures: parameter spread, fits, and correlations 

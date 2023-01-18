The folder contains the standalone code to replicate the IRF reported in figures 3,4,5,6 
of the paper. Figures 9,10,11,12 can be also replicated by including a different set of 
variables in the code, as well as some of the robustness exercises. The data also suffice 
to replicate Figure 1.

Files in the folder:

simulationLP: Writes vectors for the simulations and their differences
simulateLP: Computes the effects of plans from the betas on variables and Theta
estimatePLP: Estimate the interacted regressions and save beta, erros, and varcov matrices
estimatePLPboot: estimate regressions in the bootsrap (without cleaning variables for missing obs.
IPLPboot: bootstrap
dopicsmain: figures routine
impcompareoutout: figures subroutine
panellag: lagging function for panel (Towbin and Weber IPVAR model)
paneldiff: takes first differences for panel (Towbin and Weber IPVAR model)
minpercentile: eturns the upper and lower value of the values in X which define the (1-pctile) perecent of the values (Towbin and Weber IPVAR model)
data.mat: data in the matlab format
data.csv: data in the csv format
data.txt: data in the text format

The codes are deeply indebted to Towbin and Weber (2013) IPVAR model codes.
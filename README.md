# Modelling Drug Delivery
This MATLAB code is designed to modelling the drug concentration distribution in human tissue. This helps scientists determine how a drug applied on the skin diffuses inside the tissue and what dose is required to reach a certain depth with a given concentration.
This are key parameters when developping a drug.

## Overview

Transdermal drug delivery offers a non invasive drug administration together with reduced side effects given by oral assumptions. 
Investigation of drug penetration across the skin is important in topical pharmaceutical formulations.
This code uses the Finite Element Method (FEM) to resolve the transient diffusion-reaction equation in a 1D mesh.
When a block of human skin is considered, the problem can be divided in discrete layers represented by a linear mesh. 
The diffusion reaction equation can be used to model the concentration of the drug in each part of the skin at any point in time.
Having resolved how a drug propagates in the skin, the code finds the minimum dose needed to have the desired pharmaceutical effect.

More information and results are available in the report `Transient_FEM_modelling.pdf`

### Code

The code is structured in two parts further eplained in the report.
The first part, with the main function in `TransientFEM.m`, simply solves the diffusion-reaction equation using the FEM.
The second part, using the same main file, applies the diffusion-reaction equation to human tissue.

The code makes extensive use of functions. The main one is `TransientFEM` all other functions are called by the main or other sub functions.

### Execution

To run the code in part one, execute the `Part1Plotter.m` or the `Part1L2Norm.m` MATLAB scripts.
To run the code in part two, execute the `Part2Plotter.m` MATLAB script.
Check the functions correct working by running the Unit Tests in the files `*UT.m`

Refer to the report for further explanation on the results.

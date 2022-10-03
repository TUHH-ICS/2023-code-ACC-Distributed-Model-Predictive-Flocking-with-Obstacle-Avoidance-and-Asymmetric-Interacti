# 2023-code-ACC-Distributed Model Predictive Flocking with Obstacle Avoidance and Asymmetric Interaction Forces

## General

This repository contains an implementation of the algorithms and simulations described in 

> P.Hastedt and H. Werner, "Distributed Model Predictive Flocking with Obstacle Avoidance and Asymmetric Interaction Forces"

submitted to ACC, 2023.

It may be used to recreate and validate the simulation results and figures from the paper. To do so, run either of the two scripts `simulation.m` and `evaluation.m`.

Additionally, videos for the scenarios described in the paper are provided in the `videos` directory.

Running the simulations can take up to 10 minutes depending on the computer hardware.



## Simulation

For the simulations, an open source MAS library which can be found [on GitHub](https://github.com/TUHH-ICS/MAS-Simulation) is utilized.

At the top of `simulation.m`, the algorithm and scenario to be simulated can be selected by changing the `algorithmIndex` and `scenarioIndex` variables. The simulation results will be saved in the `simulation/out` directory and can then be used for evaluation.

## Evaluation

At the top of `evaluation.m`, the scenarios to be compared can be selected by adding the corresponding data index to the `dataSelection` array.  To evaluate additional data generated by the simulation, copy the `.mat` files from the `simulation/out` directory to the `data` directory and add the name of the data file to the `simData` array at the top of `evaluation.m`.

## Prerequisites

When downloading the code from Zenodo, the MAS-simulation submodule directory `simulation/MAS-simulation` will be empty. This can be resolved by either directly downloading the code for the paper from GitHub or by copying the source code of the [MAS library](https://github.com/TUHH-ICS/MAS-Simulation) to the corresponding directory.

The code in this repository was tested in the following environment:

* *Windows 10* Version 21H2
* *Matlab* 2021a

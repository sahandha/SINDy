# Sparse Identification of Nonlinear Dynamical Systems

## Overview:

This is a python implementation of the SYNDy algorithm introduced by Brunton et al. of University of Washington, Seattle.
See [paper.pdf](https://github.com/sahandha/SINDy/blob/master/paper.pdf) for the details of the method. You can also watch a presentation of it [here](https://www.youtube.com/watch?v=gSCa78TIldg).

The algorithm is used to identify governing equations for nonlinear systems for which we have data, but no mathematical model.

## Files:

#### demo.ipynb
A jupyter notebook showing an example of deducing a math model purely from data.

#### SINDy.py
This file includes all the implementation details of the algorithm as well as a demo example in the main method.

#### SIR.py
Implementation of the classical SIR model which is used as an example to demonstrate how SYNDy works.

#### paper.py
The PNAS publication introducing the algorithm.

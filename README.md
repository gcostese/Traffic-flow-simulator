# Traffic-flow-simulator

This repository contains some Matlab code to solve a scalar conservation law in 1D, known in the traffic flow literature as the LWR model. It uses the Godunov numerical scheme which is a first order finite volume scheme. It can accomodate any kind of flux function, not necessarily the "triangular fundamental diagram" as in the classical Daganzo's Cell Transmission Model.
There are also some Matlab code to solve Riemann problems on simple junctions, namely 1-to-2 diverge and 2-to-1 merge.

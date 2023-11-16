# Tutorial

This page explains how to run the code after download


## Instal Pre-requisites

Install auxiliar packages 

```
$ cd cerisse/Submodules/
$ ./install safe
```
It will connect to Github and download the required packages.
`$ ./install git`, will intall latest release commit in the **development** branch
of AMReX


### AMReX

This is the AMR library that controls grdi generation/movement/IO/parallelization, etc.
Is an approx 30 MB download, it wil expand to a folder 27 M , the install safe
option will install version **23.11**


### PelePhysics

This is the library that control chemistry
Is an approx 30 MB download, it wil expand to a folder 146 M, the install safe
option will install version **23.03**


### Sundials

SUNDIALS - a SUite of Nonlinear and DIfferential/ALgebraic equation Solvers. 
Approx 30 Mb zip download, expand to 
the install safe option will install version **6.2.2**

## Tutorial 1

Go to Exec folder and pick one example. In this Tutorial  we we will work with
Sod1D

```
$ cd cerisse/EB_CNS/Exec/Sod1D
```

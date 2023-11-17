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


## Tutorial 1

### 1) Go to Problem Folder

Go to Exec folder and pick one example. In this Tutorial  we we will work with
Sod1D

```
$ cd cerisse/EB_CNS/Exec/Sod1D
```

### 2) Install SUNDIALS

SUNDIALS - a SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.

controlled by AMREX Options so it can be targeted. Do it once unless toying with chemistry.
It is a 30 M install done within PelePhysics.

```
$ make SUNDIALS
```

### 3) Compile code

```
$ make
```


### 4) Run


Screenshot


### 5) See Results



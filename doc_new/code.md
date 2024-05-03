# Code Structure


# Index.h


Lots of usefuel stuff, includes teh definition of global **indicies**
For the conserved varaibles
```cpp
static constexpr int UMX=0;
static constexpr int UMY=1;
...
```
Fore the primtive (derived) variables

```cpp
static constexpr int QRHO=0;
static constexpr int QU=1;
static constexpr int QV=2;
static constexpr int QW=3
```
As well as the number of equations to solve and number of halo points

```cpp
static constexpr int NCONS=UFS + NUM_SPECIES;
static constexpr int NPRIM=QFS + NUM_SPECIES;
static constexpr int NGHOST=3; // TODO: make it an automatic parameter
```

The names of the variables will also be defined here

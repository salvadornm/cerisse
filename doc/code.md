# Code Structure


# Index.h


Lots of usefuel stuff, includes the definition of global **indicies**
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
...
```
As well as the number of equations to solve and number of halo points

```cpp
static constexpr int NCONS=UFS + NUM_SPECIES;
static constexpr int NPRIM=QFS + NUM_SPECIES;
static constexpr int NGHOST=3; // TODO: make it an automatic parameter
```

## Global Variables description

### Primitve Variables Index

| Index                      |  Default | Description                                                  |
| --------------------------- | ------------- | ------------------------------------------------------------ |
|   `QRHO`             | 0       | density    |
|   `QU`               | 1       | x-velocity    |
|   `QV`               | 2       | y-velocity    |
|   `QW`               | 3       | z-velocity    |
|   `QT`               | 4       | temperature    |
|   `QPRES`            | 5       | pressure   |
|   `QC`            | 6          | ?  |
|   `QG`            | 7          | ?   |
|   `QEINT`            | 8       | internal energy  |
|   `QFS`            | 9         | index first specie   |

### Conservative


| Index                      |  Default | Description                                                  |
| --------------------------- | ------------- | ------------------------------------------------------------ |
|   `UMX`             | 0       | density    |
|   `UMY`             | 1       | x-velocity    |
|   `UMZ`             | 2       | y-velocity    |
|   `UET`             | 3       | z-velocity    |
|   `URHO`            | 4       | temperature    |
|   `UFS`             | 4       | last index, replace density when solving species  |


### Other
| Index                      |  Default | Description                                                  |
| --------------------------- | ------------- | ------------------------------------------------------------ |
|   `NCONS`               | 5       | number of conservative variables    |
|   `NPRIM`               | 6       | number of primitive variables    |

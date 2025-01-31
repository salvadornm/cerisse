# Variables

## Code Structure

## Index.h

Lots of useful stuff, includes the definition of global **indicies.** For the conserved variables they are:

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
static constexpr int QT=4;
static constexpr int QPRES=5;
static constexpr int QC=6;
static constexpr int QG=7;
static constexpr int QEINT=8;
static constexpr int QFS=9;
```

As well as the number of equations to solve and number of halo points

```cpp
static constexpr int NCONS=UFS + NUM_SPECIES;
static constexpr int NPRIM=QFS + NUM_SPECIES;
static constexpr int NGHOST=3; // TODO: make it an automatic parameter
```

### Global Variables description

#### Primitve Variables Index

| Index   | Default | Description        |
| ------- | ------- | ------------------ |
| `QRHO`  | 0       | density            |
| `QU`    | 1       | x-velocity         |
| `QV`    | 2       | y-velocity         |
| `QW`    | 3       | z-velocity         |
| `QT`    | 4       | temperature        |
| `QPRES` | 5       | pressure           |
| `QC`    | 6       | sound speed        |
| `QG`    | 7       | gamma              |
| `QEINT` | 8       | internal energy    |
| `QFS`   | 9       | index first specie |

#### Conservative

| Index  | Default | Description                                            |
| ------ | ------- | ------------------------------------------------------ |
| `UMX`  | 0       | x-momentum _rho u_                                     |
| `UMY`  | 1       | y-momentum _rho v_                                     |
| `UMZ`  | 2       | z-momentum _rho w_                                     |
| `UET`  | 3       | total energy _rho e_                                   |
| `URHO` | 4       | density                                                |
| `UFS`  | 4       | first spec index, replace density when solving species |

#### Transport Properties

| Index   | Default | Description                       |
| ------- | ------- | --------------------------------- |
| `CMU`   | 0       | dynamic viscosity                 |
| `CXI`   | 1       | bulk viscosity                    |
| `CLAM`  | 2       | heat conductivity                 |
| `CRHOD` | 3       | first species diffusivity _rho D_ |

#### Other

| Index   | Default | Description                      |
| ------- | ------- | -------------------------------- |
| `NCONS` | 5       | number of conservative variables |
| `NPRIM` | 6       | number of primitive variables    |

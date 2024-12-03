---
icon: flask-vial
cover: >-
  https://images.unsplash.com/photo-1693919653649-27492e78899d?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxDaGVtaXN0cnl8ZW58MHx8fHwxNzMzMDY4Mzc4fDA&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Chemistry

Complex chemistry in the code is handled through the [PelePhysics](https://pelephysics.readthedocs.io/en/latest/index.html) library. To install PelePhysics, type`$ ./install.sh pelephys` in the `lib` directory. It may be alrady installed, if the code was installed using _safe_ (check  [QuickStart](quickstart.md))

## set-up

In the `GNUmakefile` the user needs to set-up the following options

```bash
# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := BurkeDryer
```

and also set the following lines exists, pointing to PelePhysics.

```bash
USE_PELEPHYSICS = TRUE
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)
```

The first time running PelePhysics, third-party libraries have to be compiled

```bash
make TPL
```

## Chemical Mechanisms

The chemical mechanism are located in `PelePhysics/Support/Mechanism/Models`

```bash
$ ls lib/PelePhysics/Support/Mechanism/Models
Aromatic_KrNara LuDME               decane_3sp          heptane_lu_88sk
BurkeDryer      LuEthylene          dodecane_lu         heptane_lu_qss
C1-C2-NO        LuEthylene_qss      dodecane_lu_qss     isooctane_lu
C1-C2-NO_qss    Make.package        dodecane_wang       list_mech
CH4_lean        NUIGalway	    dodmethair_4sp      list_qss_mech
CH4_lean_qss    Null		    drm19		methaneIons_diRenzo
Davis           PPHYS_CONSTANTS.H.  ethylene_af		ndodecane_35
FFCM1_Red       SootReaction	    grimech12		nitrogens
HP_DME          air		    grimech30		propane_fc
IonizedAir      alzeta		    grimech30-noArN	sCO2
JL4	        chem-CH4-2step	    header
Kolla           chem-H		    heptane_3sp
LiDryer         converter.sh	    heptane_fc
```

All of them are in _yaml_ fromat used in [Cantera](https://cantera.org). The list of chemical mechanisms directlly available can be obtained direclty by  `$cat PelePhysics/Support/Mechanism/Models/list_mech` with a few more that use QSSA (see `list_qss_mech` file). Opening the `mechanism.yaml` file within a directory, will give an indication of the species involved and chemical reactions

For example, in the Jones and Lindstedt mechanism (a 4-step process for hydrocarbon combustion), seven species are used. The label **phases** indicates which chemical components will be included.

```
phases:
- name: gas
  thermo: ideal-gas
  elements: [C, O, H, N]
  species: [CH4, O2, H2O, N2, CO, CO2, H2]
  kinetics: gas
  transport: mixture-averaged
  state: {T: 300.0, P: 1 atm}
```

In the above case, the 7 species will be CH4, O2, H2O, N2, CO, CO2, H2. Note, that the order is important to set the corresponding species in the PROB file. For example, CH4 will be solved. The label **species**, indicate thermodynamic and transport properties of the chemical species.

```
species:
- name: CH4
  composition: {C: 1, H: 4}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11,
      -1.02466476e+04, -4.64130376]
    - [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.0181523e-13,
      -9468.34459, 18.437318]
    note: L 8/88
  transport:
    model: gas
    geometry: nonlinear
    well-depth: 141.4
    diameter: 3.746
    polarizability: 2.6
    rotational-relaxation: 13.0
- name: O2
  ...
```

Similarly the **reactions** label describes the chemical reactiosn used.

```
reactions:
- equation: 2 CH4 + O2 => 2 CO + 4 H2  # Reaction 1
  rate-constant: {A: 3.91e+13, b: 0.0, Ea: 3.0e+04} 
  orders: {CH4: 0.5, O2: 1.25}
- equation: CH4 + H2O <=> CO + 3 H2  # Reaction 2
  rate-constant: {A: 3.0e+11, b: 0.0, Ea: 3.0e+04}
  ...
```

See details of yaml format in [YAML](https://cantera.org/tutorials/yaml/defining-phases.html)

### Generate a new mechanism

For all the available mechanisms, a Cantera yaml format is provided. If CHEMKIN files are present Pelephysics rely on Cantera’s _**ck2yaml**_ utility to convert CHEMKIN files to the Cantera yaml format. They are converter scripts to faciliatte this process. Check [PelePhysics Tutorial](https://pelephysics.readthedocs.io/en/latest/EOS.html)

## Equation of State

Using PelePhysics there are three EOS models:

* A simple GammaLaw model for a single component perfect gas. This can also be used directly without PelePhysics (see PROB)
* An multi-component ideal gas labelled **Fuego** Used for multi-component (reacting or not) calculations.
* The Soave-Redlich-Kwong cubic equation of state for a general mixture of species

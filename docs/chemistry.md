---
icon: flask-vial
cover: >-
  https://images.unsplash.com/photo-1693919653649-27492e78899d?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw2fHxDaGVtaXN0cnl8ZW58MHx8fHwxNzMzMDY4Mzc4fDA&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Chemistry

Complex chemistry in the code is handled through the [PelePhysics](https://pelephysics.readthedocs.io/en/latest/index.html) library. To install PelePhysics, type`$ ./install.sh pelephys` in the `lib` directory. It may be already installed if the code was installed using _safe_ (check [QuickStart](quickstart.md))

## set-up

In the `GNUmakefile` the user needs to set-up the following options

```bash
# PelePhysics
EOS_MODEL := FUEGO
TRANSPORT_MODEL := SIMPLE
CHEMISTRY_MODEL := BurkeDryer
```

and also  the following lines exists, pointing to PelePhysics library.

```bash
USE_PELEPHYSICS = TRUE
PELE_PHYSICS_HOME = $(abspath ../../lib/PelePhysics)
```

The first time running PelePhysics, third-party libraries have to be downloaded compiled

```bash
make TPL
```

## Chemical Mechanisms

The chemical mechanism are located in `PelePhysics/Support/Mechanism/Models`\
(PelePhysics **v23**) or `PelePhysics/Mechanisms` (PelePhysics **v25**). The mechanisms are organised in directories, that must have the same name as `CHEMISTRY_MODEL` option in the `GNUMakefile`

```bash
$ ls lib/PelePhysics/Support/Mechanism/Models
Aromatic_KrNara LuDME               decane_3sp          heptane_lu_88sk
BurkeDryer      LuEthylene          dodecane_lu         heptane_lu_qss
C1-C2-NO        LuEthylene_qss      dodecane_lu_qss     isooctane_lu
C1-C2-NO_qss    Make.package        dodecane_wang       list_mech
CH4_lean        NUIGalway           dodmethair_4sp      list_qss_mech
CH4_lean_qss    Null                drm19               methaneIons_diRenzo
Davis           PPHYS_CONSTANTS.H.  ethylene_af         ndodecane_35
FFCM1_Red       SootReaction        grimech12           nitrogens
HP_DME          air                 grimech30           propane_fc
IonizedAir      alzeta              grimech30-noArN     sCO2
JL4             chem-CH4-2step      header
Kolla           chem-H              heptane_3sp
LiDryer         converter.sh        heptane_fc
```

A local mechanism can be used by setting `USE_LOCALCHEM = TRUE` and specify a `LOCALCHEM_PATH` (absolute path) in `GNUMakefile`. To build the code, the mechanism folder needs several files: `mechanism.H`,`mechanims.cpp` and  `Make.package`  that hard-coded the chemistry, see an example of local chemistry set-up in `exm/ibm/srp`.

All mechanisms are derived using _yaml_ format from [Cantera](https://cantera.org). The list of chemical mechanisms directlly available can be obtained directly by `$cat PelePhysics/Support/Mechanism/Models/list_mech` with a few more that use QSSA (see `list_qss_mech` file). Opening the `mechanism.yaml` file within a directory, will give an indication of the species involved and chemical reactions.

For example, in the Jones and Lindstedt mechanism (a 4-step process for hydrocarbon combustion), seven species are used. In the `mechanism.yaml`  file, the label **phases** indicates which chemical components will be included.

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

In the above case, the 7 species will be $$\text{CH}_4$$, $$\text{O}_2$$, $$\text{H}_2\text{O}$$, $$\text{N}_2$$, $$\text{CO}$$, $$\text{CO}_2$$, $$\text{H}_2$$. Note that the order  will correspond to the species order in the `prob.h` file. For example, $$\text{CH}_4$$ will be solved first in the species. The label **species**, indicate thermodynamic and transport properties of the chemical species.

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

Similarly the **reactions** label describes the chemical reactions used.

```
reactions:
- equation: 2 CH4 + O2 => 2 CO + 4 H2  # Reaction 1
  rate-constant: {A: 3.91e+13, b: 0.0, Ea: 3.0e+04} 
  orders: {CH4: 0.5, O2: 1.25}
- equation: CH4 + H2O <=> CO + 3 H2  # Reaction 2
  rate-constant: {A: 3.0e+11, b: 0.0, Ea: 3.0e+04}
  ...
```

See details of _yaml_ format in [YAML](https://cantera.org/tutorials/yaml/defining-phases.html)

### Generate a new mechanism

For all the available mechanisms in PelePhysics, a Cantera _yaml_ format is provided. If CHEMKIN files are present, Pelephysics rely on Cantera’s _**ck2yaml**_ utility to convert CHEMKIN files to the Cantera _yaml_ format. They are converter scripts to facilitate this process. Check [PelePhysics Tutorial](https://pelephysics.readthedocs.io/en/latest/EOS.html)

Once the `mechanism.yaml` is generated, is good idea to check that is readable by Cantera, there are scripts in `tools/combustion` can be helpful). The _yaml_ format needs to generate two chemistry-specific files: `mechanism.cpp` and `mechanism.H`. These are the files that Cerisse/Pele codes require to run (Cantera is not needed to run cerisse). To convert use the CLI utility `ceptr` (located in `Pelephysics/Support/ceptr`), which requires the [poetry](https://python-poetry.org/docs/) package manager.

{% hint style="warning" %}
Beware that`ceptr` requires Python version between (>=3.8 and <3.11). Is best to work with environments: **pyenv**, **conda** or similar. This may require careful python install (see [Tips](tips.md))
{% endhint %}

The script command is then

```bash
$ cd ${PELE_PHYSICS_HOME}/Support/ceptr
$ poetry run convert -f ${PATH_TO_CHEMISTRY}/mechanism.yaml
```

This will create the required files in `${PATH_TO_CHEMISTRY}`A  similar script can be used to convert CHEMKIN files to _yaml_ (and then use the above script to generate the C++ files).

```bash
$ cd ${PELE_PHYSICS_HOME}/Support/ceptr
$ poetry run ck2yaml --input ${PATH_TO_CHEMISTRY}/mechanism.inp 
--thermo ${PATH_TO_CHEMISTRY}/therm.dat --transport ${PATH_TO_CHEMISTRY}/tran.dat --permissive 
```

More details in Cantera [CK2YAML](https://cantera.org/3.1/userguide/ck2yaml-tutorial.html) documentation

## Equations of State

Using PelePhysics there are three Equations of State (EOS) models:

* A simple GammaLaw model for a single component perfect gas. This can also be used directly without PelePhysics (see PROB)
* A multi-component ideal gas labelled **Fuego** Used for multi-component (reacting or not) calculations.
* The [Soave-Redlich-Kwong](https://doi.org/10.1016/0009-2509\(72\)80096-4) cubic equation of state for a general mixture of non-ideal gases

## Chemistry Integration

PelePhysics has different options to integrate the chemistrym using implicit ODE solvers: CVODE an DVODE. CVODE is part of a software family called sundials for SUite of Nonlinear and DIfferential / ALgebraic equation Solver [SUNDIALS](https://dl.acm.org/doi/10.1145/1089014.1089020). which should be installed (see [install SUNDIALS](quickstart.md#installation-sundials)). The methods implemented in CVODE are _variable-order_, _variable-step_ multistep methods, based on formulas broadly following

$$
\alpha_n y_k + \beta_n h_n \dot{y}_k = 0
$$

that require a linear solver, PelePhysics has different options, see [PelePhsyics documentation](https://amrex-combustion.github.io/PelePhysics/IntroductionToCvode.html). To activate CVODE, in the `GNUMakefile`

```
USE_SUNDIALS_PP = TRUE
```

### input options

The input file consists of specific blocks containing keywords that apply to different aspects of the problem's integration. Each block is associated with a suffix that helps determine which keywords are required, depending on the options set in the `GNUmakefile`. For example, if CVODE is enabled via the `GNUmakefile`, keywords starting with `cvode.*` are relevant. The general `ode.*` keywords are shared across all ODE integrators and are also applicable to CVODE.

### Key Parameters and Options

| Keyword                   | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ode.reactor_type`        | Switches between a **CV reactor** and a **CVH reactor**.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `cvode.solve_type`        | <p>Controls the CVODE linear integration method:<br><strong>1</strong> = Dense direct linear solver<br><strong>5</strong> = Sparse direct linear solver (requires KLU library)<br><strong>99</strong> = Krylov iterative solver</p>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| `ode.analytical_jacobian` | <p>Determines the Jacobian solver method, with different behaviors based on <code>cvode.solve_type</code>:</p><ul><li><strong>If <code>cvode.solve_type = 1</code></strong> <code>ode.analytical_jacobian = 1</code> enables the <strong>Analytical Jacobian</strong>.</li><li><strong>If <code>cvode.solve_type = 99</code></strong> <code>ode.analytical_jacobian = 1</code> activates the <strong>preconditioned GMRES solver</strong>, while setting <code>ode.analytical_jacobian = 0</code> enables the <strong>non-preconditioned GMRES solver</strong>.</li><li><strong>If <code>cvode.solve_type = 99</code></strong> and the <strong>KLU library is linked</strong>, then the preconditioned solver operates in a <strong>sparse format</strong>.</li><li><strong>If <code>cvode.solve_type = 5</code></strong>, the only valid option is <code>ode.analytical_jacobian = 1</code>.</li></ul> |

This structure allows users to configure the solver behavior efficiently based on their requirements.

### Premixed Flame Initialisation

Pre-computed profiles from 1D freely propagating premixed flames can be used to initialise a solutions.\
To create a 1D profile, a python script is used `tools/combustion/1d-flame-run.py` that uses\
Cantera and requires a chemistry files in _yaml_ format. The script can be modified to create different profiles based on different conditions. The flame profiles can be seen by

```bash
$ python 1d-flame-plot.py
```

A `pmf-Y.txt` file (or `pmf-X.txt` depending on the options) will be created. PelePhysics requires a specific \*_.dat_ format that is used in Cerisse as well. To convert the files, you can use the Python script `pmf_f2dat.py` in `lib/PelePhysics/Utility` (PelePhysics v23)

```bash
$ python pmf_f2dat.py inputfile.txt pmf.dat
```

The `pmf.dat` file will look like (for a 9-species hydrogen mechanism)

```
VARIABLES  =  "X"  "temp"  "u"  "rho"  "Y_H2" "Y_O2" "Y_H2O" "Y_H" "Y_O" "Y_OH" "Y_HO2" "Y_H2O2" "Y_N2"
 ZONE  I=454  FORMAT=POINT  SPECFORMAT=MASS 
0.0 298.0 0.534292484155303 0.9893277276003428 0.014467517837899192 0.22962878760502656 1.1722583290927352e-19 -5.887036225686357e-21 6.511836937427712e-19 -3.6814739215435833e-19 -3.4206663655715174e-18 -1.7741231617078512e-18 0.7559036945570743
```

And is the file that the code needs, the place can be selected in the input file (see [example](twodim.md#planar-flame) in `exm/planar flame`).\
This set-up does not reuire to set-up poetry, but additional options can be done following[Pelephysics documentation](https://pelephysics.readthedocs.io/en/latest/Utility.html)

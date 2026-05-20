---
icon: list-radio
---

# Options

This page explains the arguments of each class that can be defined in `prob.h`.\
Most of these options are passed as an argument in set-up through a light structure. All arguments are defined as `constexpr` so they are known and evaluated at compile time. This helps optimize the code (at compile time expense).\
A quick tip on checking the available options is to look in **defaultparm\_t** in `NumParam.h`\
where all options are defined to avoid  name conflicts.

## closures\_dt

### transport\_const\_t

Example use is

```cpp
typedef closures_dt<indicies_t, transport_const_t<viscparm_t>, ...  > ProbClosures;
```

where specific values of conductivity and viscosity are selected

```cpp
struct viscparm_t {

  public:

  static constexpr Real conductivity = 0.001;      // conductivity (for constant value)
  static constexpr Real viscosity    = 1e-5;       // viscosity    (for constant value)
};
```

### Smagorinsky\_t / WALE\_t

Both **Smagorinsky** and **WALE** are sub-classes of **LES** and have the same parameters.

Example use is

```cpp
typedef closures_dt< ... , Smagorinsky_t<LESparm_t,indicies_t> > ProbClosures;
```

with `LESparm_t` containing all parameters required for LES models. Some may not be used\
if specific models are used, but ALL have to be defined if one od these templates are used.

```cpp
struct LESparm_t {

  public:

  static constexpr int  order = 2;         
  static constexpr Real Pr_o_Prsgs = 0.8;  
  static constexpr Real Scsgs = 0.7;       
  static constexpr Real Cs = 0.18;          
  static constexpr Real CI = 0.08;        
  static constexpr bool fixDelta = false;  
  static constexpr Real Delta = 0.02;      
};
```

| Variable Name | Type   | Value   | Description                                                                                |
| ------------- | ------ | ------- | ------------------------------------------------------------------------------------------ |
| `order`       | `int`  | `2`     | Order of the numerical scheme used for gradient estimation.                                |
| `Pr_o_Prsgs`  | `Real` | `0.8`   | Ratio of molecular Prandtl number (Pr) to subgrid-scale Prandtl number (Pr<sub>sgs</sub>). |
| `Scsgs`       | `Real` | `0.7`   | Subgrid-scale Schmidt number (Sc<sub>sgs</sub>).                                           |
| `Cs`          | `Real` | `0.18`  | Smagorinsky constant; used in sgs models.                                                  |
| `CI`          | `Real` | `0.08`  | Yoshizawa constant; associated with isotropic sgs contibution.                             |
| `fixDelta`    | `bool` | `false` | Flag to indicate whether the filter width is fixed.                                        |
| `Delta`       | `Real` | `0.02`  | Filter width (only used if `fixDelta` is true)                                             |

### TFM\_t

The `TFM_t` class reads its parameters from the `atf` namespace in the input file. These values can be used to override the default thickened-flame model settings.&#x20;

| Input option                  |     Default | Description                                                                     |
| ----------------------------- | ----------: | ------------------------------------------------------------------------------- |
| `atf.thickening_factor`       |       `1.0` | Maximum artificial thickening factor, `F0`. A value of `1` disables thickening. |
| `atf.laminar_flame_speed`     |       `2.0` | Laminar flame speed, `SL0`, used by the wrinkling model.                        |
| `atf.laminar_flame_thickness` |    `1.0e-4` | Laminar flame thickness, `deltaL0`, used to normalize the filter width.         |
| `atf.Ret`                     |        `10` | Turbulent Reynolds number used in the wrinkling-efficiency model.               |
| `atf.Cburn`                   |    `2200.0` | Burned-state value used to normalize the flame sensor.                          |
| `atf.Cunburn`                 |     `298.0` | Unburned-state value used to normalize the flame sensor.                        |
| `atf.Cindex`                  | `idx_t::QT` | Primitive-variable index used by the flame sensor, defaulting to temperature.   |

Example input:

```ini
# ATF / TFM model
atf.thickening_factor = 5
atf.laminar_flame_speed = 1.2
atf.laminar_flame_thickness = 1.0e-4
atf.Ret = 10
atf.Cburn = 2200.0
atf.Cunburn = 298.0
```

TFM requires **two options** to pass to `react_sourceLES`  and  `viscous_LES` often to be defined in `user_source_t` (or other structure)

```cpp
public :  
  // ATF options
  bool static constexpr use_ATF = true; // use adaptive thickening factor
  static constexpr int ATF_model = 1;   // 1: Classic  Colin/Charlette , 2: Rathore transformation  
```

These compile-time options control whether the **Artificially Thickened Flame (ATF)** model is applied and which formulation is used to modify the chemical source terms.

```cpp
static constexpr bool use_ATF = true;
```

Enables the ATF combustion model. When enabled, reaction rates are modified to account for flame thickening and unresolved flame wrinkling. If set to `false`, no ATF correction is applied.

```cpp
static constexpr int ATF_model = 1;
```

Selects the ATF formulation:

* **`1` — Classic Colin/Charlette model**\
  Applies both flame thickening and a wrinkling efficiency correction based on sub-grid turbulence, giving a source-term scaling of **E/F** and diffusion **F**
* **`2` — Rathore transformation**\
  Applies only the thickening correction without the wrinkling efficiency model, providing a simpler transformed-source formulation (check [Theory](theory/equations/turbcomb.md#atf))

In practice, these options define how unresolved flame–turbulence interaction is represented in LES combustion simulations.

### PaSR

The `PaSR_t` class reads its parameters from the `pasr` namespace in the input file. These values can be used to override the default Partially Stirred Reactor model settings.

{% hint style="warning" %}
The mixing constant is not used yet (Jan 2026)
{% endhint %}

| Input option           | Default | Description                                                                                                                                                      |
| ---------------------- | ------: | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `pasr.mixing_constant` |   `1.0` | Mixing-model constant controlling the turbulent mixing timescale. Smaller values increase the effective reaction rate; larger values increase mixing limitation. |

Example input:

```ini
# PaSR combustion model
pasr.mixing_constant = 0.7
```

## rhs\_dt

### skew\_t

Example use is

```cpp
typedef rhs_dt<skew_t<methodparm_t, ProbClosures>,  .... > ProbRHS;
```

with

```cpp
struct methodparm_t {

  public:

  static constexpr bool dissipation = true;         
  static constexpr int  order = 4;                  
  static constexpr Real C2skew=0.1,C4skew=0.016;   

};
```

<table><thead><tr><th>Variable Name</th><th>Type</th><th width="143.140625">Value</th><th>Description</th></tr></thead><tbody><tr><td><code>dissipation</code></td><td><code>bool</code></td><td><code>true</code></td><td>Flag indicating whether dissipation is applied (true = yes).</td></tr><tr><td><code>order</code></td><td><code>int</code></td><td><code>4</code></td><td>Order of the numerical scheme used 2/4/6.</td></tr><tr><td><code>C2skew</code></td><td><code>Real</code></td><td><code>0.1</code></td><td>Coefficient for 2nd-order skew-symmetric dissipation.</td></tr><tr><td><code>C4skew</code></td><td><code>Real</code></td><td><code>0.016</code></td><td>Coefficient for 4th-order skew-symmetric dissipation.</td></tr></tbody></table>

### viscous\_t

Example use is

```cpp
typedef rhs_dt<skew_t<... viscous_t<methodparm_t, ProbClosures> ..  > ProbRHS;
```

```cpp
struct methodparm_t {

  public:

  static constexpr int  order = 2;                  // order numerical scheme
  static constexpr bool use_LES= false;             // LES model 
};
```

| Variable Name | Type   | Value   | Description                                 |
| ------------- | ------ | ------- | ------------------------------------------- |
| `order`       | `int`  | `2`     | Order of the numerical scheme used 2/4/6.   |
| `use_LES`     | `bool` | `false` | Flag indicating sgs model used (true =yes). |

If `use_LES = true`, sub-grid viscosty/conductivity/etc. will be added to the molecular viscosity.

## IBM

Example use is

```cpp
typedef eib_t<TypeWall,ibmparm_t,ProbClosures> ProbIB;
```

```cpp
struct ibmparm_t {

  public:

  static constexpr int  interp_order = 1;
  static constexpr int  extrap_order = 1;
  static constexpr Real alpha= 0.6;      
};
```

This controls IBM options, and it is needed to define stencils for extrapolation.\
Order is **interp\_order+1**, in the above example the interpolation order is 2

| Variable Name  | Type   | Value | Description                                          |
| -------------- | ------ | ----- | ---------------------------------------------------- |
| `interp_order` | `int`  | `1`   | Order of the interpolation scheme -1 used in IBM.    |
| `extrap_order` | `int`  | `1`   | Order of the extrapolation scheme -1 used in IBM.    |
| `alpha`        | `Real` | `0.6` | Ratio between first IP and IB points (in mesh units) |

## Default values

The default options are (any struct can be replace by `defaultparm_t`)

```cpp
// default values 
struct defaultparm_t {

public:


static constexpr int order = 2;                   // order of scheme (used in many options)

static constexpr bool dissipation = false;        // use in some schemes to activate dissipation

#ifdef AMREX_USE_GPIBM
static constexpr bool ibm = true;
#else
static constexpr bool ibm = false;
#endif

//  Skew 
static constexpr Real C2skew=0.1,C4skew=0.0016;   // Skew symmetric default

// Transport properties
static constexpr Real conductivity = 0.0262;      // conductivity (for constant value)
static constexpr Real viscosity   = 1.85e-5;      // viscosity    (for constant value)

// ebm/walls options
static constexpr Real Twall = 300;              // wall temperature (for isothermal wall type)
static constexpr bool solve_diffwall = false;   // do not solve viscous fluxes at teh wall (for ebm)

// viscous options 
static constexpr bool use_LES = false;

// LES options
static constexpr Real Pr_o_Prsgs = 1.0;          // Pr/Prsgs
static constexpr Real Scsgs = 0.7;              // sgs Schmidt number
static constexpr Real Cs = 0.1;                 // Smagorinsky constant
static constexpr Real CI = 0.1;                 // Yoshizawa constant
static constexpr bool fixDelta = false;         // (it will use local mesh size as filter width)
static constexpr Real Delta = 0.01;             // Filter width  

};
```

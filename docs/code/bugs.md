---
icon: bug
---

# Bugs

## Common errors

This section outlines common errors, their possible causes, and suggested solutions.

### installation errors

```
examples do not compile
```

**Fix:**

* Check requirements by executing the `checkreq.sh` script in `cerisse/bin`

```
sundials not installed
```

**Fix:**

* Ensure to do `make TPL` the first time you prepare a simulation involving reacting chemistry.

### input file

```
SIGABRT
Domain size not divisible by blocking_factor
```

**Cause:**\
The mesh domain size is not divisible by the specified `blocking_factor`.

**Fix:**

* Ensure that each domain dimension is divisible by the `blocking_factor`.
* If embedded boundaries (EB) are used, ensure the mesh is isotropic, i.e., `dx = dy = dz`.

```
SIGABRT
AMR::checkInput bad_ref_ratios
```

**Cause:**\
Probably a typo in refinement `ref_ratio`, maybe is set to 1?

**Fix:**

* Ensure that ref\_ratio is 2 or 4 (or similar number greater than 1)

```
SIGABRT
amrex::Abort::0::Exiting because either max_step and/or stop_time is less than or equal to 0. !!!
```

**Cause:**\
Probably input file not specified

**Fix:**

* Remember to run `./main3d.gnu.ex input` or similar

### IBM related

```
ERROR via ASSERT
Interpolation point weights do not sum to 1
```

**Cause:**\
Presence of extremely thin surfaces prevents valid interpolation grid points from being placed in solid regions.

**Fix:**

* Refine the mesh, or
* Coarsen the geometry to ensure proper interpolation.

```
ERROR via ASSERT
Grid point on IB surface
```

**Cause:**\
A grid point lies exactly on the immersed boundary (IB) surface, which is ambigous.

**Fix:**

* Slightly shift the domain to prevent grid lines from aligning exactly with the surface.

```
WARNING
Less than 2 interpolation points are fluid points !
```

**Cause:**\
The geometry is too detailed relative to the mesh resolution.

**Fix:**

* Clean or simplify the STL geometry (recommended), or
* Refine the mesh in that region

{% hint style="warning" %}
This is a warning, but the simulation will probably crash in the affected areas.
{% endhint %}

### Chemistry related

While compiling

```
fatal error: mechanism.H
```

**Cause:**

* The selected chemical mechanism in the `GNUMakefile` is not available (e.g., _JL4_ in PelePhysics v23), or
* The mechanism does not exist.

**Fix:**

* Check whether the mechanism is supported in your PelePhysics version (e.g., v25).
* Alternatively, convert your chemistry input from CHEMKIN or YAML format using available tools

## Warnings

Running using **clang** compiler

```
clang++: warning: argument unused during compilation: '-rdynamic' [-Wunused-command-line-argument]
```

This flag (-rdynamic) is used mostly on Linux to export all symbols for use with runtime libraries. On **macOS**, it is unused and ignored by clang++. You can safely ignore this warning.

{% hint style="danger" %}
**DEBUG** mode can generate many warnings, particularly from libraries like CGAL, which are usually harmless. Still, warnings in files like `prob.h` may point to real issues or inefficiencies and should be reviewed. It is good practice to keep the code warning-free.
{% endhint %}

## Known bugs

### Skew

{% hint style="danger" %}
**Issue:**\
6<sup>th</sup> order Skew-symmetric scheme does not work.  The cause is currently unknown, probably a bug in the coefficients.
{% endhint %}

**Workaround:**\
Use a different scheme, it will be fixed soon (2025)

### IBM-Specific

{% hint style="danger" %}
**Issue:**\
Mesh refinement fails near the immersed body at simulation step 0.
{% endhint %}

**Workaround:**\
Apply solid marker-based refinement only after time step 0, or consider using a refinement condition based on geometric criteria — for example, refine regions where `x < 1`, or where the distance to a specific point exceeds a threshold..&#x20;

{% hint style="danger" %}
**Issue:**\
Restarts using checkpoints and IBM must have the same number of levels.
{% endhint %}

**Workaround:**\
This is a bug to be fixed, hopefully soon (2025)

## Features missing

KEEP and central schemes not ready for IBM or EB boundary method

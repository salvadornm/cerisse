# Bugs

## Common errors

This section outlines common errors, their possible causes, and suggested solutions.

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

### IBM related

```
ERROR via ASSERT
Interpolation point weights do not sum to 1
```

**Cause:**\
Presence of extremely thin surfaces prevents valid grid points (GP) from being placed in solid regions.

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
This is a warning, but the simulation will probably crash in theaffected areas.
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
* Alternatively, convert your chemistry input from CHEMKIN or YAML format using available tools&#x20;

## Known bugs

### IBM-Specific

{% hint style="danger" %}
**Issue:**\
Mesh refinement fails near the immersed body at simulation step 0.
{% endhint %}

**Workaround:**\
Apply solid marker-based refinement only after time step 0 using a condition like:

{% hint style="danger" %}
**Issue:**\
Restarts using checkpoints and IBM must have the same number of levels.
{% endhint %}

**Workaround:**\
NRY

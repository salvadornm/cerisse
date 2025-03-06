# Examples

In the `exm/` folder, you will find several tests designed to compare numerical schemes. These tests can also serve as prototypes for setting up more complex cases. While they generally run at a slower pace, they can be used for verifying the accuracy of the schemes and for benchmarking against test cases from the literature. Some of the examples use chemistry or in-line diagnostics (such as probes or statistics), or geometry (IBM or EB).

## One Dimensional Examples

<table><thead><tr><th>Example</th><th width="213">Folder </th><th>Numerics</th><th>Notes</th></tr></thead><tbody><tr><td><a href="onedim.md#sods-shock-tube">Shock-Tube</a></td><td><code>numerics/riemann</code></td><td><a href="theory/equations/numerical-methods.md#riemann-solver-with-muscl">Riemann</a> <a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td></td></tr><tr><td><a href="onedim.md#shu-osher-problem">Shu-Osher</a></td><td><code>numerics/shu</code></td><td><a href="theory/equations/numerical-methods.md#riemann-solver-with-muscl">Riemann</a> <a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td></td></tr><tr><td><em>Reactive Shock</em></td><td><code>reactive_sod</code></td><td><a href="theory/equations/numerical-methods.md#teno">TENO</a></td><td><a href="chemistry.md">Chemistry</a></td></tr><tr><td><a href="onedim.md#constant-volume-reactor">Constant volume reactor</a></td><td><code>auto_ignition</code></td><td>N/A</td><td><a href="chemistry.md">Chemistry</a></td></tr></tbody></table>

## Two Dimensional Examples

<table><thead><tr><th width="276">Example</th><th width="193">Folder</th><th width="108">Numerics</th><th>Notes</th></tr></thead><tbody><tr><td><a href="twodim.md#convection-vortex">Convection vortex</a></td><td><code>numerics/covo</code></td><td>Riemann Skew</td><td></td></tr><tr><td><a href="twodim.md#shock-reflection">Shock Reflection</a></td><td><code>shock_reflect</code></td><td><a href="theory/equations/numerical-methods.md#hllc-riemann-solver">Riemann</a></td><td></td></tr><tr><td><em>Mixing Layer</em></td><td><code>mixing_layer</code></td><td><a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td>Stats</td></tr><tr><td><a href="twodim.md#supersonic-flow-forward-step">Supersonic flow over a forward step</a></td><td><code>forward_step</code></td><td><a href="theory/equations/numerical-methods.md#hllc-riemann-solver">Riemann </a><a href="theory/equations/numerical-methods.md#rusanov-scheme">Rusanov</a><br><a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td><a href="theory/ibmeb.md#embedded-boundaries">EB</a></td></tr><tr><td><a href="twodim.md#periodic-channel-laminar-flow">Periodic channel laminar flow</a></td><td><code>viscwall</code></td><td>Central</td><td>Viscosity</td></tr><tr><td><a href="twodim.md#planar-flame">Planar Flame</a></td><td><code>planar_flame</code></td><td><a href="theory/equations/numerical-methods.md#hllc-riemann-solver">Riemann</a></td><td><a href="chemistry.md">Chemistry</a><br><a href="code/viscous.md">Diffusion</a></td></tr></tbody></table>

## Three Dimensional Examples

| Example                        | Folder         | Notes        | Notes |
| ------------------------------ | -------------- | ------------ | ----- |
| _Taylor-Green Vortex_          | `numerics/tgv` | Skew         | Probe |
| _Turbulent Channel_            | `turbulent`    | Skew         |       |
| _Sphere under supersonic flow_ | `ibm`          | Skew/Rusanov | IBM   |

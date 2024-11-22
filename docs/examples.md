# Examples

In the `exm/` folder, you will find several tests designed to compare numerical schemes. These tests can also serve as prototypes for setting up more complex cases. While they generally run at a slower pace, they can be used for verifying the accuracy of the schemes and for benchmarking against test cases from the literature. Some of the exampels use chemistry or in-line diagnostics (such as probes or statistics), or geometry (IBM/EB)

## One dimensional Examples

<table><thead><tr><th>Example</th><th width="213">Folder</th><th>Numerics</th><th>Notes</th></tr></thead><tbody><tr><td><a href="onedim.md#sods-shock-tube">Shock-Tube</a></td><td><code>numerics/riemann</code></td><td><a href="theory/equations/numerical-methods.md#riemann-solver-with-muscl">Riemann</a>, <a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td></td></tr><tr><td><a href="onedim.md#shu-osher-problem">Shu-Osher</a></td><td><code>numerics/shu</code></td><td><a href="theory/equations/numerical-methods.md#riemann-solver-with-muscl">Riemann</a>, <a href="theory/equations/numerical-methods.md#skew-symmetric">Skew</a></td><td></td></tr><tr><td>Reactive Shock</td><td><code>reactive_sod</code></td><td><a href="theory/equations/numerical-methods.md#teno">TENO</a></td><td>Chemistry</td></tr></tbody></table>

## Two dimensional Examples

| Example           | Folder          | Numerics      | Notes      |
| ----------------- | --------------- | ------------- | ---------- |
| Convection vortex | `numerics/covo` | Riemann, Skew |            |
| Shock Reflection  | `shock_reflect` | Riemann       |            |
| Mixing Layer      | `mixing_layer`  | Skew          | Chem/Stats |
| _Forward Step_    | `forward_step`  | Riemann       | EB         |

## Three dimensional Examples

| Example                        | Folder         | Notes        | Notes |
| ------------------------------ | -------------- | ------------ | ----- |
| _Taylor-Green Vortex_          | `numerics/tgv` | Skew         | Probe |
| _Turbulent Channel_            | `turbulent`    | Skew         |       |
| _Sphere under supersonic flow_ | `ibm`          | Skew/Rusanov | IBM   |

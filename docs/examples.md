# Examples


In the `exm/` folder, you will find several tests designed to compare numerical schemes. These tests can also serve as prototypes for setting up more complex cases. While they generally run at a slower pace, they can be used for verifying the accuracy of the schemes and for benchmarking against test cases from the literature. Some of the exampels use chemistry or in-line diagnostics (such as probes or statistics), or geometry (IBM/EB)


## One dimensional 

| Example                   | Folder              |  Numerics          |  Notes        |
| --------------------------| ------------------- | ------------------ | ------------- |
| Shock-Tube                | `numerics/riemann`  | Riemann, Skew      |               |
| Shu-Osher                 | `numerics/shu`      | Riemann, Skew      |               |
| Reactive Shock            | `reactive_sod`      | TENO               |  Chemistry    |


## Two dimensional 

| Example                   | Folder              |   Numerics         |  Notes        |
| --------------------------| ------------------- | ------------------ | ------------- |
| Convection vortex         | `numerics/covo`     | Riemann, Skew      |               |
| Shock Reflection          | `shock_reflect`     | Riemann            |               |
| Mixing Layer              | `mixing_layer`      | Skew               | Chem/Stats    |
| *Forward Step*            | `forward_step`      | Riemann            | EB            |


## Three dimensional 

| Example                   | Folder              |  Notes             |  Notes        | 
| --------------------------| ------------------- | ------------------ | ------------- |
| *Taylor-Green Vortex*     | `numerics/tgv`      | Skew               | Probe         |
| *Turbulent Channel*       | `turbulent`         | Skew               |               |  
| *Sphere under supersonic* | `ibm`               | Skew/Rusov         | IBM           |  



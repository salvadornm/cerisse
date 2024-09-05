# Test

This set of classic/regression test-cases are located under `tst/`
The code will be tested periodically against these regression tests to validate merges, versions, etc. They are also good examples of basic functionality


| Test                        |  Type          | Description                                                  |
| --------------------------- | ---------------|    --------------------------------------------------------- |
| Test1                       | 1D             |   **Sod's tube**                              |
| Test2                       | 1D             |   **Reactive Sod's**                              |
| Test3                       | 3D             |   **TGV**                              |
| Test4                       | 3D             |   **Hypersonic Sphere**                              |



## Test 1

--8<-- "tst/test1/README.md"

## Test 2

--8<-- "tst/test2/README.md"


# Other examples

These examples are classical tests designed to compare numerical schemes and are located in the ```exm``` folder. While they typically run more slowly, they are valuable for verifying the accuracy of the schemes and for comparison with test cases from the literature.

## Convection of a Vortex

--8<-- "exm/covo/README.md"

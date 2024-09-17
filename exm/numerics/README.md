
These tests (located in ```exm/numerics```) consist of classical tests.
These tests are good to validate the behaviour 
of numerical schemes, regarding shock capturing and dissipation of a vortex, which represents how numerical schemes dissipate turbulence.


The problems are defined in sub-folders


### Convective-Vortex

Same configuration a sin ```exm/covo```, located in 
```exm/numerics/covo```

**Skew-symmetric** 4th order no shock capturing nor damping  N=128 x 128
![covnum](images/num_covo_skew_nodamp.png)


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.5** and **C=0.0016**  N=128 x 128
![covnum2](images/num_covo_skew_damp.png)


**MUSCL** Riemann solver 2th order  N=128 x 128
![covonum3](images/num_covo_muscl.png)

Table with speed

### Riemann Problem


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.5** and **C=0.0016**  N=96
![rienum1](images/num_rie_skew_damp.png)


**Skew-symmetric** 4th order,Shock Capturing and high freq damping with constants **C=0.1** and **C=0.00016**  N=96
![rienum2](images/num_rie_skew_damp2.png)


**MUSCL** Riemann solver 2th order  N=96
![rienum3](images/num_rie_muscl.png)



### Shu-Osher Problem

Thre probelm correspinbds to as M=3 shock moving into a field with small density disturbance, which can be interpreted as a entropy disturbance.
The solutyion is compared with a eference solution with a fine mesh.


Table cost gcc version 11, 1 core, all CFL fix 0.3, times are in secs approximated, all in same CPU

| Run                     | Mesh         | Run Time | 
| ------------------------| -------------|----------| 
| Reference MUSCL         | 8912         |  49.337  | 
| MUSCL                   | 256          |  0.098   | 
| Skew 4 *Cs=0.1 Cd=0.016*  | 256          |  0.08    | 


**MUSCL** Riemann solver 2th order  **N=256**
![shunum3](images/num_shu_muscl.png)

**Skew**  solver 4th order  **N=256**, *Cs=0.1* *Cd=0.016*
![shunum2](images/num_shu_skew.png)

**Skew**  solver 4th order  **N=512**, *Cs=0.1* *Cd=0.016*
![shunum2](images/num_shu_skew2.png)
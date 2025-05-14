# Performance

{% hint style="danger" %}
Work in progress (results from **cerisse1**)
{% endhint %}

## Strong single-node scaling

Strong scaling refers to how the solution time decreases as more core\
are used for a fixed total problem size.

Single node AMD Rome EPYC 7742 up to 128 cores\
Tested with gcc 13.3.0 and iimpi 2024a.

The test case is a TGV\
128x 128 x 1128, Mach=1.25, Re=1600, with fix CFL=0.3 using TENO5

<figure><img src="../.gitbook/assets/MPI.png" alt="" width="375"><figcaption><p>Strong scaling test on a single node up to 128 cores</p></figcaption></figure>

### Performance Table

| # CPU | Runtime Advance \[s] | Speedup | Efficiency |   Min |   Avg |   Max | % Advance | % Load-Imbalance |
| :---: | -------------------: | ------: | ---------: | ----: | ----: | ----: | --------: | ---------------: |
|   1   |             20297.25 |   1.000 |      1.000 | 19310 | 19310 | 19310 |    95.14% |            0.00% |
|   2   |             10571.57 |   1.920 |      0.960 |  9776 |  9779 |  9782 |    92.53% |            0.06% |
|   4   |              5562.19 |   3.649 |      0.912 |  4900 |  4922 |  4947 |    88.94% |            0.95% |
|   8   |              2946.20 |   6.889 |      0.861 |  2468 |  2470 |  2474 |    83.97% |            0.24% |
|   16  |              1521.85 |  13.337 |      0.834 |  1223 |  1231 |  1234 |    81.09% |            0.89% |
|   32  |               792.86 |  25.600 |      0.800 | 629.3 | 638.0 | 641.7 |    80.94% |            1.94% |
|   64  |               416.84 |  48.693 |      0.761 | 316.4 | 328.8 | 334.7 |    80.29% |            5.57% |
|  128  |               216.65 |  93.688 |      0.732 | 165.7 | 167.0 | 175.7 |    81.10% |            5.99% |

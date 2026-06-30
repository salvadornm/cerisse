---
icon: chart-mixed
cover: >-
  https://images.unsplash.com/photo-1488229297570-58520851e868?crop=entropy&cs=srgb&fm=jpg&ixid=M3wxOTcwMjR8MHwxfHNlYXJjaHw3fHxEYXRhfGVufDB8fHx8MTczMzA2ODU4NXww&ixlib=rb-4.0.3&q=85
coverY: 0
---

# Data Analysis

The code includes several capabilities for analyzing simulation data on the fly.\
The main features include:

* **Time statistics**: Calculation of averages and correlations.
* **Probes**: Tools to study the evolution of specific points or regions over time.

## Time statistics

The code can store the evolution of certain properties over time, storing mean and mean square

$$
\left< \phi \right> = \sum \frac{\phi}{N_t} \;\; \mbox{ and } \;\; \left< \phi^2 \right> = \sum \frac{\phi^2}{N_t}
$$

Root mean squre (rms) can be easily obtained a posteriori by:

$$
\phi_{rms} = \sqrt{  \left< \phi \right>^2 - \left< \phi^2 \right>   }
$$

This can be directly extracted in Visit, by loading the `visitexpressions.xml` file in `./Analysis`

### Set-up

To alocate space to store the statistics, in `prob.h` a new index template has to be used, `indicies_stat_t`

```cpp
using ProbClosures = closures_dt<indicies_stat_t, visc_suth_t, cond_suth_t,
                                 calorifically_perfect_gas_t<indicies_t>>;
```

This increases the storage data to allocate space for statistical variables. Two input to the template are required, **record\_velocity** and **record\_PTrho** , which determine if statistics are to be stored for velocities (includes all correlations) and pressure, temperatrue and density

```cpp
static constexpr int record_velocity   = 1;
static constexpr int record_PTrho      = 1;
```

Additionaly, the flag `cns.record_stats = 1` in the input must be set to 1 to enable the recording of statistics. If not, space for the statistics will be allocated, but no data will be written.&#x20;

{% hint style="info" %}
New variables will appear in the plot files. Note that  storing statistical data will increase greatly the size of plotting and checkpoint files.
{% endhint %}

{% hint style="danger" %}
user-specific statistics in progress
{% endhint %}

## Probes

Probes record time-evolution of seelcted quantities ,  they are handled through the **input** file.&#x20;

A probe is defined by a **field name** and a **bounding box** (`box_lo` and `box_hi`). At every sampling interval (`cns.time_probe_int`), Cerisse extracts the requested quantity over that region and stores a single representative value in the output file. By default, this value is the **volume average** over the selected box, making probes useful for monitoring the evolution of global or local flow quantities. Alternatively, the reduction operation can be configured to record the **maximum** (or **minimum**) value within the box, which is particularly useful for tracking shock strength, peak temperature, maximum vorticity, or other localized extrema. Setting `box_lo` equal (or nearly equal) to `box_hi` effectively defines a point probe, allowing time histories to be recorded at a specific location.

```
# DIAGNOSTICS
amr.data_log = time_probe.log
cns.record_probe = 1
cns.time_probe_lev = 0
cns.time_probe_int = 1
cns.time_probes = KineticEnergy Enstrophy Density PressurePoint
KineticEnergy.field_name = kinetic_energy
KineticEnergy.box_lo = 0.0 0.0 0.0
KineticEnergy.box_hi = 6.28318 6.28318 6.28318
Enstrophy.field_name = enstrophy
Enstrophy.box_lo = 0.0 0.0 0.0
Enstrophy.box_hi = 6.28318 6.28318 6.28318
Density.field_name = Density
Density.box_lo = 0.0 0.0 0.0
Density.box_hi = 6.28318 6.28318 6.28318
PressurePoint.field_name = pressure
PressurePoint.box_lo = 3.14159 3.14159 3.14159
PressurePoint.box_hi = 3.19068 3.19068 3.19068
```

In the example above, a file named `time_probe.log` is defined to store data, including four quantities: two derived ones—_energy_ and _enstrophy_—and two direct ones—_density_ and _pressure_. The first two are averaged over a box of size (0, 6.28), while the latter correspond to a much smaller box (a "point")

The `cns.record_probe = 1` is needed to store the data. An example can be found in the [constant volume reactor ](onedim.md#constant-volume-reactor)where temperature is tracked over time. Results can then be plot with conventional python scripts.


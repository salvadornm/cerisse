# Spherical pressure wave (Lodato section 4.4)

Based on G. Lodato, P. Domingo and L. Vervisch, *Three-dimensional boundary
conditions for direct and large-eddy simulation of compressible viscous flows*,
J. Comput. Phys. 227 (2008), 5105-5143, section 4.4, Eqs. (96)-(97), Figs. 12/15.
The global error is Eq. (91). This example prepares comparable small/wider runs;
it is not a claim to reproduce the paper's complete 3D-NSCBC implementation.

## Physical and numerical setup

- Cube centered at the origin, side L = 0.013 m, 64^3 cells.
- p = 101325 [1 + 0.001 exp(-r^2/(2 Rp^2))] Pa, Rp = 0.00065 m.
- Uniform initial T = 300 K; rho = p/(R T), not uniform density or an
  isentropic density perturbation. All three initial velocities are zero.
- Perfect air: gamma = 1.4, molecular weight 28.9 g/mol. Table 2 Sutherland
  viscosity uses mu_ref = 1.827e-5 Pa s, T_ref = 291.15 K, S = 120 K;
  conductivity is mu cp / Pr, with Pr = 0.72.
- Fourth-order skew-symmetric convection without artificial dissipation,
  fourth-order centered diffusion, three-stage third-order RK.
- Type-3 pressure-relaxed outflow on all six faces: p_target = 101325 Pa,
  sigma = 0.28 (T1), Mmax = 0 (rest-state reference), Lchar = domain side.
  `inputs` uses the standard LODI option (`use_transverse = 0`).
- Fixed dt = 1e-7 s in both domains; 840 steps cover t/t_ref = 0.3,
  where the paper uses t_ref = 2.8e-4 s. Reduce dt when refining the mesh.
- `input.wider` doubles each side to 0.026 m and uses 128^3 cells: **same dx**,
  initial pulse radius/amplitude/location, time step, and physical probe boxes.
  Its Lchar is doubled. It needs approximately eight times the interior storage.

The 64^3 resolution is a chosen starting point, not a mesh size attributed to
section 4.4. Perform grid/time convergence before interpreting small differences.

## Build and run

From this directory:

```sh
make -j4
mpirun -np 4 ./main3d.gnu.MPI.ex inputs > screen.small 2>&1
mpirun -np 4 ./main3d.gnu.MPI.ex input.wider > screen.wider 2>&1
```

Outputs are separated: `plot_small/`, `plot_wider/`, `probes_small.csv`, and
`probes_wider.csv`. Start fresh runs for comparison; preserve outputs before
rerunning a different boundary treatment. Plotfiles every 20 steps contain
pressure for the full-domain error; probes are sampled every step with 17
significant digits (`cns.time_probe_precision`). The default elsewhere remains 6.

## Sampling locations

The paper samples exactly at a face center L1, an edge center L2 and a corner L3.
CERISSE's probes sample **cell centers**, so these are nearby interior averages:

| Probe | Small-domain indices (zero based) | Representative location [m] |
|---|---|---|
| L1 | (63, 31:32, 31:32), 4 cells | (0.0063984375, 0, 0) |
| L2 | (63, 63, 31:32), 2 cells | (0.0063984375, 0.0063984375, 0) |
| L3 | (63, 63, 63), 1 cell | (0.0063984375, 0.0063984375, 0.0063984375) |

Wider indices are shifted by 32 in each direction. Keep the same physical
sample volumes in both runs. For a mesh-convergence study, update the boxes
in both input files to track the last interior layer and the center plane.
The probe CSV header uses array indices rather than probe names; numeric columns
are always time, pressure_L1, pressure_L2, pressure_L3.

## Compare with the wider reference

Install numpy and matplotlib for local errors; additionally install yt for
plotfile/global errors. Use the same Python environment for these dependencies.

```sh
python3 compare_pressure.py probes_small.csv probes_wider.csv --output error_T1.png
python3 compare_pressure.py probes_small.csv probes_wider.csv \
  --small-plots plot_small --wider-plots plot_wider --output error_T1_full.png
```

The first command plots the three local errors |p-p_wider|/p_wider (Eq. 97).
The second adds the global relative L2 norm over **only the small domain**
(Eq. 91), producing the four-panel layout of each Fig. 15 column. Error CSVs
are saved beside the figure. The script requires matched times and, for global
errors, matching cell spacing/alignment and level-0 grids; it does not silently
interpolate histories. Local metrics are errors of probe-volume averages.

The vertical dashed line is t/t_ref = 0.2. The paper warns that the doubled-domain
reference becomes affected by its own boundary around this time; later differences
are not a clean estimate of small-domain boundary error. Approximate pulse arrival
times are L/(2c), sqrt(2)L/(2c), sqrt(3)L/(2c), c approximately 347.6 m/s.

## T2 and transverse comparisons

Fig. 15 uses sigma = 0.28 (T1) and sigma = 3.00 (T2). To run T2, override **all six**
face coefficients; run the wider reference with the same overrides and distinct
output names:

```sh
sigma_args=""
for face in xlo xhi ylo yhi zlo zhi; do
  sigma_args="$sigma_args cns.nscbc_${face}_sigma=3.0"
done
mpirun -np 4 ./main3d.gnu.MPI.ex inputs $sigma_args \
  amr.data_log=probes_T2.csv amr.plot_file=plot_T2/plt > screen.T2 2>&1
mpirun -np 4 ./main3d.gnu.MPI.ex input.wider $sigma_args \
  amr.data_log=probes_T2_wider.csv amr.plot_file=plot_T2_wider/plt > screen.T2_wider 2>&1
```

For an exploratory comparison with CERISSE's transverse correction, similarly
set `cns.nscbc_<face>_use_transverse=1` on all six faces in both runs, using separate
output names. The supplied beta = 5e-5 is a **fixed acoustic-Mach-scale proxy**,
not a dynamically measured maximum Mach number. The paper uses beta = Ma_max for
T1/T2; CERISSE currently accepts a fixed scalar. Its T3 experiment instead uses
beta = 0.5 and sigma = 0.28 (Fig. 16).

Known solver differences matter especially for this test:

1. `src/set/nscbc.h` disables transverse terms in multi-boundary ghost cells;
   the paper's coupled edge/corner conditions are not implemented.
2. The current outflow implementation does not explicitly zero incoming convected
   waves during local backflow, as the paper describes on p. 5128.
3. CERISSE advances cell-centered ghost states, while the paper uses auxiliary
   boundary-face nodes. Our probes are near-boundary averages, not face values.

Consequently, use this case to quantify the present solver and future boundary
improvements. Do not label its transverse curves as an exact reproduction of
Lodato's 3D-NSCBC results.

## Verification performed

- Release build and two-rank MPI startup/advance/shutdown checks on the supplied
  64^3 and 128^3 configurations (two steps each, including plotfiles).
- Checked Eq. (96), uniform 300 K, zero initial momentum, and matching pressure
  fields in the shared region using the generated plotfiles.
- Checked matched probe volumes, error normalization, rejection of mismatched
  time series, and the local/global comparison script on actual output.
- A separate 16^3 smoke test with adjusted probe boxes ran all 840 steps through
  t/t_ref = 0.3 with finite positive probe pressures. This is a stability smoke
  test, not a resolution study or validation of Fig. 15's error magnitudes.

The full-duration 64^3/128^3 comparison and T2/transverse sweeps remain to be run.

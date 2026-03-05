# CODEBUDDY.md This file provides guidance to CodeBuddy when working with code in this repository.

## Common Development Commands

### Running the Main Simulation
- `main_GDW` - Execute the primary GDW aerodynamic simulation with BEM initialization. Uses built-in parameter set in the file.
- `results = BEM(params)` - Run the unified BEM+GDW solver with custom parameters. Pass a `params` struct to configure simulation settings.

### Testing Individual Components
- `test_bem_initialize` - Test BEM steady-state initialization convergence.
- `test_gdw_rhs` - Validate GDW right-hand side computation against analytical cases.
- `test_aero_interp` - Check airfoil interpolation accuracy across thickness range.

### Data Processing
- `resolve_blade_geometry(params)` - Load blade geometry from beamModel.mat or override with params.
- `plot_aero_comparison(results)` - Generate performance comparison plots between GDW and steady BEM.

## Code Architecture

This is a MATLAB-based wind turbine aerodynamics simulator implementing Blade Element Momentum (BEM) theory for initialization and Generalized Dynamic Wake (GDW) theory for time-domain simulation. The architecture follows a clear separation of concerns: input resolution, steady-state initialization, unsteady wake evolution, and post-processing.

### Core Workflow
1. **Input Resolution** (`resolve_blade_geometry.m`) loads blade geometry from `beamModel.mat` containing distance, chord, twist, and Thickness arrays. Falls back to manual parameter override via `params` struct. Also loads airfoil data from `cl_cd_cm.mat` through `aeroInterp`.

2. **BEM Initialization** (`bem_initialize.m`) computes steady-state axial/tangential induction factors using iterative BEM with Prandtl tip/hub losses and Glauert correction for high loading. Provides physically realistic initial conditions for GDW.

3. **GDW Time Marching** (`main_GDW.m` or `BEM.m`) evolves the wake using state-space form:
   - Pressure coefficients `tau` computed from blade loading via `compute_tau`
   - State derivatives evaluated by `gdw_rhs` implementing $[M]\dot{\alpha} + [L]^{-1}[V]\{\alpha\} = \frac{1}{2}\{\tau\}$
   - Time integration with ABM4 predictor-corrector via `abm4_update`

4. **Aerodynamic Forces** computed each time step using induced velocities from GDW state, with blade geometry, local flow angles, and airfoil coefficients via `aeroInterp`.

5. **Output Processing** (`plot_aero_comparison.m`) provides dynamic monitoring during simulation and final comparison plots between GDW results and BEM steady baseline.

### Key Data Structures
- `params` struct: Configuration container for rotor parameters, solver settings, and options
- `results` struct: Contains time-series of sectional forces (L,D,M), total forces, and meta information
- `mode` struct: GDW mode definitions (harmonics, shape functions, matrix precomputations)
- `beamModel` struct: Blade geometry with distance, chord, twist, Thickness fields

### Critical Dependencies
- `cl_cd_cm.mat`: Airfoil coefficient tables for thickness-based interpolation
- `beamModel.mat`: Blade structural model with geometry definition
- All utility functions in `./utils/` directory must be on MATLAB path

### Coordinate System Convention
- x-axis: Rotor axis (positive downstream)
- z-axis: Blade spanwise direction 
- y-axis: Blade chordwise direction
- GDW uses ellipsoidal coordinates internally but transforms to this Cartesian system for force calculations

## Important Implementation Notes

### Solver Configuration
- Use `params.N_harmonics` to control GDW fidelity (typical range 0-4)
- Set `params.useParallel = true` to enable blade-loop parallelization
- Enable `params.useSparse = true` for memory efficiency with large numbers of elements
- Monitor convergence with `params.monitor = true` for real-time force plots

### Inflow Angle Parameters

#### Pitch Angle (`params.pitch`)
- **Parameter**: `params.pitch` — collective blade pitch angle in **degrees**, default `0`
- **Sign convention**: Positive pitch rotates the blade leading edge towards feathering (rotation about negative z-axis / blade spanwise axis)
- **Effect**: Applied as an additive offset to the geometric twist angle when computing angle of attack: `Alpha = Phi - deg2rad(theta + pitch_deg)`
- **Output**: `results.meta.pitch_deg`, `results.meta.pitch_rad`

#### Yaw Angle (`params.gamma`)
- **Parameter**: `params.gamma` — rotor yaw angle in **degrees**, default `0`
- **Sign convention**: Positive yaw means the rotor axis rotates clockwise (viewed from above) away from the wind direction, i.e. wind approaches the rotor plane from the left side
- **Effect**: Produces a lateral in-plane velocity component `V0_lat = V0 * cos(yita) * sin(gamma)` that varies the effective tangential velocity with azimuth angle
- **Output**: `results.meta.gamma_deg`, `results.meta.gamma_rad`

#### Tilt Angle (`params.yita`)
- **Parameter**: `params.yita` — rotor tilt (shaft uptilt) angle in **degrees**, default `0`
- **Sign convention**: Positive tilt means the rotor axis tilts upward (nose up), so wind approaches the rotor plane from below
- **Effect**: Produces a constant in-plane vertical velocity component `V0_elev = V0 * sin(yita)` added to the tangential velocity at every blade element
- **Output**: `results.meta.yita_deg`, `results.meta.yita_rad`

#### Freestream Velocity Decomposition
Given `V0` along the global x-axis (aligned with undisturbed wind), the three parameters together decompose into:
```
V0_ax   = V0 · cos(yita) · cos(gamma)     % axial (perpendicular to rotor plane)
V0_lat  = V0 · cos(yita) · sin(gamma)     % lateral in-plane (yaw-induced)
V0_elev = V0 · sin(yita)                  % vertical in-plane (tilt-induced)
```
The resultant in-plane speed `V0_inplane = sqrt(V0_lat² + V0_elev²)` feeds the GDW advance ratio `mu`.

#### Skewed Wake Correction (AeroDyn Theory Manual Eq. 17–19)
When yaw or tilt is non-zero, an equivalent skew angle is computed and applied:
1. `gamma_eff = atan2(V0_inplane, V0_ax)` — effective inflow skew angle
2. `chi_eff = (0.6 · a_avg + 1) · gamma_eff` — wake skew angle (Eq. 19, Burton et al. 2001)
3. Per-element azimuthal induction correction (Eq. 17, Pitt & Peters 1981):
   `a_skew = a_BEM · [1 + (15π/32) · tan(chi_eff/2) · (r/R) · cos(ψ)]`
4. GDW matrices (`precompute_gdw_matrices`) and `lambda_f` (Eq. 89) both use `chi_eff`
- **Output**: `results.meta.gamma_eff_rad`, `results.meta.chi_eff`

#### Hub Offset Length (`params.hub_offset`)
- **Parameter**: `params.hub_offset` — extra distance `L` (metres) from rotation centre to blade root, default `0`
- **Sign convention**: Non-negative value; increases the effective rotor diameter to `2·(R + L)`
- **Effect**: Each blade element's actual rotation radius becomes `r_rot = r + L`, so the base tangential velocity is `Omega · (r + L)`. The additional `Omega·L` term uniformly increases the tangential velocity at every element, raising lift and torque along the full span
- **Scope**: Applied to both BEM initialization and GDW time-marching tangential velocity; GDW inflow states and `rhat = r/R` remain unchanged
- **Output**: `results.meta.hub_L`

#### Cone Angle (`params.cone`)
- **Parameter**: `params.cone` — blade coning angle in **degrees**, default `0`
- **Sign convention**: Positive cone angle means blades deflect downwind (towards the lee side), forming a downwind cone
- **Effect**: Decomposes the axial freestream `V0_ax` (after yaw/tilt correction) into two components in the blade's local frame:
  - Normal to the local rotation plane: `V0_ax_cone = V0_ax · cos(cone)` — effective axial inflow driving induction
  - Along the blade spanwise direction (tangential): `V0_ax_tan = −V0_ax · sin(cone)` — positive cone reduces effective tangential speed
- **GDW constraint**: The GDW inflow states and wake matrices assume a flat rotor disc and are **not** modified by cone angle; only the velocity components at each blade element are corrected
- **Output**: `results.meta.cone_deg`, `results.meta.cone_rad`

#### Combined Velocity Assembly (per blade element)
```
Uax = V0_ax·cos(cone) − u_ind_GDW
Ut  = Omega·(r+L)·(1+a'_BEM) + V0_lat·cos(ψ) + V0_elev − V0_ax·sin(cone)
```
> **Note**: The axial induction is supplied entirely by the GDW induced velocity `u_ind_GDW`.
> The BEM axial induction factor `a_BEM` is **not** subtracted separately — doing so would
> double-count the axial induction, since the GDW state equations already account for it
> implicitly through the iterative wake evolution.
> The tangential induction `a'_BEM` from BEM is still used for `Ut` because GDW does not
> solve for tangential induction.

### Initial Conditions
GDW requires BEM initialization for stability. The `bem_initialize` function provides induction factors that serve as physically meaningful starting conditions for the wake state evolution.

### Time Integration
The ABM4 implementation includes startup phase with explicit Euler steps for the first 3 iterations before switching to predictor-corrector mode, ensuring numerical stability.

### GDW Inflow Gain Matrix — Cross-Harmonic Coupling (Ltilde)

> **[BUG FIX — Core Algorithm]** Previous versions of `precompute_gdw_matrices.m` incorrectly assumed that the row harmonic index `r` (induced-velocity harmonics) and the column harmonic index `m` (pressure harmonics) must be equal, resulting in a block-diagonal `Ltilde` that was always symmetric. This violates Eq. 74–76 of the AeroDyn Theory Manual, where `r` and `m` are independent.

#### Correct Implementation (current)

`precompute_gdw_matrices` now builds two **global** matrices:

| Matrix | Size | Governing equation |
|--------|------|--------------------|
| `LtildeC` | `N_c × N_c` | Eq. 74 (`r=0`) and Eq. 75 (`r>0`, cosine) |
| `LtildeS` | `N_s × N_s` | Eq. 76 (`r>0`, sine) |

where `N_c = Σ S_h` (h=0…N), `N_s = Σ S_h` (h=1…N).

Row index `(r, j)` and column index `(m, n)` are **fully independent**, so every pair `(r, m)` can contribute a non-zero entry:

```
Eq.74  r=0:  Ltilde_{jn}^{0m,c}  = X^m · Γ_{jn}^{0m}          (all m)
Eq.75  r>0:  Ltilde_{jn}^{rm,c}  = (X^|m-r| + (−1)^l · X^|m+r|) · Γ_{jn}^{rm}
Eq.76  r>0:  Ltilde_{jn}^{rm,s}  = (X^|m-r| − (−1)^l · X^|m+r|) · Γ_{jn}^{rm}
             l = min(r, m),  X^0 ≡ 1  (manual note after Eq.78)
```

#### Impact of the fix

- `LtildeC` and `LtildeS` are in general **non-symmetric** (as physically expected)
- All cross-harmonic coupling paths are now included; previously they were entirely absent
- The GDW state equation is now solved as a **single global system** rather than `N+1` independent per-harmonic sub-problems
- `main_GDW.m` now uses global state vectors `alphaVec` / `betaVec` and the helper arrays `offC` / `offS` to index into them

#### API change

`precompute_gdw_matrices` signature changed from:
```matlab
[Mdiag, LtildeC, LtildeS] = precompute_gdw_matrices(mode, chi)
% returned cell arrays: Mdiag{h+1}, LtildeC{h+1}, LtildeS{h+1}
```
to:
```matlab
[Mdiag_c, Mdiag_s, LtildeC, LtildeS, offC, offS] = precompute_gdw_matrices(mode, chi)
% returns full matrices + offset arrays
```

### Performance Optimization
- Precomputes GDW matrices (`precompute_gdw_matrices.m`) to avoid repeated calculations
- Uses active index filtering to skip blade root sections below `RootOffset * R`
- Vectorized airfoil interpolation and force calculations where possible

### Output Interpretation
- `results.total.Lift/Drag/Moment`: Time series of overall rotor forces (all blades summed)
- `results.total.Thrust/Torque`: Time series of rotor thrust and shaft torque
- `results.bemTotals`: Steady BEM reference values for comparison
- `results.section.*`: Detailed sectional data `[Nstations × B × Nt]` for each blade element over time
- `results.blade.*`: Per-blade span-integrated loads `[B × Nt]` — each row is one blade, NOT averaged across blades
  - `results.blade.L/D/M`: Lift, drag, and moment per blade per time step (N, N, N·m)
  - `results.blade.Thrust/Torque`: Axial thrust and shaft torque per blade (N, N·m)
  - `results.blade.psi`: Azimuth angle of each blade at each time step (rad)
- `results.meta`: Simulation metadata including all inflow angle parameters

### Visualization (`plot_aero_comparison`)
Generates two figures:

**Figure 1 — Total rotor loads vs. time** (unchanged):
- Subplot 1/2/3: total lift / drag / moment time series vs. BEM steady baseline

**Figure 2 — Per-blade dynamic loads** (2×3 layout):
| | Lift | Drag | Torque |
|---|---|---|---|
| **vs. Azimuth angle** | each blade's lift vs. ψ (°) | each blade's drag vs. ψ (°) | each blade's torque vs. ψ (°) |
| **vs. Time** | each blade's lift vs. t (s) | each blade's drag vs. t (s) | each blade's torque vs. t (s) |

Azimuth subplots use data from the last full rotation period; blades are shown with separate colored lines (no averaging). Under yaw/tilt conditions, the three blade curves will be offset in phase, clearly revealing azimuth-dependent asymmetric loading.

### Dynamic Monitor (`params.monitor = true`)
Real-time 2×3 monitoring window during time integration:
- **Top-left (2 columns)**: Each blade's x-direction force (thrust) — one colored line per blade
- **Top-right**: Total x-direction force (thrust)
- **Bottom row**: Total thrust | Total lift | Total moment
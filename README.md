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

### Pitch Angle Parameter (`params.pitch`)
- **Parameter**: `params.pitch` — collective blade pitch angle in **degrees**, default `0`
- **Sign convention**: Positive pitch rotates the blade leading edge towards feathering (rotation about negative z-axis / blade spanwise axis)
- **Coordinate transform**: The freestream velocity `V0` (along rotor x-axis) is decomposed via a rotation matrix into:
  - Axial component: `V0_ax = V0 * cos(pitch_rad)` — perpendicular to rotor plane, drives BEM induction and GDW wake
  - In-plane component: `V0_ip = -V0 * sin(pitch_rad)` — tangential component added to the local blade tangential velocity
- **Effect on aerodynamics**:
  - Axial inflow angle `Phi = atan2(V0_ax - u_ind, Omega*r*(1+a') + V0_ip)`
  - Angle of attack `Alpha = Phi - theta` (where `theta` is the geometric twist)
  - Increasing `pitch` (positive) reduces axial component, decreases loading
- **Usage example**:
  ```matlab
  params.pitch = 5;   % 5-degree pitch towards feathering
  main_GDW;
  ```
- **Output**: `results.meta.pitch_deg` and `results.meta.pitch_rad` record the applied pitch angle

### Initial Conditions
GDW requires BEM initialization for stability. The `bem_initialize` function provides induction factors that serve as physically meaningful starting conditions for the wake state evolution.

### Time Integration
The ABM4 implementation includes startup phase with explicit Euler steps for the first 3 iterations before switching to predictor-corrector mode, ensuring numerical stability.

### Performance Optimization
- Precomputes GDW matrices (`precompute_gdw_matrices.m`) to avoid repeated calculations
- Uses active index filtering to skip blade root sections below `RootOffset * R`
- Vectorized airfoil interpolation and force calculations where possible

### Output Interpretation
- `results.total.Lift/Drag/Moment`: Time series of overall rotor forces
- `results.bemTotals`: Steady BEM reference values for comparison
- `results.section.*`: Detailed sectional data for each blade element over time
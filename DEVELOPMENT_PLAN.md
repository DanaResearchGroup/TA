# Strategic Development Plan

Development roadmap for the TA (Thermal Analysis) simulator, organized into phases from the current gas-phase prototype toward a production-ready polymer thermal analysis tool.

---

## Phase 1: Multi-Stage Temperature Programs

**Goal**: Support arbitrary temperature profiles (ramp-hold-ramp, cooling segments, isothermal holds).

**Why**: Real TGA/DTA experiments rarely use a single linear ramp. Standard methods (e.g., ASTM E1131) specify multi-segment programs, and kinetic analysis techniques like ICTAC require multiple heating rates.

**Work items**:

- Generalize the `Func1` wall heat-flux to accept a list of `(rate_C_per_min, T_target_C)` segments, where `rate=0` means isothermal hold.
- Build a `TemperatureProgram` helper class that computes `T_furnace(t)` and the corresponding `Q(t)` for the wall, handling segment boundaries and sign changes (cooling).
- Add a `cooling_rate` path: when `rate < 0`, the furnace cools and the reactor follows with inverted thermal lag.
- Update `SimulationResult` to store the program metadata so signals and plots can annotate segment boundaries.

**Validation**: Test that a ramp-hold-ramp program produces a flat reactor temperature during the hold (within thermal-lag tolerance), and that cooling segments invert the DTA sign.

---

## Phase 2: RMG Polymer Mechanism Integration

**Goal**: Replace `gri30.yaml` with RMG-generated polymer pyrolysis mechanisms as the default use case.

**Why**: This is the project's core motivation. RMG outputs pseudo-homogeneous Cantera phases where polymer chains, radical intermediates, and volatile products coexist in a single `Solution`. The current phase-filtering architecture was designed for exactly this.

**Work items**:

- Add a mechanism loader/validator that reads an RMG-generated YAML, extracts species metadata (molecular weight, phase label if present), and auto-classifies species into condensed vs. volatile using MW thresholds or RMG labels.
- Handle large mechanisms (hundreds of species, thousands of reactions) without excessive memory: store only condensed + top-N volatile species in `SimulationResult`, with an option to record all.
- Add a `MechanismInfo` dataclass exposing species count, reaction count, temperature validity range, and element composition for logging and validation.
- Test with at least one real RMG polystyrene or polyethylene mechanism to validate that TGA curves show physically reasonable mass-loss onset temperatures and residue fractions.

**Validation**: Compare simulated TGA onset temperature and DTG peak temperature against published experimental data for a well-characterized polymer (e.g., polystyrene in N2).

---

## Phase 3: Cantera Interface / Multiphase Reactor

**Goal**: Move beyond the single-phase approximation by using Cantera's `Interface` and bulk-phase objects for true condensed-phase kinetics.

**Why**: The current approach treats all species as gas-phase, which misrepresents thermodynamic properties (density, heat capacity, vapor pressure) of polymer chains. A proper multiphase model would capture evaporation kinetics and condensed-phase heat capacity correctly.

**Work items**:

- Implement a `MultiphaseTASimulator` that creates a Cantera bulk `SolutionArray` (condensed) coupled to a gas-phase `Solution` via an `Interface`.
- Model mass transfer: volatile products generated in the condensed phase evaporate into the gas phase based on vapor-liquid equilibrium or a rate expression.
- Track condensed-phase mass directly (no phase-filtering needed) — TGA becomes a true mass measurement.
- Keep the existing `TASimulator` as the lightweight/fast option; `MultiphaseTASimulator` is the physics-accurate path.

**Validation**: For a species with known vapor pressure (e.g., naphthalene), verify that the evaporation rate and TGA curve match analytical or published experimental results.

---

## Phase 4: DSC Signal and Calibration

**Goal**: Add Differential Scanning Calorimetry (DSC) signal computation and instrument calibration framework.

**Why**: DSC is the modern replacement for DTA and is more commonly reported in literature. Adding proper calibration infrastructure makes simulated signals directly comparable to experimental data.

**Work items**:

- Implement `compute_dsc()` in `signals.py`: heat flow rate [mW] rather than temperature difference, using the wall heat flux and a proper baseline subtraction (empty-crucible reference run).
- Add a `CalibrationProfile` class that stores instrument-specific parameters (thermocouple sensitivity, heat-flow calibration factor, temperature offset) loaded from a JSON/YAML file.
- Support baseline subtraction: run a blank (empty crucible) simulation, store it, and subtract from sample runs to isolate reaction-related heat flow.
- Add enthalpy-of-transition integration: integrate the DSC peak area to compute reaction enthalpy [J/g].

**Validation**: Simulate melting of indium (known enthalpy 28.6 J/g) and verify that peak integration recovers the correct value within calibration tolerance.

---

## Phase 5: Parameter Estimation and Kinetic Analysis

**Goal**: Enable extraction of kinetic parameters (activation energy, pre-exponential factor) from simulated or experimental TGA/DTG data.

**Why**: The primary use of TGA in polymer science is kinetic analysis. Model-free methods (Kissinger, Ozawa-Flynn-Wall, Friedman) require multiple heating rates and are standard practice.

**Work items**:

- Add a `MultiRateExperiment` runner that executes simulations at several heating rates (e.g., 5, 10, 20, 40 C/min) and collects the results.
- Implement model-free kinetic methods:
  - **Kissinger**: activation energy from DTG peak temperature vs. heating rate.
  - **Ozawa-Flynn-Wall (OFW)**: isoconversional activation energy from TGA at multiple rates.
  - **Friedman**: differential isoconversional method.
- Return kinetic parameters as a `KineticResult` dataclass with confidence intervals.
- Add an Arrhenius plot generator to `plotter.py`.

**Validation**: Generate synthetic TGA data from a known single-step Arrhenius reaction, recover the input Ea and A within 5% using each method.

---

## Phase 6: Experimental Data Import and Comparison

**Goal**: Load real instrument data and overlay it with simulations for model validation.

**Why**: The simulator is only useful if its predictions can be compared against experiments. TA instrument vendors (TA Instruments, Netzsch, Mettler-Toledo) export data in vendor-specific formats.

**Work items**:

- Build parsers for common export formats: CSV (generic), NGBL (Netzsch), and TA Instruments `.txt` exports.
- Define a `ExperimentalData` dataclass mirroring `SimulationResult` fields (temperature, TGA, DTG, DTA/DSC, MS channels).
- Add `plot_comparison()` to `plotter.py`: overlay simulated and experimental curves on the same axes with residual subplot.
- Compute goodness-of-fit metrics: RMSE, R^2, and onset/peak temperature error for automated comparison.

**Validation**: Load a published polystyrene TGA dataset, overlay with a simulation using an RMG mechanism, and report quantitative agreement metrics.

---

## Phase 7: Performance and Scalability

**Goal**: Handle large RMG mechanisms (500+ species) and long multi-segment experiments efficiently.

**Work items**:

- Profile `TASimulator.run()` to identify bottlenecks (Cantera advance vs. array recording vs. object creation).
- Add selective species recording: only store mass/mole fractions for species of interest (condensed + user-specified volatiles), not all 500+.
- Implement checkpointing: save intermediate `SimulationResult` to HDF5 at configurable intervals so long runs can be resumed.
- Add progress reporting via a callback or `tqdm` integration for interactive use.
- Consider Cantera's `SolutionArray` for batch storage instead of per-species dictionaries.

**Validation**: Benchmark a 500-species mechanism at 10,000 time steps; target < 60 s wall-clock time and < 500 MB memory.

---

## Infrastructure (Cross-Cutting)

These items support all phases and should be addressed incrementally.

### CI/CD
- Add GitHub Actions workflow: run `pytest` on push/PR against Python 3.12+ with Cantera from conda-forge.
- Add a slow-test marker (`@pytest.mark.slow`) for full-mechanism runs; run only on `main` merges.

### Documentation
- Set up Sphinx with autodoc to generate API reference from the existing docstrings.
- Add a tutorial notebook (`examples/tutorial.ipynb`) walking through simulation, signal computation, and plotting.

### Packaging
- Publish to PyPI (pure-Python portion) with Cantera as an optional dependency documented in install instructions.
- Add `py.typed` marker and verify `mypy --strict` passes for the `ta` package.

### Testing
- Add integration tests that run a full simulate-signal-plot pipeline and assert output shapes and value ranges.
- Add property-based tests (Hypothesis) for signal functions: e.g., TGA is always in [0, 100], DTG of monotonic TGA is always non-positive.
- Add a `conftest.py` with shared fixtures for common simulator configurations.

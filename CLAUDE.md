# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Thermal analysis (TGA/DTA/DTG/MS) simulator built on [Cantera](https://cantera.org). Models a two-crucible furnace experiment with a linear temperature ramp (or a sequence of ramps) using operator-split gas-phase kinetics and polymer melt-phase kinetics.

## Commands

```bash
# Setup (conda required; Cantera is only available via conda-forge)
conda env create -f environment.yml
conda activate ta_env
pip install -e .

# Run all tests
pytest

# Run a single test file or test
pytest tests/test_simulator.py
pytest tests/test_simulator.py::TestThermalLag::test_reactor_lags_furnace_during_heating

# Run tests with coverage
pytest --cov=ta --cov-report=term-missing

# Run tests in parallel
pytest -n auto

# Run the example benchmark
python examples/run_benchmark.py
```

## Architecture

The package has three modules forming a pipeline: **simulate -> compute signals -> plot**.

- **`ta/temperature_program.py`** — `TemperatureProgram` and `TemperatureSegment` define multi-segment furnace profiles (ramps, isothermal holds, cooling). The program precomputes segment boundaries and provides fast `T_furnace_K(t)` / `T_furnace_C(t)` lookup via `bisect`.

- **`ta/simulator.py`** — `TASimulator` sets up a Cantera `ConstPressureReactor` (sample pan) connected to a `Reservoir` (furnace) via a `Wall` with heat-transfer coefficient and area. The furnace temperature profile (single ramp or multi-segment `TemperatureProgram`) is implemented as a `Func1` on the wall, so the energy equation is fully coupled (no operator splitting). Returns a `SimulationResult` dataclass with furnace/reactor temperature, mass fractions, mole fractions, wall heat flux, and density.

- **`ta/signals.py`** — Pure NumPy post-processing functions that derive instrument signals from `SimulationResult`:
  - `compute_tga` — sum of condensed-species mass fractions × 100
  - `compute_dtg` — time derivative of TGA (units: %/min)
  - `compute_dta` — wall heat flux normalised by current condensed mass [W/kg]
  - `compute_ms` — mole fractions of volatile (target gas) species, with optional `top_n` ranking

- **`ta/plotter.py`** — `plot_ta_figure` produces a two-panel matplotlib figure: Panel A has DTA/TGA/DTG on shared and twin axes; Panel B shows MS traces on a log scale, shared X axis for both panel of the experiment's temperature.

## Key Concepts

- **Condensed vs. volatile species**: Species in `condensed_species` stay in the crucible (tracked for TGA mass loss); all others are volatile (tracked for MS).
- **Wall-based heat transfer**: The furnace `Reservoir` stays at the initial temperature; a `Func1` heat-flux term on the `Wall` adds `U·(T_program(t) - T_initial)` so the effective furnace temperature follows the program. This avoids the Cantera limitation that `Reservoir` state cannot be updated mid-integration. Works for ramps, holds, and cooling segments alike.
- **Thermal lag**: The reactor temperature lags behind the furnace due to the finite wall `U·A`. Steady-state lag = `m·cp·rate / (U·A)`. Exo/endothermic reactions modulate this lag, producing the DTA signal.

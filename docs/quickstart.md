# Quick Start

## Single-ramp simulation

The simplest experiment is a linear temperature ramp.  This example uses
Cantera's built-in GRI-Mech 3.0 mechanism with ethane as a stand-in for
a condensed-phase polymer.

```python
from ta import (
    TASimulator,
    compute_tga, compute_dtg, compute_dta, compute_ms,
)
from ta.plotter import plot_ta_figure

# 1. Configure and run the simulation
sim = TASimulator(
    mechanism="gri30.yaml",
    sample_composition="C2H6:1",
    condensed_species=["C2H6"],
    T_initial_C=50.0,
    T_final_C=600.0,
    heating_rate_C_per_min=5.0,
)
result = sim.run()

# 2. Compute signals
tga = compute_tga(result)
dtg = compute_dtg(tga, result.time_s)
dta = compute_dta(result)
ms  = compute_ms(result, top_n=5)

# 3. Plot
fig, axes = plot_ta_figure(result.furnace_temperature_C, tga, dtg, dta, ms)
fig.savefig("ta_result.png", dpi=150)
```

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `mechanism` | `"gri30.yaml"` | Cantera YAML mechanism file |
| `sample_composition` | `"CH4:1"` | Initial composition (string or dict) |
| `sample_mass_mg` | `10.0` | Sample mass [mg] -- sets reactor volume |
| `bath_gas` | `"AR"` | Inert carrier gas species |
| `T_initial_C` | `50.0` | Start temperature [deg C] |
| `T_final_C` | `600.0` | End temperature [deg C] |
| `heating_rate_C_per_min` | `5.0` | Heating rate [deg C / min] |
| `condensed_species` | auto | Species tracked for TGA mass loss |
| `dt` | `1.0` | Recording interval [s] |
| `wall_area` | `1e-4` | Furnace-sample wall area [m^2] |
| `htc` | `100.0` | Wall heat-transfer coefficient [W/m^2/K] |

## Understanding the output

`result` is a {class}`~ta.simulator.SimulationResult` dataclass containing:

- **`furnace_temperature_C`** / **`reactor_temperature_C`** -- the prescribed
  furnace temperature and the actual sample temperature (which lags due to
  finite wall heat transfer).
- **`mass_fractions`** / **`mole_fractions`** -- per-species dictionaries of
  time-series arrays.
- **`wall_heat_flux`** -- heat transfer rate [W] through the wall (positive
  means heat flowing into the sample).

## Signal definitions

| Signal | Formula | Units |
|---|---|---|
| **TGA** | $\sum Y_\text{condensed}(t) \times 100$ | % |
| **DTG** | $\mathrm{d}(\text{TGA}) / \mathrm{d}t$ | %/min |
| **DTA** | $\dot{Q}_\text{wall} / m_\text{condensed}$ | W/kg |
| **MS**  | Mole fractions of volatile species | -- |

## Using a dict for composition

```python
sim = TASimulator(
    mechanism="gri30.yaml",
    sample_composition={"C2H6": 0.7, "CH4": 0.3},
    condensed_species=["C2H6", "CH4"],
)
```

# Temperature Programs

Real thermal analysis experiments use multi-segment temperature profiles:
ramp to a target, hold isothermally, ramp again, or cool back down.
The {class}`~ta.temperature_program.TemperatureProgram` class lets you
define arbitrary piecewise-linear furnace programs.

## Concepts

A temperature program is a sequence of **segments**, each described by a
{class}`~ta.temperature_program.TemperatureSegment`:

| Segment type | `rate_C_per_min` | `T_target_C` | `hold_min` |
|---|---|---|---|
| Heating ramp | `> 0` | target temperature | ignored |
| Cooling ramp | `< 0` | target temperature | ignored |
| Isothermal hold | `0` | current temperature | hold duration |

The program starts at `T_initial_C` and executes each segment in order.
Ramp segments heat or cool at a constant rate until the target is reached;
hold segments keep the temperature constant for the specified duration.

## Creating a program

### Single ramp (legacy shortcut)

If you only need a single linear ramp, you can either use the legacy
`TASimulator` parameters or the factory method:

```python
from ta import TemperatureProgram

prog = TemperatureProgram.single_ramp(
    T_initial_C=50.0,
    T_final_C=600.0,
    rate_C_per_min=10.0,
)
```

### Ramp -- hold -- ramp

A common pattern in TGA methods (e.g., ASTM E1131): ramp to an
intermediate temperature, hold to allow equilibration, then ramp to the
final temperature.

```python
from ta import TemperatureProgram, TemperatureSegment

prog = TemperatureProgram(50.0, [
    TemperatureSegment(rate_C_per_min=10.0, T_target_C=200.0),
    TemperatureSegment(rate_C_per_min=0.0,  T_target_C=200.0, hold_min=30.0),
    TemperatureSegment(rate_C_per_min=5.0,  T_target_C=600.0),
])

print(f"Total time: {prog.total_time_s / 60:.0f} min")  # 15 + 30 + 80 = 125 min
```

### Ramp -- hold -- ramp -- cool

Add a cooling segment with a negative rate:

```python
prog = TemperatureProgram(50.0, [
    TemperatureSegment(10.0, 300.0),                  # heat to 300
    TemperatureSegment(0.0,  300.0, hold_min=5.0),    # hold 5 min
    TemperatureSegment(5.0,  500.0),                   # heat to 500
    TemperatureSegment(-10.0, 200.0),                  # cool to 200
])
```

## Inspecting a program

### Segment boundaries

The {attr}`~ta.temperature_program.TemperatureProgram.segment_boundaries`
property returns `(time_s, temperature_C)` tuples at each transition point
(including the initial state):

```python
for i, (t, T) in enumerate(prog.segment_boundaries):
    print(f"  Boundary {i}: t = {t/60:.1f} min, T = {T:.0f} deg C")
```

### Temperature at a given time

```python
# Temperature at t = 120 s
T = prog.T_furnace_C(120.0)
```

The lookup is O(log n) in the number of segments, so it is fast even
for complex programs.

## Passing to the simulator

Pass the program via the `temperature_program` parameter.  When provided,
it overrides `T_initial_C`, `T_final_C`, and `heating_rate_C_per_min`:

```python
from ta import TASimulator

sim = TASimulator(
    mechanism="gri30.yaml",
    sample_composition="C2H6:0.7, CH4:0.3",
    condensed_species=["C2H6", "CH4"],
    temperature_program=prog,
    dt=2.0,
)
result = sim.run()
```

The {class}`~ta.simulator.SimulationResult` stores the program in its
`temperature_program` attribute so downstream code can access segment
boundaries for annotation or analysis.

## Validation rules

The constructor validates each segment:

- **At least one segment** is required.
- **Ramp direction**: the sign of `rate_C_per_min` must match the
  direction from the current temperature to `T_target_C` (positive rate
  for heating, negative for cooling).
- **No zero-distance ramps**: if `T_target_C` equals the current
  temperature, use `rate_C_per_min=0` (a hold) instead.
- **Hold duration**: isothermal segments (`rate_C_per_min=0`) require
  `hold_min > 0`.

```python
# This raises ValueError: rate is positive but T_target < T_current
TemperatureProgram(200.0, [
    TemperatureSegment(10.0, 100.0),
])
```

## Full example

```python
from ta import (
    TASimulator,
    TemperatureProgram,
    TemperatureSegment,
    compute_tga, compute_dtg, compute_dta, compute_ms,
)
from ta.plotter import plot_ta_figure

prog = TemperatureProgram(50.0, [
    TemperatureSegment(10.0, 300.0),
    TemperatureSegment(0.0,  300.0, hold_min=5.0),
    TemperatureSegment(5.0,  500.0),
    TemperatureSegment(-10.0, 200.0),
])

result = TASimulator(
    mechanism="gri30.yaml",
    sample_composition="C2H6:0.7, CH4:0.3",
    condensed_species=["C2H6", "CH4"],
    temperature_program=prog,
    dt=2.0,
).run()

tga = compute_tga(result)
dtg = compute_dtg(tga, result.time_s)
dta = compute_dta(result)
ms  = compute_ms(result, top_n=5)

fig, _ = plot_ta_figure(result.furnace_temperature_C, tga, dtg, dta, ms)
fig.savefig("multistage.png", dpi=150)
```

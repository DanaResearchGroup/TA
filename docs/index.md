# TA -- Thermal Analysis Simulator

Simulate simultaneous TGA / DTA / DTG / MS experiments using [Cantera](https://cantera.org).

TA models a sample crucible inside a furnace connected by a heat-transfer wall.
The furnace follows a user-defined temperature program (ramps, isothermal holds,
cooling segments) while the sample heats through the wall, producing realistic
thermal lag.  Chemical reactions in the sample generate the DTA and MS signals;
phase filtering separates condensed-species mass loss (TGA) from volatile
products (MS).

## How it works

```text
Temperature Program  -->  TASimulator  -->  SimulationResult
                                                  |
                            +---------------------+-------------------+
                            |                     |                   |
                       compute_tga()        compute_dta()       compute_ms()
                            |                     |                   |
                       compute_dtg()              |                   |
                            |                     |                   |
                            +---------------------+-------------------+
                                                  |
                                          plot_ta_figure()
```

```{toctree}
:maxdepth: 2
:caption: Contents

installation
quickstart
temperature_programs
api
```

# TA - Thermal Analysis Simulator

Simulate simultaneous TGA/DTA/DTG/MS experiments using [Cantera](https://cantera.org).

## Installation

```bash
conda env create -f environment.yml
conda activate ta_env
pip install -e .
```

## Quick start

```python
from ta import TASimulator, compute_tga, compute_dtg, compute_dta, compute_ms
from ta.plotter import plot_ta_figure

sim = TASimulator(
    mechanism="gri30.yaml",
    sample_composition="C2H6:1",
    condensed_species=["C2H6"],
)
result = sim.run()

tga = compute_tga(result)
dtg = compute_dtg(tga, result.time_s)
dta = compute_dta(result)
ms  = compute_ms(result, top_n=5)

plot_ta_figure(result.temperature_C, tga, dtg, dta, ms)
```

## Running tests

```bash
pytest
```

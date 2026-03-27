# Installation

## Requirements

- Python >= 3.12
- [Cantera](https://cantera.org) >= 3.0 (available via conda-forge)
- NumPy, SciPy, Matplotlib

Cantera is distributed through conda-forge, so a conda (or mamba/micromamba)
environment is recommended.

## Setup

1. Create the conda environment:

   ```bash
   conda env create -f environment.yml
   conda activate ta_env
   ```

2. Install the package in editable mode:

   ```bash
   pip install -e .
   ```

3. Verify the installation:

   ```bash
   pytest
   ```

## Quick check

```python
from ta import TASimulator

sim = TASimulator(mechanism="gri30.yaml", sample_composition="CH4:1")
result = sim.run()
print(f"Recorded {len(result.time_s)} time steps")
```

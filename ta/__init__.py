"""TA: Thermal Analysis simulation package using Cantera."""

from ta.simulator import TASimulator, SimulationResult
from ta.signals import compute_tga, compute_dtg, compute_dta, compute_ms
from ta.plotter import plot_ta_figure

__version__ = "0.1.0"
__all__ = [
    "TASimulator",
    "SimulationResult",
    "compute_tga",
    "compute_dtg",
    "compute_dta",
    "compute_ms",
    "plot_ta_figure",
]

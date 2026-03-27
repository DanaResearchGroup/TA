"""Calculation of TGA, DTG, DTA, and MS signals from simulation data.

Each function takes a :class:`~ta.simulator.SimulationResult` (or the
relevant arrays directly) and returns a NumPy array suitable for
plotting.

Signal definitions
------------------
* **TGA** — sum of condensed-species mass fractions × 100.
* **DTG** — time derivative of TGA [%/min].
* **DTA** — wall heat flux normalised by current condensed mass [W/kg].
* **MS**  — mole fractions of volatile (target gas) species.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ta.simulator import SimulationResult


# ------------------------------------------------------------------
# TGA
# ------------------------------------------------------------------
def compute_tga(
    result: SimulationResult,
    condensed_species: list[str] | None = None,
) -> np.ndarray:
    r"""Mass remaining as a percentage of the condensed-species mass fraction.

    .. math::

        \text{TGA}(t) = \left(\sum_{k \in \text{condensed}} Y_k(t)\right)
                         \times 100

    Parameters
    ----------
    result : SimulationResult
    condensed_species : list[str] or None
        Override the condensed-species list stored in *result*.

    Returns
    -------
    np.ndarray
        TGA signal [%].
    """
    species = condensed_species if condensed_species is not None else result.condensed_species
    n = len(result.time_s)
    Y_cond = np.zeros(n)
    for sp in species:
        if sp in result.mass_fractions:
            Y_cond += result.mass_fractions[sp]
    return Y_cond * 100.0


# ------------------------------------------------------------------
# DTG
# ------------------------------------------------------------------
def compute_dtg(tga: np.ndarray, time_s: np.ndarray) -> np.ndarray:
    """Time derivative of the TGA signal.

    Uses :func:`numpy.gradient` with time converted to **minutes** so
    the result is in %/min.

    Parameters
    ----------
    tga : np.ndarray
        TGA signal [%].
    time_s : np.ndarray
        Time array [s].

    Returns
    -------
    np.ndarray
        DTG signal [%/min].
    """
    time_min = time_s / 60.0
    return np.gradient(tga, time_min)


# ------------------------------------------------------------------
# DTA
# ------------------------------------------------------------------
def compute_dta(
    result: SimulationResult,
    condensed_species: list[str] | None = None,
    sensitivity: float = 1.0,
) -> np.ndarray:
    r"""Differential-thermal-analysis signal.

    The DTA signal is the wall heat-transfer rate normalised by the
    current condensed mass:

    .. math::

        \text{DTA}(t) = \frac{\dot{Q}_\text{wall}(t)}
                             {m_\text{reactor}(t) \cdot \text{TGA}(t)/100}
                         \times S

    where *S* is an optional calibration *sensitivity* factor that
    converts W/kg to µV/mg (instrument output).

    Sign convention: positive DTA = net heat flowing **into** the
    reactor (endothermic event or heating lag).

    Parameters
    ----------
    result : SimulationResult
    condensed_species : list[str] or None
        Override the condensed-species list stored in *result*.
    sensitivity : float
        Calibration constant [µV·kg / (W·mg)] (default 1.0, i.e.
        the output is in W/kg).

    Returns
    -------
    np.ndarray
        DTA signal.
    """
    tga = compute_tga(result, condensed_species)
    condensed_fraction = tga / 100.0
    condensed_mass_kg = result.reactor_mass * condensed_fraction

    with np.errstate(divide="ignore", invalid="ignore"):
        dta = np.where(
            condensed_mass_kg > 1e-30,
            result.wall_heat_flux / condensed_mass_kg,
            0.0,
        )
    return dta * sensitivity


# ------------------------------------------------------------------
# MS
# ------------------------------------------------------------------
def compute_ms(
    result: SimulationResult,
    target_gases: list[str] | None = None,
    top_n: int | None = None,
) -> dict[str, np.ndarray]:
    """Gas-phase mole fractions as a proxy for MS ion intensity.

    Parameters
    ----------
    result : SimulationResult
    target_gases : list[str] or None
        Explicit list of gas species to return.  If *None*, all
        non-condensed, non-bath-gas species with any non-zero mole
        fraction are included.
    top_n : int or None
        If set, keep only the *top_n* species ranked by their
        **peak** mole fraction over the simulation.

    Returns
    -------
    dict[str, np.ndarray]
        ``{species_name: mole_fraction_array}``.
    """
    excluded = set(result.condensed_species) | {result.bath_gas}

    if target_gases is not None:
        candidates = [s for s in target_gases if s in result.mole_fractions]
    else:
        candidates = [s for s in result.species_names if s not in excluded]

    ms: dict[str, np.ndarray] = {}
    for sp in candidates:
        arr = result.mole_fractions[sp]
        if np.any(arr > 0):
            ms[sp] = arr

    if top_n is not None and len(ms) > top_n:
        ranked = sorted(ms, key=lambda s: float(np.max(ms[s])), reverse=True)
        ms = {s: ms[s] for s in ranked[:top_n]}

    return ms

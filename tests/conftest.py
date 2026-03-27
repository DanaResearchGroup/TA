"""Shared fixtures for the TA test suite."""

import numpy as np
import pytest

from ta.simulator import SimulationResult


@pytest.fixture()
def mock_result() -> SimulationResult:
    """A minimal mock SimulationResult with linear condensed-species loss.

    - 100 time steps over 600 s (50-600 deg C furnace ramp).
    - Two species: "COND" (condensed, Y drops 1.0 -> 0.5) and "VOL" (volatile).
    - Constant wall heat flux of 0.1 W.
    - Reactor mass 1e-5 kg (10 mg).
    """
    n = 100
    time_s = np.linspace(0, 600, n)
    condensed_Y = np.linspace(1.0, 0.5, n)
    volatile_Y = 1.0 - condensed_Y

    return SimulationResult(
        time_s=time_s,
        furnace_temperature_C=np.linspace(50, 600, n),
        reactor_temperature_C=np.linspace(50, 599, n),
        reactor_mass=np.full(n, 1e-5),
        mass_fractions={"COND": condensed_Y, "VOL": volatile_Y},
        mole_fractions={"COND": condensed_Y.copy(), "VOL": volatile_Y.copy()},
        wall_heat_flux=np.full(n, 0.1),
        density=np.ones(n),
        species_names=["COND", "VOL"],
        condensed_species=["COND"],
        sample_mass_mg=10.0,
        bath_gas="VOL",
        bath_gas_flow_sccm=50.0,
        pressure_Pa=101325.0,
    )


@pytest.fixture()
def three_species_result() -> SimulationResult:
    """Mock result with three species: two condensed and one volatile."""
    n = 80
    time_s = np.linspace(0, 480, n)
    y_a = np.linspace(0.6, 0.2, n)
    y_b = np.linspace(0.3, 0.1, n)
    y_gas = 1.0 - y_a - y_b

    return SimulationResult(
        time_s=time_s,
        furnace_temperature_C=np.linspace(50, 400, n),
        reactor_temperature_C=np.linspace(50, 399, n),
        reactor_mass=np.full(n, 1e-5),
        mass_fractions={"A": y_a, "B": y_b, "GAS": y_gas},
        mole_fractions={"A": y_a.copy(), "B": y_b.copy(), "GAS": y_gas.copy()},
        wall_heat_flux=np.full(n, 0.05),
        density=np.ones(n),
        species_names=["A", "B", "GAS"],
        condensed_species=["A", "B"],
        sample_mass_mg=10.0,
        bath_gas="GAS",
        bath_gas_flow_sccm=50.0,
        pressure_Pa=101325.0,
    )

"""Tests for signal computation functions."""

import numpy as np
import pytest

from ta.signals import compute_tga, compute_dtg, compute_dta, compute_ms
from ta.simulator import SimulationResult


# ------------------------------------------------------------------
# TGA
# ------------------------------------------------------------------
class TestTGA:
    """Verify TGA signal = sum(Y_condensed) * 100."""

    def test_tga_drops_with_condensed_loss(self, mock_result):
        """If condensed Y drops from 1.0 to 0.5, TGA drops from 100 to 50."""
        tga = compute_tga(mock_result)

        assert tga[0] == pytest.approx(100.0)
        assert tga[-1] == pytest.approx(50.0)

    def test_tga_constant_when_no_loss(self):
        """If condensed Y stays at 1.0, TGA stays at 100."""
        n = 50
        result = SimulationResult(
            time_s=np.linspace(0, 300, n),
            furnace_temperature_C=np.linspace(50, 300, n),
            reactor_temperature_C=np.linspace(50, 299, n),
            reactor_mass=np.full(n, 1e-5),
            mass_fractions={"COND": np.ones(n), "VOL": np.zeros(n)},
            mole_fractions={"COND": np.ones(n), "VOL": np.zeros(n)},
            wall_heat_flux=np.full(n, 0.1),
            density=np.ones(n),
            species_names=["COND", "VOL"],
            condensed_species=["COND"],
            sample_mass_mg=10.0,
            bath_gas="VOL",
            bath_gas_flow_sccm=50.0,
            pressure_Pa=101325.0,
        )
        tga = compute_tga(result)
        np.testing.assert_allclose(tga, 100.0)

    def test_tga_custom_condensed_species(self, mock_result):
        """Passing an explicit condensed_species list overrides the default."""
        tga = compute_tga(mock_result, condensed_species=["VOL"])

        # VOL starts at 0.0 and rises to 0.5
        assert tga[0] == pytest.approx(0.0)
        assert tga[-1] == pytest.approx(50.0)

    def test_tga_multiple_condensed(self, three_species_result):
        """TGA should sum mass fractions from all condensed species."""
        tga = compute_tga(three_species_result)

        # A: 0.6 -> 0.2,  B: 0.3 -> 0.1  => total: 0.9 -> 0.3
        assert tga[0] == pytest.approx(90.0)
        assert tga[-1] == pytest.approx(30.0)

    def test_tga_unknown_species_ignored(self, mock_result):
        """Species not in mass_fractions should be silently skipped."""
        tga = compute_tga(mock_result, condensed_species=["COND", "NONEXISTENT"])

        # Only COND contributes
        assert tga[0] == pytest.approx(100.0)
        assert tga[-1] == pytest.approx(50.0)

    def test_tga_empty_condensed_species(self, mock_result):
        """An empty condensed list should give TGA = 0 everywhere."""
        tga = compute_tga(mock_result, condensed_species=[])
        np.testing.assert_allclose(tga, 0.0)

    def test_tga_range(self, mock_result):
        """TGA should always be between 0 and 100 for valid mass fractions."""
        tga = compute_tga(mock_result)
        assert np.all(tga >= 0.0)
        assert np.all(tga <= 100.0)


# ------------------------------------------------------------------
# DTG
# ------------------------------------------------------------------
class TestDTG:
    """Verify DTG accurately computes the derivative of a TGA array."""

    def test_linear_tga_gives_constant_dtg(self):
        """A linearly decreasing TGA should give a constant DTG."""
        time_s = np.linspace(0, 600, 200)  # 10 minutes
        tga = 100.0 - 5.0 * (time_s / 60.0)  # -5 %/min

        dtg = compute_dtg(tga, time_s)

        np.testing.assert_allclose(
            dtg[5:-5],
            -5.0,
            atol=0.01,
            err_msg="DTG should be constant for a linear TGA",
        )

    def test_constant_tga_gives_zero_dtg(self):
        """No mass change -> DTG = 0."""
        time_s = np.linspace(0, 600, 100)
        tga = np.full(100, 100.0)

        dtg = compute_dtg(tga, time_s)

        np.testing.assert_allclose(dtg, 0.0, atol=1e-10)

    def test_dtg_units_per_minute(self):
        """DTG must be expressed in %/min, not %/s."""
        # TGA drops from 100 -> 0 over 120 s = 2 min  =>  -50 %/min
        time_s = np.linspace(0, 120, 80)
        tga = 100.0 * (1.0 - time_s / 120.0)

        dtg = compute_dtg(tga, time_s)

        np.testing.assert_allclose(
            dtg[5:-5],
            -50.0,
            atol=1.0,
            err_msg="DTG should be in %/min",
        )

    def test_exponential_tga(self):
        """Derivative of an exponential TGA should match the analytical
        result  d/dt [100 * exp(-k*t)] = -100*k*exp(-k*t)  (in %/min)."""
        k = 0.005  # 1/s
        time_s = np.linspace(0, 200, 500)
        tga = 100.0 * np.exp(-k * time_s)

        dtg = compute_dtg(tga, time_s)
        expected = -100.0 * k * 60.0 * np.exp(-k * time_s)  # %/min

        np.testing.assert_allclose(
            dtg[10:-10],
            expected[10:-10],
            rtol=0.02,
            err_msg="DTG of exponential TGA does not match analytical derivative",
        )

    def test_dtg_output_length(self):
        """DTG array must be the same length as input TGA."""
        time_s = np.linspace(0, 100, 50)
        tga = np.linspace(100, 80, 50)
        dtg = compute_dtg(tga, time_s)
        assert len(dtg) == len(tga)


# ------------------------------------------------------------------
# DTA
# ------------------------------------------------------------------
class TestDTA:
    """Verify DTA = wall_heat_flux / condensed_mass."""

    def test_dta_scales_with_heat_flux(self):
        """Doubling wall heat flux should double DTA."""
        n = 50

        def _make(q_val):
            return SimulationResult(
                time_s=np.linspace(0, 300, n),
                furnace_temperature_C=np.linspace(50, 300, n),
                reactor_temperature_C=np.linspace(50, 299, n),
                reactor_mass=np.full(n, 1e-5),
                mass_fractions={"C": np.ones(n)},
                mole_fractions={"C": np.ones(n)},
                wall_heat_flux=np.full(n, q_val),
                density=np.ones(n),
                species_names=["C"],
                condensed_species=["C"],
                sample_mass_mg=10.0,
                bath_gas="NONE",
                bath_gas_flow_sccm=50.0,
                pressure_Pa=101325.0,
            )

        dta_1x = compute_dta(_make(1.0))
        dta_2x = compute_dta(_make(2.0))

        np.testing.assert_allclose(dta_2x, 2.0 * dta_1x, rtol=1e-10)

    def test_dta_zero_when_no_heat_flux(self, mock_result):
        """Zero wall heat flux should give zero DTA."""
        mock_result.wall_heat_flux = np.zeros_like(mock_result.wall_heat_flux)
        dta = compute_dta(mock_result)
        np.testing.assert_allclose(dta, 0.0)

    def test_dta_handles_zero_condensed_mass(self):
        """DTA should not blow up when condensed mass reaches zero."""
        n = 50
        result = SimulationResult(
            time_s=np.linspace(0, 300, n),
            furnace_temperature_C=np.linspace(50, 300, n),
            reactor_temperature_C=np.linspace(50, 299, n),
            reactor_mass=np.full(n, 1e-5),
            mass_fractions={"C": np.zeros(n)},
            mole_fractions={"C": np.zeros(n)},
            wall_heat_flux=np.ones(n),
            density=np.ones(n),
            species_names=["C"],
            condensed_species=["C"],
            sample_mass_mg=10.0,
            bath_gas="NONE",
            bath_gas_flow_sccm=50.0,
            pressure_Pa=101325.0,
        )
        dta = compute_dta(result)
        assert np.all(np.isfinite(dta))

    def test_dta_sensitivity_factor(self, mock_result):
        """The sensitivity parameter should scale DTA linearly."""
        dta_base = compute_dta(mock_result, sensitivity=1.0)
        dta_scaled = compute_dta(mock_result, sensitivity=3.0)

        np.testing.assert_allclose(dta_scaled, 3.0 * dta_base, rtol=1e-10)

    def test_dta_inversely_proportional_to_mass(self):
        """Halving condensed mass should double DTA for same heat flux."""
        n = 50

        def _make(cond_y):
            return SimulationResult(
                time_s=np.linspace(0, 300, n),
                furnace_temperature_C=np.linspace(50, 300, n),
                reactor_temperature_C=np.linspace(50, 299, n),
                reactor_mass=np.full(n, 1e-5),
                mass_fractions={"C": np.full(n, cond_y)},
                mole_fractions={"C": np.full(n, cond_y)},
                wall_heat_flux=np.ones(n),
                density=np.ones(n),
                species_names=["C"],
                condensed_species=["C"],
                sample_mass_mg=10.0,
                bath_gas="NONE",
                bath_gas_flow_sccm=50.0,
                pressure_Pa=101325.0,
            )

        dta_full = compute_dta(_make(1.0))
        dta_half = compute_dta(_make(0.5))

        np.testing.assert_allclose(dta_half, 2.0 * dta_full, rtol=1e-10)


# ------------------------------------------------------------------
# MS
# ------------------------------------------------------------------
class TestMS:
    """Verify MS returns mole fractions for volatile species."""

    def test_ms_excludes_condensed_and_bath(self, mock_result):
        """MS should not include condensed species or bath gas."""
        ms = compute_ms(mock_result)
        assert "COND" not in ms
        assert "VOL" not in ms  # VOL is the bath gas in the mock

    def test_ms_explicit_target_gases(self, mock_result):
        """Explicit target_gases should be returned even if they are condensed."""
        ms = compute_ms(mock_result, target_gases=["COND"])
        assert "COND" in ms

    def test_ms_top_n(self):
        """top_n should limit the number of returned species."""
        n = 50
        fracs = {}
        mole_fracs = {}
        for i in range(10):
            name = f"SP{i}"
            fracs[name] = np.full(n, 0.01 * (i + 1))
            mole_fracs[name] = np.full(n, 0.01 * (i + 1))

        result = SimulationResult(
            time_s=np.linspace(0, 100, n),
            furnace_temperature_C=np.linspace(50, 200, n),
            reactor_temperature_C=np.linspace(50, 199, n),
            reactor_mass=np.full(n, 1e-5),
            mass_fractions=fracs,
            mole_fractions=mole_fracs,
            wall_heat_flux=np.zeros(n),
            density=np.ones(n),
            species_names=[f"SP{i}" for i in range(10)],
            condensed_species=[],
            sample_mass_mg=10.0,
            bath_gas="NONE",
            bath_gas_flow_sccm=50.0,
            pressure_Pa=101325.0,
        )
        ms = compute_ms(result, top_n=3)
        assert len(ms) == 3

    def test_ms_top_n_selects_highest_peak(self):
        """top_n should keep the species with the largest peak mole fraction."""
        n = 50
        mole_fracs = {
            "LOW": np.full(n, 0.01),
            "MID": np.full(n, 0.05),
            "HIGH": np.full(n, 0.10),
        }
        result = SimulationResult(
            time_s=np.linspace(0, 100, n),
            furnace_temperature_C=np.linspace(50, 200, n),
            reactor_temperature_C=np.linspace(50, 199, n),
            reactor_mass=np.full(n, 1e-5),
            mass_fractions=mole_fracs,
            mole_fractions=mole_fracs,
            wall_heat_flux=np.zeros(n),
            density=np.ones(n),
            species_names=["LOW", "MID", "HIGH"],
            condensed_species=[],
            sample_mass_mg=10.0,
            bath_gas="NONE",
            bath_gas_flow_sccm=50.0,
            pressure_Pa=101325.0,
        )
        ms = compute_ms(result, top_n=1)
        assert list(ms.keys()) == ["HIGH"]

    def test_ms_excludes_zero_mole_fraction(self):
        """Species with zero mole fraction everywhere should be excluded."""
        n = 30
        result = SimulationResult(
            time_s=np.linspace(0, 100, n),
            furnace_temperature_C=np.linspace(50, 200, n),
            reactor_temperature_C=np.linspace(50, 199, n),
            reactor_mass=np.full(n, 1e-5),
            mass_fractions={"A": np.zeros(n), "B": np.full(n, 0.01)},
            mole_fractions={"A": np.zeros(n), "B": np.full(n, 0.01)},
            wall_heat_flux=np.zeros(n),
            density=np.ones(n),
            species_names=["A", "B"],
            condensed_species=[],
            sample_mass_mg=10.0,
            bath_gas="NONE",
            bath_gas_flow_sccm=50.0,
            pressure_Pa=101325.0,
        )
        ms = compute_ms(result)
        assert "A" not in ms
        assert "B" in ms

    def test_ms_empty_when_all_excluded(self, mock_result):
        """If all species are condensed or bath gas, MS should be empty."""
        ms = compute_ms(mock_result)
        assert ms == {}

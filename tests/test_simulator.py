"""Tests for the TASimulator class."""

import numpy as np
import pytest

ct = pytest.importorskip("cantera")

from ta.simulator import TASimulator, SimulationResult


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
_FAST_COMMON = dict(
    mechanism="gri30.yaml",
    T_initial_C=50.0,
    T_final_C=200.0,
    heating_rate_C_per_min=10.0,
    dt=2.0,
)


# ------------------------------------------------------------------
# Thermal-lag tests
# ------------------------------------------------------------------
class TestThermalLag:
    """Verify that the reactor temperature lags behind the furnace."""

    def test_reactor_lags_furnace_during_heating(self):
        """With finite wall U·A the sample must lag the furnace ramp."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
            wall_area=1e-4,
            htc=100.0,
        )
        result = sim.run()

        mid = len(result.time_s) // 2
        lag = result.furnace_temperature_C[mid] - result.reactor_temperature_C[mid]
        assert lag > 0, (
            f"Expected reactor to lag furnace, but "
            f"T_furnace={result.furnace_temperature_C[mid]:.4f}, "
            f"T_reactor={result.reactor_temperature_C[mid]:.4f}"
        )

    def test_larger_htc_reduces_lag(self):
        """Higher heat-transfer coefficient should reduce thermal lag."""
        common = dict(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
            wall_area=1e-4,
        )

        result_low = TASimulator(**common, htc=50.0).run()
        result_high = TASimulator(**common, htc=500.0).run()

        mid = min(len(result_low.time_s), len(result_high.time_s)) // 2
        lag_low = (
            result_low.furnace_temperature_C[mid]
            - result_low.reactor_temperature_C[mid]
        )
        lag_high = (
            result_high.furnace_temperature_C[mid]
            - result_high.reactor_temperature_C[mid]
        )
        assert lag_high < lag_low, (
            f"Higher U should reduce lag: lag_high={lag_high:.4f}, "
            f"lag_low={lag_low:.4f}"
        )

    def test_wall_heat_flux_positive_during_heating(self):
        """During heating, net heat should flow into the reactor (positive)."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
        )
        result = sim.run()

        # Skip t=0 (furnace and reactor at same temperature)
        assert np.all(result.wall_heat_flux[1:] > 0), (
            "Wall heat flux should be positive when furnace leads reactor"
        )


# ------------------------------------------------------------------
# Output-consistency tests
# ------------------------------------------------------------------
class TestOutputConsistency:
    """Verify that the simulator produces arrays of consistent length."""

    def test_array_lengths_match(self):
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="CH4:1",
            condensed_species=["CH4"],
        )
        result = sim.run()

        n = len(result.time_s)
        assert n > 0, "No data points produced"
        assert len(result.furnace_temperature_C) == n
        assert len(result.reactor_temperature_C) == n
        assert len(result.reactor_mass) == n
        assert len(result.wall_heat_flux) == n
        assert len(result.density) == n

        for sp in result.species_names:
            assert len(result.mass_fractions[sp]) == n, (
                f"mass_fractions['{sp}'] length mismatch"
            )
            assert len(result.mole_fractions[sp]) == n, (
                f"mole_fractions['{sp}'] length mismatch"
            )

    def test_result_is_dataclass(self):
        sim = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="AR:1",
            T_initial_C=50.0,
            T_final_C=60.0,
            heating_rate_C_per_min=30.0,
            bath_gas="AR",
            dt=5.0,
        )
        result = sim.run()
        assert isinstance(result, SimulationResult)

    def test_mass_is_conserved(self):
        """In a closed ConstPressureReactor the total mass must not change."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="CH4:1",
            condensed_species=["CH4"],
        )
        result = sim.run()

        np.testing.assert_allclose(
            result.reactor_mass,
            result.reactor_mass[0],
            rtol=1e-6,
            err_msg="Reactor mass should be conserved in a closed system",
        )

    def test_mass_fractions_sum_to_one(self):
        """Mass fractions across all species must sum to 1 at every step."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="C2H6:0.5, CH4:0.5",
            condensed_species=["C2H6", "CH4"],
        )
        result = sim.run()

        Y_sum = sum(result.mass_fractions[sp] for sp in result.species_names)
        np.testing.assert_allclose(Y_sum, 1.0, atol=1e-10)

    def test_mole_fractions_sum_to_one(self):
        """Mole fractions across all species must sum to 1 at every step."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="CH4:1",
            condensed_species=["CH4"],
        )
        result = sim.run()

        X_sum = sum(result.mole_fractions[sp] for sp in result.species_names)
        np.testing.assert_allclose(X_sum, 1.0, atol=1e-10)

    def test_furnace_temperature_monotonic(self):
        """Furnace temperature must increase monotonically for a positive ramp."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
        )
        result = sim.run()

        diffs = np.diff(result.furnace_temperature_C)
        assert np.all(diffs >= 0), "Furnace temperature should be monotonically increasing"

    def test_initial_temperatures_match(self):
        """At t=0, furnace and reactor should both be at T_initial."""
        T_init = 100.0
        sim = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="AR:1",
            bath_gas="AR",
            T_initial_C=T_init,
            T_final_C=200.0,
            heating_rate_C_per_min=20.0,
            dt=5.0,
        )
        result = sim.run()

        assert result.furnace_temperature_C[0] == pytest.approx(T_init, abs=0.01)
        assert result.reactor_temperature_C[0] == pytest.approx(T_init, abs=0.01)

    def test_time_starts_at_zero(self):
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
        )
        result = sim.run()

        assert result.time_s[0] == pytest.approx(0.0)


# ------------------------------------------------------------------
# Composition input formats
# ------------------------------------------------------------------
class TestCompositionInput:
    """Verify that different composition formats produce valid results."""

    def test_string_composition(self):
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="CH4:0.8, C2H6:0.2",
            condensed_species=["CH4"],
        )
        result = sim.run()
        assert len(result.time_s) > 1

    def test_dict_composition(self):
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition={"CH4": 0.8, "C2H6": 0.2},
            condensed_species=["CH4"],
        )
        result = sim.run()
        assert len(result.time_s) > 1


# ------------------------------------------------------------------
# Auto-detection of condensed species
# ------------------------------------------------------------------
class TestCondensedSpeciesAutoDetect:
    """Verify auto-detection when condensed_species=None."""

    def test_auto_detects_nonzero_non_bath(self):
        """All species with Y>0 except bath gas should be condensed."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="CH4:0.7, C2H6:0.3",
            bath_gas="AR",
            condensed_species=None,
        )
        assert "CH4" in sim.condensed_species
        assert "C2H6" in sim.condensed_species
        assert "AR" not in sim.condensed_species

    def test_auto_detect_excludes_bath_gas(self):
        """Bath gas should never appear in auto-detected condensed list."""
        sim = TASimulator(
            **_FAST_COMMON,
            sample_composition="AR:1",
            bath_gas="AR",
            condensed_species=None,
        )
        assert "AR" not in sim.condensed_species
        assert len(sim.condensed_species) == 0


# ------------------------------------------------------------------
# Metadata propagation
# ------------------------------------------------------------------
class TestMetadata:
    """Verify that metadata fields are correctly propagated to SimulationResult."""

    def test_metadata_fields(self):
        sim = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="CH4:1",
            sample_mass_mg=5.0,
            bath_gas="N2",
            bath_gas_flow_sccm=25.0,
            T_initial_C=50.0,
            T_final_C=100.0,
            heating_rate_C_per_min=20.0,
            pressure_atm=2.0,
            condensed_species=["CH4"],
            dt=5.0,
        )
        result = sim.run()

        assert result.sample_mass_mg == 5.0
        assert result.bath_gas == "N2"
        assert result.bath_gas_flow_sccm == 25.0
        assert result.pressure_Pa == pytest.approx(2.0 * ct.one_atm)
        assert result.condensed_species == ["CH4"]
        assert "CH4" in result.species_names


# ------------------------------------------------------------------
# No-crash / reactive tests
# ------------------------------------------------------------------
class TestNoCrash:
    """Verify that the simulator completes for reactive mixtures."""

    def test_reactive_mixture(self):
        sim = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="C2H6:0.5, CH4:0.5",
            T_initial_C=100.0,
            T_final_C=300.0,
            heating_rate_C_per_min=20.0,
            condensed_species=["C2H6", "CH4"],
            dt=5.0,
        )
        result = sim.run()

        assert len(result.time_s) > 0
        assert result.furnace_temperature_C[0] == pytest.approx(100.0, abs=1.0)
        assert result.furnace_temperature_C[-1] <= 300.0 + 1.0


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------
class TestValidation:
    """Input validation."""

    def test_bad_bath_gas_raises(self):
        with pytest.raises(ValueError, match="Bath gas"):
            TASimulator(
                mechanism="gri30.yaml",
                sample_composition="CH4:1",
                bath_gas="KRYPTONITE",
            )

    def test_bad_condensed_species_raises(self):
        with pytest.raises(ValueError, match="Condensed species"):
            TASimulator(
                mechanism="gri30.yaml",
                sample_composition="CH4:1",
                condensed_species=["UNOBTANIUM"],
            )

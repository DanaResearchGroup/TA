"""Integration tests: full simulate -> signals -> plot pipeline."""

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ct = pytest.importorskip("cantera")

from ta.simulator import TASimulator
from ta.signals import compute_tga, compute_dtg, compute_dta, compute_ms
from ta.plotter import plot_ta_figure
from ta.temperature_program import TemperatureProgram, TemperatureSegment


class TestFullPipeline:
    """Run the full workflow and verify outputs are self-consistent."""

    def test_inert_gas_pipeline(self):
        """Inert gas (AR) should show no mass loss and no DTA events."""
        result = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="AR:1",
            bath_gas="AR",
            T_initial_C=50.0,
            T_final_C=200.0,
            heating_rate_C_per_min=20.0,
            condensed_species=[],
            dt=5.0,
        ).run()

        tga = compute_tga(result)
        dtg = compute_dtg(tga, result.time_s)
        dta = compute_dta(result)
        ms = compute_ms(result, top_n=3)

        # TGA should be 0 (no condensed species)
        np.testing.assert_allclose(tga, 0.0)
        # DTG should be 0
        np.testing.assert_allclose(dtg, 0.0, atol=1e-10)
        # DTA should be 0 (no condensed mass to normalise by)
        np.testing.assert_allclose(dta, 0.0)

        # Plot should not raise
        fig, axes = plot_ta_figure(
            result.furnace_temperature_C, tga, dtg, dta, ms,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_reactive_pipeline(self):
        """Reactive mixture should produce non-trivial signals."""
        result = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="C2H6:0.7, CH4:0.3",
            T_initial_C=50.0,
            T_final_C=400.0,
            heating_rate_C_per_min=20.0,
            condensed_species=["C2H6", "CH4"],
            dt=5.0,
        ).run()

        tga = compute_tga(result)
        dtg = compute_dtg(tga, result.time_s)
        dta = compute_dta(result)
        ms = compute_ms(result, top_n=5)

        # Basic shape checks
        assert len(tga) == len(result.time_s)
        assert len(dtg) == len(result.time_s)
        assert len(dta) == len(result.time_s)

        # TGA starts at 100% (all condensed)
        assert tga[0] == pytest.approx(100.0, abs=0.1)

        # All signals finite
        assert np.all(np.isfinite(tga))
        assert np.all(np.isfinite(dtg))
        assert np.all(np.isfinite(dta))
        for arr in ms.values():
            assert np.all(np.isfinite(arr))

        # Plot should not raise
        fig, axes = plot_ta_figure(
            result.furnace_temperature_C, tga, dtg, dta, ms,
        )
        assert isinstance(fig, plt.Figure)
        assert len(axes) == 4
        plt.close(fig)

    def test_signals_use_correct_result_fields(self):
        """Verify that signals consume the expected SimulationResult fields."""
        result = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="CH4:1",
            T_initial_C=50.0,
            T_final_C=150.0,
            heating_rate_C_per_min=30.0,
            condensed_species=["CH4"],
            dt=5.0,
        ).run()

        tga = compute_tga(result)
        dtg = compute_dtg(tga, result.time_s)
        dta = compute_dta(result)
        ms = compute_ms(result, target_gases=["H2O", "CO2"])

        # TGA should reference mass_fractions
        assert tga[0] == pytest.approx(
            result.mass_fractions["CH4"][0] * 100.0
        )

        # MS should return dict keyed by species name
        assert isinstance(ms, dict)
        for sp, arr in ms.items():
            assert sp in result.species_names
            np.testing.assert_array_equal(arr, result.mole_fractions[sp])

    def test_multistage_pipeline(self):
        """Ramp-hold-ramp program with reactive mixture produces valid signals."""
        prog = TemperatureProgram(50.0, [
            TemperatureSegment(20.0, 200.0),
            TemperatureSegment(0.0, 200.0, hold_min=2.0),
            TemperatureSegment(10.0, 400.0),
        ])
        result = TASimulator(
            mechanism="gri30.yaml",
            sample_composition="C2H6:0.7, CH4:0.3",
            condensed_species=["C2H6", "CH4"],
            dt=5.0,
            temperature_program=prog,
        ).run()

        tga = compute_tga(result)
        dtg = compute_dtg(tga, result.time_s)
        dta = compute_dta(result)
        ms = compute_ms(result, top_n=5)

        assert len(tga) == len(result.time_s)
        assert tga[0] == pytest.approx(100.0, abs=0.1)
        assert np.all(np.isfinite(tga))
        assert np.all(np.isfinite(dtg))
        assert np.all(np.isfinite(dta))

        fig, axes = plot_ta_figure(
            result.furnace_temperature_C, tga, dtg, dta, ms,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

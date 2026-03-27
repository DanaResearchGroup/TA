"""Tests for the plotting module."""

import numpy as np
import pytest

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for CI / headless tests
import matplotlib.pyplot as plt

from ta.plotter import plot_ta_figure


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def _synthetic_data(n=100):
    """Return (temperature_C, tga, dtg, dta, ms) for plotting tests."""
    T = np.linspace(50, 600, n)
    tga = 100.0 - 50.0 * (T - 50.0) / 550.0
    dtg = np.gradient(tga, T)
    dta = np.sin((T - 50.0) / 550.0 * np.pi) * 10.0
    rng = np.random.default_rng(42)
    ms = {"H2O": np.abs(rng.normal(1.0, 0.3, n))}
    return T, tga, dtg, dta, ms


# ------------------------------------------------------------------
# Basic creation
# ------------------------------------------------------------------
class TestPlotTAFigure:
    """Verify that plot_ta_figure creates a valid figure."""

    def test_basic_plot(self):
        """The function should return a Figure with four Axes."""
        T, tga, dtg, dta, ms = _synthetic_data()

        fig, axes = plot_ta_figure(T, tga, dtg, dta, ms)

        assert isinstance(fig, plt.Figure)
        assert len(axes) == 4
        plt.close(fig)

    def test_empty_ms_dict(self):
        """An empty MS dict should not raise."""
        T, tga, dtg, dta, _ = _synthetic_data()

        fig, axes = plot_ta_figure(T, tga, dtg, dta, {})

        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_multiple_ms_species(self):
        """Should handle several MS traces without error."""
        T, tga, dtg, dta, _ = _synthetic_data(80)
        rng = np.random.default_rng(0)
        ms = {f"SP{i}": np.abs(rng.normal(1, 0.5, 80)) for i in range(5)}

        fig, _ = plot_ta_figure(T, tga, dtg, dta, ms)
        plt.close(fig)


# ------------------------------------------------------------------
# Axis labels and title
# ------------------------------------------------------------------
class TestPlotLabels:
    """Verify that labels and title are applied correctly."""

    def test_default_title(self):
        T, tga, dtg, dta, ms = _synthetic_data(50)
        fig, (ax_dta, *_) = plot_ta_figure(T, tga, dtg, dta, ms)

        assert ax_dta.get_title() == "Thermal Analysis"
        plt.close(fig)

    def test_custom_title(self):
        T, tga, dtg, dta, ms = _synthetic_data(50)
        fig, (ax_dta, *_) = plot_ta_figure(
            T, tga, dtg, dta, ms, title="Custom Title"
        )

        assert ax_dta.get_title() == "Custom Title"
        plt.close(fig)

    def test_custom_axis_labels(self):
        T, tga, dtg, dta, ms = _synthetic_data(50)
        fig, (ax_dta, ax_tga, ax_dtg, ax_ms) = plot_ta_figure(
            T, tga, dtg, dta, ms,
            dta_label="My DTA",
            tga_label="My TGA",
            dtg_label="My DTG",
            ms_label="My MS",
        )

        assert ax_dta.get_ylabel() == "My DTA"
        assert ax_tga.get_ylabel() == "My TGA"
        assert ax_dtg.get_ylabel() == "My DTG"
        assert ax_ms.get_ylabel() == "My MS"
        plt.close(fig)


# ------------------------------------------------------------------
# Custom figsize
# ------------------------------------------------------------------
class TestPlotFigsize:
    def test_custom_figsize(self):
        T, tga, dtg, dta, ms = _synthetic_data(50)
        fig, _ = plot_ta_figure(T, tga, dtg, dta, ms, figsize=(12, 6))

        w, h = fig.get_size_inches()
        assert w == pytest.approx(12.0)
        assert h == pytest.approx(6.0)
        plt.close(fig)

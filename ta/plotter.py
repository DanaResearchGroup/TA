"""Matplotlib routines for thermal-analysis figure generation.

The main entry point is :func:`plot_ta_figure`, which produces a
two-panel plot:

* **Panel A** — DTA (left Y), TGA and DTG (twin right Y-axes) vs.
  temperature.
* **Panel B** — MS ion intensities (log scale) vs. temperature.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np


def plot_ta_figure(
    temperature_C: np.ndarray,
    tga: np.ndarray,
    dtg: np.ndarray,
    dta: np.ndarray,
    ms_dict: dict[str, np.ndarray],
    *,
    figsize: tuple[float, float] = (10, 8),
    dta_label: str = r"DTA ($\mu$V/mg)",
    tga_label: str = "TGA (mass %)",
    dtg_label: str = "DTG (%/min)",
    ms_label: str = "Ion intensity (a.u.)",
    title: str = "Thermal Analysis",
) -> tuple[plt.Figure, tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]]:
    """Create a two-panel thermal-analysis figure.

    Parameters
    ----------
    temperature_C : np.ndarray
        Temperature [deg C] (shared X-axis for both panels).
    tga, dtg, dta : np.ndarray
        Signals computed by :mod:`ta.signals`.
    ms_dict : dict[str, np.ndarray]
        Species name -> mole-fraction array for the MS panel.
    figsize : tuple[float, float]
        Figure dimensions in inches.
    dta_label, tga_label, dtg_label, ms_label : str
        Y-axis labels.
    title : str
        Figure super-title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    axes : tuple[Axes, Axes, Axes, Axes]
        ``(ax_dta, ax_tga, ax_dtg, ax_ms)``
    """
    fig, (ax_top, ax_ms) = plt.subplots(
        2, 1, figsize=figsize, sharex=True,
    )

    # ---- colours ----
    c_dta = "tab:blue"
    c_tga = "tab:red"
    c_dtg = "tab:green"

    # ==== Panel A: DTA / TGA / DTG ====
    ax_dta: plt.Axes = ax_top

    # DTA — left axis
    ax_dta.plot(temperature_C, dta, color=c_dta, label="DTA")
    ax_dta.set_ylabel(dta_label, color=c_dta)
    ax_dta.tick_params(axis="y", labelcolor=c_dta)

    # TGA — first right axis
    ax_tga: plt.Axes = ax_dta.twinx()
    ax_tga.plot(temperature_C, tga, color=c_tga, label="TGA")
    ax_tga.set_ylabel(tga_label, color=c_tga)
    ax_tga.tick_params(axis="y", labelcolor=c_tga)

    # DTG — second right axis (offset outward)
    ax_dtg: plt.Axes = ax_dta.twinx()
    ax_dtg.spines["right"].set_position(("outward", 60))
    ax_dtg.plot(temperature_C, dtg, color=c_dtg, label="DTG")
    ax_dtg.set_ylabel(dtg_label, color=c_dtg)
    ax_dtg.tick_params(axis="y", labelcolor=c_dtg)

    # combined legend
    all_lines = (
        ax_dta.get_lines() + ax_tga.get_lines() + ax_dtg.get_lines()
    )
    all_labels = [ln.get_label() for ln in all_lines]
    ax_dta.legend(all_lines, all_labels, loc="best")

    ax_dta.set_title(title)

    # ==== Panel B: MS ====
    if ms_dict:
        for name, intensity in ms_dict.items():
            ax_ms.plot(temperature_C, intensity, label=name)
        ax_ms.set_yscale("log")
        ax_ms.set_ylim(1e-2, 1e2)
        ax_ms.legend(loc="best", fontsize="small", ncol=2)

    ax_ms.set_xlabel(r"Temperature ($\degree$C)")
    ax_ms.set_ylabel(ms_label)

    fig.tight_layout()
    return fig, (ax_dta, ax_tga, ax_dtg, ax_ms)
